"""
scripts/run_rdp5.py
===================
Wrapper around RDP5CL.exe for both:

  (A) Snakemake batch mode  — called by the `rdp5_recombination` rule
      Input:  a per-genotype aligned FASTA  (results/{virus}/alignments/{genotype}_aligned.fasta)
      Output: results/{virus}/recombination/{genotype}_recombinants.csv   (raw RDP5 CSV)
              results/{virus}/recombination/{genotype}_recombinants.tsv   (normalised, pipeline-ready)

  (B) Dashboard user-sequence mode — called by pipeline_runner.py
      Input:  jobs/<job_id>/input.fasta
      Output: <--outdir>/recombinants.tsv  (same normalised format)

CLI
---
  # Batch / Snakemake mode
  python scripts/run_rdp5.py \
      --fasta  results/hbv/alignments/A_aligned.fasta \
      --virus  hbv \
      --outdir results/hbv/recombination \
      --genotype A

  # Dashboard / user-submission mode
  python scripts/run_rdp5.py \
      --fasta  jobs/<id>/input.fasta \
      --virus  hbv \
      --outdir jobs/<id>/rdp5_out

Output TSV columns (normalised)
--------------------------------
  sequence_id       Accession / sequence header from the input FASTA
  is_recombinant    true | false
  breakpoint_start  Integer nt position (0 if none)
  breakpoint_end    Integer nt position (0 if none)
  p_value           Best p-value across all RDP5 methods (empty if none)
  methods           Pipe-separated list of methods that detected recombination
  parent_1          Inferred major parent (empty if not recombinant)
  parent_2          Inferred minor parent (empty if not recombinant)

Configuration (environment variables)
--------------------------------------
  RDP5_EXE      Full path to RDP5CL.exe
                Default: looks for RDP5CL.exe next to this script,
                         then in PATH, then in ~/RDP5/
  RDP5_WINE     Wine executable to use on Linux/macOS
                Default: wine
  RDP5_TIMEOUT  Seconds to wait for RDP5 to finish  (default: 3600)
  RDP5_MIN_SEQS Minimum sequences needed to run RDP5 (default: 4)
                Alignments smaller than this are skipped gracefully.
"""

import argparse
import csv
import logging
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("run_rdp5")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

def _env(key, default):
    return os.environ.get(key, default)

RDP5_TIMEOUT  = int(_env("RDP5_TIMEOUT",  "3600"))
RDP5_MIN_SEQS = int(_env("RDP5_MIN_SEQS", "4"))
WINE          = _env("RDP5_WINE", "wine")


def _find_rdp5_exe() -> Path:
    """
    Locate RDP5CL.exe in order of preference:
      1. RDP5_EXE environment variable
      2. PlayOnMac / Wine installation paths (contains all helper exe/dll files)
      3. Alongside this script  (scripts/RDP5CL.exe)
      4. In PATH
      5. ~/RDP5/RDP5CL.exe
    """
    env_path = os.environ.get("RDP5_EXE")
    if env_path:
        p = Path(env_path)
        if p.exists():
            return p
        raise FileNotFoundError(f"RDP5_EXE is set but not found: {p}")

    candidates = []

    # On macOS/Linux, check PlayOnMac or Wine prefixes first to use the full installation
    if not _on_windows():
        playonmac_prefix = Path.home() / "Library" / "PlayOnMac" / "wineprefix" / "RDP5.84_"
        if playonmac_prefix.exists():
            candidates.extend([
                playonmac_prefix / "drive_c" / "Program Files" / "RDP5" / "RDP5CL.exe",
                playonmac_prefix / "drive_c" / "Program Files (x86)" / "RDP5" / "RDP5CL.exe",
            ])
        wine_prefix = Path.home() / ".wine"
        candidates.extend([
            wine_prefix / "drive_c" / "Program Files" / "RDP5" / "RDP5CL.exe",
            wine_prefix / "drive_c" / "Program Files (x86)" / "RDP5" / "RDP5CL.exe",
        ])

    candidates.extend([
        Path(__file__).parent / "RDP5CL.exe",
        Path(__file__).parent / "RDP5" / "RDP5CL.exe",
        Path.home() / "RDP5" / "RDP5CL.exe",
        Path.home() / "rdp5" / "RDP5CL.exe",
    ])

    for c in candidates:
        if c.exists():
            return c

    # Try PATH (on Windows, .exe is fine; on Linux via Wine, we call wine <path>)
    in_path = shutil.which("RDP5CL.exe") or shutil.which("RDP5CL")
    if in_path:
        return Path(in_path)

    raise FileNotFoundError(
        "Cannot find RDP5CL.exe.\n"
        "Set the RDP5_EXE environment variable to its full path:\n"
        "  export RDP5_EXE=/path/to/RDP5CL.exe\n"
        "Or place RDP5CL.exe in the same directory as this script."
    )


def _on_windows() -> bool:
    return platform.system() == "Windows"


def _ensure_local_wine_image() -> str:
    image_name = "local-wine:latest"
    # Check if image exists
    res = subprocess.run(["docker", "images", "-q", image_name], capture_output=True, text=True)
    if res.stdout.strip():
        return image_name
    
    log.info("Building local-wine Docker image (this runs once and is very lightweight)...")
    dockerfile_content = """FROM debian:stable-slim
RUN dpkg --add-architecture i386 && \\
    apt-get update && \\
    apt-get install -y --no-install-recommends wine wine32 && \\
    rm -rf /var/lib/apt/lists/*
"""
    # Create the tmp directory inside the project root (which starts with /home/...)
    project_root = Path(__file__).parent.parent.resolve()
    tmp_base = project_root / ".tmp_docker_build"
    tmp_base.mkdir(exist_ok=True)
    
    try:
        with tempfile.TemporaryDirectory(dir=str(tmp_base)) as tmpdir:
            df_path = Path(tmpdir) / "Dockerfile"
            df_path.write_text(dockerfile_content)
            build_res = subprocess.run(
                ["docker", "build", "-t", image_name, tmpdir],
                capture_output=True,
                text=True
            )
            if build_res.returncode != 0:
                log.error("Failed to build local-wine image: %s", build_res.stderr)
                raise RuntimeError(f"Failed to build local-wine image: {build_res.stderr}")
    finally:
        shutil.rmtree(tmp_base, ignore_errors=True)
    log.info("Successfully built local-wine Docker image.")
    return image_name


def _build_cmd(rdp5_exe: Path, fasta_path: Path, out_prefix: Path) -> list[str]:
    """
    Build the RDP5CL.exe command line.

    RDP5 writes its output files using the input filename as the base,
    so we give it the FASTA in the output dir to keep everything together.

    Flags used:
      -nor   Skip writing the large .rdp5 project file (we only need the CSV)
    """
    exe_str   = str(rdp5_exe)
    fasta_str = str(fasta_path)

    if _on_windows():
        return [exe_str, f"-f{fasta_str}", "-nor"]
    else:
        # Linux / macOS — run through Wine
        wine_exe = None
        
        # On macOS, check if PlayOnMac optimized wine32on64 is available
        if platform.system() == "Darwin":
            playonmac_wine = Path("/Applications/PlayOnMac.app/Contents/Resources/unix/wine/bin/wine32on64")
            if playonmac_wine.exists():
                wine_exe = str(playonmac_wine)
                # Set WINEPREFIX to PlayOnMac prefix to avoid compatibility issues
                playonmac_prefix = Path.home() / "Library" / "PlayOnMac" / "wineprefix" / "RDP5.84_"
                if playonmac_prefix.exists():
                    os.environ["WINEPREFIX"] = str(playonmac_prefix)
                    log.info("Using PlayOnMac optimized WINEPREFIX: %s", playonmac_prefix)
                os.environ["WINEDEBUG"] = "-all"
                log.info("Using PlayOnMac optimized wine32on64: %s", wine_exe)

        if not wine_exe:
            wine_exe = shutil.which(WINE)
            
        if wine_exe:
            return [wine_exe, exe_str, f"-f{fasta_str}", "-nor"]
            
        # Fallback to Docker
        docker_exe = shutil.which("docker")
        if docker_exe:
            log.info("Wine not found on host. Attempting to run via Docker container...")
            image_name = _ensure_local_wine_image()
            
            rdp5_dir = rdp5_exe.parent.resolve()
            exe_rel = rdp5_exe.name
            fasta_rel = fasta_path.name
            
            uid = os.getuid() if hasattr(os, "getuid") else 0
            gid = os.getgid() if hasattr(os, "getgid") else 0
            
            return [
                docker_exe, "run", "--rm",
                "-v", f"{rdp5_dir}:/work",
                "-w", "/work",
                "-e", "WINEDEBUG=-all",
                "-u", f"{uid}:{gid}",
                image_name,
                "sh", "-c", f"mkdir -p '/tmp/wineprefix/drive_c/Program Files/RDP5' && cp RDP.ini PairsScores BinProbs '/tmp/wineprefix/drive_c/Program Files/RDP5/' && WINEPREFIX=/tmp/wineprefix wine {exe_rel} -f{fasta_rel} -nor"
            ]
            
        raise EnvironmentError(
            f"Wine not found (looked for '{WINE}') and Docker is not available. "
            "Install Wine or Docker to run RDP5CL.exe on Linux:\n"
            "  sudo apt install wine  OR  install docker"
        )



# ---------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------

def _count_sequences(fasta_path: Path) -> int:
    count = 0
    with fasta_path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
    return count


def _read_sequence_ids(fasta_path: Path) -> list[str]:
    ids = []
    with fasta_path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0].strip())
    return ids


# ---------------------------------------------------------------------------
# RDP5 CSV parser
# ---------------------------------------------------------------------------

# RDP5 CSV column names vary slightly between versions.
# We normalise them here.
_RECOMB_SEQ_COLS   = ["Recombinant", "RecombinantSequence", "Sequence", "Recombinant Sequence(s)"]
_PARENT1_COLS      = ["MajorParent",  "Parent1", "ParentA", "Major Parental Sequence(s)"]
_PARENT2_COLS      = ["MinorParent",  "Parent2", "ParentB", "Minor Parental Sequence(s)"]
_BP_START_COLS     = ["BreakpointBegin", "BP_Start", "Breakpoint1Start", "Start1", "Begin"]
_BP_END_COLS       = ["BreakpointEnd",   "BP_End",   "Breakpoint1End",   "End1", "End"]
_PVALUE_COLS       = ["Pvalue", "P-value", "BestPvalue", "Best_pvalue"]
_METHODS_COLS      = ["Methods", "DetectingMethods", "Method"]


def _first_col(row: dict, candidates: list[str]) -> str:
    for c in candidates:
        if c in row and row[c] not in ("", None):
            return str(row[c]).strip()
    return ""


def _parse_rdp5_csv(csv_path: Path, all_seq_ids: list[str]) -> list[dict]:
    """
    Parse RDP5's summary CSV and return a list of normalised dicts,
    one per sequence in the original input.

    RDP5 only writes rows for sequences it identifies as recombinants;
    all others are inferred to be non-recombinant.
    """
    recombinants: dict[str, list[dict]] = {}   # seq_id → list of recomb events

    if csv_path.exists():
        with csv_path.open(newline="", encoding="utf-8-sig") as fh:
            all_lines = fh.readlines()

        # Find the header line (starts with Recombination Event Number or contains Recombinant Sequence(s))
        start_idx = -1
        for idx, line in enumerate(all_lines):
            if "Recombinant Sequence(s)" in line or "Recombination Event Number" in line:
                start_idx = idx
                break

        if start_idx != -1:
            lines = all_lines[start_idx:]
            reader = csv.DictReader(lines)
            for row in reader:
                # Clean row keys and values (strip whitespace, keep first of stripped duplicate keys)
                clean_row = {}
                for k, v in row.items():
                    if k is not None:
                        k_clean = k.strip()
                        if k_clean not in clean_row:
                            clean_row[k_clean] = v

                seq_id = _first_col(clean_row, _RECOMB_SEQ_COLS)
                if not seq_id:
                    continue
                # Strip leading caret (^), which RDP5 adds to flag potential misidentifications
                seq_id = seq_id.lstrip("^").strip()

                bp_start = _first_col(clean_row, _BP_START_COLS)
                bp_end   = _first_col(clean_row, _BP_END_COLS)

                # Clean any formatting suffixes/prefixes in breakpoints (like *, $, ~, #)
                if bp_start:
                    bp_start = str(bp_start).rstrip("*$~#").strip()
                if bp_end:
                    bp_end = str(bp_end).rstrip("*$~#").strip()

                try:
                    bp_start = int(float(bp_start)) if bp_start else 0
                    bp_end   = int(float(bp_end))   if bp_end   else 0
                except ValueError:
                    bp_start = bp_end = 0

                # Clean parent names of any carets as well
                parent_1 = _first_col(clean_row, _PARENT1_COLS).lstrip("^").strip()
                parent_2 = _first_col(clean_row, _PARENT2_COLS).lstrip("^").strip()

                # Find best p-value across methods
                p_val = _first_col(clean_row, _PVALUE_COLS)
                if not p_val:
                    method_cols = ["RDP", "GENECONV", "Bootscan", "Maxchi", "Chimaera", "SiSscan", "PhylPro", "LARD", "3Seq"]
                    best_p = 1.0
                    for m in method_cols:
                        val = clean_row.get(m)
                        if val and val.strip().upper() != "NS":
                            try:
                                p_num = float(val.strip())
                                if p_num < best_p:
                                    best_p = p_num
                            except ValueError:
                                pass
                    if best_p < 1.0:
                        p_val = str(best_p)

                # Find detecting methods list
                methods = _first_col(clean_row, _METHODS_COLS)
                if not methods:
                    methods_list = []
                    method_cols = ["RDP", "GENECONV", "Bootscan", "Maxchi", "Chimaera", "SiSscan", "PhylPro", "LARD", "3Seq"]
                    for m in method_cols:
                        val = clean_row.get(m)
                        if val and val.strip().upper() != "NS" and val.strip() != "":
                            methods_list.append(m)
                    methods = "|".join(methods_list)

                event = {
                    "bp_start": bp_start,
                    "bp_end":   bp_end,
                    "p_value":  p_val,
                    "methods":  methods,
                    "parent_1": parent_1,
                    "parent_2": parent_2,
                }
                recombinants.setdefault(seq_id, []).append(event)
    else:
        log.warning("RDP5 CSV not found at %s — treating all as non-recombinant.", csv_path)

    # Build one output row per input sequence
    results = []
    for sid in all_seq_ids:
        events = recombinants.get(sid, [])
        if events:
            # Take the breakpoint with the smallest p-value
            def pval_key(e):
                try:
                    return float(e["p_value"])
                except (ValueError, TypeError):
                    return 1.0

            best       = min(events, key=pval_key)
            all_bps    = [[e["bp_start"], e["bp_end"]] for e in events if e["bp_start"]]
            all_methods = "|".join(sorted({e["methods"] for e in events if e["methods"]}))

            results.append({
                "sequence_id":      sid,
                "is_recombinant":   "true",
                "breakpoint_start": best["bp_start"],
                "breakpoint_end":   best["bp_end"],
                "p_value":          best["p_value"],
                "methods":          all_methods,
                "parent_1":         best["parent_1"],
                "parent_2":         best["parent_2"],
            })
        else:
            results.append({
                "sequence_id":      sid,
                "is_recombinant":   "false",
                "breakpoint_start": 0,
                "breakpoint_end":   0,
                "p_value":          "",
                "methods":          "",
                "parent_1":         "",
                "parent_2":         "",
            })

    return results


def _write_tsv(rows: list[dict], out_path: Path):
    fieldnames = [
        "sequence_id", "is_recombinant",
        "breakpoint_start", "breakpoint_end",
        "p_value", "methods", "parent_1", "parent_2",
    ]
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    log.info("Wrote normalised TSV → %s", out_path)


# ---------------------------------------------------------------------------
# Core runner
# ---------------------------------------------------------------------------

def run_rdp5(fasta_path: Path, outdir: Path, virus: str = None, label: str = "sequences") -> Path:
    """
    Run RDP5CL.exe on `fasta_path`, write outputs to `outdir`,
    return the path to the normalised TSV.

    This function is safe to call from both Snakemake rules and
    pipeline_runner.py.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Preflight checks ──────────────────────────────────────────────────
    n_seqs = _count_sequences(fasta_path)
    seq_ids = _read_sequence_ids(fasta_path)

    if n_seqs < RDP5_MIN_SEQS:
        log.warning(
            "%s: only %d sequence(s) — RDP5 needs at least %d. "
            "Writing non-recombinant TSV and skipping.",
            label, n_seqs, RDP5_MIN_SEQS,
        )
        rows = [
            {
                "sequence_id": sid, "is_recombinant": "false",
                "breakpoint_start": 0, "breakpoint_end": 0,
                "p_value": "", "methods": "", "parent_1": "", "parent_2": "",
            }
            for sid in seq_ids
        ]
        out_tsv = outdir / "recombinants.tsv"
        _write_tsv(rows, out_tsv)
        return out_tsv

    rdp5_exe = _find_rdp5_exe()
    log.info("Using RDP5CL.exe: %s", rdp5_exe)
    rdp5_dir = rdp5_exe.parent

    # ── Copy RDP.ini from root directory depending on virus ────────────────
    if virus:
        ini_filename = f"{virus.lower()}_RDP.ini"
        src_ini = Path(__file__).parent.parent / ini_filename
        if src_ini.exists():
            dest_ini = rdp5_dir / "RDP.ini"
            log.info("Copying configuration: %s -> %s", src_ini, dest_ini)
            shutil.copy(src_ini, dest_ini)
        else:
            log.warning("Configuration file %s not found. RDP5 will run with defaults/existing RDP.ini", src_ini)


    # ── Copy FASTA into rdp5_dir so RDP5 writes its outputs there ──────────
    # RDP5 places output files next to the input file, named after it.
    work_fasta = rdp5_dir / f"{label}.fasta"
    shutil.copy(fasta_path, work_fasta)

    # ── Build and run command ─────────────────────────────────────────────
    # Pass relative path to avoid space and backslash translation bugs in Wine
    cmd = _build_cmd(rdp5_exe, Path(f"{label}.fasta"), outdir)
    log.info("Running: %s", " ".join(cmd))
    log.info("Timeout: %d s | Sequences: %d | Label: %s", RDP5_TIMEOUT, n_seqs, label)

    start = time.time()
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(rdp5_dir),       # RDP5 looks for RDP.ini and helper exes in cwd
            capture_output=True,
            text=True,
            timeout=RDP5_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        if work_fasta.exists():
            work_fasta.unlink()
        raise RuntimeError(
            f"RDP5 timed out after {RDP5_TIMEOUT} s on {label}.\n"
            f"Increase RDP5_TIMEOUT (current: {RDP5_TIMEOUT})."
        )

    elapsed = time.time() - start
    log.info("RDP5 finished in %.1f s (exit code %d)", elapsed, proc.returncode)

    # RDP5 returns non-zero even on success in some versions — we check
    # for output rather than exit code.
    if proc.returncode not in (0, 1):
        log.warning(
            "RDP5 exited with code %d. Will attempt to parse output anyway.",
            proc.returncode,
        )
        if proc.stdout:
            log.info("RDP5 STDOUT:\n%s", proc.stdout[-3000:])
        if proc.stderr:
            log.warning("RDP5 STDERR:\n%s", proc.stderr[-3000:])
    else:
        if proc.stdout:
            log.debug("STDOUT:\n%s", proc.stdout[-3000:])
        if proc.stderr:
            log.debug("STDERR:\n%s", proc.stderr[-3000:])

    # ── Locate and Move the CSV RDP5 wrote ─────────────────────────────────
    csv_src = rdp5_dir / f"{label}.fasta.csv"
    csv_dest = outdir / f"{label}.csv"

    if csv_src.exists():
        shutil.move(str(csv_src), str(csv_dest))
        log.info("Moved CSV: %s -> %s", csv_src, csv_dest)
    else:
        # Fallback: find any CSV ending with .fasta.csv or containing label
        csvs = list(rdp5_dir.glob(f"*{label}*.fasta.csv"))
        if not csvs:
            csvs = list(rdp5_dir.glob(f"*{label}*.csv"))
        
        # Exclude RecIDTests.csv to avoid matching wrong files
        csvs = [c for c in csvs if "RecIDTests" not in c.name]

        if csvs:
            csv_src = csvs[0]
            shutil.move(str(csv_src), str(csv_dest))
            log.info("Moved fallback CSV: %s -> %s", csv_src, csv_dest)
        else:
            # Cleanup FASTA before raising error
            if work_fasta.exists():
                work_fasta.unlink()
            if proc.stdout:
                log.error("RDP5 STDOUT on failure:\n%s", proc.stdout)
            if proc.stderr:
                log.error("RDP5 STDERR on failure:\n%s", proc.stderr)
            raise RuntimeError(
                f"RDP5 failed to generate any CSV output for {label} in {rdp5_dir}.\n"
                f"Exit code: {proc.returncode}\n"
                f"Check command arguments or Wine execution logs."
            )

    # ── Parse and write normalised TSV ────────────────────────────────────
    rows    = _parse_rdp5_csv(csv_dest, seq_ids)
    out_tsv = outdir / "recombinants.tsv"
    _write_tsv(rows, out_tsv)

    # ── Cleanup intermediate files in rdp5_dir ─────────────────────────────
    # Delete the copied FASTA and any other generated output files in rdp5_dir
    for p in rdp5_dir.glob(f"{label}*"):
        try:
            if p.exists() and p.is_file():
                p.unlink()
                log.debug("Cleaned up temp file: %s", p)
        except Exception as e:
            log.warning("Failed to clean up temp file %s: %s", p, e)

    n_recomb = sum(1 for r in rows if r["is_recombinant"] == "true")
    log.info("Result: %d / %d sequences flagged as recombinant.", n_recomb, len(rows))

    return out_tsv


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _parse_args():
    p = argparse.ArgumentParser(
        description="Wrapper around RDP5CL.exe — produces a normalised recombinants.tsv"
    )
    p.add_argument("--fasta",    required=True, help="Input FASTA alignment")
    p.add_argument("--virus",    required=True, choices=["hbv", "hcv", "hev"],
                   help="Virus type (used for logging only)")
    p.add_argument("--outdir",   required=True, help="Directory for all output files")
    p.add_argument("--genotype", default=None,
                   help="Genotype label (used as output file prefix in batch mode)")
    p.add_argument("--min-seqs", type=int, default=RDP5_MIN_SEQS,
                   help=f"Minimum sequences to attempt RDP5 (default: {RDP5_MIN_SEQS})")
    p.add_argument("--timeout",  type=int, default=RDP5_TIMEOUT,
                   help=f"RDP5 timeout in seconds (default: {RDP5_TIMEOUT})")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()

    # Allow CLI overrides of globals
    RDP5_MIN_SEQS = args.min_seqs
    RDP5_TIMEOUT  = args.timeout

    label = args.genotype or Path(args.fasta).stem

    out_tsv = run_rdp5(
        fasta_path = Path(args.fasta),
        outdir     = Path(args.outdir),
        virus      = args.virus,
        label      = label,
    )
    print(f"Done → {out_tsv}")
