"""
scripts/run_rdp5.py
===================
Wrapper around 3seq for both:

  (A) Snakemake batch mode  — called by the `rdp5_recombination` rule
      Input:  a per-genotype aligned FASTA  (results/{virus}/alignments/{genotype}_aligned.fasta)
      Output: results/{virus}/recombination/{genotype}/recombinants.tsv   (normalised, pipeline-ready)

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
  p_value           Best p-value across all 3Seq runs
  methods           3Seq
  parent_1          Inferred major parent (empty if not recombinant)
  parent_2          Inferred minor parent (empty if not recombinant)
"""

import argparse
import csv
import logging
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("run_3seq")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

RDP5_TIMEOUT  = int(os.environ.get("RDP5_TIMEOUT",  "3600"))
RDP5_MIN_SEQS = int(os.environ.get("RDP5_MIN_SEQS", "4"))


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
# 3seq CSV parser
# ---------------------------------------------------------------------------

def _parse_3seq_csv(csv_path: Path, all_seq_ids: list[str]) -> list[dict]:
    """
    Parse the 3seq output CSV file and return a list of normalised dicts,
    one per sequence in all_seq_ids (preserving order).
    """
    recombinants = {}

    with csv_path.open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header:
            log.warning("Empty 3seq CSV file: %s", csv_path)
            return []

        for row in reader:
            if not row or len(row) < 13:
                continue

            parent_1_raw = row[0].strip()
            parent_2_raw = row[1].strip()
            child_raw = row[2].strip()
            p_val_raw = row[6].strip()
            bp_str = row[12].strip()

            # Map to first token of the raw headers to match all_seq_ids
            child_id = child_raw.split()[0] if child_raw else ""
            parent_1 = parent_1_raw.split()[0] if parent_1_raw else ""
            parent_2 = parent_2_raw.split()[0] if parent_2_raw else ""

            if not child_id:
                continue

            # Parse breakpoints (e.g. "1742-1751 & 2585-2590")
            bp_start = 0
            bp_end = 0
            if "&" in bp_str:
                parts = bp_str.split("&")
                if len(parts) == 2:
                    start_part = parts[0].strip()
                    end_part = parts[1].strip()

                    def get_midpoint(r_str):
                        r_str = r_str.split()[0]
                        if "-" in r_str:
                            sub_parts = r_str.split("-")
                            try:
                                return (int(sub_parts[0]) + int(sub_parts[1])) // 2
                            except ValueError:
                                return 0
                        else:
                            try:
                                return int(r_str)
                            except ValueError:
                                return 0

                    bp_start = get_midpoint(start_part)
                    bp_end = get_midpoint(end_part)

            event = {
                "bp_start": bp_start,
                "bp_end": bp_end,
                "p_value": p_val_raw,
                "methods": "3Seq",
                "parent_1": parent_1,
                "parent_2": parent_2,
            }
            recombinants.setdefault(child_id, []).append(event)

    # Build one output row per input sequence
    results = []
    for sid in all_seq_ids:
        events = recombinants.get(sid, [])
        if not events:
            # Fallback prefix/suffix matching
            for cid, evs in recombinants.items():
                if cid.startswith(sid) or sid.startswith(cid):
                    events = evs
                    break

        if events:
            def pval_key(e):
                try:
                    return float(e["p_value"])
                except (ValueError, TypeError):
                    return 1.0

            best = min(events, key=pval_key)
            results.append({
                "sequence_id": sid,
                "is_recombinant": "true",
                "breakpoint_start": best["bp_start"],
                "breakpoint_end": best["bp_end"],
                "p_value": best["p_value"],
                "methods": "3Seq",
                "parent_1": best["parent_1"],
                "parent_2": best["parent_2"],
            })
        else:
            results.append({
                "sequence_id": sid,
                "is_recombinant": "false",
                "breakpoint_start": 0,
                "breakpoint_end": 0,
                "p_value": "",
                "methods": "",
                "parent_1": "",
                "parent_2": "",
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
    Run 3seq on `fasta_path`, write outputs to `outdir`,
    return the path to the normalised TSV.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Preflight checks ──────────────────────────────────────────────────
    n_seqs = _count_sequences(fasta_path)
    seq_ids = _read_sequence_ids(fasta_path)

    if n_seqs < RDP5_MIN_SEQS:
        log.warning(
            "%s: only %d sequence(s) — 3seq needs at least %d. "
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

    # Find 3seq binary
    threeseq_exe = shutil.which("3seq")
    if not threeseq_exe:
        # Check standard miniconda location
        candidates = [
            Path("/Users/derektshiabuila/miniconda3/bin/3seq"),
            Path.home() / "miniconda3" / "bin" / "3seq",
            Path.home() / "anaconda3" / "bin" / "3seq",
        ]
        for c in candidates:
            if c.exists():
                threeseq_exe = str(c)
                break
    if not threeseq_exe:
        threeseq_exe = "3seq"

    log.info("Using 3seq executable: %s", threeseq_exe)

    run_id = f"run_{label}"
    # Clean up old run files
    for p in outdir.glob(f"{run_id}.3s.*"):
        try:
            p.unlink()
        except Exception as e:
            log.warning("Failed to clean up old 3seq output file %s: %s", p, e)

    # 3seq command
    cmd = [
        str(threeseq_exe),
        "-full",
        str(fasta_path.resolve()),
        "-id", run_id,
        "-p"
    ]
    log.info("Running: %s", " ".join(cmd))
    log.info("Cwd: %s", outdir)

    start = time.time()
    try:
        # Input "Y\n" to standard input to auto-confirm memory limits
        proc = subprocess.run(
            cmd,
            cwd=str(outdir),
            input="Y\n",
            capture_output=True,
            text=True,
            timeout=RDP5_TIMEOUT,
        )
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"3seq timed out after {RDP5_TIMEOUT} s on {label}."
        )

    elapsed = time.time() - start
    log.info("3seq finished in %.1f s (exit code %d)", elapsed, proc.returncode)

    if proc.returncode != 0:
        log.error("3seq STDOUT:\n%s", proc.stdout)
        log.error("3seq STDERR:\n%s", proc.stderr)
        raise RuntimeError(f"3seq exited with code {proc.returncode}")

    # Check for CSV output
    csv_src = outdir / f"{run_id}.3s.rec.csv"
    if not csv_src.exists():
        log.warning("3seq CSV not found at %s. No recombination events detected.", csv_src)
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

    # Parse and write normalised TSV
    rows = _parse_3seq_csv(csv_src, seq_ids)
    out_tsv = outdir / "recombinants.tsv"
    _write_tsv(rows, out_tsv)

    # Cleanup intermediate 3seq files
    for p in outdir.glob(f"{run_id}.3s.*"):
        if p.exists() and p.name != f"{run_id}.3s.rec.csv":
            try:
                p.unlink()
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
        description="Wrapper around 3seq — produces a normalised recombinants.tsv"
    )
    p.add_argument("--fasta",    required=True, help="Input FASTA alignment")
    p.add_argument("--virus",    required=True, choices=["hbv", "hcv", "hev"],
                   help="Virus type (used for logging only)")
    p.add_argument("--outdir",   required=True, help="Directory for all output files")
    p.add_argument("--genotype", default=None,
                   help="Genotype label (used as output file prefix in batch mode)")
    p.add_argument("--min-seqs", type=int, default=RDP5_MIN_SEQS,
                   help=f"Minimum sequences to attempt 3seq (default: {RDP5_MIN_SEQS})")
    p.add_argument("--timeout",  type=int, default=RDP5_TIMEOUT,
                   help=f"3seq timeout in seconds (default: {RDP5_TIMEOUT})")
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
