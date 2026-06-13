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
    import concurrent.futures
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Preflight checks ──────────────────────────────────────────────────
    n_seqs = _count_sequences(fasta_path)
    seq_ids = _read_sequence_ids(fasta_path)

    # Parse and split into parents (references starting with ref_) and children (query seqs)
    parents = []
    children = []
    
    current_header = None
    current_seq = []
    
    with fasta_path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if current_header:
                    seq_entry = (current_header, "".join(current_seq))
                    if current_header[1:].startswith("ref_"):
                        parents.append(seq_entry)
                    else:
                        children.append(seq_entry)
                current_header = line.strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_header:
            seq_entry = (current_header, "".join(current_seq))
            if current_header[1:].startswith("ref_"):
                parents.append(seq_entry)
            else:
                children.append(seq_entry)

    run_id = f"run_{label}"
    parents_path = outdir / f"{run_id}_parents.fasta"

    use_two_file_mode = len(parents) > 0 and len(children) > 0

    if use_two_file_mode:
        target_seq_ids = [c[0][1:].split()[0].strip() for c in children]
        # Write parents file
        with parents_path.open("w") as fh:
            for h, s in parents:
                fh.write(f"{h}\n{s}\n")
        log.info("Parsed %d references (parents) and %d query sequences (children).", len(parents), len(children))
    else:
        target_seq_ids = seq_ids
        log.info("Running in single-file mode: %d sequences", n_seqs)

    if not use_two_file_mode and n_seqs < RDP5_MIN_SEQS:
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
            for sid in target_seq_ids
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

    # Clean up old run files
    for p in outdir.glob(f"{run_id}.3s.*"):
        try:
            p.unlink()
        except Exception as e:
            log.warning("Failed to clean up old 3seq output file %s: %s", p, e)

    # Find P-value table file
    script_dir = Path(__file__).resolve().parent
    ptable_path = script_dir / "PVT.3SEQ.2017.700"
    if not ptable_path.exists():
        ptable_path = Path.cwd() / "PVT.3SEQ.2017.700"

    start = time.time()

    if use_two_file_mode:
        # Determine parallel workers and chunk size
        max_workers = min(8, os.cpu_count() or 4)
        chunk_size = (len(children) + max_workers - 1) // max_workers
        chunks = [children[i:i + chunk_size] for i in range(0, len(children), chunk_size)]
        log.info("Splitting %d children into %d chunks of size %d for parallel 3seq runs (workers: %d)...",
                 len(children), len(chunks), chunk_size, max_workers)

        def _run_chunk(chunk_seqs, idx):
            chunk_children_path = outdir / f"{run_id}_children_chunk{idx}.fasta"
            with chunk_children_path.open("w") as fh:
                for h, s in chunk_seqs:
                    fh.write(f"{h}\n{s}\n")
            
            chunk_run_id = f"{run_id}_chunk{idx}"
            cmd = [
                str(threeseq_exe),
                "-full",
                str(parents_path.resolve()),
                str(chunk_children_path.resolve()),
            ]
            if ptable_path.exists():
                cmd.extend(["-ptable", str(ptable_path.resolve())])
            cmd.extend([
                "-id", chunk_run_id,
                "-p"
            ])

            proc = subprocess.run(
                cmd,
                cwd=str(outdir),
                input="Y\n",
                capture_output=True,
                text=True,
                timeout=RDP5_TIMEOUT,
            )
            if proc.returncode != 0:
                log.error("Chunk %d failed. STDOUT:\n%s", idx, proc.stdout)
                log.error("Chunk %d failed. STDERR:\n%s", idx, proc.stderr)
                raise RuntimeError(f"3seq chunk {idx} exited with code {proc.returncode}")
            
            # Clean up chunk children fasta file
            try:
                chunk_children_path.unlink()
            except Exception:
                pass
            return idx

        # Run chunks in parallel
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(_run_chunk, chunk, idx) for idx, chunk in enumerate(chunks)]
            for fut in concurrent.futures.as_completed(futures):
                fut.result() # raises exception if chunk failed

        # Merge chunk outputs
        csv_src = outdir / f"{run_id}.3s.rec.csv"
        longrec_src = outdir / f"{run_id}.3s.longRec"
        
        unique_longrec_lines = set()
        header_written = False
        
        with csv_src.open("w", encoding="utf-8") as out_fh:
            for idx in range(len(chunks)):
                chunk_csv = outdir / f"{run_id}_chunk{idx}.3s.rec.csv"
                if chunk_csv.exists():
                    with chunk_csv.open("r", encoding="utf-8") as in_fh:
                        header = in_fh.readline()
                        if not header_written and header:
                            out_fh.write(header)
                            header_written = True
                        for line in in_fh:
                            out_fh.write(line)
                    try:
                        chunk_csv.unlink()
                    except Exception:
                        pass
                
                # Load and accumulate longRec chunk files
                chunk_longrec = outdir / f"{run_id}_chunk{idx}.3s.longRec"
                if chunk_longrec.exists():
                    with chunk_longrec.open("r", encoding="utf-8") as in_fh:
                        for line in in_fh:
                            unique_longrec_lines.add(line.strip())
                    try:
                        chunk_longrec.unlink()
                    except Exception:
                        pass
                
                # Cleanup other 3seq chunk files
                for p in outdir.glob(f"{run_id}_chunk{idx}.3s.*"):
                    try:
                        p.unlink()
                    except Exception:
                        pass
                        
        if unique_longrec_lines:
            with longrec_src.open("w", encoding="utf-8") as out_fh:
                for line in sorted(unique_longrec_lines):
                    out_fh.write(f"{line}\n")
    else:
        # Single file mode
        cmd = [
            str(threeseq_exe),
            "-full",
            str(fasta_path.resolve()),
        ]
        if ptable_path.exists():
            cmd.extend(["-ptable", str(ptable_path.resolve())])
        cmd.extend([
            "-id", run_id,
            "-p"
        ])
        log.info("Running standard 3seq: %s", " ".join(cmd))
        
        try:
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

        if proc.returncode != 0:
            log.error("3seq STDOUT:\n%s", proc.stdout)
            log.error("3seq STDERR:\n%s", proc.stderr)
            raise RuntimeError(f"3seq exited with code {proc.returncode}")

    elapsed = time.time() - start
    log.info("3seq finished in %.1f s", elapsed)

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
            for sid in target_seq_ids
        ]
        out_tsv = outdir / "recombinants.tsv"
        _write_tsv(rows, out_tsv)
        return out_tsv

    # Parse and write normalised TSV
    rows = _parse_3seq_csv(csv_src, target_seq_ids)
    out_tsv = outdir / "recombinants.tsv"
    _write_tsv(rows, out_tsv)

    # Cleanup temporary parents file
    if parents_path.exists():
        try:
            parents_path.unlink()
        except Exception:
            pass

    # Cleanup intermediate 3seq files
    for p in outdir.glob(f"{run_id}.3s.*"):
        if p.exists() and p.name not in [f"{run_id}.3s.rec.csv", f"{run_id}.3s.longRec"]:
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
