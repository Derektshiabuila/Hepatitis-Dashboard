#!/usr/bin/env python3
"""
scripts/validate_recombination.py
==================================
Automated validation of 3seq recombination results using:
  1. 3s.longRec length filtering (>= 100bp).
  2. Cross-method alignment distance validation (RDP/OpenRDP proxy).
  3. Regional ML trees (biological clustering validation via IQ-TREE).

Usage:
  python scripts/validate_recombination.py --virus hbv
"""

import argparse
import csv
import logging
import os
import shutil
import subprocess
from pathlib import Path
from Bio import SeqIO, Phylo

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S"
)
log = logging.getLogger("validate_recomb")

# Try to import openrdp, or clone it first if not present
import sys
try:
    import openrdp
except ImportError:
    project_root = Path(__file__).resolve().parent.parent
    openrdp_path = project_root / "scratch" / "OpenRDP"
    if not openrdp_path.exists():
        openrdp_path = Path.cwd() / "scratch" / "OpenRDP"
        
    if not openrdp_path.exists():
        openrdp_path = project_root / "scratch" / "OpenRDP"
        log.info("OpenRDP not found locally. Cloning PoonLab/OpenRDP from GitHub...")
        openrdp_path.parent.mkdir(parents=True, exist_ok=True)
        cmd = ["git", "clone", "https://github.com/PoonLab/OpenRDP.git", str(openrdp_path)]
        subprocess.run(cmd, check=True)
        
    sys.path.insert(0, str(openrdp_path.resolve()))
    try:
        import openrdp
    except ImportError as e:
        log.error("Failed to import openrdp after cloning: %s", e)
        raise e

def load_long_rec_candidates(longrec_path: Path) -> set:
    """Load sequence IDs from the .3s.longRec file."""
    candidates = set()
    if not longrec_path.exists():
        return candidates
    with longrec_path.open() as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("ref_"):
                # Clean header to get accession ID (first token)
                candidates.add(line.split()[0])
    return candidates

def get_alignment_sequences(aligned_path: Path, ids_to_extract: list) -> dict:
    """Extract sequences from an aligned FASTA file for specific IDs."""
    seqs = {}
    if not aligned_path.exists():
        return seqs
    for rec in SeqIO.parse(aligned_path, "fasta"):
        # Match by first token of ID
        acc = rec.id.split()[0]
        if acc in ids_to_extract:
            seqs[acc] = str(rec.seq).upper()
        # Also match fallback prefix/suffix
        else:
            for target_id in ids_to_extract:
                if acc.startswith(target_id) or target_id.startswith(acc):
                    seqs[target_id] = str(rec.seq).upper()
                    break
    return seqs

def calculate_identity(seq1: str, seq2: str, start: int, end: int, invert: bool = False) -> float:
    """Calculate nucleotide identity between two aligned sequences in a region (1-indexed)."""
    # Convert to 0-indexed coordinates
    length = len(seq1)
    matches = 0
    total = 0
    
    for i in range(length):
        in_region = (start - 1) <= i < end
        if invert:
            in_region = not in_region
            
        if in_region:
            char1 = seq1[i]
            char2 = seq2[i]
            # Ignore gaps and ambiguous nucleotides
            if char1 not in ('-', 'N') and char2 not in ('-', 'N'):
                total += 1
                if char1 == char2:
                    matches += 1
                    
    return (matches / total) if total > 0 else 0.0

def sanitize_sequence_for_openrdp(seq_str: str) -> str:
    """Replace any non-standard IUPAC characters (anything except A, C, G, T, -) with N."""
    sanitized = []
    for char in seq_str.upper():
        if char in ('A', 'C', 'G', 'T', '-'):
            sanitized.append(char)
        else:
            sanitized.append('N')
    return "".join(sanitized)

def run_openrdp_validation(
    aligned_path: Path, candidate_id: str, p1_id: str, p2_id: str, temp_dir: Path
) -> bool:
    """Run OpenRDP cross-validation on the triplet alignment."""
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    # Load all records
    all_recs = list(SeqIO.parse(aligned_path, "fasta"))
    extracted = {}
    
    # Locate candidate and parent sequences
    for rec in all_recs:
        acc = rec.id.split()[0]
        if acc == candidate_id:
            extracted[candidate_id] = rec
        elif acc == p1_id:
            extracted[p1_id] = rec
        elif acc == p2_id:
            extracted[p2_id] = rec
            
    # Fallback to fuzzy prefix/suffix matching if not found exactly
    for acc in [candidate_id, p1_id, p2_id]:
        if acc not in extracted:
            for rec in all_recs:
                rec_acc = rec.id.split()[0]
                if rec_acc.startswith(acc) or acc.startswith(rec_acc):
                    extracted[acc] = rec
                    break
                    
    if len(extracted) < 3:
        log.warning("Could not extract all sequences for OpenRDP validation: %s, %s, %s", candidate_id, p1_id, p2_id)
        return False
        
    temp_fasta = temp_dir / f"temp_{candidate_id}.fasta"
    
    # Write sanitized sequences to temporary FASTA
    sanitized_recs = []
    for acc in [candidate_id, p1_id, p2_id]:
        rec = extracted[acc]
        rec_copy = rec.__class__(
            id=rec.id,
            name=rec.name,
            description=rec.description,
            seq=rec.seq.__class__(sanitize_sequence_for_openrdp(str(rec.seq)))
        )
        sanitized_recs.append(rec_copy)
        
    SeqIO.write(sanitized_recs, temp_fasta, "fasta")
    
    try:
        # Run Scanner with standard methods (rdp, maxchi, chimaera)
        scanner = openrdp.Scanner(methods=("rdp", "maxchi", "chimaera"), verbose=False)
        results = scanner.run_scans(str(temp_fasta))
        
        # Verify if any method detected a significant event (p-value < 0.05)
        for method in results.keys():
            events = results[method]
            for e in events:
                p_val = e.get("pvalue", 1.0)
                if p_val < 0.05:
                    log.info("OpenRDP: validated %s using method '%s' (p-value: %.2E)", candidate_id, method, p_val)
                    return True
        return False
    except Exception as e:
        log.warning("OpenRDP validation failed for %s: %s", candidate_id, e)
        return False
    finally:
        if temp_fasta.exists():
            try:
                temp_fasta.unlink()
            except Exception:
                pass

def run_regional_tree_validation(
    aligned_path: Path, candidate_id: str, p1_id: str, p2_id: str,
    start: int, end: int, scratch_dir: Path
) -> bool:
    """Build regional trees using IQ-TREE and verify parent clustering."""
    scratch_dir.mkdir(parents=True, exist_ok=True)
    
    # Load all records to pick outgroups
    all_recs = list(SeqIO.parse(aligned_path, "fasta"))
    outgroups = []
    for rec in all_recs:
        acc = rec.id.split()[0]
        if acc not in (candidate_id, p1_id, p2_id) and acc.startswith("ref_"):
            outgroups.append(rec)
            if len(outgroups) >= 2:
                break
                
    # Combine candidate, parents, and outgroups
    outgroup_ids = {r.id for r in outgroups}
    selected_recs = []
    for rec in all_recs:
        acc = rec.id.split()[0]
        if acc in (candidate_id, p1_id, p2_id) or rec.id in outgroup_ids:
            selected_recs.append(rec)
            
    # Paths for sliced alignments
    recomb_aln_path = scratch_dir / "recomb_slice.fasta"
    nonrecomb_aln_path = scratch_dir / "nonrecomb_slice.fasta"
    
    # Slice alignments
    with recomb_aln_path.open("w") as recomb_fh, nonrecomb_aln_path.open("w") as nonrecomb_fh:
        for rec in selected_recs:
            seq_str = str(rec.seq)
            # Recombinant region slice
            recomb_seq = seq_str[start - 1:end]
            # Non-recombinant region slice
            nonrecomb_seq = seq_str[:start - 1] + seq_str[end:]
            
            recomb_fh.write(f">{rec.id}\n{recomb_seq}\n")
            nonrecomb_fh.write(f">{rec.id}\n{nonrecomb_seq}\n")
            
    # Find IQ-TREE binary
    iqtree_bin = shutil.which("iqtree2") or shutil.which("iqtree")
    if not iqtree_bin:
        log.warning("IQ-TREE executable not found (neither iqtree2 nor iqtree). Skipping tree validation.")
        return False

    # Run IQ-TREE on recomb slice
    recomb_prefix = scratch_dir / "recomb_tree"
    cmd_recomb = [iqtree_bin, "-s", str(recomb_aln_path), "--prefix", str(recomb_prefix), "-m", "GTR", "-nt", "1", "-redo"]
    proc = subprocess.run(cmd_recomb, capture_output=True, text=True)
    if proc.returncode != 0:
        log.warning("IQ-TREE recomb slice failed for %s. STDERR: %s", candidate_id, proc.stderr)
    
    # Run IQ-TREE on non-recomb slice
    nonrecomb_prefix = scratch_dir / "nonrecomb_tree"
    cmd_nonrecomb = [iqtree_bin, "-s", str(nonrecomb_aln_path), "--prefix", str(nonrecomb_prefix), "-m", "GTR", "-nt", "1", "-redo"]
    proc = subprocess.run(cmd_nonrecomb, capture_output=True, text=True)
    if proc.returncode != 0:
        log.warning("IQ-TREE non-recomb slice failed for %s. STDERR: %s", candidate_id, proc.stderr)
    
    recomb_tree_path = scratch_dir / "recomb_tree.treefile"
    nonrecomb_tree_path = scratch_dir / "nonrecomb_tree.treefile"
    
    if not recomb_tree_path.exists() or not nonrecomb_tree_path.exists():
        log.warning("IQ-TREE failed to generate treefiles for %s.", candidate_id)
        return False
        
    try:
        # Load trees
        recomb_tree = Phylo.read(recomb_tree_path, "newick")
        nonrecomb_tree = Phylo.read(nonrecomb_tree_path, "newick")
        
        # Helper to find node matching ID
        def find_clade_name(tree, target_id):
            for leaf in tree.get_terminals():
                if leaf.name.split()[0] == target_id or target_id in leaf.name:
                    return leaf.name
            return None
            
        cand_name_r = find_clade_name(recomb_tree, candidate_id)
        p1_name_r = find_clade_name(recomb_tree, p1_id)
        p2_name_r = find_clade_name(recomb_tree, p2_id)
        
        cand_name_nr = find_clade_name(nonrecomb_tree, candidate_id)
        p1_name_nr = find_clade_name(nonrecomb_tree, p1_id)
        p2_name_nr = find_clade_name(nonrecomb_tree, p2_id)
        
        if not all((cand_name_r, p1_name_r, p2_name_r, cand_name_nr, p1_name_nr, p2_name_nr)):
            log.warning("Could not map all taxa in regional trees for %s.", candidate_id)
            return False
            
        # Calculate cophenetic tree distances
        dist_p1_r = recomb_tree.distance(cand_name_r, p1_name_r)
        dist_p2_r = recomb_tree.distance(cand_name_r, p2_name_r)
        
        dist_p1_nr = nonrecomb_tree.distance(cand_name_nr, p1_name_nr)
        dist_p2_nr = nonrecomb_tree.distance(cand_name_nr, p2_name_nr)
        
        # Validation condition:
        # Recomb tree: Closer to Parent 2 (minor parent)
        # Non-recomb tree: Closer to Parent 1 (major parent)
        passed_recomb = dist_p2_r < dist_p1_r
        passed_nonrecomb = dist_p1_nr < dist_p2_nr
        
        return passed_recomb and passed_nonrecomb
        
    except Exception as e:
        log.warning("Error parsing regional trees for %s: %s", candidate_id, e)
        return False
    finally:
        # Cleanup temp iqtree files
        for f in scratch_dir.glob("recomb_*"):
            try: f.unlink()
            except Exception: pass
        for f in scratch_dir.glob("nonrecomb_*"):
            try: f.unlink()
            except Exception: pass

def main():
    parser = argparse.ArgumentParser(description="Validate 3seq recombination events")
    parser.add_argument("--virus", required=True, choices=["hbv", "hcv", "hev"], help="Virus type")
    parser.add_argument("--results-dir", default=None, help="Root results directory")
    args = parser.parse_args()
    
    virus = args.virus
    base_dir = Path(args.results_dir) if args.results_dir else Path(f"results/{virus}")
    recomb_summary_file = base_dir / "all_recombination.tsv"
    
    if not recomb_summary_file.exists():
        log.error("Recombination summary file %s not found. Run the pipeline first.", recomb_summary_file)
        sys.exit(1)
        
    log.info("Starting validation for %s recombination events...", virus.upper())
    
    # Load all candidate events
    candidates = []
    with recomb_summary_file.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("is_recombinant") == "true":
                candidates.append(row)
                
    log.info("Loaded %d recombinant candidate sequences from %s.", len(candidates), recomb_summary_file.name)
    
    validated_rows = []
    stats = {
        "total": len(candidates),
        "long_rec_passed": 0,
        "dist_passed": 0,
        "tree_passed": 0
    }
    
    scratch_dir = Path("scratch/recomb_val")
    scratch_dir.mkdir(parents=True, exist_ok=True)
    
    warned_missing_paths = set()
    
    for cand in candidates:
        candidate_id = cand["sequence_id"]
        # Find genotype of this sequence
        # Locate recombinants.tsv files to map sequence genotype
        genotype = None
        for p in base_dir.glob("recombination/*/recombinants.tsv"):
            with p.open() as r_fh:
                r_reader = csv.DictReader(r_fh, delimiter="\t")
                for r_row in r_reader:
                    if r_row["sequence_id"] == candidate_id and r_row["is_recombinant"] == "true":
                        genotype = p.parent.name
                        break
            if genotype:
                break
                
        if not genotype:
            log.warning("Genotype not found for candidate %s, skipping.", candidate_id)
            continue
            
        # ── Step 1: LongRec Filter ─────────────────────────────────────────
        longrec_path = base_dir / "recombination" / genotype / f"run_{genotype}.3s.longRec"
        if not longrec_path.exists():
            if longrec_path not in warned_missing_paths:
                log.warning("LongRec file %s not found. Re-run with '-R rdp5_recombination' to generate/preserve it.", longrec_path)
                warned_missing_paths.add(longrec_path)
            continue
            
        long_rec_set = load_long_rec_candidates(longrec_path)
        
        if candidate_id not in long_rec_set:
            log.debug("Candidate %s rejected: not in .3s.longRec (segment too short).", candidate_id)
            continue
        stats["long_rec_passed"] += 1
        
        # ── Step 2: Distance Identity Cross-Validation ──────────────────────
        aligned_path = base_dir / "alignments" / f"{genotype}_aligned.fasta"
        p1_id = cand["parent_1"]
        p2_id = cand["parent_2"]
        start = int(cand["breakpoint_start"])
        end = int(cand["breakpoint_end"])
        
        extracted_seqs = get_alignment_sequences(aligned_path, [candidate_id, p1_id, p2_id])
        if candidate_id not in extracted_seqs or p1_id not in extracted_seqs or p2_id not in extracted_seqs:
            log.warning("Could not extract alignments for triplet: %s (P1: %s, P2: %s)", candidate_id, p1_id, p2_id)
            continue
            
        cand_seq = extracted_seqs[candidate_id]
        p1_seq = extracted_seqs[p1_id]
        p2_seq = extracted_seqs[p2_id]
        
        # Identity in recombinant region (expect candidate closer to P2)
        ident_p1_r = calculate_identity(cand_seq, p1_seq, start, end, invert=False)
        ident_p2_r = calculate_identity(cand_seq, p2_seq, start, end, invert=False)
        
        # Identity in non-recombinant region (expect candidate closer to P1)
        ident_p1_nr = calculate_identity(cand_seq, p1_seq, start, end, invert=True)
        ident_p2_nr = calculate_identity(cand_seq, p2_seq, start, end, invert=True)
        
        # Check parent switching
        passed_dist = (ident_p2_r > ident_p1_r) and (ident_p1_nr > ident_p2_nr)
        
        if not passed_dist:
            log.debug("Candidate %s rejected: failed distance parent-switching check.", candidate_id)
            continue
            
        # ── Step 2b: OpenRDP Cross-Validation ────────────────────────────────
        passed_openrdp = run_openrdp_validation(aligned_path, candidate_id, p1_id, p2_id, scratch_dir)
        if not passed_openrdp:
            log.debug("Candidate %s rejected: failed OpenRDP cross-method validation.", candidate_id)
            continue
            
        stats["dist_passed"] += 1
        
        # ── Step 3: Regional ML Tree Validation ─────────────────────────────
        passed_tree = run_regional_tree_validation(aligned_path, candidate_id, p1_id, p2_id, start, end, scratch_dir)
        
        if not passed_tree:
            log.debug("Candidate %s rejected: failed regional tree clustering check.", candidate_id)
            continue
        stats["tree_passed"] += 1
        
        # If all steps pass, it is a validated recombinant sequence!
        validated_row = cand.copy()
        validated_row["genotype"] = genotype
        validated_row["p1_identity_recomb"] = f"{ident_p1_r:.4f}"
        validated_row["p2_identity_recomb"] = f"{ident_p2_r:.4f}"
        validated_row["p1_identity_nonrecomb"] = f"{ident_p1_nr:.4f}"
        validated_row["p2_identity_nonrecomb"] = f"{ident_p2_nr:.4f}"
        validated_rows.append(validated_row)
        log.info("Recombinant sequence VALIDATED: %s (Genotype: %s, Parents: %s & %s)", candidate_id, genotype, p1_id, p2_id)

    # Write output tsv
    out_file = base_dir / "validated_recombinants.tsv"
    fieldnames = [
        "sequence_id", "genotype", "is_recombinant", "breakpoint_start", "breakpoint_end", 
        "p_value", "methods", "parent_1", "parent_2",
        "p1_identity_recomb", "p2_identity_recomb", "p1_identity_nonrecomb", "p2_identity_nonrecomb"
    ]
    with out_file.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(validated_rows)
        
    # Clean up scratch dir
    try:
        shutil.rmtree(scratch_dir)
    except Exception:
        pass
        
    print("\n" + "="*80)
    print(f"RECOMBINATION VALIDATION SUMMARY FOR {virus.upper()}")
    print("="*80)
    print(f"Total Candidate Sequences Tested : {stats['total']}")
    print(f"Passed Step 1 (LongRec Filter)   : {stats['long_rec_passed']} / {stats['total']} ({stats['long_rec_passed']/stats['total']*100:.1f}%)" if stats['total'] > 0 else "N/A")
    print(f"Passed Step 2 (Dist + OpenRDP)   : {stats['dist_passed']} / {stats['long_rec_passed']} ({stats['dist_passed']/stats['long_rec_passed']*100:.1f}%)" if stats['long_rec_passed'] > 0 else "N/A")
    print(f"Passed Step 3 (Regional Trees)   : {stats['tree_passed']} / {stats['dist_passed']} ({stats['tree_passed']/stats['dist_passed']*100:.1f}%)" if stats['dist_passed'] > 0 else "N/A")
    print(f"Total Validated Recombinants     : {len(validated_rows)}")
    print(f"Output TSV file written to       : {out_file.resolve()}")
    print("="*80 + "\n")

if __name__ == "__main__":
    main()
