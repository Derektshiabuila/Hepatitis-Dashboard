import sys
import csv
import logging
from pathlib import Path
from Bio import SeqIO
import argparse

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-8s  %(message)s")
log = logging.getLogger("diagnose_val")

# Dynamic import logic for OpenRDP
try:
    import openrdp
except ImportError:
    project_root = Path(__file__).resolve().parent.parent
    openrdp_path = project_root / "scratch" / "OpenRDP"
    if not openrdp_path.exists():
        openrdp_path = Path.cwd() / "scratch" / "OpenRDP"
    if openrdp_path.exists():
        sys.path.insert(0, str(openrdp_path.resolve()))
        import openrdp
    else:
        log.error("OpenRDP not found. Diagnostic might skip OpenRDP validation.")

def load_long_rec_candidates(longrec_path: Path) -> set:
    candidates = set()
    if not longrec_path.exists():
        return candidates
    with longrec_path.open() as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("ref_"):
                candidates.add(line.split()[0])
    return candidates

def get_alignment_sequences(aligned_path: Path, sequence_ids: list) -> dict:
    seqs = {}
    for rec in SeqIO.parse(aligned_path, "fasta"):
        acc = rec.id.split()[0]
        if acc in sequence_ids:
            seqs[acc] = str(rec.seq)
    return seqs

def calculate_identity(seq1: str, seq2: str, start: int, end: int, invert: bool = False) -> float:
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
            if char1 not in ('-', 'N') and char2 not in ('-', 'N'):
                total += 1
                if char1 == char2:
                    matches += 1
    return (matches / total) if total > 0 else 0.0

def sanitize_sequence_for_openrdp(seq_str: str) -> str:
    sanitized = []
    for char in seq_str.upper():
        if char in ('A', 'C', 'G', 'T', '-'):
            sanitized.append(char)
        else:
            sanitized.append('N')
    return "".join(sanitized)

def run_openrdp_validation(aligned_path: Path, candidate_id: str, p1_id: str, p2_id: str, temp_dir: Path) -> bool:
    if "openrdp" not in sys.modules:
        return False
    temp_dir.mkdir(parents=True, exist_ok=True)
    all_recs = list(SeqIO.parse(aligned_path, "fasta"))
    extracted = {}
    for rec in all_recs:
        acc = rec.id.split()[0]
        if acc == candidate_id: extracted[candidate_id] = rec
        elif acc == p1_id: extracted[p1_id] = rec
        elif acc == p2_id: extracted[p2_id] = rec
    for acc in [candidate_id, p1_id, p2_id]:
        if acc not in extracted:
            for rec in all_recs:
                rec_acc = rec.id.split()[0]
                if rec_acc.startswith(acc) or acc.startswith(rec_acc):
                    extracted[acc] = rec
                    break
    if len(extracted) < 3:
        return False
    temp_fasta = temp_dir / f"temp_diag_{candidate_id}.fasta"
    sanitized_recs = []
    for acc in [candidate_id, p1_id, p2_id]:
        rec = extracted[acc]
        rec_copy = rec.__class__(
            id=rec.id, name=rec.name, description=rec.description,
            seq=rec.seq.__class__(sanitize_sequence_for_openrdp(str(rec.seq)))
        )
        sanitized_recs.append(rec_copy)
    SeqIO.write(sanitized_recs, temp_fasta, "fasta")
    try:
        scanner = openrdp.Scanner(methods=("rdp", "maxchi", "chimaera"), verbose=False)
        results = scanner.run_scans(str(temp_fasta))
        for method in results.keys():
            events = results[method]
            for e in events:
                if e.get("pvalue", 1.0) < 0.05:
                    return True
        return False
    except Exception:
        return False
    finally:
        if temp_fasta.exists():
            try: temp_fasta.unlink()
            except Exception: pass

def main():
    parser = argparse.ArgumentParser(description="Diagnose 3seq recombination validation")
    parser.add_argument("--virus", default="hbv", help="Virus type")
    parser.add_argument("--results-dir", default="results/hbv", help="Results directory")
    args = parser.parse_args()
    
    base_dir = Path(args.results_dir)
    recomb_summary_file = base_dir / "all_recombination.tsv"
    
    if not recomb_summary_file.exists():
        log.error("Recombination summary file %s not found.", recomb_summary_file)
        sys.exit(1)
        
    candidates = []
    with recomb_summary_file.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("is_recombinant") == "true":
                candidates.append(row)
                
    log.info("Loaded %d candidate recombinant sequences.", len(candidates))
    
    stats = {
        "total": len(candidates),
        "long_rec_passed": 0,
        "dist_orientation_a": 0,
        "dist_orientation_b": 0,
        "dist_passed_either": 0,
        "openrdp_passed_after_either": 0,
        "openrdp_passed_after_a": 0,
        "total_failures_detailed": {
            "not_in_longrec": 0,
            "failed_both_dist": 0,
            "failed_openrdp": 0
        }
    }
    
    scratch_dir = Path("scratch/recomb_val_diag")
    scratch_dir.mkdir(parents=True, exist_ok=True)
    
    limit = 500  # Diagnose up to 500 candidates
    tested = 0
    
    for cand in candidates[:limit]:
        tested += 1
        candidate_id = cand["sequence_id"]
        
        # Locate genotype from filename
        genotype = None
        for p in base_dir.glob("recombination/*/recombinants.tsv"):
            with p.open() as r_fh:
                r_reader = csv.DictReader(r_fh, delimiter="\t")
                for r_row in r_reader:
                    if r_row["sequence_id"] == candidate_id and r_row["is_recombinant"] == "true":
                        genotype = p.parent.name
                        break
            if genotype: break
            
        if not genotype:
            continue
            
        longrec_path = base_dir / "recombination" / genotype / f"run_{genotype}.3s.longRec"
        if not longrec_path.exists():
            continue
            
        long_rec_set = load_long_rec_candidates(longrec_path)
        if candidate_id not in long_rec_set:
            stats["total_failures_detailed"]["not_in_longrec"] += 1
            continue
            
        stats["long_rec_passed"] += 1
        
        aligned_path = base_dir / "alignments" / f"{genotype}_aligned.fasta"
        p1_id = cand["parent_1"]
        p2_id = cand["parent_2"]
        start = int(cand["breakpoint_start"])
        end = int(cand["breakpoint_end"])
        
        extracted_seqs = get_alignment_sequences(aligned_path, [candidate_id, p1_id, p2_id])
        if candidate_id not in extracted_seqs or p1_id not in extracted_seqs or p2_id not in extracted_seqs:
            continue
            
        cand_seq = extracted_seqs[candidate_id]
        p1_seq = extracted_seqs[p1_id]
        p2_seq = extracted_seqs[p2_id]
        
        ident_p1_r = calculate_identity(cand_seq, p1_seq, start, end, invert=False)
        ident_p2_r = calculate_identity(cand_seq, p2_seq, start, end, invert=False)
        ident_p1_nr = calculate_identity(cand_seq, p1_seq, start, end, invert=True)
        ident_p2_nr = calculate_identity(cand_seq, p2_seq, start, end, invert=True)
        
        # Orientation A: P2 closer in recombinant region, P1 closer in non-recombinant region
        passed_a = (ident_p2_r > ident_p1_r) and (ident_p1_nr > ident_p2_nr)
        
        # Orientation B: P1 closer in recombinant region, P2 closer in non-recombinant region
        passed_b = (ident_p1_r > ident_p2_r) and (ident_p2_nr > ident_p1_nr)
        
        if passed_a:
            stats["dist_orientation_a"] += 1
        if passed_b:
            stats["dist_orientation_b"] += 1
            
        if passed_a or passed_b:
            stats["dist_passed_either"] += 1
            
            # Run OpenRDP
            passed_openrdp = run_openrdp_validation(aligned_path, candidate_id, p1_id, p2_id, scratch_dir)
            if passed_openrdp:
                stats["openrdp_passed_after_either"] += 1
                if passed_a:
                    stats["openrdp_passed_after_a"] += 1
            else:
                stats["total_failures_detailed"]["failed_openrdp"] += 1
        else:
            stats["total_failures_detailed"]["failed_both_dist"] += 1
            
    # Print results
    print("\n" + "="*80)
    print(f"DIAGNOSTIC REPORT FOR HBV RECOMBINATION VALIDATION (Sampled first {tested} candidates)")
    print("="*80)
    print(f"Total Candidates Sampled         : {tested}")
    print(f"Step 1 Passed (LongRec filter)   : {stats['long_rec_passed']} / {tested} (failed: {stats['total_failures_detailed']['not_in_longrec']})")
    print(f"Passed Distance Orientation A    : {stats['dist_orientation_a']} / {stats['long_rec_passed']}")
    print(f"Passed Distance Orientation B    : {stats['dist_orientation_b']} / {stats['long_rec_passed']}")
    print(f"Passed Either Distance Check     : {stats['dist_passed_either']} / {stats['long_rec_passed']} (failed both: {stats['total_failures_detailed']['failed_both_dist']})")
    print(f"OpenRDP Passed (from Either)     : {stats['openrdp_passed_after_either']} / {stats['dist_passed_either']} (failed: {stats['total_failures_detailed']['failed_openrdp']})")
    print(f"OpenRDP Passed (from A only)     : {stats['openrdp_passed_after_a']} / {stats['dist_orientation_a']}")
    print("="*80 + "\n")

if __name__ == "__main__":
    main()
