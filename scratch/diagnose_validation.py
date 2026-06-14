import sys
import csv
import logging
from pathlib import Path
from Bio import SeqIO
import argparse

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-8s  %(message)s")
log = logging.getLogger("diagnose_val")

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

def calculate_identity_detailed(seq1: str, seq2: str, start: int, end: int, invert: bool = False):
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
    pct = (matches / total) if total > 0 else 0.0
    return pct, total

def main():
    parser = argparse.ArgumentParser(description="Diagnose 3seq recombination validation detail")
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
    
    print("\n" + "="*100)
    print("DETAILED DATA FOR FIRST 10 CANDIDATES:")
    print("="*100)
    
    printed_count = 0
    for cand in candidates:
        if printed_count >= 15:
            break
            
        candidate_id = cand["sequence_id"]
        
        # Locate genotype
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
            continue
            
        aligned_path = base_dir / "alignments" / f"{genotype}_aligned.fasta"
        p1_id = cand["parent_1"]
        p2_id = cand["parent_2"]
        start = int(cand["breakpoint_start"])
        end = int(cand["breakpoint_end"])
        
        extracted_seqs = get_alignment_sequences(aligned_path, [candidate_id, p1_id, p2_id])
        if candidate_id not in extracted_seqs or p1_id not in extracted_seqs or p2_id not in extracted_seqs:
            print(f"Candidate {candidate_id}: Alignment sequences not found!")
            continue
            
        cand_seq = extracted_seqs[candidate_id]
        p1_seq = extracted_seqs[p1_id]
        p2_seq = extracted_seqs[p2_id]
        
        ident_p1_r, tot_p1_r = calculate_identity_detailed(cand_seq, p1_seq, start, end, invert=False)
        ident_p2_r, tot_p2_r = calculate_identity_detailed(cand_seq, p2_seq, start, end, invert=False)
        ident_p1_nr, tot_p1_nr = calculate_identity_detailed(cand_seq, p1_seq, start, end, invert=True)
        ident_p2_nr, tot_p2_nr = calculate_identity_detailed(cand_seq, p2_seq, start, end, invert=True)
        
        passed_a = (ident_p2_r > ident_p1_r) and (ident_p1_nr > ident_p2_nr)
        passed_b = (ident_p1_r > ident_p2_r) and (ident_p2_nr > ident_p1_nr)
        
        print(f"Candidate  : {candidate_id} (Genotype: {genotype})")
        print(f"Parents    : P1={p1_id}, P2={p2_id}")
        print(f"Breakpoints: Start={start}, End={end} (Length: {end - start})")
        print(f"Alignment L: {len(cand_seq)}")
        print(f"Recomb Region Identity:")
        print(f"  - to P1: {ident_p1_r:.4f} (over {tot_p1_r} bases)")
        print(f"  - to P2: {ident_p2_r:.4f} (over {tot_p2_r} bases)")
        print(f"Non-Recomb Region Identity:")
        print(f"  - to P1: {ident_p1_nr:.4f} (over {tot_p1_nr} bases)")
        print(f"  - to P2: {ident_p2_nr:.4f} (over {tot_p2_nr} bases)")
        print(f"Validation Status:")
        print(f"  - Passed A (P2 in recomb, P1 outside): {passed_a}")
        print(f"  - Passed B (P1 in recomb, P2 outside): {passed_b}")
        print("-" * 50)
        
        printed_count += 1

if __name__ == "__main__":
    main()
