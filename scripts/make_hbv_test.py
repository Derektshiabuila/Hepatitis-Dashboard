import os
import subprocess
from pathlib import Path

def parse_fasta(filepath):
    sequences = []
    current_header = None
    current_seq = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences.append((current_header, ''.join(current_seq)))
    return sequences

def main():
    project_root = Path(__file__).parent.parent.resolve()
    
    # Paths
    ref_path = project_root / "refs" / "hbv" / "HBV_refs_GLUE.fasta"
    seq_path = project_root / "results" / "hbv" / "hbv_sequences.fasta"
    
    out_dir = project_root / "test_data" / "hbv"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    combined_path = out_dir / "test_hbv_combined.fasta"
    aligned_path = out_dir / "test_hbv_aligned.fasta"
    
    # 1. Read HBV refs (all of them)
    print(f"Reading references from {ref_path}...")
    refs = parse_fasta(ref_path)
    print(f"Loaded {len(refs)} reference sequences.")
    
    # 2. Read first 10 HBV samples
    print(f"Reading samples from {seq_path}...")
    samples = []
    count = 0
    with open(seq_path, 'r') as f:
        current_header = None
        current_seq = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header:
                    samples.append((current_header, ''.join(current_seq)))
                    count += 1
                    if count >= 10:
                        break
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header and count < 10:
            samples.append((current_header, ''.join(current_seq)))
            
    print(f"Loaded {len(samples)} sample sequences.")
    
    # 3. Combine them
    print(f"Writing combined FASTA to {combined_path}...")
    with open(combined_path, 'w') as f:
        for header, seq in refs:
            f.write(f"{header}\n{seq}\n")
        for header, seq in samples:
            f.write(f"{header}\n{seq}\n")
            
    # 4. Run MAFFT alignment
    print("Running MAFFT to align sequences...")
    cmd = ["mafft", "--auto", str(combined_path)]
    with open(aligned_path, 'w') as out_f:
        res = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
    
    if res.returncode == 0:
        print(f"Successfully created aligned test file: {aligned_path}")
    else:
        print(f"Error running MAFFT: {res.stderr}")
        
    # Copy to true_runner test_data directory if it exists
    true_runner_out = project_root / "true_runner" / "test_data" / "hbv"
    true_runner_out.mkdir(parents=True, exist_ok=True)
    
    # copy files
    import shutil
    shutil.copy2(combined_path, true_runner_out / "test_hbv_combined.fasta")
    shutil.copy2(aligned_path, true_runner_out / "test_hbv_aligned.fasta")
    print(f"Copied test files to true_runner: {true_runner_out}")

if __name__ == '__main__':
    main()
