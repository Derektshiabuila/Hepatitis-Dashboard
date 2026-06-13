#!/usr/bin/env python3
# scripts/debug_glue.py
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
print(f"Project root: {PROJECT_ROOT}")

# Create a dummy HCV FASTA file
dummy_fasta = PROJECT_ROOT / "temp_dummy_hcv.fasta"
dummy_fasta.write_text(""">AB049090_dummy
ATGGCGTGGGGGATGGGGATGATGGTGGTGGGGGGGGTGGTGGGGGGGGTGGTGGGGGGG
""")

output_dir = PROJECT_ROOT / "results" / "debug_job"
output_dir.mkdir(parents=True, exist_ok=True)

glue_wrapper = PROJECT_ROOT / "scripts" / "run_glue.sh"
print(f"Glue wrapper path: {glue_wrapper}")
print(f"Exists: {glue_wrapper.exists()}")

cmd = ["bash", str(glue_wrapper), str(dummy_fasta), "hcv", str(output_dir)]
print(f"Running command: {' '.join(cmd)}")

try:
    result = subprocess.run(
        cmd,
        cwd=str(PROJECT_ROOT),
        capture_output=True,
        text=True,
        timeout=300,
    )
    print("\n--- RETURN CODE ---")
    print(result.returncode)
    print("\n--- STDOUT ---")
    print(result.stdout)
    print("\n--- STDERR ---")
    print(result.stderr)
finally:
    # Clean up
    if dummy_fasta.exists():
        dummy_fasta.unlink()
    import shutil
    shutil.rmtree(output_dir, ignore_errors=True)
