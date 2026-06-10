#!/bin/bash
# scripts/run_glue.sh
# Usage: ./run_glue.sh <fasta_file> <virus> <output_dir>

set -euo pipefail

FASTA_FILE=$1
VIRUS=$2
OUTPUT_DIR=$3

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Copy FASTA to work directory (create split_fasta subdir as expected by GLUE)
mkdir -p "$OUTPUT_DIR/split_fasta"

# Split user FASTA into individual sequences and build glue_cmd.glue
python3 -c "
import os
from pathlib import Path

fasta_file = Path('$FASTA_FILE')
out_dir = Path('$OUTPUT_DIR')
split_dir = out_dir / 'split_fasta'
split_dir.mkdir(parents=True, exist_ok=True)

# Read and split sequences
seqs = {}
with open(fasta_file) as f:
    curr_id = None
    curr_seq = []
    for line in f:
        if line.startswith('>'):
            if curr_id:
                seqs[curr_id] = ''.join(curr_seq)
            curr_id = line[1:].strip().split()[0]
            # Replace characters that might be problematic in filenames or GLUE
            curr_id = ''.join(c for c in curr_id if c.isalnum() or c in '._-')
            curr_seq = []
        else:
            curr_seq.append(line.strip())
    if curr_id:
        seqs[curr_id] = ''.join(curr_seq)

# Determine GLUE module based on virus
if '$VIRUS' == 'hbv':
    glue_module = 'hdrReportingController invoke-function reportFastaCli'
elif '$VIRUS' == 'hcv':
    glue_module = 'phdrReportingController invoke-function reportFasta'
else:
    glue_module = 'hevMaxLikelihoodGenotyper genotype file -f'

# Write split fasta files and build GLUE command list
glue_cmds = [
    'project $VIRUS',
    'console set log-level INFO'
]

for sid, seq in seqs.items():
    seq_path = split_dir / f'{sid}.fasta'
    with open(seq_path, 'w') as out_f:
        out_f.write(f'>{sid}\n{seq}\n')
    
    # Add GLUE commands for this sequence
    glue_cmds.append(f'console set next-cmd-output-file /work/{sid}.xml')
    glue_cmds.append(f'module {glue_module} /work/split_fasta/{sid}.fasta')

# Write glue_cmd.glue
with open(out_dir / 'glue_cmd.glue', 'w') as cmd_f:
    cmd_f.write('\n'.join(glue_cmds) + '\nexit\n')
"

# Ensure the MySQL container for this virus is running
MYSQL_CONTAINER="gluetools-mysql-${VIRUS}"
if ! docker ps --format '{{.Names}}' | grep -q "^${MYSQL_CONTAINER}$"; then
    echo "MySQL container '${MYSQL_CONTAINER}' is not running. Starting it..." >&2
    docker start "$MYSQL_CONTAINER" 2>/dev/null || {
        echo "ERROR: Could not start ${MYSQL_CONTAINER}. Is it created?" >&2
        echo "Try: docker start ${MYSQL_CONTAINER}" >&2
        exit 1
    }
fi

# Wait for MySQL to accept connections (up to 30 seconds)
echo "Waiting for MySQL in ${MYSQL_CONTAINER} to be ready..."
for i in $(seq 1 30); do
    docker exec "$MYSQL_CONTAINER" mysqladmin ping -h localhost --silent 2>/dev/null && break
    sleep 1
done

# Run GLUE using Docker, linked to the virus-specific MySQL container
docker run --rm \
    --platform linux/amd64 \
    --link "${MYSQL_CONTAINER}:gluetools-mysql" \
    -v "$OUTPUT_DIR:/work" \
    cvrbioinformatics/gluetools:latest \
    gluetools.sh -f /work/glue_cmd.glue -n

# Check if any XML output was created
if ls "$OUTPUT_DIR"/*.xml >/dev/null 2>&1; then
    echo "GLUE analysis completed successfully"
else
    echo "WARNING: GLUE did not produce any XML reports" >&2
fi