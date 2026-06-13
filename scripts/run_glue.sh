#!/bin/bash
# scripts/run_glue.sh
# Usage: ./run_glue.sh <fasta_file> <virus> <output_dir>

set -euo pipefail

FASTA_FILE=$1
VIRUS=$2
OUTPUT_DIR=$3

# Create output directory
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")

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

# Generate a temporary config file with allowPublicKeyRetrieval=true to bypass MySQL 8.0/driver connection issues
cat > "$OUTPUT_DIR/gluetools-config.xml" << 'EOF'
<gluetools>
	<database>
		<username>gluetools</username>
		<password>glue12345</password>
		<jdbcUrl>jdbc:mysql://gluetools-mysql:3306/GLUE_TOOLS?characterEncoding=UTF-8&amp;useSSL=false&amp;allowPublicKeyRetrieval=true</jdbcUrl>
	</database>
	<properties>
		<!-- BLAST related config -->
	    <property>
			<name>gluetools.core.programs.blast.blastn.executable</name>
			<value>/opt/gluetools/blast/ncbi-blast-2.2.31+/bin/blastn</value>
		</property>
		<property>
			<name>gluetools.core.programs.blast.tblastn.executable</name>
			<value>/opt/gluetools/blast/ncbi-blast-2.2.31+/bin/tblastn</value>
		</property>
		<property>
			<name>gluetools.core.programs.blast.makeblastdb.executable</name>
			<value>/opt/gluetools/blast/ncbi-blast-2.2.31+/bin/makeblastdb</value>
		</property>
		<property>
			<name>gluetools.core.programs.blast.search.threads</name>
			<value>4</value>
		</property>
		<property>
			<name>gluetools.core.programs.blast.temp.dir</name>
			<value>/opt/gluetools/tmp/blastfiles</value>
		</property>
		<property>
			<name>gluetools.core.programs.blast.db.dir</name>
			<value>/opt/gluetools/tmp/blastdbs</value>
		</property>
		<!-- RAxML-specific config -->
		<property>
			<name>gluetools.core.programs.raxml.raxmlhpc.executable</name>
			<value>/opt/gluetools/raxml/bin/raxmlHPC-PTHREADS-SSE3</value>
		</property>
		<property>
			<name>gluetools.core.programs.raxml.raxmlhpc.cpus</name>
			<value>4</value>
		</property>
		<property>
			<name>gluetools.core.programs.raxml.temp.dir</name>
			<value>/opt/gluetools/tmp/raxmlfiles</value>
		</property>
		<!-- MAFFT-specific config -->
		<property>
			<name>gluetools.core.programs.mafft.executable</name>
			<value>/usr/local/bin/mafft</value>
		</property>
		<property>
			<name>gluetools.core.programs.mafft.cpus</name>
			<value>4</value>
		</property>
		<property>
			<name>gluetools.core.programs.mafft.temp.dir</name>
			<value>/opt/gluetools/tmp/mafftfiles</value>
		</property>
		<!-- JModelTest -->
		<property>
			<name>gluetools.core.programs.jmodeltester.jar</name>
			<value>/opt/gluetools/jModelTest/lib/jModelTest.jar</value>
		</property>
		<property>
			<name>gluetools.core.programs.jmodeltester.temp.dir</name>
			<value>/opt/gluetools/tmp/jmodeltest</value>
		</property>
		<property>
			<name>gluetools.core.programs.jmodeltester.cpus</name>
			<value>4</value>
		</property>
		<!-- tbl2asn-->
		<property>
			<name>gluetools.core.programs.tbl2asn.executable</name>
			<value>/opt/gluetools/tbl2asn/bin/tbl2asn</value>
		</property>
		<property>
			<name>gluetools.core.programs.tbl2asn.temp.dir</name>
			<value>/opt/gluetools/tmp/tbl2asn</value>
		</property>
		<!-- ClusterPicker -->
		<property>
			<name>gluetools.core.programs.clusterPicker.jarPath</name>
			<value>/opt/gluetools/clusterPicker/lib/ClusterPicker_1.2.5.jar</value>
		</property>
		<property>
			<name>gluetools.core.programs.clusterPicker.temp.dir</name>
			<value>/opt/gluetools/tmp/clusterPicker</value>
		</property>
		<!-- SAM/BAM file processing -->
		<property>
			<name>gluetools.core.sam.temp.dir</name>
			<value>/opt/gluetools/tmp/sam</value>
		</property>
		<property>
			<name>gluetools.core.sam.cpus</name>
			<value>4</value>
		</property>
		<!-- Cayenne -->
		<property>
			<name>cayenne.querycache.size</name>
			<value>30000</value>
		</property>
	</properties>
</gluetools>
EOF

# Translate path if running inside a DooD container where /app maps to a host directory
HOST_OUTPUT_DIR="$OUTPUT_DIR"
if [ -n "${HEP_HOST_PROJECT_ROOT:-}" ]; then
    if [[ "$OUTPUT_DIR" == /app* ]]; then
        HOST_OUTPUT_DIR="${HEP_HOST_PROJECT_ROOT}${OUTPUT_DIR#/app}"
    fi
fi

# Run GLUE using Docker, linked to the virus-specific MySQL container
docker run --rm \
    --platform linux/amd64 \
    --link "${MYSQL_CONTAINER}:gluetools-mysql" \
    -v "$HOST_OUTPUT_DIR:/work" \
    cvrbioinformatics/gluetools:latest \
    java -jar /opt/gluetools/lib/gluetools-core.jar -c /work/gluetools-config.xml -f /work/glue_cmd.glue -n

rm -f "$OUTPUT_DIR/gluetools-config.xml"

# Check if any XML output was created
if ls "$OUTPUT_DIR"/*.xml >/dev/null 2>&1; then
    echo "GLUE analysis completed successfully"
else
    echo "WARNING: GLUE did not produce any XML reports" >&2
fi