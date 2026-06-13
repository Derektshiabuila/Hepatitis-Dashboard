import yaml
import glob
import os

configfile: "config.yaml"

cfg = config["viruses"]
VIRUSES = list(cfg.keys())

# Define functions to get samples dynamically
def get_samples(virus):
    fasta_files = glob.glob(f"results/{virus}/split_fasta/*.fasta")
    return [os.path.splitext(os.path.basename(f))[0] for f in fasta_files]

def samples_from_checkpoint(wildcards):
    ckpt = checkpoints.split_fasta_per_sequence.get(virus=wildcards.virus)
    fasta_dir = ckpt.output[0]
    return glob_wildcards(f"{fasta_dir}/{{sample}}.fasta").sample

# ============================================
# HELPER FUNCTION: Get genotypes from checkpoint
# ============================================
def get_recombination_inputs(wildcards):
    """Dynamically determines required RDP5 outputs based on checkpoint results"""
    ckpt = checkpoints.get_recombination_files.get(virus=wildcards.virus)
    checkpoint_dir = ckpt.output[0]
    list_file = os.path.join(checkpoint_dir, "_all_genotypes.txt")
    
    if not os.path.exists(list_file):
        return []
    
    with open(list_file, 'r') as f:
        genotypes = [line.strip() for line in f if line.strip()]
        
    # This properly forces Snakemake to track and execute rdp5_recombination per genotype
    return expand("results/{virus}/recombination/{genotype}/recombinants.tsv", 
                  virus=wildcards.virus, genotype=genotypes)


# Final outputs rule
rule all:
    input:
        expand("results/{virus}/final_resistance.tsv", virus=["hbv", "hcv", "hev"]),
        expand("results/{virus}/all_recombination.tsv", virus=["hbv", "hcv", "hev"]),
        expand("results/{virus}/metadata_with_recomb.tsv", virus=VIRUSES),


###########################################
# 1 — DOWNLOAD (NCBI)
###########################################
rule download_sequences:
    output:
        fasta="results/{virus}/{virus}_sequences.fasta",
        meta="results/{virus}/{virus}_metadata.tsv"
    params:
        query=lambda w: cfg[w.virus]["ncbi_query"],
        email=config["downloads"]["email"],
        virus="{virus}",
        minlength=lambda w: 2900 if w.virus == "hbv" else 8000,
        maxlength=lambda w: 3500 if w.virus == "hbv" else 11000
    conda: "envs/python.yaml"
    script:
        "scripts/download_sequences.py"
        
###########################################
# Population retriever
###########################################
rule retrieve_population:
    output:
        "data/population_by_country_year.csv"
    conda: "envs/python.yaml"
    script:
        "scripts/population_retriever.py"

###########################################
# Split fasta for GLUE
###########################################
checkpoint split_fasta_per_sequence:
    input:
        "results/{virus}/{virus}_sequences.fasta"
    output:
        directory("results/{virus}/split_fasta")
    conda: "envs/python.yaml"
    shell:
        r"""
        mkdir -p {output}
        awk '
          /^>/ {{
            if (out) close(out);
            name = substr($1, 2);
            out = "{output}/" name ".fasta"
          }}
          {{ print > out }}
        ' {input}
        """

###########################################
# 5 — HBV GLUE Mutation Analysis
###########################################
rule glue_analysis:
    input:
        fasta="results/{virus}/split_fasta/{sample}.fasta"
    output:
        xml="results/{virus}/glue/{sample}.xml"
    threads: 4
    resources:
        gluemysql=1
    params:
        project=lambda w: w.virus,
        glue_module=lambda w:
            "hdrReportingController invoke-function reportFastaCli" if w.virus == "hbv"
            else "phdrReportingController invoke-function reportFasta" if w.virus == "hcv"
            else "hevMaxLikelihoodGenotyper genotype file -f",
        retries=1,
        sleep_time=10
    shell:
        r"""
        mkdir -p results/{wildcards.virus}/glue
        SAMPLE="{wildcards.sample}"
        
        MAX_RETRIES={params.retries}
        SLEEP_TIME={params.sleep_time}
        RETRY_COUNT=0
        
        for i in {{1..30}}; do
            docker exec gluetools-mysql-{wildcards.virus} \
                mysqladmin ping -h localhost --silent && break
            sleep 1
        done
        
        cat > results/{wildcards.virus}/gluetools-config_$SAMPLE.xml << EOF
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

        while [ $RETRY_COUNT -le $MAX_RETRIES ]; do
        cat > results/{wildcards.virus}/glue_cmd_$SAMPLE.glue << EOF
project {params.project}
console set log-level INFO
console set next-cmd-output-file /work/result.xml
module {params.glue_module} /work/split_fasta/$SAMPLE.fasta
EOF
        
            timeout --signal=SIGKILL 600s docker run --rm \
              --platform linux/amd64 \
              --link gluetools-mysql-{wildcards.virus}:gluetools-mysql \
              -e MAFFT_NTHREAD=1 \
              -e OMP_NUM_THREADS=1 \
              -e OPENBLAS_NUM_THREADS=1 \
              -v $(realpath results/{wildcards.virus}):/work \
              cvrbioinformatics/gluetools:latest \
              java -jar /opt/gluetools/lib/gluetools-core.jar -c /work/gluetools-config_$SAMPLE.xml -f /work/glue_cmd_$SAMPLE.glue -n || true
        
            if [ -f "results/{wildcards.virus}/result.xml" ]; then
                cp results/{wildcards.virus}/result.xml {output.xml}
                rm -f results/{wildcards.virus}/result.xml
                rm -f results/{wildcards.virus}/gluetools-config_$SAMPLE.xml results/{wildcards.virus}/glue_cmd_$SAMPLE.glue
                exit 0
            fi
        
            RETRY_COUNT=$((RETRY_COUNT+1))
            sleep $SLEEP_TIME
        done
        rm -f results/{wildcards.virus}/gluetools-config_$SAMPLE.xml results/{wildcards.virus}/glue_cmd_$SAMPLE.glue

        VIRUS="{wildcards.virus}"
        echo "GLUE failed for $SAMPLE — writing empty XML" >&2

        cat > {output.xml} << EOF
<glueFailure>
  <sample>$SAMPLE</sample>
  <virus>$VIRUS</virus>
  <reason>hdrReportingController null input</reason>
</glueFailure>
EOF
        exit 0
        """

###########################################
# Generate mutation tsv
###########################################
rule parse_glue_xml:
    input:
        xml="results/{virus}/glue/{sample}.xml"
    output:
        tsv="results/{virus}/resistance/{sample}.tsv"
    conda: "envs/python.yaml"
    shell:
        r"""
        mkdir -p results/{wildcards.virus}/resistance

        if [ "{wildcards.virus}" = "hbv" ]; then
            python scripts/parse_hbv_glue_xml.py {input.xml} {output.tsv}
        elif [ "{wildcards.virus}" = "hcv" ]; then
            python scripts/parse_hcv_glue_xml.py {input.xml} {output.tsv}
        else
            python scripts/parse_hev_glue_xml.py {input.xml} {output.tsv}
        fi

        python scripts/validate_schema.py {output.tsv}
        """

############################################
# Merge tables
############################################        
rule merge_resistance_tables:
    input:
        lambda w: expand(
            "results/{virus}/resistance/{sample}.tsv",
            virus=w.virus,
            sample=samples_from_checkpoint(w)
        )
    output:
        "results/{virus}/final_resistance.tsv"
    shell:
        r"""
        set -euo pipefail
        header_file=$(awk 'FNR==1 && NF>1 {{print FILENAME; exit}}' {input})

        if [ -z "$header_file" ]; then
            echo "[ERROR] No valid TSV files found for {wildcards.virus}" >&2
            exit 1
        fi

        head -n 1 "$header_file" > {output}
        for f in {input}; do
            if [ -s "$f" ]; then
                tail -n +2 "$f" >> {output} || true
            fi
        done
        """

# ============================================
# CHECKPOINT: Get recombination files list
# ============================================
checkpoint get_recombination_files:
    input:
        "results/{virus}/final_resistance.tsv"
    output:
        directory("results/{virus}/recombination_files_list")
    conda: "envs/python.yaml"
    run:
        import pandas as pd
        import os
        
        print(f"Reading resistance file: {input[0]}")
        df = pd.read_csv(input[0], sep='\t', dtype=str)
        
        genotype_col = None
        for col in ['genotype', 'Genotype', 'GT', 'gt', 'genotype_assignment']:
            if col in df.columns:
                genotype_col = col
                break
        
        if genotype_col is None:
            os.makedirs(output[0], exist_ok=True)
            with open(os.path.join(output[0], "_all_genotypes.txt"), 'w') as f:
                pass
            return
        
        genotypes = df[genotype_col].astype(str).str.strip().dropna().unique().tolist()
        valid_genotypes = [str(g).strip() for g in genotypes if str(g).strip().lower() not in ['nan', 'none', 'unknown', '']]
        
        os.makedirs(output[0], exist_ok=True)
        with open(os.path.join(output[0], "_all_genotypes.txt"), 'w') as f:
            for g in valid_genotypes:
                f.write(f"{g}\n")


# ============================================
# GENOTYPE ALIGNMENT (MAFFT)
# Triggered per-genotype via the checkpoint
# ============================================
rule genotype_alignment:
    input:
        resistance="results/{virus}/final_resistance.tsv",
        sequences="results/{virus}/{virus}_sequences.fasta",
        refs=lambda w: (f"refs/{w.virus}/ref_msa.fasta"
                        if os.path.exists(f"refs/{w.virus}/ref_msa.fasta")
                        else [])
    output:
        aligned="results/{virus}/alignments/{genotype}_aligned.fasta"
    params:
        genotype="{genotype}",
        virus="{virus}"
    threads: 4
    conda: "envs/mafft.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/{wildcards.virus}/alignments

        GENO_CLEAN=$(echo "{params.genotype}" | tr ' ' '_')

        python -c "
import pandas as pd, sys
df = pd.read_csv('{input.resistance}', sep='\t', dtype=str)
genotype_col = next((c for c in ['genotype','Genotype','GT','gt','genotype_assignment'] if c in df.columns), None)
sample_col   = next((c for c in ['sample','Sample','sequence_id','accession','sequence_name'] if c in df.columns), None)
if not genotype_col or not sample_col:
    print('ERROR: missing genotype or sample column', file=sys.stderr); sys.exit(1)
ids = df[df[genotype_col].astype(str).str.strip() == '{params.genotype}'][sample_col].dropna().astype(str).str.strip().tolist()
open('/tmp/{params.virus}_' + sys.argv[1] + '_samples.txt', 'w').write('\n'.join(ids) + '\n')
print('Found %d samples for genotype %s' % (len(ids), '{params.genotype}'))
" "$GENO_CLEAN"

        if [ -s "/tmp/{params.virus}_${{GENO_CLEAN}}_samples.txt" ]; then
            seqkit grep -f "/tmp/{params.virus}_${{GENO_CLEAN}}_samples.txt" \
                {input.sequences} > "/tmp/{params.virus}_${{GENO_CLEAN}}_seqs.fasta"
        else
            touch "/tmp/{params.virus}_${{GENO_CLEAN}}_seqs.fasta"
            echo "Warning: no samples for genotype {params.genotype}" >&2
        fi

        REF="{input.refs}"
        SEQS="/tmp/{params.virus}_${{GENO_CLEAN}}_seqs.fasta"
        COMBINED="/tmp/{params.virus}_${{GENO_CLEAN}}_combined.fasta"

        HAS_REF=0
        if [ -f "$REF" ] && [ -s "$REF" ]; then
            HAS_REF=1
        fi

        HAS_SEQS=0
        if [ -s "$SEQS" ]; then
            HAS_SEQS=1
        fi

        ALIGN_TARGET=""
        if [ "$HAS_REF" -eq 1 ] && [ "$HAS_SEQS" -eq 1 ]; then
            cat "$REF" "$SEQS" > "$COMBINED"
            ALIGN_TARGET="$COMBINED"
        elif [ "$HAS_REF" -eq 1 ]; then
            ALIGN_TARGET="$REF"
        elif [ "$HAS_SEQS" -eq 1 ]; then
            ALIGN_TARGET="$SEQS"
        fi

        if [ -n "$ALIGN_TARGET" ] && [ -s "$ALIGN_TARGET" ]; then
            n_seqs=$(grep -c "^>" "$ALIGN_TARGET" || true)
            if [ "$n_seqs" -ge 2 ]; then
                echo "Aligning $n_seqs sequences for genotype {params.genotype} using MAFFT..." >&2
                mafft --thread {threads} --auto "$ALIGN_TARGET" > "{output.aligned}"
            elif [ "$n_seqs" -eq 1 ]; then
                echo "Only 1 sequence found for genotype {params.genotype}. Copying directly without alignment." >&2
                cat "$ALIGN_TARGET" > "{output.aligned}"
            else
                printf "" > "{output.aligned}"
            fi
        else
            printf "" > "{output.aligned}"
            echo "Warning: no sequences found for genotype {params.genotype}" >&2
        fi

        rm -f "/tmp/{params.virus}_${{GENO_CLEAN}}_"*.fasta \
              "/tmp/{params.virus}_${{GENO_CLEAN}}_samples.txt" 2>/dev/null || true
        """


# ============================================
# RDP5 RECOMBINATION
# ============================================
rule rdp5_recombination:
    input:
        aligned="results/{virus}/alignments/{genotype}_aligned.fasta"
    output:
        tsv="results/{virus}/recombination/{genotype}/recombinants.tsv"
    params:
        outdir  = "results/{virus}/recombination/{genotype}",
        genotype= "{genotype}",
        virus   = "{virus}",
    threads: 1
    resources:
        rdp5=1
    log:
        "logs/{virus}/rdp5_{genotype}.log"
    conda: "envs/recombination.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}" "logs/{wildcards.virus}"

        n_seqs=$(grep -c "^>" "{input.aligned}" || true)
        if [ "$n_seqs" -lt 4 ]; then
            echo "$(date) Skipping {params.genotype}: only $n_seqs sequence(s) (need ≥4)" | tee -a "{log}"
            printf "sequence_id\tis_recombinant\tbreakpoint_start\tbreakpoint_end\tp_value\tmethods\tparent_1\tparent_2\n" > "{output.tsv}"
            exit 0
        fi

        echo "$(date) Running RDP5 on {params.genotype} ($n_seqs sequences)" | tee -a "{log}"
        python scripts/run_rdp5.py \
            --fasta    "{input.aligned}" \
            --virus    "{params.virus}" \
            --outdir   "{params.outdir}" \
            --genotype "{params.genotype}" >> "{log}" 2>&1
        """


# ============================================
# MERGE RECOMBINATION (Aggregates Checkpoint Data)
# ============================================
rule merge_recombination:
    input:
        # Using the lambda function intercepts the execution graph correctly
        recomb_files = lambda wildcards: get_recombination_inputs(wildcards)
    output:
        "results/{virus}/all_recombination.tsv"
    log:
        "logs/{virus}/merge_recombination.log"
    run:
        import os
        
        # Write default header if no files exist
        header_str = "sequence_id\tis_recombinant\tbreakpoint_start\tbreakpoint_end\tp_value\tmethods\tparent_1\tparent_2\n"
        
        if not input.recomb_files:
            with open(output[0], 'w') as f:
                f.write(header_str)
            return

        # Write header from first available file
        with open(input.recomb_files[0], 'r') as src:
            header = src.readline()
            
        with open(output[0], 'w') as dst:
            dst.write(header)
            
        # Append all details
        for recomb_file in input.recomb_files:
            if os.path.exists(recomb_file):
                with open(recomb_file, 'r') as src:
                    lines = src.readlines()
                    if len(lines) > 1:
                        with open(output[0], 'a') as dst:
                            dst.writelines(lines[1:])


# ============================================
# ANNOTATE RECOMBINATION
# ============================================
rule annotate_recombination:
    input:
        resistance = "results/{virus}/final_resistance.tsv",
        recomb     = "results/{virus}/all_recombination.tsv"
    output:
        "results/{virus}/metadata_with_recomb.tsv"
    log:
        "logs/{virus}/annotate_recombination.log"
    conda: "envs/python.yaml"
    shell:
        r"""
        set -euo pipefail
        python - <<'PYEOF' 2>&1 | tee -a {log}
import pandas as pd
import sys

res_df = pd.read_csv("{input.resistance}", sep="\t", dtype=str)
rec_df = pd.read_csv("{input.recomb}", sep="\t", dtype=str)

id_col = None
for col in ["sample", "Sample", "sequence_id", "accession", "sequence_name"]:
    if col in res_df.columns:
        id_col = col
        break

if id_col is None:
    sys.exit(1)

if "sequence_id" in rec_df.columns:
    rec_df = rec_df.rename(columns={{"sequence_id": id_col}})

available_cols = [col for col in ["is_recombinant", "breakpoint_start", "breakpoint_end", "p_value", "methods", "parent_1", "parent_2"] if col in rec_df.columns]
if id_col not in available_cols:
    available_cols.insert(0, id_col)

merged = res_df.merge(rec_df[available_cols], on=id_col, how="left")
merged["is_recombinant"] = merged["is_recombinant"].fillna("false")
merged.to_csv("{output}", sep="\t", index=False)
PYEOF
        """

# ============================================
# VALIDATE RECOMBINATION
# ============================================
rule validate_recombination:
    input:
        recomb = "results/{virus}/all_recombination.tsv"
    output:
        val = "results/{virus}/validated_recombinants.tsv"
    log:
        "logs/{virus}/validate_recombination.log"
    conda: "envs/validation.yaml"
    shell:
        r"""
        set -euo pipefail
        python scripts/validate_recombination.py --virus {wildcards.virus} 2>&1 | tee {log}
        """

###########################################
# 7 — FINAL MERGE
###########################################
rule merge_everything:
    input:
        tree="results/{virus}/aligned.fasta.treefile",
        recomb="results/{virus}/metadata_with_recomb.tsv",
    output:
        "results/{virus}/final_report.tsv"
    conda: "envs/python.yaml"
    script:
        "scripts/merge_everything.py"