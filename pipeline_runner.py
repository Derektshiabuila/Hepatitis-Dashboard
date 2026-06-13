"""
pipeline_runner.py
==================
Bridges the Dash user-sequence submission tab to the analysis backend.

Architecture
------------

  Dash callback  ──►  dispatch_user_pipeline()
                           │
                           ├─ writes FASTA + job metadata to jobs/<job_id>/
                           ├─ submits a background thread (or Celery task)
                           └─ returns job_id immediately

  dcc.Interval polls  ──►  get_job_status(job_id)
                               │
                               └─ reads jobs/<job_id>/state.json
                                    { state, progress, label, results? }

  On completion, results dict is:
  {
    "sequences": [
      { "id", "virus", "genotype", "is_recombinant",
        "breakpoints", "nearest_ref", "epa_score" }, ...
    ],
    "newick":       str,          # pruned subtree Newick
    "sequence_map": { seq_id: { orfs, mutations, breakpoints } }
  }

Pipeline steps (run inside a background thread)
------------------------------------------------
  1. Write user FASTA  →  jobs/<id>/input.fasta
  2. MAFFT --addfragments  (add user seqs to existing reference MSA)
  3. EPA-ng placement  (place onto pre-built reference tree)
  4. Genotype assignment  (from nearest reference neighbour annotations)
  5. IQ-TREE2 subtree rebuild  (user seqs + nearest N refs)
  6. GLUE mutation analysis  (reuse existing Snakemake rules via subprocess)
  7. RDP5 recombination  (via existing Snakefile target)
  8. Serialise results  →  jobs/<id>/state.json { state: "done", results: ... }

Each step writes progress updates so the poll callback can show a live bar.

Configuration
-------------
Set these either as environment variables or by editing DEFAULTS below.

  HEP_JOBS_DIR        directory for per-job working files    (default: /tmp/hep_dash_jobs)
  HEP_PROJECT_ROOT    root of the Snakemake project          (default: cwd)
  HEP_REF_MSA_DIR     directory with per-virus reference MSAs  e.g. refs/hbv/ref_msa.fasta
  HEP_REF_TREE_DIR    directory with pre-built reference trees  e.g. refs/hbv/ref_tree.treefile
  HEP_REF_ANNOT_DIR   per-virus reference annotation TSVs (id → genotype)
  HEP_SNAKEMAKE_CORES cores to pass to snakemake --cores     (default: 4)
  HEP_EPA_NG          path to epa-ng binary                  (default: epa-ng)
  HEP_IQTREE          path to iqtree2 binary                 (default: iqtree2)
  HEP_MAFFT           path to mafft binary                   (default: mafft)

Usage in user_sequence_analysis.py
-----------------------------------
  from pipeline_runner import dispatch_user_pipeline, get_job_status
"""

import json
import logging
import os
import re
import shutil
import subprocess
import tempfile
import threading
import uuid
from datetime import datetime
from pathlib import Path

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

def _env(key, default):
    return os.environ.get(key, default)


def _resolve_conda_bin(binary: str, env_name: str = "phylo") -> str:
    """
    Resolve a binary to its full path, checking in order:
      1. Active conda environment  ($CONDA_PREFIX/bin/<binary>)
      2. Named conda environment   (conda info --base)/envs/<env_name>/bin/<binary>)
      3. Common conda base paths   (~/miniconda3, ~/anaconda3, /opt/conda)
      4. Plain PATH lookup         (shutil.which)

    Raises FileNotFoundError with a helpful install message if not found.
    """
    # 1. Active conda env
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidate = Path(conda_prefix) / "bin" / binary
        if candidate.exists():
            return str(candidate)

    # 2. Named conda env via `conda info --base`
    conda_base = None
    try:
        result = subprocess.run(
            ["conda", "info", "--base"],
            capture_output=True, text=True, timeout=10,
        )
        if result.returncode == 0:
            conda_base = result.stdout.strip()
    except Exception:
        pass

    # 3. Common fallback base paths
    if not conda_base:
        for p in [Path.home() / "miniconda3", Path.home() / "anaconda3",
                  Path("/opt/conda"), Path("/usr/local/conda")]:
            if p.exists():
                conda_base = str(p)
                break

    if conda_base:
        candidate = Path(conda_base) / "envs" / env_name / "bin" / binary
        if candidate.exists():
            return str(candidate)

    # 4. PATH
    in_path = shutil.which(binary)
    if in_path:
        return in_path

    raise FileNotFoundError(
        f"'{binary}' not found.\n"
        f"Install it by running:\n"
        f"  conda env create -f envs/phylo.yaml\n"
        f"Or set the environment variable to its full path, e.g.:\n"
        f"  export HEP_{binary.upper().replace('-', '_')}=/path/to/{binary}\n"
        f"Then restart the dashboard."
    )


JOBS_DIR        = Path(_env("HEP_JOBS_DIR",        "/tmp/hep_dash_jobs"))
PROJECT_ROOT    = Path(_env("HEP_PROJECT_ROOT",    os.getcwd())).resolve()
REF_MSA_DIR     = Path(_env("HEP_REF_MSA_DIR",     PROJECT_ROOT / "refs"))
REF_TREE_DIR    = Path(_env("HEP_REF_TREE_DIR",    PROJECT_ROOT / "refs"))
REF_ANNOT_DIR   = Path(_env("HEP_REF_ANNOT_DIR",   PROJECT_ROOT / "refs"))
SNAKEMAKE_CORES = int(_env("HEP_SNAKEMAKE_CORES",  "4"))

# Binaries: use env var override if set, otherwise resolve from conda env
EPA_NG = _env("HEP_EPA_NG",  None) or _resolve_conda_bin("epa-ng",  "phylo")
IQTREE = _env("HEP_IQTREE",  None) or _resolve_conda_bin("iqtree2", "phylo")
MAFFT  = _env("HEP_MAFFT",   None) or _resolve_conda_bin("mafft",   "phylo")

NEAREST_REF_N = 20   # number of nearest reference sequences to include in subtree

# ---------------------------------------------------------------------------
# Job state helpers
# ---------------------------------------------------------------------------

def _job_dir(job_id: str) -> Path:
    d = JOBS_DIR / job_id
    d.mkdir(parents=True, exist_ok=True)
    return d


def _write_state(job_id: str, state: str, progress: int,
                 label: str, results: dict = None, error: str = None):
    payload = {
        "state":    state,
        "progress": progress,
        "label":    label,
        "updated":  datetime.utcnow().isoformat(),
    }
    if results is not None:
        payload["results"] = results
    if error is not None:
        payload["error"] = error

    path = _job_dir(job_id) / "state.json"
    tmp  = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(payload))
    tmp.replace(path)   # atomic on POSIX


def get_job_status(job_id: str) -> dict:
    """
    Read the current state for a job.
    Returns a dict with keys: state, progress, label, results?, error?
    Returns {"state": "not_found"} if the job directory does not exist.
    """
    path = _job_dir(job_id) / "state.json"
    try:
        return json.loads(path.read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return {"state": "not_found", "progress": 0, "label": "Job not found"}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def dispatch_user_pipeline(validated_data: dict, run_recombination: bool = False) -> str:
    """
    Write the input FASTA, spin up a background thread, return the job_id.

    Parameters
    ----------
    validated_data : dict
        The payload stored in useq-validated-store:
        {
          "sequences":      [{"id": str, "seq": str, "length": int}, ...],
          "detected_virus": "HBV" | "HCV" | "HEV" | None,
          "count":          int,
        }
    run_recombination : bool
        Whether to run RDP5 recombination testing as part of this run.

    Returns
    -------
    str  — UUID job identifier
    """
    job_id = str(uuid.uuid4())
    jdir   = _job_dir(job_id)

    # Write input FASTA
    fasta_path = jdir / "input.fasta"
    with fasta_path.open("w") as fh:
        for rec in validated_data["sequences"]:
            fh.write(f">{rec['id']}\n{rec['seq']}\n")

    # Persist metadata including full validated_data
    (jdir / "meta.json").write_text(json.dumps({
        "job_id":        job_id,
        "submitted_at":  datetime.utcnow().isoformat(),
        "detected_virus": validated_data.get("detected_virus"),
        "count":         validated_data.get("count", 0),
        "validated_data": validated_data,
    }))

    _write_state(job_id, "pending", 2, "Job queued…")

    # Launch background thread
    t = threading.Thread(
        target=_run_pipeline,
        args=(job_id, fasta_path, validated_data, run_recombination),
        daemon=True,
    )
    t.start()

    return job_id


# ---------------------------------------------------------------------------
# Pipeline orchestration (runs in a background thread)
# ---------------------------------------------------------------------------

def _run_pipeline(job_id: str, fasta_path: Path, validated_data: dict, run_recombination: bool):
    jdir = _job_dir(job_id)
    try:
        virus = (validated_data.get("detected_virus") or "HBV").lower()

        # ── Step 1: Detect virus if not already known ─────────────────────
        _write_state(job_id, "running", 8, "Detecting virus type…")
        if not validated_data.get("detected_virus"):
            virus = _blast_detect_virus(fasta_path, jdir) or "hbv"
        virus_lc = virus.lower()

        # ── Step 2: MAFFT add-fragment alignment ──────────────────────────
        _write_state(job_id, "running", 18, "Aligning with MAFFT (add-fragment mode)…")
        ref_msa   = REF_MSA_DIR / virus_lc / "ref_msa.fasta"
        aligned   = jdir / "user_aligned.fasta"
        _mafft_add_fragments(fasta_path, ref_msa, aligned, jdir)

        # ── Step 3: EPA-ng placement ───────────────────────────────────────
        _write_state(job_id, "running", 32, "Placing sequences on reference tree (EPA-ng)…")
        ref_tree  = REF_TREE_DIR / virus_lc / "ref_tree.treefile"
        epa_out   = jdir / "epa_out"
        placements = _epa_ng_place(aligned, ref_msa, ref_tree, epa_out, jdir)

        # ── Step 4: Genotype assignment ────────────────────────────────────
        _write_state(job_id, "running", 46, "Assigning genotypes…")
        ref_annot  = REF_ANNOT_DIR / virus_lc / "ref_annotations.tsv"
        genotype_map = _assign_genotypes(placements, ref_annot)

        # ── Step 5: Subtree rebuild with IQ-TREE ──────────────────────────
        _write_state(job_id, "running", 58, "Building phylogenetic subtree (IQ-TREE2)…")
        subtree_seqs = jdir / "subtree_input.fasta"
        newick       = _build_subtree(
            fasta_path, placements, ref_msa, subtree_seqs, jdir
        )

        # ── Step 6: GLUE mutation analysis ────────────────────────────────
        _write_state(job_id, "running", 70, "Running GLUE mutation analysis…")
        mutation_results = _run_glue_mutation(fasta_path, virus_lc, jdir)

        # ── Step 7: RDP5 recombination via Snakemake ──────────────────────
        recomb_results = {}
        if run_recombination:
            _write_state(job_id, "running", 82, "Running recombination analysis (RDP5)…")
            aligned_fasta = jdir / "user_aligned.fasta"
            rdp_input = aligned_fasta if aligned_fasta.exists() else fasta_path
            recomb_results = _run_recombination(rdp_input, virus_lc, jdir)
        else:
            logger.info("Recombination analysis skipped as per job settings.")

        # ── Step 8: Assemble results ───────────────────────────────────────
        _write_state(job_id, "running", 94, "Assembling results…")
        results = _assemble_results(
            validated_data["sequences"],
            virus_lc,
            genotype_map,
            recomb_results,
            mutation_results,
            newick,
            recombination_run=run_recombination,
        )

        _write_state(job_id, "done", 100, "Analysis complete.", results=results)
        logger.info("Job %s completed successfully.", job_id)

    except Exception as exc:
        logger.exception("Job %s failed.", job_id)
        _write_state(job_id, "error", 0, "Pipeline error.", error=str(exc))


# ---------------------------------------------------------------------------
# Step implementations
# ---------------------------------------------------------------------------

def _run(cmd: list[str], cwd: Path, label: str, timeout: int = 600) -> subprocess.CompletedProcess:
    """Run a shell command, raise on non-zero exit."""
    logger.info("[%s] %s", label, " ".join(str(c) for c in cmd))
    result = subprocess.run(
        cmd,
        cwd=str(cwd),
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"{label} failed (exit {result.returncode}):\n"
            f"STDOUT: {result.stdout[-2000:]}\n"
            f"STDERR: {result.stderr[-2000:]}"
        )
    return result


def _blast_detect_virus(fasta_path: Path, jdir: Path) -> str:
    """
    BLASTn the first sequence against per-virus reference DBs.
    Returns lowercase virus name: "hbv", "hcv", or "hev".
    Falls back to "hbv" if BLAST fails or DBs are absent.
    """
    for virus in ["hbv", "hcv", "hev"]:
        db = REF_MSA_DIR / virus / "blast_db" / f"{virus}_ref"
        if not db.with_suffix(".nin").exists() and not db.with_suffix(".nhr").exists():
            continue
        out = jdir / f"blast_{virus}.txt"
        try:
            _run(
                ["blastn", "-query", str(fasta_path), "-db", str(db),
                 "-outfmt", "6 qseqid sseqid pident length",
                 "-max_target_seqs", "1", "-out", str(out)],
                cwd=jdir, label=f"BLAST/{virus}", timeout=120,
            )
            if out.exists() and out.stat().st_size > 0:
                return virus
        except Exception:
            continue
    return "hbv"


def _mafft_add_fragments(user_fasta: Path, ref_msa: Path,
                         output: Path, jdir: Path):
    """
    Add user sequences to the existing reference MSA without re-aligning it.
    Uses MAFFT --addfragments for short/partial sequences, --add for full genomes.
    """
    if not ref_msa.exists():
        raise FileNotFoundError(
            f"Reference MSA not found: {ref_msa}\n"
            "Run the full pipeline first to generate reference alignments, "
            "or point HEP_REF_MSA_DIR at your existing MSA directory."
        )
    _run(
        [MAFFT,
         "--addfragments", str(user_fasta),
         "--thread", "4",
         "--auto",
         "--keeplength",        # don't change existing alignment length
         str(ref_msa)],
        cwd=jdir, label="MAFFT",
        timeout=300,
    )
    # MAFFT writes to stdout; capture and redirect
    result = subprocess.run(
        [MAFFT,
         "--addfragments", str(user_fasta),
         "--thread", "4",
         "--auto",
         "--keeplength",
         str(ref_msa)],
        cwd=str(jdir),
        capture_output=True,
        text=True,
        timeout=300,
    )
    if result.returncode != 0:
        raise RuntimeError(f"MAFFT failed:\n{result.stderr[-2000:]}")
    output.write_text(result.stdout)


def _epa_ng_place(aligned_query: Path, ref_msa: Path,
                  ref_tree: Path, epa_out: Path, jdir: Path) -> list[dict]:
    """
    Run EPA-ng to place user sequences onto the fixed reference tree.
    Returns a list of placement dicts:
        [{"query_id": str, "ref_id": str, "like_weight": float, "distal_length": float}, ...]
    """
    if not ref_tree.exists():
        raise FileNotFoundError(
            f"Reference tree not found: {ref_tree}\n"
            "Build it with: iqtree2 -s <ref_msa> -m GTR+G -bb 1000 --prefix refs/<virus>/ref_tree"
        )

    epa_out.mkdir(exist_ok=True)

    # EPA-ng requires an evaluated model file from RAxML, not just a model string.
    # Expected at: refs/{virus}/RAxML_info.{virus}_ref
    # Generate once per virus with:
    #   raxmlHPC -f e -s <ref_msa> -t <ref_tree> -n {virus}_ref -m GTRGAMMAX -w refs/{virus}/
    virus_lc  = ref_tree.parent.name          # infer virus from refs/<virus>/
    info_file = ref_tree.parent / f"RAxML_info.{virus_lc}_ref"

    cmd = [EPA_NG,
           "--tree",    str(ref_tree),
           "--ref-msa", str(ref_msa),
           "--query",   str(aligned_query),
           "--outdir",  str(epa_out),
           "--redo"]

    if info_file.exists():
        cmd += ["--model", str(info_file)]
        logger.info("Using RAxML model file: %s", info_file)
    else:
        # Fallback: generic GTR+G — EPA-ng will warn but still run
        cmd += ["--model", "GTR+G"]
        logger.warning(
            "RAxML info file not found at %s. "
            "EPA-ng will use generic GTR+G without optimised parameters.\n"
            "Generate it once per virus with:\n"
            "  conda run -n phylo raxmlHPC -f e -s %s -t %s -n %s_ref -m GTRGAMMAX -w %s",
            info_file, ref_msa, ref_tree, virus_lc, ref_tree.parent,
        )

    _run(cmd, cwd=str(jdir), label="EPA-ng", timeout=600)

    # Parse the output jplace file
    jplace_files = list(epa_out.glob("*.jplace"))
    if not jplace_files:
        raise FileNotFoundError("EPA-ng did not produce a .jplace file")

    return _parse_jplace(jplace_files[0])


def _find_descendants(tree_str: str, edge_num: int, valid_labels: set = None) -> list[str]:
    """Find all valid leaf descendant accessions under the internal branch matching edge_num inside the Newick tree string."""
    import re
    target = f"{{{edge_num}}}"
    idx = tree_str.find(target)
    if idx == -1:
        return []
    
    # Find the closing parenthesis before this edge support/distance
    p_close = -1
    for i in range(idx - 1, -1, -1):
        if tree_str[i] == ")":
            p_close = i
            break
        elif tree_str[i] == ",":
            break
            
    if p_close == -1:
        return []
        
    # Walk backward to find the matching opening parenthesis
    depth = 0
    p_open = -1
    for i in range(p_close, -1, -1):
        if tree_str[i] == ")":
            depth += 1
        elif tree_str[i] == "(":
            depth -= 1
            if depth == 0:
                p_open = i
                break
                
    if p_open == -1:
        return []
        
    substring = tree_str[p_open:p_close]
    leaves = re.findall(r'([\'"]?[A-Za-z0-9_.-]+[\'"]?):', substring)
    cleaned = [l.replace("'", "").replace('"', "").strip() for l in leaves if l.strip()]
    
    if valid_labels:
        cleaned = [l for l in cleaned if l in valid_labels]
    return cleaned


def _parse_jplace(jplace_path: Path) -> list[dict]:
    """Parse EPA-ng .jplace file with robust label mapping."""
    import re
    
    data = json.loads(jplace_path.read_text())
    fields = data.get("fields", [])
    
    # Find field indices
    try:
        edge_idx = fields.index("edge_num")
        like_idx = fields.index("like_weight_ratio")
        distal_idx = fields.index("distal_length")
    except ValueError:
        edge_idx, like_idx, distal_idx = 0, 3, 4
    
    # Extract mapping from edge numbers to original labels
    tree_str = data.get("tree", "")
    edge_to_ref = {}
    
    # Method 1: Look for pattern "label:distance{number}"
    # e.g., ref_AB056513:0.1029124771{0}
    pattern = r'([A-Za-z_][A-Za-z0-9_.]*):[\d.eE-]+\{(\d+)\}'
    for match in re.finditer(pattern, tree_str):
        label = match.group(1)
        edge_num = float(match.group(2))
        edge_to_ref[edge_num] = label
    
    # Method 3: Fallback - load reference MSA to get order
    if not edge_to_ref:
        # Try to load reference MSA
        ref_msa_path = jplace_path.parent.parent / "ref_msa.fasta"
        if ref_msa_path.exists():
            labels = []
            with ref_msa_path.open() as f:
                for line in f:
                    if line.startswith(">"):
                        labels.append(line[1:].strip())
            
            # Extract numeric edges in order
            numeric_edges = re.findall(r'\{([0-9.]+)\}', tree_str)
            for i, edge_str in enumerate(numeric_edges):
                if i < len(labels):
                    edge_to_ref[float(edge_str)] = labels[i]
    
    logger.info(f"Mapped {len(edge_to_ref)} edges to sequence labels")
    
    # Parse placements
    placements = []
    for pquery in data.get("placements", []):
        # Get query name
        if "n" in pquery and pquery["n"]:
            query_id = pquery["n"][0]
        else:
            query_id = "unknown"
        
        # Get best placement
        best = max(pquery["p"], key=lambda p: p[like_idx])
        edge_num = best[edge_idx]
        
        # Convert to float for lookup
        try:
            edge_key = float(edge_num) if isinstance(edge_num, str) else edge_num
        except (ValueError, TypeError):
            edge_key = edge_num
        
        # Get reference label
        ref_label = edge_to_ref.get(edge_key)
        original_edge = None
        if not ref_label:
            try:
                edge_int = int(edge_key)
            except (ValueError, TypeError):
                edge_int = None
            
            if edge_int is not None:
                descendant_refs = _find_descendants(tree_str, edge_int, set(edge_to_ref.values()))
                if descendant_refs:
                    ref_label = descendant_refs[0]
                    original_edge = f"edge_{edge_num}"
            
            if not ref_label:
                ref_label = f"edge_{edge_num}"
        
        placements.append({
            "query_id": query_id,
            "edge_num": edge_num,
            "ref_label": ref_label,
            "original_edge": original_edge,
            "like_weight": best[like_idx],
            "distal_length": best[distal_idx],
        })
    
    return placements


def _assign_genotypes(placements: list[dict], annot_path: Path) -> dict[str, dict]:
    """
    Map query sequences to genotypes using nearest reference annotations.

    annot_path TSV format (tab-separated, with header):
        ref_id  genotype  [subtype]

    Returns { query_id: {"genotype": str, "nearest_ref": str, "epa_score": float} }
    """
    # Load annotations
    annot = {}
    if annot_path.exists():
        for line in annot_path.read_text().splitlines()[1:]:   # skip header
            parts = line.split("\t")
            if len(parts) >= 2:
                annot[parts[0].strip()] = parts[1].strip()
    else:
        logger.warning("Annotation file not found: %s — genotype will be 'Unknown'", annot_path)

    result = {}
    for p in placements:
        qid      = p["query_id"]
        ref_id   = p.get("ref_label", str(p.get("edge_num", "?")))
        original_edge = p.get("original_edge")
        
        genotype = annot.get(ref_id, "Unknown")
        # If exact match fails, try prefix match (accession without version suffix)
        if genotype == "Unknown":
            ref_prefix = ref_id.split(".")[0]
            genotype = next(
                (v for k, v in annot.items() if k.split(".")[0] == ref_prefix),
                "Unknown"
            )
            
        display_ref = f"{original_edge} (resolved to {ref_id})" if original_edge else ref_id
        
        result[qid] = {
            "genotype":    genotype,
            "nearest_ref": display_ref,
            "epa_score":   round(p["like_weight"], 4),
        }
    return result


def _build_subtree(user_fasta: Path, placements: list[dict],
                   ref_msa: Path, subtree_seqs: Path, jdir: Path) -> str:
    """
    Concatenate user sequences + nearest NEAREST_REF_N reference sequences,
    realign with MAFFT, and build a small ML tree with IQ-TREE2.
    Returns the Newick string (or "" on failure).
    """
    # Collect nearest reference labels (real accessions from jplace tree)
    ref_ids = {p.get("ref_label", str(p["edge_num"])) for p in placements}
    logger.info("Subtree: extracting %d reference IDs from ref MSA", len(ref_ids))

    # Extract those references from the ref MSA
    ref_seqs = _extract_sequences(ref_msa, ref_ids, max_seqs=NEAREST_REF_N)
    ref_count = ref_seqs.count(">")
    logger.info("Subtree: found %d matching reference sequences", ref_count)

    # If we still have too few references, fall back to taking the first N from the MSA
    if ref_count < 2:
        logger.warning(
            "Too few reference sequences found by ID (%d). "
            "Falling back to first %d sequences from ref MSA.",
            ref_count, NEAREST_REF_N,
        )
        ref_seqs = _extract_first_n(ref_msa, NEAREST_REF_N)
        ref_count = ref_seqs.count(">")

    # Write combined FASTA
    with subtree_seqs.open("w") as fh:
        fh.write(user_fasta.read_text())
        fh.write(ref_seqs)

    total_seqs = subtree_seqs.read_text().count(">")
    if total_seqs < 3:
        logger.warning("Only %d sequences for subtree — skipping IQ-TREE2 (need ≥3)", total_seqs)
        return ""

    # Realign with MAFFT
    aligned = jdir / "subtree_aligned.fasta"
    result = subprocess.run(
        [MAFFT, "--thread", "4", "--auto", str(subtree_seqs)],
        capture_output=True, text=True, timeout=300, cwd=str(jdir),
    )
    if result.returncode == 0 and result.stdout.count(">") >= 3:
        aligned.write_text(result.stdout)
    else:
        logger.warning("MAFFT subtree alignment failed; using unaligned input")
        shutil.copy(subtree_seqs, aligned)

    # Final sequence count check before IQ-TREE2
    aligned_count = aligned.read_text().count(">")
    if aligned_count < 3:
        logger.warning("Aligned subtree has %d sequences — skipping IQ-TREE2", aligned_count)
        return ""

    # IQ-TREE2
    prefix = jdir / "subtree"
    try:
        _run(
            [IQTREE, "-s", str(aligned), "-m", "GTR+G",
             "--prefix", str(prefix), "-nt", "AUTO", "--redo", "-fast"],
            cwd=str(jdir), label="IQ-TREE2 subtree", timeout=900,
        )
        treefile = prefix.with_suffix(".treefile")
        if treefile.exists():
            return treefile.read_text().strip()
    except Exception as e:
        logger.warning("IQ-TREE2 subtree failed: %s", e)

    return ""


def _extract_sequences(fasta_path: Path, ids: set, max_seqs: int = 20) -> str:
    """
    Extract up to max_seqs sequences from a FASTA file matching any of the given IDs.
    Simple grep-style parser — no Biopython dependency.
    """
    out_lines = []
    current_id = None
    current_keep = False
    count = 0
    for line in fasta_path.read_text().splitlines():
        if line.startswith(">"):
            current_id = line[1:].split()[0]
            current_keep = (current_id in ids) and (count < max_seqs)
            if current_keep:
                count += 1
        if current_keep:
            out_lines.append(line)
    return "\n".join(out_lines) + ("\n" if out_lines else "")


def _extract_first_n(fasta_path: Path, n: int) -> str:
    """Extract the first n sequences from a FASTA file (fallback when ID matching fails)."""
    out_lines = []
    count = 0
    active = False
    for line in fasta_path.read_text().splitlines():
        if line.startswith(">"):
            count += 1
            if count > n:
                break
            active = True
        if active:
            out_lines.append(line)
    return "\n".join(out_lines) + ("\n" if out_lines else "")


def _run_glue_mutation(user_fasta: Path, virus: str, jdir: Path) -> dict:
    """
    Run the existing GLUE mutation analysis for user-submitted sequences
    by invoking the relevant Snakemake rules directly.

    Steps:
      1. Copy user FASTA into the Snakemake results tree under a special
         job-scoped directory so it doesn't collide with production data.
      2. Call snakemake with --config to override the input path.
      3. Parse the resulting resistance TSV.

    Returns { seq_id: { "mutations": [...], "drug_resistance": [...] } }
    """
    job_results_dir = PROJECT_ROOT / "results" / f"user_{jdir.name}" / virus

    try:
        # Write user FASTA in the expected Snakemake location
        fasta_dest = job_results_dir / f"{virus}_sequences.fasta"
        fasta_dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(user_fasta, fasta_dest)

        # Call GLUE directly via its groovy script rather than routing
        # through Snakemake (avoids DAG resolution issues with dynamic targets).
        # Expects: scripts/run_glue.sh <fasta> <virus> <outdir>
        # which internally calls: gluetools.sh -p batchfile scripts/glue/<virus>_resistance.glue
        glue_wrapper = PROJECT_ROOT / "scripts" / "run_glue.sh"
        if glue_wrapper.exists():
            _run(
                ["bash", str(glue_wrapper),
                 str(fasta_dest), virus, str(job_results_dir)],
                cwd=str(PROJECT_ROOT),
                label="GLUE",
                timeout=3600,
            )
        else:
            # Fallback: invoke Snakemake targeting just the per-sequence GLUE rule.
            # Use --allowed-rules to avoid DAG ambiguity with dynamic checkpoint targets.
            _run(
                [
                    "snakemake",
                    "--cores", str(SNAKEMAKE_CORES),
                    "--use-conda",
                    "--snakefile", str(PROJECT_ROOT / "Snakefile"),
                    "--config", f"user_job_id={jdir.name}",
                    "--allowed-rules", "glue_analysis", "parse_glue_xml", "merge_resistance_tables",
                    "--", str(job_results_dir / "final_resistance.tsv"),
                ],
                cwd=str(PROJECT_ROOT),
                label="Snakemake/GLUE",
                timeout=3600,
            )

        # Parse all generated XML files
        all_results = {}
        for xml_file in job_results_dir.glob("*.xml"):
            if xml_file.name != "glue_cmd.glue":
                res = _parse_user_glue_xml(xml_file, virus)
                all_results.update(res)
        
        # If no XML was parsed, fallback to TSV if present
        if not all_results:
            resistance_tsv = job_results_dir / "final_resistance.tsv"
            if resistance_tsv.exists():
                return _parse_resistance_tsv(resistance_tsv)
        return all_results

    except Exception as e:
        logger.warning("GLUE mutation step failed: %s", e)
        return {}
    finally:
        # Clean up the job-scoped Snakemake results dir
        try:
            shutil.rmtree(job_results_dir.parent, ignore_errors=True)
        except Exception:
            pass


def _parse_resistance_tsv(tsv_path: Path) -> dict:
    """
    Parse the merged resistance TSV produced by merge_resistance_tables.
    Returns { sample_id: { mutations: list, drugs: list } }
    """
    if not tsv_path.exists():
        return {}
    results = {}
    try:
        import csv
        with tsv_path.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                # Identify sample column
                sample = (row.get("sample") or row.get("Sample") or
                          row.get("sequence_id") or row.get("accession") or "?")
                if sample not in results:
                    results[sample] = {"mutations": [], "drugs": []}
                mutation = row.get("mutation") or row.get("Mutation")
                drug     = row.get("drug") or row.get("Drug") or row.get("gene")
                if mutation:
                    results[sample]["mutations"].append(mutation)
                if drug and drug not in results[sample]["drugs"]:
                    results[sample]["drugs"].append(drug)
    except Exception as e:
        logger.warning("Failed to parse resistance TSV: %s", e)
    return results


def _parse_user_glue_xml(xml_path: Path, virus: str) -> dict:
    """
    Parse the raw XML report generated by GLUE (result.xml) and extract
    resistance mutations and drugs for each sequence.
    Returns { sample_id: { "mutations": list, "drugs": list } }
    """
    import xml.etree.ElementTree as ET
    if not xml_path.exists():
        return {}

    results = {}
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        virus_lc = virus.lower()

        if virus_lc == "hbv":
            reports = root.findall(".//hbvReport")
            for rep in reports:
                seq_res = rep.find("sequenceResult")
                if seq_res is None:
                    continue
                sample_id = seq_res.findtext("id")
                if not sample_id:
                    continue

                if sample_id not in results:
                    results[sample_id] = {"mutations": [], "drugs": []}

                # Antiviral resistance
                for ar in seq_res.findall(".//antiviralResistance"):
                    mut = ar.findtext("description")
                    drug = ar.findtext("drug/name")
                    if mut and mut not in ["NA", "none", "None", "unknown", "?", "NA_mutation"]:
                        if mut not in results[sample_id]["mutations"]:
                            results[sample_id]["mutations"].append(mut)
                    if drug and drug not in ["NA", "none", "None", "unknown", "?"]:
                        if drug not in results[sample_id]["drugs"]:
                            results[sample_id]["drugs"].append(drug)

                # Vaccine escape
                for ve in seq_res.findall(".//vaccineEscape"):
                    mut = ve.findtext("description")
                    if mut and mut not in ["NA", "none", "None", "unknown", "?", "NA_mutation"]:
                        if mut not in results[sample_id]["mutations"]:
                            results[sample_id]["mutations"].append(mut)
                        if "Vaccine" not in results[sample_id]["drugs"]:
                            results[sample_id]["drugs"].append("Vaccine")

        elif virus_lc == "hcv":
            reports = []
            if root.tag == "phdrReport":
                reports = [root]
            else:
                reports = root.findall(".//phdrReport")

            for rep in reports:
                seq_res = rep.find("sequenceResult")
                if seq_res is None:
                    continue
                sample_id = seq_res.findtext("id")
                if not sample_id:
                    continue

                if sample_id not in results:
                    results[sample_id] = {"mutations": [], "drugs": []}

                # 1. rasScanResults
                for ras in seq_res.findall(".//rasScanResults"):
                    present = ras.findtext("present")
                    if present != "true":
                        continue

                    var_name = ras.findtext("variationName")
                    gene = ras.findtext("featureName")
                    mut = var_name.split(":")[-1] if var_name and ":" in var_name else var_name
                    
                    if mut and mut not in ["NA", "none", "None", "unknown", "?", "NA_mutation"]:
                        mut_label = f"{gene}_{mut}" if gene and not mut.startswith(gene) else mut
                        if mut_label not in results[sample_id]["mutations"]:
                            results[sample_id]["mutations"].append(mut_label)

                    for ard in ras.findall(".//alignmentRasDrug"):
                        drug = ard.findtext("drug")
                        res_cat = ard.findtext("resistanceCategory")
                        if drug and res_cat not in ["insignificant", "-"]:
                            if drug not in results[sample_id]["drugs"]:
                                results[sample_id]["drugs"].append(drug)

                # 2. drugScores
                for ds in seq_res.findall(".//drugScores"):
                    for da in ds.findall(".//drugAssessments"):
                        drug = da.findtext(".//drug/id")
                        score = da.findtext("drugScore")
                        if drug and score not in ["no_significant_resistance_detected", "insignificant"]:
                            if drug not in results[sample_id]["drugs"]:
                                results[sample_id]["drugs"].append(drug)

                            for ras_cat in ["rasScores_category_I", "rasScores_category_II", "rasScores_category_III"]:
                                for ras in da.findall(f".//{ras_cat}"):
                                    mut = ras.findtext("structure")
                                    gene = ras.findtext("gene")
                                    if mut and mut not in ["NA", "none", "None", "unknown", "?", "NA_mutation"]:
                                        mut_label = f"{gene}_{mut}" if gene and not mut.startswith(gene) else mut
                                        if mut_label not in results[sample_id]["mutations"]:
                                            results[sample_id]["mutations"].append(mut_label)

                # 3. substitutionsOfInterest
                for sub in seq_res.findall(".//substitutionsOfInterest"):
                    mut = sub.findtext("displayStructure")
                    gene = sub.findtext("virusProtein")
                    if mut and mut not in ["NA", "none", "None", "unknown", "?", "NA_mutation"]:
                        mut_label = f"{gene}_{mut}" if gene and not mut.startswith(gene) else mut
                        if mut_label not in results[sample_id]["mutations"]:
                            results[sample_id]["mutations"].append(mut_label)

        elif virus_lc == "hev":
            reports = []
            if root.tag == "hevReport":
                reports = [root]
            else:
                reports = root.findall(".//hevReport")

            for rep in reports:
                seq_res = rep.find("sequenceResult")
                if seq_res is None:
                    continue
                sample_id = seq_res.findtext("id")
                if not sample_id:
                    continue
                if sample_id not in results:
                    results[sample_id] = {"mutations": [], "drugs": []}

    except Exception as e:
        logger.warning("Failed to parse GLUE XML: %s", e)

    return results


def _run_recombination(user_fasta: Path, virus: str, jdir: Path) -> dict:
    """
    Run 3seq on user sequence(s) and perform validation steps (OpenRDP and regional trees).
    Classifies results into "high_confidence", "needs_review", or "none" (discarded).
    """
    wrapper = PROJECT_ROOT / "scripts" / "run_rdp5.py"
    rdp_out = jdir / "rdp5_out"
    rdp_out.mkdir(exist_ok=True)

    if not wrapper.exists():
        logger.warning("3seq wrapper not found at %s — recombination step skipped.", wrapper)
        return {}

    try:
        # Run 3seq
        _run(
            ["python", str(wrapper),
             "--fasta",  str(user_fasta),
             "--virus",  virus,
             "--outdir", str(rdp_out)],
            cwd=str(PROJECT_ROOT),
            label="3seq",
            timeout=1800,
        )
        
        tsv_path = rdp_out / "recombinants.tsv"
        if not tsv_path.exists():
            return {}
            
        # Import validation script functions dynamically
        import sys
        sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
        import validate_recombination
        
        aligned_path = jdir / "user_aligned.fasta"
        if not aligned_path.exists():
            aligned_path = user_fasta
            
        results = {}
        import csv
        with tsv_path.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sid = row.get("sequence_id", "?")
                is_r = row.get("is_recombinant", "false").lower() in ("true", "1", "yes")
                p1 = row.get("parent_1")
                p2 = row.get("parent_2")
                bp_s = row.get("breakpoint_start")
                bp_e = row.get("breakpoint_end")
                
                results[sid] = {
                    "is_recombinant": False,
                    "validation_status": "none",
                    "breakpoints": []
                }
                
                if is_r:
                    # Parse breakpoints
                    start, end = 0, 0
                    try:
                        start = int(bp_s)
                        end = int(bp_e)
                    except (ValueError, TypeError):
                        pass
                        
                    # Run validation checks
                    passed_len = False
                    passed_dist = False
                    passed_openrdp = False
                    passed_tree = False
                    
                    if p1 and p2 and start > 0 and end > 0:
                        # 1. Length check: is it in the longRec file?
                        longrec_path = next(rdp_out.glob("*.3s.longRec"), None)
                        long_rec_set = set()
                        if longrec_path and longrec_path.exists():
                            long_rec_set = validate_recombination.load_long_rec_candidates(longrec_path)
                        if sid in long_rec_set:
                            passed_len = True
                        
                        # 2. Distance check: parent switching
                        extracted_seqs = validate_recombination.get_alignment_sequences(aligned_path, [sid, p1, p2])
                        if sid in extracted_seqs and p1 in extracted_seqs and p2 in extracted_seqs:
                            cand_seq = extracted_seqs[sid]
                            p1_seq = extracted_seqs[p1]
                            p2_seq = extracted_seqs[p2]
                            
                            ident_p1_r = validate_recombination.calculate_identity(cand_seq, p1_seq, start, end, invert=False)
                            ident_p2_r = validate_recombination.calculate_identity(cand_seq, p2_seq, start, end, invert=False)
                            ident_p1_nr = validate_recombination.calculate_identity(cand_seq, p1_seq, start, end, invert=True)
                            ident_p2_nr = validate_recombination.calculate_identity(cand_seq, p2_seq, start, end, invert=True)
                            
                            if (ident_p2_r > ident_p1_r) and (ident_p1_nr > ident_p2_nr):
                                passed_dist = True
                        
                        # 3. OpenRDP cross-method validation
                        passed_openrdp = validate_recombination.run_openrdp_validation(
                            aligned_path, sid, p1, p2, rdp_out
                        )
                        
                        # 4. Regional ML tree validation
                        passed_tree = validate_recombination.run_regional_tree_validation(
                            aligned_path, sid, p1, p2, start, end, rdp_out
                        )

                    # Now categorize according to the rules:
                    # High confidence = 3SEQ + OpenRDP/RDP5 support + regional tree support
                    # (Note: also must pass basic validation: length and distance checks)
                    # Candidate = 3SEQ only (passed length/distance basic checks, but fails one or both of OpenRDP/tree support)
                    # Discard/hidden = fails validation (fails length check or fails distance check)
                    
                    if passed_len and passed_dist:
                        if passed_openrdp and passed_tree:
                            val_status = "high_confidence"
                        else:
                            val_status = "needs_review"
                    else:
                        # Fails basic validation -> Discarded
                        val_status = "none"
                        is_r = False
                            
                    if is_r:
                        results[sid]["is_recombinant"] = True
                        results[sid]["validation_status"] = val_status
                        if start > 0 and end > 0:
                            results[sid]["breakpoints"].append([start, end])
                            
        return results
    except Exception as e:
        logger.warning("3seq/validation step failed: %s", e)
        return {}


def _parse_recombinants_tsv(tsv_path: Path) -> dict:
    """
    Parse the recombinants TSV produced by the RDP5 wrapper.
    Expected columns: sequence_id, is_recombinant, breakpoint_start, breakpoint_end
    Returns { seq_id: { is_recombinant: bool, breakpoints: [[s,e], ...] } }
    """
    if not tsv_path.exists():
        return {}
    results = {}
    try:
        import csv
        with tsv_path.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sid = row.get("sequence_id", "?")
                is_r = row.get("is_recombinant", "false").lower() in ("true", "1", "yes")
                bp_s = row.get("breakpoint_start")
                bp_e = row.get("breakpoint_end")
                if sid not in results:
                    results[sid] = {"is_recombinant": False, "breakpoints": []}
                if is_r:
                    results[sid]["is_recombinant"] = True
                if bp_s and bp_e and str(bp_s) != "0" and str(bp_e) != "0":
                    try:
                        results[sid]["breakpoints"].append([int(bp_s), int(bp_e)])
                    except ValueError:
                        pass
    except Exception as e:
        logger.warning("Failed to parse recombinants TSV: %s", e)
    return results


# ---------------------------------------------------------------------------
# Result assembly
# ---------------------------------------------------------------------------

# Approximate ORF coordinates (reference genome positions).
# These are used to annotate the sequence map; adjust to match your references.
ORF_COORDS = {
    "hbv": [
        (0,    832,  "PreS1/S2"),
        (155,  832,  "PreS2"),
        (388,  832,  "S"),
        (1374, 1838, "P(RT)"),
        (1814, 2454, "P(RNaseH)"),
        (1901, 2452, "X"),
        (2307, 3182, "Core"),
    ],
    "hcv": [
        (341,  914,  "Core"),
        (914,  1490, "E1"),
        (1490, 2580, "E2"),
        (2580, 2769, "p7"),
        (2769, 3419, "NS2"),
        (3419, 5312, "NS3"),
        (5312, 5474, "NS4A"),
        (5474, 6257, "NS4B"),
        (6257, 7601, "NS5A"),
        (7601, 9374, "NS5B"),
    ],
    "hev": [
        (26,   5137, "ORF1"),
        (5138, 7127, "ORF2"),
        (4934, 5477, "ORF3"),
    ],
}


def _assemble_results(sequences, virus_lc, genotype_map, recomb_results,
                      mutation_results, newick, recombination_run: bool = False) -> dict:
    """
    Combine all step outputs into the results dict consumed by render_results().
    """
    seq_list   = []
    seq_map    = {}
    orfs       = ORF_COORDS.get(virus_lc, [])

    for rec in sequences:
        sid     = rec["id"]
        geno    = genotype_map.get(sid, {})
        recomb  = recomb_results.get(sid, {"is_recombinant": False, "validation_status": "none", "breakpoints": []})
        muts    = mutation_results.get(sid, {"mutations": [], "drugs": []})

        seq_list.append({
            "id":             sid,
            "virus":          virus_lc.upper(),
            "genotype":       geno.get("genotype", "Unknown"),
            "is_recombinant": recomb["is_recombinant"],
            "validation_status": recomb.get("validation_status", "none"),
            "breakpoints":    recomb["breakpoints"],
            "nearest_ref":    geno.get("nearest_ref", "—"),
            "epa_score":      geno.get("epa_score", 0.0),
            "mutations":      muts.get("mutations", []),
            "drugs":          muts.get("drugs", []),
        })

        # Build per-sequence map annotations
        mutation_annotations = []
        for mut in muts.get("mutations", []):
            # Try to extract a genomic position from the mutation label
            match = re.search(r"(\d+)", mut)
            pos   = int(match.group(1)) if match else 0
            
            offset = 0
            mut_lower = mut.lower()
            if virus_lc == "hbv":
                if mut_lower.startswith("rt"): offset = 1374
                elif mut_lower.startswith("s"): offset = 388
                elif mut_lower.startswith("c"): offset = 2307
            elif virus_lc == "hcv":
                if "ns3" in mut_lower: offset = 3419
                elif "ns5a" in mut_lower: offset = 6257
                elif "ns5b" in mut_lower: offset = 7601
                
            nt_pos = offset + (pos * 3) if pos > 0 else 0
            mutation_annotations.append([nt_pos, mut, "#e00603"])

        seq_map[sid] = {
            "orfs":        orfs,
            "mutations":   mutation_annotations,
            "breakpoints": recomb["breakpoints"],
        }

    return {
        "sequences":    seq_list,
        "newick":       newick,
        "sequence_map": seq_map,
        "recombination_run": recombination_run,
    }


def _parse_fasta_simple(text: str) -> list[dict]:
    seqs = []
    curr_id = None
    curr_seq = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if curr_id:
                seqs.append({"id": curr_id, "seq": "".join(curr_seq), "length": len("".join(curr_seq))})
            curr_id = line[1:].split()[0]
            curr_seq = []
        else:
            curr_seq.append(line)
    if curr_id:
        seqs.append({"id": curr_id, "seq": "".join(curr_seq), "length": len("".join(curr_seq))})
    return seqs


def dispatch_recombination_only(job_id: str) -> bool:
    """
    Submit a background thread to run recombination analysis on an existing, finished basic job.
    """
    jdir = _job_dir(job_id)
    meta_path = jdir / "meta.json"
    state_path = jdir / "state.json"
    
    if not meta_path.exists() or not state_path.exists():
        logger.error("Cannot run recombination: job directory or state not found for %s", job_id)
        return False
        
    try:
        meta = json.loads(meta_path.read_text())
        state = json.loads(state_path.read_text())
    except Exception as e:
        logger.error("Failed to parse metadata/state for job %s: %s", job_id, e)
        return False
        
    validated_data = meta.get("validated_data")
    if not validated_data:
        # Reconstruct from input.fasta
        fasta_path = jdir / "input.fasta"
        if fasta_path.exists():
            seqs = _parse_fasta_simple(fasta_path.read_text())
            validated_data = {
                "sequences": seqs,
                "detected_virus": meta.get("detected_virus"),
                "count": len(seqs)
            }
        else:
            logger.error("Cannot reconstruct validated data: input.fasta missing for job %s", job_id)
            return False
            
    _write_state(job_id, "running", 80, "Queuing recombination analysis…")
    
    # Launch background thread for recombination only
    t = threading.Thread(
        target=_run_recombination_only_pipeline,
        args=(job_id, validated_data, state.get("results", {})),
        daemon=True,
    )
    t.start()
    return True


def _run_recombination_only_pipeline(job_id: str, validated_data: dict, old_results: dict):
    jdir = _job_dir(job_id)
    try:
        virus = (validated_data.get("detected_virus") or "HBV").lower()
        virus_lc = virus.lower()
        
        _write_state(job_id, "running", 85, "Running recombination analysis (RDP5)…")
        aligned_fasta = jdir / "user_aligned.fasta"
        fasta_path = jdir / "input.fasta"
        rdp_input = aligned_fasta if aligned_fasta.exists() else fasta_path
        
        recomb_results = _run_recombination(rdp_input, virus_lc, jdir)
        
        _write_state(job_id, "running", 95, "Assembling updated results…")
        
        sequences = validated_data.get("sequences", [])
        genotype_map = {}
        mutation_results = {}
        newick = old_results.get("newick", "")
        
        for s in old_results.get("sequences", []):
            sid = s["id"]
            genotype_map[sid] = {
                "genotype": s.get("genotype", "Unknown"),
                "nearest_ref": s.get("nearest_ref", "—"),
                "epa_score": s.get("epa_score", 0.0),
            }
            mutation_results[sid] = {
                "mutations": s.get("mutations", []),
                "drugs": s.get("drugs", []),
            }
            
        results = _assemble_results(
            sequences,
            virus_lc,
            genotype_map,
            recomb_results,
            mutation_results,
            newick,
            recombination_run=True
        )
        
        _write_state(job_id, "done", 100, "Analysis complete.", results=results)
        logger.info("Job %s recombination update completed successfully.", job_id)
    except Exception as exc:
        logger.exception("Job %s recombination update failed.", job_id)
        _write_state(job_id, "error", 0, "Recombination pipeline error.", error=str(exc))
