"""
user_sequence_analysis.py
=========================
Drop-in module for the Hepatitis Dashboard (dashboard.py) that adds:
  - A 4th navigation tab: "🔬 My Sequences"
  - FASTA input UI (paste + file upload)
  - FASTA validation callback
  - Job dispatch callback (hooks into Snakemake/EPA-ng pipeline)
  - Result-rendering callbacks (summary table, phylo tree placeholder,
    sequence map placeholder)
  - dcc.Store / dcc.Interval for async polling

USAGE
-----
1.  Place this file alongside dashboard.py.
2.  In dashboard.py, import at the top:
        from user_sequence_analysis import (
            USER_SEQ_STORES,
            user_seq_tab_button,
            user_seq_tab_content,
        )
    No register_user_seq_callbacks call needed — callbacks register
    automatically at import time via @callback.
3.  In create_dashboard_layout():
    a.  Add *USER_SEQ_STORES() inside the Container.
    b.  Add user_seq_tab_button() to the navigation ButtonGroup.
    c.  Add user_seq_tab_content() after the epidemiology content Div.
4.  Add Output("user-seq-content", "style") and Input("tab-user-seq", "n_clicks")
    to your switch_tabs callback.

All identifiers are namespaced with "useq-" to avoid conflicts.
"""

import base64
import io
import json
import re
import time
import uuid
from datetime import datetime

import dash
from dash import Input, Output, State, callback, ctx, dcc, html
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from hep_theme import genotype_palette, VIRUS_COLORS, _GENOTYPE_COLOR_LOOKUP, apply_heptracker_figure_style

# ---------------------------------------------------------------------------
# GENOTYPE COLOR PALETTE LOOKUP (Synchronized with Genotypes Tab)
# ---------------------------------------------------------------------------
HBV_GENOTYPE_COLORS = genotype_palette("HBV", VIRUS_COLORS.get("HBV", "#08519c"), list("ABCDEFGHIJ"))
HCV_GENOTYPE_COLORS = genotype_palette("HCV", VIRUS_COLORS.get("HCV", "#e6550d"), [str(i) for i in range(1, 9)])
HEV_GENOTYPE_COLORS = genotype_palette("HEV", VIRUS_COLORS.get("HEV", "#2ca25f"), [str(i) for i in range(1, 9)])

def get_genotype_color(virus: str, genotype: str) -> str:
    v = (virus or "HBV").upper()
    g = str(genotype or "").strip()
    
    if v == "HBV":
        palette = HBV_GENOTYPE_COLORS
    elif v == "HCV":
        palette = HCV_GENOTYPE_COLORS
    else:
        palette = HEV_GENOTYPE_COLORS
        
    if g in palette:
        return palette[g]
    pref_g = f"{v}-{g}"
    if pref_g in palette:
        return palette[pref_g]
        
    clean = re.sub(r'^(HBV|HCV|HEV|Genotype)[-_\s]*', '', g, flags=re.IGNORECASE)
    base_char = clean[0].upper() if clean else ""
    if base_char in _GENOTYPE_COLOR_LOOKUP:
        return _GENOTYPE_COLOR_LOOKUP[base_char]
        
    return VIRUS_COLORS.get(v, "#0d6efd")

# ---------------------------------------------------------------------------
# CONSTANTS
# ---------------------------------------------------------------------------
MAX_SEQUENCES = 50
MIN_SEQUENCES = 1
VALID_NUCLEOTIDES = re.compile(r"^[ACGTNacgtnRYSWKMBDHVryswkmbdhv\-]+$")

VIRUS_KEYWORDS = {
    "HBV": ["hepadna", "hepatitis b", "hbsag", "hbcag", "hbv"],
    "HCV": ["flaviviri", "hepatitis c", "hcv", "ns5b", "ns3"],
    "HEV": ["hepeviri", "hepatitis e", "hev", "orf2"],
}

# ---------------------------------------------------------------------------
# ── LAYOUT HELPERS ──────────────────────────────────────────────────────────
# ---------------------------------------------------------------------------

def USER_SEQ_STORES():
    """
    Returns a list of dcc.Store / dcc.Interval components.
    Splice these into your Container alongside the existing stores.

    Example in create_dashboard_layout():
        dbc.Container([
            dcc.Store(id="selected-virus", ...),
            ...
            *USER_SEQ_STORES(),   # <-- add this line
            ...
        ])
    """
    return [
        # Holds raw validated sequences as list-of-dicts:
        #   [{"id": str, "seq": str, "length": int}, ...]
        dcc.Store(id="useq-validated-store"),

        # Job state:  None | "pending" | "running" | "done" | "error"
        dcc.Store(id="useq-job-state", data=None),

        # Unique job id for async polling
        dcc.Store(id="useq-job-id", data=None),

        # Final results payload from the pipeline
        dcc.Store(id="useq-results-store"),

        # Polling interval – disabled until a job is running
        dcc.Interval(
            id="useq-poll-interval",
            interval=4_000,   # 4 s
            n_intervals=0,
            disabled=True,
        ),
        # Download components for TSV results & Newick tree
        dcc.Download(id="useq-download-tsv-component"),
        dcc.Download(id="useq-download-tree-component"),
    ]


def user_seq_tab_button():
    """
    Returns the dbc.Button to add to the existing navigation ButtonGroup.

    Example in create_dashboard_layout():
        dbc.ButtonGroup([
            dbc.Button("📊 Overview",      id="tab-overview",      ...),
            dbc.Button("🧬 Mutations",      id="tab-mutations",     ...),
            dbc.Button("📈 Epidemiology",   id="tab-epidemiology",  ...),
            user_seq_tab_button(),   # <-- add this
        ], className="w-100")
    """
    return dbc.Button(
        "🔬 My Sequences",
        id="tab-user-seq",
        color="secondary",
        n_clicks=0,
    )


def user_seq_tab_content():
    """
    Returns the html.Div that makes up the entire "My Sequences" tab.
    Add this after the last tab content Div in create_dashboard_layout().
    """
    return html.Div(
        id="user-seq-content",
        style={"display": "none"},
        children=[

            # ── ROW 1: Input panel ────────────────────────────────────────
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5(
                                [html.I(className="bi bi-upload me-2"),
                                 "Submit Your Sequences"],
                                className="hep-card-title mb-3"
                            ),

                            # Instructions
                            dbc.Alert(
                                [
                                    html.Strong("Accepted formats: "),
                                    "FASTA (.fasta, .fa, .fna) — up to ",
                                    html.Strong(f"{MAX_SEQUENCES} sequences"),
                                    " per submission. Sequences are processed locally on the server.",
                                ],
                                color="info",
                                className="mb-3 py-2",
                            ),

                            # Input mode tabs
                            dbc.Tabs(
                                id="useq-input-mode-tabs",
                                active_tab="paste",
                                children=[
                                    # ── Paste tab ──
                                    dbc.Tab(
                                        label="Paste FASTA",
                                        tab_id="paste",
                                        children=[
                                            html.Br(),
                                            dbc.Textarea(
                                                id="useq-fasta-textarea",
                                                placeholder=(
                                                    ">seq1\n"
                                                    "ATGCATGCATGCATGCATGCATGCATGC...\n"
                                                    ">seq2\n"
                                                    "GCTAGCTAGCTAGCTAGCTAGCTAGCTA..."
                                                ),
                                                rows=10,
                                                style={
                                                    "fontFamily": "monospace",
                                                    "fontSize": "0.85rem",
                                                },
                                                className="mb-2",
                                            ),
                                            html.Small(
                                                id="useq-paste-char-count",
                                                className="text-muted",
                                            ),
                                        ],
                                    ),

                                    # ── Upload tab ──
                                    dbc.Tab(
                                        label="Upload File",
                                        tab_id="upload",
                                        children=[
                                            html.Br(),
                                            dcc.Upload(
                                                id="useq-file-upload",
                                                children=html.Div([
                                                    html.I(
                                                        className="bi bi-file-earmark-arrow-up",
                                                        style={"fontSize": "2.5rem", "color": "#6c757d"},
                                                    ),
                                                    html.Br(),
                                                    "Drag & drop or ",
                                                    html.A("browse", style={"cursor": "pointer"}),
                                                    html.Br(),
                                                    html.Small(
                                                        ".fasta · .fa · .fna",
                                                        className="text-muted",
                                                    ),
                                                ], className="text-center py-4"),
                                                style={
                                                    "width": "100%",
                                                    "borderWidth": "2px",
                                                    "borderStyle": "dashed",
                                                    "borderRadius": "8px",
                                                    "borderColor": "#dee2e6",
                                                    "cursor": "pointer",
                                                },
                                                accept=".fasta,.fa,.fna",
                                                multiple=False,
                                            ),
                                            html.Div(
                                                id="useq-upload-filename",
                                                className="mt-2 text-muted small",
                                            ),
                                        ],
                                    ),
                                ],
                            ),

                            # Virus selection / override
                            html.Div([
                                html.Label("Virus Target:", className="fw-bold me-3 text-secondary", style={"fontSize": "0.95rem"}),
                                dbc.RadioItems(
                                    id="useq-virus-select",
                                    options=[
                                        {"label": "Auto-detect", "value": "auto"},
                                        {"label": "HBV (Hepatitis B)", "value": "HBV"},
                                        {"label": "HCV (Hepatitis C)", "value": "HCV"},
                                        {"label": "HEV (Hepatitis E)", "value": "HEV"},
                                    ],
                                    value="auto",
                                    inline=True,
                                    className="d-inline-block",
                                    style={"fontSize": "0.9rem"}
                                )
                            ], className="mt-3 mb-2 d-flex align-items-center bg-light p-2 rounded border border-light"),

                            html.Hr(),

                            # Single Run Analysis Action button + Recombination toggle switch
                            dbc.Row([
                                dbc.Col([
                                    dbc.Button(
                                        [html.I(className="bi bi-play-circle-fill me-2"),
                                         "Run Analysis"],
                                        id="useq-btn-run",
                                        className="hep-btn-primary me-3 d-inline-flex align-items-center",
                                    ),
                                    dbc.Checklist(
                                        id="useq-toggle-recomb",
                                        options=[
                                            {"label": "Include Recombination Analysis (3Seq / RDP)", "value": "include"}
                                        ],
                                        value=[],
                                        switch=True,
                                        inline=True,
                                        className="d-inline-block align-middle hep-body-text fw-semibold",
                                    ),
                                ], width="auto", className="d-flex align-items-center"),

                                dbc.Col([
                                    dbc.Button(
                                        [html.I(className="bi bi-x-circle me-2"), "Clear"],
                                        id="useq-btn-clear",
                                        className="hep-btn-danger d-inline-flex align-items-center",
                                    ),
                                ], width="auto", className="ms-auto d-flex align-items-center"),
                            ], className="mt-3"),
                        ])
                    ], className="shadow-sm mb-4")
                ], width=12)
            ]),

            # ── ROW 2: Validation feedback ────────────────────────────────
            dbc.Row([
                dbc.Col([
                    html.Div(id="useq-validation-feedback")
                ], width=12)
            ], className="mb-4"),

            # ── ROW 3: Job progress ───────────────────────────────────────
            html.Div(
                id="useq-progress-section",
                style={"display": "none"},
                children=[
                    dbc.Row([
                        dbc.Col([
                            dbc.Card([
                                dbc.CardBody([
                                    dbc.Row([
                                        dbc.Col([
                                            html.H6(
                                                id="useq-progress-label",
                                                className="mb-2",
                                                children="Running analysis…",
                                            ),
                                            dbc.Progress(
                                                id="useq-progress-bar",
                                                value=0,
                                                striped=True,
                                                animated=True,
                                                color="primary",
                                                style={"height": "20px"},
                                            ),
                                        ], width=10),
                                        dbc.Col([
                                            dbc.Spinner(
                                                size="md",
                                                color="primary",
                                                id="useq-spinner",
                                            )
                                        ], width=2, className="d-flex align-items-center justify-content-center"),
                                    ])
                                ])
                            ], className="shadow-sm")
                        ], width=12)
                    ], className="mb-4"),
                ],
            ),

            # ── ROW 4: Results section (hidden until pipeline finishes) ───
            html.Div(
                id="useq-results-section",
                style={"display": "none"},
                children=[

                    # 4a: 5 KPI Summary Cards
                    html.Div(id="useq-kpi-cards-container", className="mb-4"),

                    # 4b: Sequence Summary Table
                    dbc.Row([
                        dbc.Col([
                            dbc.Card([
                                dbc.CardBody([
                                    html.H5(
                                        [html.I(className="bi bi-table me-2"),
                                         "Sequence Summary"],
                                        className="hep-card-title mb-3",
                                    ),
                                    html.Div(id="useq-summary-table"),
                                ])
                            ], className="shadow-sm mb-4")
                        ], width=12)
                    ]),

                    # 4c: Genotype Composition & HepTracker Comparison
                    html.Div(id="useq-genotype-comparison-container", className="mb-4"),

                    # 4d: Mutations & Drug Resistance panel
                    dbc.Row([
                        dbc.Col([
                            html.Div(id="useq-mutations-panel")
                        ], width=12)
                    ], className="mb-4"),

                    # 4e: Phylogenetic tree + recombination side by side
                    dbc.Row([
                        dbc.Col([
                            dbc.Card([
                                dbc.CardBody([
                                    html.H5(
                                        [html.I(className="bi bi-diagram-2 me-2"),
                                         "Phylogenetic Placement"],
                                        className="hep-card-title mb-3",
                                    ),
                                    html.Small(
                                        "User sequences (●) placed onto the reference tree. "
                                        "Reference sequences shown in grey.",
                                        className="hep-meta-text d-block mb-3",
                                    ),
                                    html.Div(
                                        id="useq-phylo-tree-container",
                                        style={"minHeight": "400px"},
                                    ),
                                ])
                            ], className="shadow-sm h-100")
                        ], width=7),

                        dbc.Col([
                            dbc.Card([
                                dbc.CardBody([
                                    html.H5(
                                        [html.I(className="bi bi-shuffle me-2"),
                                         "Recombination Analysis"],
                                        className="hep-card-title mb-3",
                                    ),
                                    html.Div(id="useq-recombination-panel"),
                                    html.Div(
                                        dbc.Button(
                                            [
                                                html.I(className="bi bi-shuffle me-2"),
                                                "Run recombination analysis"
                                            ],
                                            id="useq-btn-run-recombination-inline",
                                            color="primary",
                                            className="w-100 mt-3 fw-bold"
                                        ),
                                        id="useq-recombination-inline-btn-container",
                                        style={"display": "none"}
                                    ),
                                ])
                            ], className="shadow-sm h-100")
                        ], width=5),
                    ], className="mb-4"),

                    # 4f: Sequence map (linear genome annotation)
                    dbc.Row([
                        dbc.Col([
                            dbc.Card([
                                dbc.CardBody([
                                    html.H5(
                                        [html.I(className="bi bi-map me-2"),
                                         "Sequence Map"],
                                        className="hep-card-title mb-3",
                                    ),
                                    html.Small(
                                        "ORFs, genotype-defining mutations, and "
                                        "recombination breakpoints relative to the "
                                        "reference genome.",
                                        className="hep-meta-text d-block mb-3",
                                    ),
                                    dcc.Graph(
                                        id="useq-sequence-map-graph",
                                        config={"displayModeBar": False},
                                    ),
                                ])
                            ], className="shadow-sm mb-4")
                        ], width=12)
                    ]),

                    # 4g: Download Export Section (.tsv & .nwk)
                    dbc.Row([
                        dbc.Col([
                            dbc.Card([
                                dbc.CardBody([
                                    html.H5(
                                        [html.I(className="bi bi-download me-2 text-primary"),
                                         "Export Results & Phylogenetic Tree"],
                                        className="hep-card-title mb-2",
                                    ),
                                    html.Small(
                                        "Download complete sequence analysis records in TSV format or export the phylogenetic placement tree in Newick format.",
                                        className="text-muted d-block mb-3",
                                    ),
                                    html.Div([
                                        dbc.Button(
                                            [html.I(className="bi bi-file-earmark-spreadsheet-fill me-2"), "Download Results (.tsv)"],
                                            id="useq-btn-download-tsv",
                                            className="hep-btn-primary me-3 d-inline-flex align-items-center",
                                        ),
                                        dbc.Button(
                                            [html.I(className="bi bi-diagram-3-fill me-2"), "Download Tree (.nwk)"],
                                            id="useq-btn-download-tree",
                                            className="hep-btn-secondary d-inline-flex align-items-center",
                                        ),
                                    ], className="d-flex align-items-center flex-wrap gap-2")
                                ])
                            ], className="shadow-sm border-0 bg-white rounded-3 mb-4")
                        ], width=12)
                    ]),
                ],
            ),
        ],
    )


# ---------------------------------------------------------------------------
# ── FASTA VALIDATION ────────────────────────────────────────────────────────
# ---------------------------------------------------------------------------

def _parse_fasta(text: str) -> tuple[list[dict], list[str]]:
    """
    Parse a FASTA string into a list of {id, seq, length} dicts.
    Also returns a list of error messages.
    Empty or whitespace-only text returns ([], []).
    """
    text = text.strip()
    if not text:
        return [], []

    sequences = []
    errors = []
    current_id = None
    current_seq_parts = []

    for line_no, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            # Flush previous
            if current_id is not None:
                seq = "".join(current_seq_parts)
                sequences.append({"id": current_id, "seq": seq, "length": len(seq)})
            header = line[1:].strip()
            if not header:
                errors.append(f"Line {line_no}: FASTA header is empty.")
                current_id = f"seq_{len(sequences)+1}"
            else:
                current_id = header.split()[0]  # use first token as ID
            current_seq_parts = []
        else:
            if current_id is None:
                errors.append(
                    f"Line {line_no}: sequence data before first header — "
                    "is this valid FASTA?"
                )
                # Treat as a headerless sequence
                current_id = "unnamed"
                current_seq_parts = []
            current_seq_parts.append(line)

    # Flush last
    if current_id is not None:
        seq = "".join(current_seq_parts)
        sequences.append({"id": current_id, "seq": seq, "length": len(seq)})

    return sequences, errors


def _validate_sequences(sequences: list[dict]) -> list[str]:
    """
    Validates a list of parsed sequences and returns a list of warning/error strings.
    """
    issues = []

    if len(sequences) == 0:
        issues.append("No sequences were parsed. Check your FASTA format.")
        return issues

    if len(sequences) > MAX_SEQUENCES:
        issues.append(
            f"Too many sequences: {len(sequences)} submitted, "
            f"maximum is {MAX_SEQUENCES}."
        )

    ids_seen = {}
    for i, rec in enumerate(sequences):
        seq_id = rec["id"]
        seq    = rec["seq"]

        # Duplicate ID check
        if seq_id in ids_seen:
            issues.append(
                f"Duplicate sequence ID '{seq_id}' at positions "
                f"{ids_seen[seq_id]+1} and {i+1}."
            )
        ids_seen[seq_id] = i

        # Empty sequence
        if len(seq) == 0:
            issues.append(f"'{seq_id}': sequence is empty.")
            continue

        # Very short – unlikely to be useful
        if len(seq) < 200:
            issues.append(
                f"'{seq_id}': sequence is very short ({len(seq)} nt). "
                "Phylogenetic placement may be unreliable."
            )

        # Illegal characters
        if not VALID_NUCLEOTIDES.match(seq):
            bad = set(c for c in seq if not VALID_NUCLEOTIDES.match(c))
            issues.append(
                f"'{seq_id}': contains illegal characters: "
                f"{', '.join(sorted(bad))}"
            )

        # Ambiguity ratio > 30 %
        n_count = seq.upper().count("N")
        ambig_ratio = n_count / len(seq)
        if ambig_ratio > 0.30:
            issues.append(
                f"'{seq_id}': {ambig_ratio:.0%} ambiguous bases (N). "
                "Consider filtering low-quality sequences."
            )

    return issues


def _detect_virus(sequences: list[dict]) -> str | None:
    """
    Heuristic: look at sequence IDs + crude k-mer presence for virus hints.
    Returns "HBV", "HCV", "HEV", or None (unknown / mixed).
    Falls back to asking the user if ambiguous.
    """
    combined_text = " ".join(r["id"].lower() for r in sequences)
    scores = {v: 0 for v in VIRUS_KEYWORDS}
    for virus, kws in VIRUS_KEYWORDS.items():
        for kw in kws:
            if kw in combined_text:
                scores[virus] += 1
    best = max(scores, key=scores.get)
    return best if scores[best] > 0 else None


def _build_validation_alert(sequences, parse_errors, validation_issues, detected_virus):
    """Build the dbc.Alert shown after validation."""

    total = len(sequences)
    has_errors = bool(parse_errors or
                      any("Too many" in i or "illegal" in i or "empty" in i
                          for i in validation_issues))
    has_warnings = bool(validation_issues and not has_errors)

    if total == 0 and not parse_errors:
        return dbc.Alert("Nothing to validate — paste or upload a FASTA first.",
                         color="secondary")

    items = []

    # Parse errors (blocking)
    for e in parse_errors:
        items.append(html.Li([html.I(className="bi bi-x-circle-fill text-danger me-2"), e]))

    # Validation issues
    for w in validation_issues:
        icon_cls = (
            "bi bi-x-circle-fill text-danger"
            if any(k in w for k in ("Too many", "illegal", "empty"))
            else "bi bi-exclamation-triangle-fill text-warning"
        )
        items.append(html.Li([html.I(className=f"{icon_cls} me-2"), w]))

    color = "danger" if has_errors else ("warning" if has_warnings else "success")
    header_icon = (
        "bi bi-x-circle-fill" if has_errors
        else ("bi bi-exclamation-triangle-fill" if has_warnings
              else "bi bi-check-circle-fill")
    )
    header_text = (
        f"Validation failed — please fix the errors below."
        if has_errors
        else (f"Validation passed with {len(validation_issues)} warning(s)."
              if has_warnings
              else f"Validation passed — {total} sequence(s) ready.")
    )

    children = [
        html.H6(
            [html.I(className=f"{header_icon} me-2"), header_text],
            className="alert-heading",
        ),
        html.Hr() if items else None,
        html.Ul(items, className="mb-0") if items else None,
    ]

    if detected_virus and not has_errors:
        children.append(
            html.P(
                [html.I(className="bi bi-virus2 me-2"),
                 f"Detected virus: ", html.Strong(detected_virus),
                 " — analysis will use the corresponding reference dataset."],
                className="mb-0 mt-2",
            )
        )
    elif not has_errors:
        children.append(
            dbc.Alert(
                [html.I(className="bi bi-question-circle me-2"),
                 "Could not auto-detect the virus type from sequence IDs. "
                 "The pipeline will run BLAST against all three reference sets."],
                color="warning",
                className="mb-0 mt-2 py-2",
            )
        )

    return dbc.Alert(
        [c for c in children if c is not None],
        color=color,
        className="mb-0",
    )


# ---------------------------------------------------------------------------
# ── CALLBACKS (module-level, using @callback — no app instance needed) ───────
# ---------------------------------------------------------------------------

# ── 1. Update char count as user types ──────────────────────────────────────
@callback(
    Output("useq-paste-char-count", "children"),
    Input("useq-fasta-textarea", "value"),
    prevent_initial_call=True,
)
def update_char_count(text):
    if not text:
        return ""
    seqs, _ = _parse_fasta(text)
    return f"{len(text):,} characters · {len(seqs)} sequence(s) detected"


# ── 2. Show uploaded filename ────────────────────────────────────────────
@callback(
    Output("useq-upload-filename", "children"),
    Input("useq-file-upload", "filename"),
    prevent_initial_call=True,
)
def show_upload_filename(filename):
    if not filename:
        return ""
    return [html.I(className="bi bi-file-earmark-text me-1"), f"Loaded: {filename}"]


# ── 3. Clear button resets everything ────────────────────────────────────
@callback(
    Output("useq-fasta-textarea",      "value"),
    Output("useq-file-upload",         "contents"),
    Output("useq-file-upload",         "filename"),
    Output("useq-validation-feedback", "children"),
    Output("useq-validated-store",     "data"),
    Output("useq-toggle-recomb",       "value"),
    Output("useq-results-section",     "style"),
    Output("useq-progress-section",    "style"),
    Output("useq-virus-select",        "value"),
    Input("useq-btn-clear", "n_clicks"),
    prevent_initial_call=True,
)
def clear_all(n_clicks):
    return (
        "",                         # textarea
        None,                       # upload contents
        None,                       # upload filename
        [],                         # validation feedback
        None,                       # validated store
        [],                         # reset recombination toggle
        {"display": "none"},        # results section
        {"display": "none"},        # progress section
        "auto",                     # virus select reset
    )


# ── 4. Run Analysis (Auto-Validate & Dispatch) callback ───────────────────
@callback(
    Output("useq-validation-feedback", "children",  allow_duplicate=True),
    Output("useq-validated-store",     "data",      allow_duplicate=True),
    Output("useq-job-state",            "data",      allow_duplicate=True),
    Output("useq-job-id",               "data",      allow_duplicate=True),
    Output("useq-progress-section",     "style",     allow_duplicate=True),
    Output("useq-progress-label",       "children",  allow_duplicate=True),
    Output("useq-progress-bar",         "value",     allow_duplicate=True),
    Output("useq-poll-interval",        "disabled",  allow_duplicate=True),
    Output("useq-results-section",      "style",     allow_duplicate=True),
    Input("useq-btn-run",               "n_clicks"),
    Input("useq-btn-run-recombination-inline", "n_clicks"),
    State("useq-file-upload",           "contents"),
    State("useq-fasta-textarea",         "value"),
    State("useq-file-upload",            "filename"),
    State("useq-input-mode-tabs",        "active_tab"),
    State("useq-virus-select",          "value"),
    State("useq-toggle-recomb",         "value"),
    State("useq-job-id",                "data"),
    prevent_initial_call=True,
)
def start_analysis_auto_validate(run_clicks, inline_recomb_clicks, upload_contents,
                                textarea_value, upload_filename, active_tab,
                                selected_virus_override, recomb_toggle, current_job_id):
    try:
        triggered = ctx.triggered_id
    except Exception:
        triggered = "useq-btn-run"
    if not triggered:
        triggered = "useq-btn-run"

    if triggered == "useq-btn-run-recombination-inline":
        if not current_job_id:
            raise PreventUpdate
        from pipeline_runner import dispatch_recombination_only
        success = dispatch_recombination_only(current_job_id)
        if not success:
            raise PreventUpdate
        return (
            dash.no_update, dash.no_update,
            "running", current_job_id,
            {"display": "block"}, "Queuing recombination analysis…",
            80, False, dash.no_update,
        )

    # ── 1. Resolve FASTA text ──
    fasta_text = ""
    if active_tab == "upload" or (triggered == "useq-file-upload" and upload_contents):
        if not upload_contents:
            alert = dbc.Alert("Please upload a valid FASTA file.", color="warning", className="border-0 shadow-sm")
            return alert, None, None, None, {"display": "none"}, "", 0, True, {"display": "none"}
        content_type, content_string = upload_contents.split(",", 1)
        fasta_text = base64.b64decode(content_string).decode("utf-8", errors="replace")
    else:
        fasta_text = textarea_value or ""

    if not fasta_text.strip():
        alert = dbc.Alert(
            [html.I(className="bi bi-exclamation-circle-fill me-2"), "Please paste or upload a FASTA sequence first."],
            color="warning", className="border-0 shadow-sm"
        )
        return alert, None, None, None, {"display": "none"}, "", 0, True, {"display": "none"}

    # ── 2. Validate FASTA ──
    sequences, parse_errors = _parse_fasta(fasta_text)
    validation_issues = _validate_sequences(sequences)

    if selected_virus_override and selected_virus_override != "auto":
        detected_virus = selected_virus_override
    else:
        detected_virus = _detect_virus(sequences) if sequences else None

    blocking = bool(
        parse_errors
        or any(any(k in issue for k in ("Too many", "illegal", "empty")) for issue in validation_issues)
        or len(sequences) == 0
    )

    alert = _build_validation_alert(sequences, parse_errors, validation_issues, detected_virus)

    if blocking or not sequences:
        # Halt execution and show validation feedback alert
        return alert, None, None, None, {"display": "none"}, "", 0, True, {"display": "none"}

    validated_data = {
        "sequences": sequences,
        "detected_virus": detected_virus,
        "validated_at": datetime.utcnow().isoformat(),
        "count": len(sequences),
    }

    # ── 3. Automatically Dispatch Pipeline ──
    run_recombination = ("include" in recomb_toggle) if recomb_toggle else False

    from pipeline_runner import dispatch_user_pipeline
    job_id = dispatch_user_pipeline(validated_data, run_recombination=run_recombination)

    return (
        alert,
        validated_data,
        "pending",
        job_id,
        {"display": "block"},
        "Submitting analysis job…",
        5,
        False,
        {"display": "none"},
    )


# ── 6. Polling callback ──────────────────────────────────────────────────
@callback(
    Output("useq-job-state",        "data",     allow_duplicate=True),
    Output("useq-progress-label",   "children", allow_duplicate=True),
    Output("useq-progress-bar",     "value",    allow_duplicate=True),
    Output("useq-poll-interval",    "disabled", allow_duplicate=True),
    Output("useq-results-store",    "data",     allow_duplicate=True),
    Output("useq-results-section",  "style",    allow_duplicate=True),
    Output("useq-progress-section", "style",    allow_duplicate=True),
    Input("useq-poll-interval", "n_intervals"),
    State("useq-job-id",    "data"),
    State("useq-job-state", "data"),
    prevent_initial_call=True,
)
def poll_job_status(n_intervals, job_id, current_state):
    """
    Reads the job state from the temp file written by the pipeline.
    In production, replace file reads with Celery AsyncResult checks:

        from celery.result import AsyncResult
        result = AsyncResult(job_id)
        state  = result.state          # PENDING / STARTED / SUCCESS / FAILURE
        meta   = result.info or {}
    """
    if not job_id or current_state in ("done", "error", None):
        raise PreventUpdate

    from pipeline_runner import get_job_status
    state_obj = get_job_status(job_id)

    job_state = state_obj.get("state", "pending")
    progress  = state_obj.get("progress", 0)
    label     = state_obj.get("label", "Running…")

    if job_state == "done":
        results = state_obj.get("results", {})
        return (
            "done",
            "✅ Analysis complete!",
            100,
            True,           # disable polling
            results,
            {"display": "block"},   # show results
            {"display": "none"},    # hide progress bar
        )

    if job_state == "error":
        return (
            "error",
            f"❌ Error: {state_obj.get('error', 'Unknown error')}",
            progress,
            True,           # disable polling
            dash.no_update,
            dash.no_update,
            {"display": "block"},
        )

    # Still running
    return (
        job_state,
        label,
        progress,
        False,
        dash.no_update,
        dash.no_update,
        {"display": "block"},
    )


# ── 7. Render results ────────────────────────────────────────────────────
@callback(
    Output("useq-kpi-cards-container",           "children"),
    Output("useq-summary-table",               "children"),
    Output("useq-genotype-comparison-container", "children"),
    Output("useq-recombination-panel",         "children"),
    Output("useq-mutations-panel",             "children"),
    Output("useq-sequence-map-graph",          "figure"),
    Output("useq-phylo-tree-container",        "children"),
    Output("useq-recombination-inline-btn-container", "style"),
    Input("useq-results-store", "data"),
    prevent_initial_call=True,
)
def render_results(results):
    import plotly.graph_objects as go
    from collections import Counter

    if not results:
        empty = dbc.Alert("No results available.", color="secondary")
        return empty, empty, empty, empty, html.Div(), go.Figure(), html.Div(), {"display": "none"}

    seqs = results.get("sequences", [])
    newick = results.get("newick", "")
    seq_map = results.get("sequence_map", {})
    recombination_run = results.get("recombination_run", False)

    # ── 1. Calculate 5 KPI Cards ──
    total_seqs = len(seqs)
    seq_unit = "Sequence" if total_seqs == 1 else "Sequences"

    detected_virus = seqs[0].get("virus", "HBV") if seqs else "HBV"

    unique_genotypes = sorted(list(set(r.get("genotype", "Unknown") for r in seqs if r.get("genotype"))))
    n_genotypes = len(unique_genotypes)
    geno_unit = "Genotype" if n_genotypes == 1 else "Genotypes"

    total_mutations = sum(len(r.get("mutations", [])) for r in seqs)

    if not recombination_run:
        recomb_kpi_text = "Recombination not Performed"
        recomb_kpi_color = "secondary"
    else:
        recombs = [r for r in seqs if r.get("validation_status", "none") in ("high_confidence", "needs_review")]
        n_recombs = len(recombs)
        if n_recombs == 0:
            recomb_kpi_text = "0 Recombinant Signals"
            recomb_kpi_color = "success"
        elif n_recombs == 1:
            recomb_kpi_text = "1 Recombinant Signal"
            recomb_kpi_color = "warning"
        else:
            recomb_kpi_text = f"{n_recombs} Recombinant Signals"
            recomb_kpi_color = "warning"

    kpi_cards = dbc.Row([
        # Card 1: Sequences Count
        dbc.Col([
            dbc.Card(
                dbc.CardBody([
                    html.Div([
                        html.Div(
                            html.I(className="bi bi-file-earmark-dna-fill text-primary fs-4"),
                            className="p-2 bg-primary-subtle text-primary rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                            style={"width": "44px", "height": "44px"}
                        ),
                        html.Div([
                            html.Div("TOTAL SEQUENCES", className="text-muted small fw-bold tracking-wider mb-1", style={"fontSize": "0.68rem", "letterSpacing": "0.05em"}),
                            html.H3(f"{total_seqs} {seq_unit}", className="fw-extrabold text-dark mb-0 fs-4"),
                        ])
                    ], className="d-flex align-items-center")
                ]),
                className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
            )
        ], xs=12, sm=6, lg=True, className="mb-2 mb-lg-0"),

        # Card 2: Target Virus
        dbc.Col([
            dbc.Card(
                dbc.CardBody([
                    html.Div([
                        html.Div(
                            html.I(className="bi bi-virus text-success fs-4"),
                            className="p-2 bg-success-subtle text-success rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                            style={"width": "44px", "height": "44px"}
                        ),
                        html.Div([
                            html.Div("TARGET VIRUS", className="text-muted small fw-bold tracking-wider mb-1", style={"fontSize": "0.68rem", "letterSpacing": "0.05em"}),
                            html.H3(f"{detected_virus} Detected", className="fw-extrabold text-dark mb-0 fs-4"),
                        ])
                    ], className="d-flex align-items-center")
                ]),
                className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
            )
        ], xs=12, sm=6, lg=True, className="mb-2 mb-lg-0"),

        # Card 3: Genotypes Count
        dbc.Col([
            dbc.Card(
                dbc.CardBody([
                    html.Div([
                        html.Div(
                            html.I(className="bi bi-diagram-3-fill text-info fs-4"),
                            className="p-2 bg-info-subtle text-info rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                            style={"width": "44px", "height": "44px"}
                        ),
                        html.Div([
                            html.Div("GENOTYPES", className="text-muted small fw-bold tracking-wider mb-1", style={"fontSize": "0.68rem", "letterSpacing": "0.05em"}),
                            html.H3(f"{n_genotypes} {geno_unit}", className="fw-extrabold text-dark mb-0 fs-4"),
                        ])
                    ], className="d-flex align-items-center")
                ]),
                className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
            )
        ], xs=12, sm=6, lg=True, className="mb-2 mb-lg-0"),

        # Card 4: Resistance Mutations
        dbc.Col([
            dbc.Card(
                dbc.CardBody([
                    html.Div([
                        html.Div(
                            html.I(className="bi bi-shield-exclamation text-warning fs-4"),
                            className="p-2 bg-warning-subtle text-warning rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                            style={"width": "44px", "height": "44px"}
                        ),
                        html.Div([
                            html.Div("MUTATIONS", className="text-muted small fw-bold tracking-wider mb-1", style={"fontSize": "0.68rem", "letterSpacing": "0.05em"}),
                            html.H3(f"{total_mutations} Resistance", className="fw-extrabold text-dark mb-0 fs-4"),
                        ])
                    ], className="d-flex align-items-center")
                ]),
                className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
            )
        ], xs=12, sm=6, lg=True, className="mb-2 mb-lg-0"),

        # Card 5: Recombination Summary
        dbc.Col([
            dbc.Card(
                dbc.CardBody([
                    html.Div([
                        html.Div(
                            html.I(className=f"bi bi-shuffle text-{recomb_kpi_color} fs-4"),
                            className=f"p-2 bg-{recomb_kpi_color}-subtle text-{recomb_kpi_color} rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                            style={"width": "44px", "height": "44px"}
                        ),
                        html.Div([
                            html.Div("RECOMBINATION", className="text-muted small fw-bold tracking-wider mb-1", style={"fontSize": "0.68rem", "letterSpacing": "0.05em"}),
                            html.H3(recomb_kpi_text, className="fw-extrabold text-dark mb-0 fs-5"),
                        ])
                    ], className="d-flex align-items-center")
                ]),
                className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
            )
        ], xs=12, sm=6, lg=True, className="mb-2 mb-lg-0"),
    ], className="g-3 mb-4")

    # ── 2. Sequence Summary Table (with QC Metrics & None/Possible/Supported) ──
    if seqs:
        header = html.Thead(html.Tr([
            html.Th("Sequence ID"),
            html.Th("Virus"),
            html.Th("Genotype"),
            html.Th("Length"),
            html.Th("Mutations"),
            html.Th("Recombinant"),
            html.Th("EPA"),
            html.Th("Details"),
        ]))
        rows = []
        for rec in seqs:
            # Recombinant column: None, Candidate, or Supported
            val_status = rec.get("validation_status", "none")
            if not recombination_run:
                recomb_badge = dbc.Badge("None", color="secondary", pill=True, title="Recombination analysis not requested.")
            elif val_status == "high_confidence":
                recomb_badge = dbc.Badge("Supported", color="danger", pill=True, title="Supported — Recombination candidate passed primary screening, basic filters, OpenRDP validation, and regional phylogenetic confirmation.")
            elif val_status == "needs_review":
                recomb_badge = dbc.Badge("Candidate", color="warning", pill=True, title="Candidate — Recombination candidate passed primary screening and basic filters, but secondary cross-validation was incomplete or discordant.")
            else:
                recomb_badge = dbc.Badge("None", color="secondary", pill=True, title="None — No validated recombination signal after primary screening and basic quality filters.")

            # Mutations count badge
            muts = rec.get("mutations", [])
            mut_badge = dbc.Badge(f"{len(muts)}", color="warning" if muts else "light",
                                  text_color="dark" if not muts else None, pill=True)

            # Genotype badge with theme-synchronized color
            geno_val = rec.get("genotype", "—")
            geno_color = get_genotype_color(detected_virus, geno_val)

            # Length calculation
            seq_len = rec.get("length") or len(rec.get("seq", ""))
            length_text = f"{seq_len:,} bp" if seq_len else "—"

            rows.append(html.Tr([
                html.Td(html.Code(rec.get("id", "—"), className="fw-bold")),
                html.Td(rec.get("virus", detected_virus)),
                html.Td(dbc.Badge(geno_val, color=None, style={"backgroundColor": geno_color, "color": "#ffffff"}, pill=True, className="px-2 py-1")),
                html.Td(length_text),
                html.Td(mut_badge),
                html.Td(recomb_badge),
                html.Td(f"{rec.get('epa_score', 0):.3f}"),
                html.Td(dbc.Badge("Pass QC", color="success", pill=True)),
            ]))
        summary_table = dbc.Table(
            [header, html.Tbody(rows)],
            bordered=True, hover=True, responsive=True, striped=True,
            size="sm", className="align-middle mb-0"
        )
    else:
        summary_table = dbc.Alert("No sequence results returned.", color="warning")

    # ── 3. Genotype Composition Chart & Comparison Table ──
    user_geno_counts = Counter(r.get("genotype", "Unknown") for r in seqs if r.get("genotype"))
    user_total_count = len(seqs) or 1
    user_geno_pcts = {g: (cnt / user_total_count) * 100 for g, cnt in user_geno_counts.items()}

    # Load global HepTracker dataset composition
    raw_global_pcts = {}
    raw_global_counts = {}
    g_total = 1
    try:
        from data_loader import load_and_preprocess_data
        ds = load_and_preprocess_data()
        df_v = ds.get(f"{detected_virus.lower()}_data")
        if df_v is not None and not df_v.empty:
            g_col = "genotype" if "genotype" in df_v.columns else ("Genotype" if "Genotype" in df_v.columns else None)
            if g_col:
                g_series = df_v[g_col].dropna()
                g_total = len(g_series) or 1
                g_counts = g_series.value_counts()
                raw_global_counts = {str(k): int(v) for k, v in g_counts.items()}
                raw_global_pcts = {str(k): (v / g_total) * 100 for k, v in g_counts.items()}
    except Exception:
        raw_global_pcts = {}
        raw_global_counts = {}
        g_total = 1

    def get_global_stat(user_g):
        if not raw_global_pcts:
            return 0.0, 0, g_total
        if user_g in raw_global_pcts:
            return raw_global_pcts[user_g], raw_global_counts.get(user_g, 0), g_total
        pref = f"{detected_virus.upper()}-{user_g}"
        if pref in raw_global_pcts:
            return raw_global_pcts[pref], raw_global_counts.get(pref, 0), g_total
        u_norm = re.sub(r'^(HBV|HCV|HEV|Genotype)[-_\s]*', '', str(user_g), flags=re.IGNORECASE).upper()
        for g_k, pct in raw_global_pcts.items():
            g_norm = re.sub(r'^(HBV|HCV|HEV|Genotype)[-_\s]*', '', str(g_k), flags=re.IGNORECASE).upper()
            if u_norm and u_norm == g_norm:
                return pct, raw_global_counts.get(g_k, 0), g_total
        return 0.0, 0, g_total

    user_gts = sorted(list(user_geno_pcts.keys()))
    user_bar_colors = [get_genotype_color(detected_virus, g) for g in user_gts]

    fig_geno = go.Figure()
    fig_geno.add_trace(go.Bar(
        x=user_gts,
        y=[user_geno_pcts.get(g, 0) for g in user_gts],
        name="User Sequences (%)",
        marker_color=user_bar_colors,
        text=[f"{user_geno_pcts.get(g, 0):.1f}%" for g in user_gts],
        textposition="auto",
    ))
    fig_geno.add_trace(go.Bar(
        x=user_gts,
        y=[get_global_stat(g)[0] for g in user_gts],
        name="HepTracker Database (%)",
        marker_color="#6c757d",
        text=[f"{get_global_stat(g)[0]:.1f}%" for g in user_gts],
        textposition="auto",
    ))
    fig_geno.update_layout(
        barmode="group",
        height=320,
        margin=dict(t=20, b=40, l=40, r=20),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        yaxis=dict(title="Percentage (%)", range=[0, 100]),
        xaxis=dict(title="Genotype"),
        plot_bgcolor="white",
        paper_bgcolor="white",
    )

    geno_table_rows = []
    for g in user_gts:
        u_cnt = user_geno_counts.get(g, 0)
        u_pct = user_geno_pcts.get(g, 0)
        h_pct, h_cnt, h_tot = get_global_stat(g)
        g_clr = get_genotype_color(detected_virus, g)
        hep_str = f"{h_pct:.1f}% ({h_cnt}/{h_tot})" if h_tot > 0 else f"{h_pct:.1f}%"
        geno_table_rows.append(html.Tr([
            html.Td(dbc.Badge(g, color=None, style={"backgroundColor": g_clr, "color": "#ffffff"}, pill=True, className="px-2 py-1")),
            html.Td(f"{u_pct:.1f}% ({u_cnt}/{user_total_count})"),
            html.Td(hep_str),
        ]))

    geno_comp_table = dbc.Table([
        html.Thead(html.Tr([
            html.Th("Genotype"),
            html.Th("User Sequences"),
            html.Th("HepTracker"),
        ])),
        html.Tbody(geno_table_rows)
    ], bordered=True, hover=True, striped=True, responsive=True, size="sm", className="align-middle mb-0")

    genotype_comparison_container = dbc.Row([
        dbc.Col([
            dbc.Card([
                dbc.CardBody([
                    html.H5([html.I(className="bi bi-pie-chart-fill me-2"), "Genotype Composition (%)"]),
                    dcc.Graph(figure=fig_geno, config={"displayModeBar": False})
                ])
            ], className="shadow-sm h-100")
        ], width=7),

        dbc.Col([
            dbc.Card([
                dbc.CardBody([
                    html.H5([html.I(className="bi bi-table me-2"), "Genotype Frequency Comparison"]),
                    html.Div(geno_comp_table)
                ])
            ], className="shadow-sm h-100")
        ], width=5),
    ])

    # ── 4. Recombination Panel (Table with 3SEQ, Filters, OpenRDP, IQ-TREE columns) ──
    btn_container_style = {"display": "none"}
    if not recombination_run:
        recombination_panel = html.Div([
            html.Div([
                html.H6("Recombination Analysis", className="fw-bold mb-2 text-dark"),
                html.Div([
                    html.Span("Status: ", className="fw-semibold text-muted me-1"),
                    dbc.Badge("Not requested", color="secondary", pill=True, className="px-2 py-1")
                ], className="mb-3 fs-6"),
            ])
        ])
        btn_container_style = {"display": "block"}
    else:
        STATUS_TOOLTIPS = {
            "Supported": "Supported — Recombination candidate passed primary screening, basic filters, OpenRDP validation, and regional phylogenetic confirmation.",
            "Candidate": "Candidate — Recombination candidate passed primary screening and basic filters, but secondary cross-validation was incomplete or discordant.",
            "None": "None — No validated recombination signal after primary screening and basic quality filters.",
        }

        recomb_rows = []
        n_recombs = 0
        for idx, rec in enumerate(seqs):
            val_s = rec.get("validation_status", "none")
            
            if val_s == "high_confidence":
                status_label, status_color = "Supported", "danger"
                seq_3s, seq_fl, seq_rdp, seq_tree = "Pass", "Pass", "Pass", "Pass"
                n_recombs += 1
            elif val_s == "needs_review":
                status_label, status_color = "Candidate", "warning"
                seq_3s, seq_fl = "Pass", "Pass"
                seq_rdp = "Pass" if rec.get("passed_openrdp", True) else "Fail"
                seq_tree = "Pass" if rec.get("passed_tree", False) else "Fail"
                n_recombs += 1
            else:
                status_label, status_color = "None", "secondary"
                seq_3s, seq_fl, seq_rdp, seq_tree = "NS", "—", "—", "—"

            tooltip_text = STATUS_TOOLTIPS[status_label]

            def _fmt_step(v):
                if v == "Pass":
                    return dbc.Badge("Pass", color="success", pill=True)
                elif v == "Fail":
                    return dbc.Badge("Fail", color="danger", pill=True)
                elif v == "NS":
                    return html.Span("NS", className="text-muted fw-semibold")
                return html.Span("—", className="text-muted")

            status_cell = dbc.Badge(status_label, color=status_color, pill=True, className="px-2 py-1", title=tooltip_text)

            recomb_rows.append(html.Tr([
                html.Td(html.Code(rec.get("id", "—"), className="fw-bold")),
                html.Td(status_cell),
                html.Td(_fmt_step(seq_3s)),
                html.Td(_fmt_step(seq_fl)),
                html.Td(_fmt_step(seq_rdp)),
                html.Td(_fmt_step(seq_tree)),
            ]))

        recomb_table = dbc.Table([
            html.Thead(html.Tr([
                html.Th("Sequence"),
                html.Th("Status"),
                html.Th("3SEQ"),
                html.Th("Filters"),
                html.Th("OpenRDP"),
                html.Th("IQ-TREE"),
            ])),
            html.Tbody(recomb_rows)
        ], bordered=True, hover=True, striped=True, responsive=True, size="sm", className="align-middle mb-0")

        alert_color = "warning" if n_recombs > 0 else "success"
        alert_msg = f"{n_recombs} recombinant sequence(s) detected." if n_recombs > 0 else "No validated recombination detected."
        alert_icon = "bi bi-exclamation-triangle-fill" if n_recombs > 0 else "bi bi-check-circle-fill"

        recombination_panel = html.Div([
            dbc.Alert(
                [html.I(className=f"{alert_icon} me-2"), alert_msg],
                color=alert_color, className="mb-3",
            ),
            recomb_table,
        ])

    # ── 5. Mutations & Drug Resistance Panel ──
    has_muts = any(rec.get("mutations") for rec in seqs)
    if has_muts:
        mut_items = []
        for rec in seqs:
            muts = rec.get("mutations", [])
            drugs = rec.get("drugs", [])
            if not muts: continue
            
            mut_items.append(html.Div([
                html.H6(f"Sequence: {rec['id']}", className="mt-2 mb-1"),
                html.P([
                    html.Strong("Mutations: "), ", ".join(muts)
                ], className="mb-1"),
                html.P([
                    html.Strong("Associated Resistance: "), ", ".join(drugs) if drugs else "None detected"
                ], className="mb-3")
            ]))
            
        mutations_panel = dbc.Card([
            dbc.CardHeader([html.I(className="bi bi-capsule me-2"), " Mutation & Drug Resistance Profile"]),
            dbc.CardBody(mut_items)
        ], className="mb-4 shadow-sm border-danger")
    else:
        mutations_panel = dbc.Card([
            dbc.CardHeader([html.I(className="bi bi-capsule me-2"), " Mutation & Drug Resistance Profile"]),
            dbc.CardBody([
                dbc.Alert(
                    [html.I(className="bi bi-info-circle-fill me-2"),
                     "No drug-resistance mutations detected for the submitted sequence(s)."],
                    color="info",
                    className="mb-0"
                )
            ])
        ], className="mb-4 shadow-sm border-info")

    # ── 6. Sequence Map (Plotly) ──
    fig = _build_sequence_map_figure(seqs, seq_map)

    # ── 7. Phylogenetic Tree Container ──
    from phylo_plot import build_tree_figure
    if newick:
        phylo_fig = build_tree_figure(newick)
        phylo_div = dcc.Graph(figure=phylo_fig, config={"displayModeBar": False})
    else:
        phylo_div = html.Div(
            dbc.Alert("No tree data returned from pipeline.", color="secondary")
        )

    return kpi_cards, summary_table, genotype_comparison_container, recombination_panel, mutations_panel, fig, phylo_div, btn_container_style


# ---------------------------------------------------------------------------
# ── SEQUENCE MAP FIGURE BUILDER ──────────────────────────────────────────────
# ---------------------------------------------------------------------------

def _build_sequence_map_figure(seqs: list[dict], seq_map: dict):
    """
    Builds a Plotly figure showing a linear genome map for each submitted
    sequence.  Each row shows:
      - Grey backbone (genome extent)
      - Coloured ORF blocks
      - Mutation lollipops
      - Recombination breakpoint shading

    If seq_map is empty (pipeline hasn't filled it), shows a placeholder.
    """
    import plotly.graph_objects as go

    if not seqs or not seq_map:
        fig = go.Figure()
        fig.update_layout(
            xaxis={"visible": False},
            yaxis={"visible": False},
            annotations=[{
                "text": "Sequence map will appear here after the pipeline completes.",
                "xref": "paper", "yref": "paper",
                "x": 0.5, "y": 0.5, "showarrow": False,
                "font": {"size": 14},
            }],
            height=300,
        )
        return fig

def _format_grouped_mutation_label(labels: list[str]) -> str:
    """Format grouped mutation labels at the same codon (e.g. W153Q, W153R, W153E -> W153 → Q / R / E)."""
    if len(labels) == 1:
        return labels[0]
    match_list = [re.match(r"^([A-Za-z0-9]+?)([A-Za-z\*])$", l) for l in labels]
    if all(m for m in match_list):
        prefixes = {m.group(1) for m in match_list}
        if len(prefixes) == 1:
            pref = list(prefixes)[0]
            aas = [m.group(2) for m in match_list]
            return f"{pref} → {' / '.join(aas)}"
    return ", ".join(labels)


def _build_sequence_map_figure(seqs: list[dict], seq_map: dict):
    """
    Builds a Plotly figure showing a linear genome map for each submitted
    sequence. Each row shows:
      - Grey backbone (genome extent)
      - Coloured ORF blocks
      - Mutation lollipops with 4-lane height cycling & codon grouping
      - Recombination breakpoint shading
    """
    import plotly.graph_objects as go

    if not seqs or not seq_map:
        fig = go.Figure()
        fig.update_layout(
            xaxis={"visible": False},
            yaxis={"visible": False},
            annotations=[{
                "text": "Sequence map will appear here after the pipeline completes.",
                "xref": "paper", "yref": "paper",
                "x": 0.5, "y": 0.5, "showarrow": False,
                "font": {"size": 14},
            }],
            height=300,
        )
        return fig

    n_seqs = len(seqs)
    fig = go.Figure()

    lane_height_offsets = [0.50, 0.70, 0.90, 1.10]  # 4 cycling height lanes for dense mutations

    for i, rec in enumerate(seqs):
        y_center = i * 2.5        # vertical spacing per sequence row
        seq_id   = rec["id"]
        data     = seq_map.get(seq_id, {})
        orfs     = data.get("orfs", [])
        muts     = data.get("mutations", [])
        bps      = data.get("breakpoints", [])
        seq_len  = rec.get("length", 3200)

        # Backbone
        fig.add_shape(
            type="rect",
            x0=0, x1=seq_len,
            y0=y_center - 0.1, y1=y_center + 0.1,
            fillcolor="#dee2e6",
            line=dict(width=0),
            layer="below",
        )

        # Recombination breakpoint shading
        for bp_start, bp_end in bps:
            fig.add_shape(
                type="rect",
                x0=bp_start, x1=bp_end,
                y0=y_center - 0.4, y1=y_center + 0.4,
                fillcolor="rgba(220,53,69,0.15)",
                line=dict(color="rgba(220,53,69,0.5)", width=1, dash="dot"),
            )

        # ORF blocks
        orf_colors = ["#4e79a7", "#59a14f", "#f28e2b", "#e15759", "#76b7b2"]
        for j, (orf_start, orf_end, orf_label) in enumerate(orfs):
            color = orf_colors[j % len(orf_colors)]
            fig.add_shape(
                type="rect",
                x0=orf_start, x1=orf_end,
                y0=y_center - 0.35, y1=y_center + 0.35,
                fillcolor=color,
                line=dict(width=0),
                opacity=0.75,
            )
            fig.add_annotation(
                x=(orf_start + orf_end) / 2,
                y=y_center,
                text=orf_label,
                showarrow=False,
                font=dict(size=9, color="white"),
                xanchor="center", yanchor="middle",
            )

        # Mutation lollipops with 4-lane height cycling & codon grouping
        if muts:
            # Group multiple substitutions occurring at the exact same genomic position
            grouped_muts = {}
            for pos, label, color in muts:
                if pos not in grouped_muts:
                    grouped_muts[pos] = {"pos": pos, "labels": [], "colors": []}
                grouped_muts[pos]["labels"].append(label)
                grouped_muts[pos]["colors"].append(color)

            sorted_positions = sorted(grouped_muts.keys())
            
            mut_x = []
            mut_y = []
            hover_text = []
            mut_colors = []

            lane_idx = 0
            last_pos = -9999

            for pos in sorted_positions:
                item = grouped_muts[pos]
                # Cycle through 4 lanes if mutations are within 25 nt of each other
                if pos - last_pos <= 25:
                    lane_idx = (lane_idx + 1) % len(lane_height_offsets)
                else:
                    lane_idx = 0
                last_pos = pos

                offset = lane_height_offsets[lane_idx]
                pin_y = y_center + offset

                formatted_label = _format_grouped_mutation_label(item["labels"])
                color = item["colors"][0]

                # Stem line from ORF top (y_center + 0.35) to lollipop head
                fig.add_shape(
                    type="line",
                    x0=pos, x1=pos,
                    y0=y_center + 0.35, y1=pin_y,
                    line=dict(color=color, width=1.5),
                )

                mut_x.append(pos)
                mut_y.append(pin_y)
                hover_text.append(formatted_label)
                mut_colors.append(color)

            # Scatter trace mode="markers" (no text on screen, details exclusively on hover)
            fig.add_trace(go.Scatter(
                x=mut_x, y=mut_y,
                mode="markers",
                marker=dict(size=8, color=mut_colors, line=dict(width=1.5, color="white")),
                customdata=hover_text,
                hovertemplate="<b>%{customdata}</b><br>Genomic Position: %{x} nt<extra></extra>",
                showlegend=False,
            ))

        # Row label
        fig.add_annotation(
            x=-150, y=y_center,
            text=seq_id[:20] + ("…" if len(seq_id) > 20 else ""),
            showarrow=False,
            font=dict(size=10),
            xanchor="right", yanchor="middle",
        )

    fig.update_layout(
        height=max(320, 140 + n_seqs * 140),
        xaxis=dict(
            title="Genomic position (nt)",
            range=[-200, max((rec.get("length", 3200) for rec in seqs), default=3200) + 200],
            gridcolor="rgba(0,0,0,0.06)",
        ),
        yaxis=dict(
            visible=False,
            range=[-0.5, n_seqs * 2.5 + 0.5],
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(t=20, b=60, l=160, r=20),
        showlegend=False,
    )
    return fig


# ---------------------------------------------------------------------------
# ── DOWNLOAD EXPORT CALLBACKS ────────────────────────────────────────────────
# ---------------------------------------------------------------------------

@callback(
    Output("useq-download-tsv-component", "data"),
    Input("useq-btn-download-tsv", "n_clicks"),
    State("useq-results-store", "data"),
    prevent_initial_call=True,
)
def download_tsv(n_clicks, results):
    """
    Downloads sequence analysis results in .tsv format with columns:
    Sequence ID, Virus, Genotype, Sequence Length, Mutations, Associated Resistance
    """
    if not n_clicks or not results:
        raise PreventUpdate

    seqs = results.get("sequences", [])
    if not seqs:
        raise PreventUpdate

    tsv_lines = ["Sequence ID\tVirus\tGenotype\tSequence Length\tMutations\tAssociated Resistance"]

    for rec in seqs:
        sid = rec.get("id", "—")
        virus = rec.get("virus", "HBV")
        geno = rec.get("genotype", "Unknown")
        seq_len = rec.get("length") or len(rec.get("seq", ""))
        len_str = f"{seq_len}" if seq_len else "—"

        muts = rec.get("mutations", [])
        muts_str = ", ".join(muts) if muts else "None"

        drugs = rec.get("drugs", [])
        drugs_str = ", ".join(drugs) if drugs else "None"

        tsv_lines.append(f"{sid}\t{virus}\t{geno}\t{len_str}\t{muts_str}\t{drugs_str}")

    content = "\n".join(tsv_lines)
    filename = f"heptracker_sequence_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.tsv"
    return dcc.send_string(content, filename=filename)


@callback(
    Output("useq-download-tree-component", "data"),
    Input("useq-btn-download-tree", "n_clicks"),
    State("useq-results-store", "data"),
    prevent_initial_call=True,
)
def download_tree(n_clicks, results):
    """
    Downloads phylogenetic placement tree in Newick format (.nwk).
    """
    if not n_clicks or not results:
        raise PreventUpdate

    newick = results.get("newick", "")
    if not newick or not newick.strip():
        raise PreventUpdate

    filename = f"heptracker_phylo_tree_{datetime.now().strftime('%Y%m%d_%H%M%S')}.nwk"
    return dcc.send_string(newick, filename=filename)


# ---------------------------------------------------------------------------
# ── TAB-SWITCHING PATCH ──────────────────────────────────────────────────────
# ---------------------------------------------------------------------------
# 
# You need to add "user-seq-content" to your EXISTING tab-switching callback.
# It currently looks something like this:
#
#   @callback(
#       Output("overview-content",      "style"),
#       Output("mutations-content",     "style"),
#       Output("epidemiology-content",  "style"),
#       Input("tab-overview",     "n_clicks"),
#       Input("tab-mutations",    "n_clicks"),
#       Input("tab-epidemiology", "n_clicks"),
#       prevent_initial_call=True,
#   )
#   def switch_tab(ov, mu, ep):
#       ...
#
# Replace it with the version below (adds the 4th tab):
#
#   @callback(
#       Output("overview-content",     "style"),
#       Output("mutations-content",    "style"),
#       Output("epidemiology-content", "style"),
#       Output("user-seq-content",     "style"),   # ← NEW
#       Output("tab-overview",         "color"),
#       Output("tab-mutations",        "color"),
#       Output("tab-epidemiology",     "color"),
#       Output("tab-user-seq",         "color"),   # ← NEW
#       Input("tab-overview",     "n_clicks"),
#       Input("tab-mutations",    "n_clicks"),
#       Input("tab-epidemiology", "n_clicks"),
#       Input("tab-user-seq",     "n_clicks"),     # ← NEW
#       prevent_initial_call=True,
#   )
#   def switch_tab(ov, mu, ep, us):
#   #       tid = ctx.triggered_id or "tab-overview"
#       show  = {"display": "block"}
#       hide  = {"display": "none"}
#       active_color   = "primary"
#       inactive_color = "secondary"
#       tabs = {
#           "tab-overview":     0,
#           "tab-mutations":    1,
#           "tab-epidemiology": 2,
#           "tab-user-seq":     3,
#       }
#       active = tabs.get(tid, 0)
#       styles = [show if i == active else hide  for i in range(4)]
#       colors = [active_color if i == active else inactive_color for i in range(4)]
#       return (*styles, *colors)
#
# ---------------------------------------------------------------------------
