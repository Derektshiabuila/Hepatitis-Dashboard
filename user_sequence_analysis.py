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
                                className="mb-3"
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

                            # Validate + Submit buttons
                            dbc.Row([
                                dbc.Col([
                                    dbc.Button(
                                        [html.I(className="bi bi-check2-circle me-2"),
                                         "Validate Sequences"],
                                        id="useq-btn-validate",
                                        color="secondary",
                                        outline=True,
                                        className="me-2",
                                    ),
                                    dbc.Button(
                                        [html.I(className="bi bi-lightning-charge-fill me-2"),
                                         "Run Basic Analysis"],
                                        id="useq-btn-run",
                                        color="success",
                                        disabled=True,   # enabled only after validation passes
                                        className="me-2",
                                    ),
                                    dbc.Button(
                                        [html.I(className="bi bi-shuffle me-2"),
                                         "Run Recombination Analysis"],
                                        id="useq-btn-recomb",
                                        color="primary",
                                        disabled=True,   # enabled only after validation passes
                                    ),
                                ], width="auto"),

                                dbc.Col([
                                    dbc.Button(
                                        [html.I(className="bi bi-x-circle me-2"), "Clear"],
                                        id="useq-btn-clear",
                                        color="danger",
                                        outline=True,
                                        size="sm",
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

                    # 4a: Summary table
                    dbc.Row([
                        dbc.Col([
                            dbc.Card([
                                dbc.CardBody([
                                    html.H5(
                                        [html.I(className="bi bi-table me-2"),
                                         "Sequence Summary"],
                                        className="mb-3",
                                    ),
                                    html.Div(id="useq-summary-table"),
                                ])
                            ], className="shadow-sm mb-4")
                        ], width=12)
                    ]),

                    # 4b: Mutations & Drug Resistance panel
                    dbc.Row([
                        dbc.Col([
                            html.Div(id="useq-mutations-panel")
                        ], width=12)
                    ], className="mb-4"),

                    # 4c: Phylogenetic tree + recombination side by side
                    dbc.Row([
                        dbc.Col([
                            dbc.Card([
                                dbc.CardBody([
                                    html.H5(
                                        [html.I(className="bi bi-diagram-2 me-2"),
                                         "Phylogenetic Placement"],
                                        className="mb-3",
                                    ),
                                    html.Small(
                                        "User sequences (●) placed onto the reference tree. "
                                        "Reference sequences shown in grey.",
                                        className="text-muted d-block mb-3",
                                    ),
                                    # The tree is rendered by your front-end viewer;
                                    # this Div receives the Newick string as a data attribute
                                    # and is hydrated by a clientside callback (see below).
                                    html.Div(
                                        id="useq-phylo-tree-container",
                                        style={"minHeight": "400px"},
                                    ),
                                ])
                            ], className="shadow-sm")
                        ], width=7),

                        dbc.Col([
                            dbc.Card([
                                dbc.CardBody([
                                    html.H5(
                                        [html.I(className="bi bi-shuffle me-2"),
                                         "Recombination Analysis"],
                                        className="mb-3",
                                    ),
                                    html.Div(id="useq-recombination-panel"),
                                    html.Div(
                                        dbc.Button(
                                            [
                                                html.I(className="bi bi-shuffle me-2"),
                                                "Run Recombination Analysis (3Seq)"
                                            ],
                                            id="useq-btn-run-recombination-inline",
                                            color="primary",
                                            className="w-100 mt-3"
                                        ),
                                        id="useq-recombination-inline-btn-container",
                                        style={"display": "none"}
                                    ),
                                ])
                            ], className="shadow-sm")
                        ], width=5),
                    ], className="mb-4"),

                    # 4c: Sequence map (linear genome annotation)
                    dbc.Row([
                        dbc.Col([
                            dbc.Card([
                                dbc.CardBody([
                                    html.H5(
                                        [html.I(className="bi bi-map me-2"),
                                         "Sequence Map"],
                                        className="mb-3",
                                    ),
                                    html.Small(
                                        "ORFs, genotype-defining mutations, and "
                                        "recombination breakpoints relative to the "
                                        "reference genome.",
                                        className="text-muted d-block mb-3",
                                    ),
                                    dcc.Graph(
                                        id="useq-sequence-map-graph",
                                        config={"displayModeBar": False},
                                    ),
                                ])
                            ], className="shadow-sm")
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
    Output("useq-btn-run",             "disabled"),
    Output("useq-btn-recomb",          "disabled"),
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
        True,                       # run button disabled
        True,                       # recomb button disabled
        {"display": "none"},        # results section
        {"display": "none"},        # progress section
        "auto",                     # virus select reset
    )


# ── 4. Validate callback ─────────────────────────────────────────────────
# Triggered by the Validate button OR when a file is uploaded.
@callback(
    Output("useq-validation-feedback", "children",  allow_duplicate=True),
    Output("useq-validated-store",     "data",      allow_duplicate=True),
    Output("useq-btn-run",             "disabled",  allow_duplicate=True),
    Output("useq-btn-recomb",          "disabled",  allow_duplicate=True),
    Input("useq-btn-validate",  "n_clicks"),
    Input("useq-file-upload",   "contents"),
    State("useq-fasta-textarea", "value"),
    State("useq-file-upload",    "filename"),
    State("useq-input-mode-tabs", "active_tab"),
    State("useq-virus-select",  "value"),
    prevent_initial_call=True,
)
def validate_sequences(n_clicks, upload_contents, textarea_value,
                       upload_filename, active_tab, selected_virus_override):
    """
    1. Determine the FASTA source (paste vs upload).
    2. Parse the FASTA.
    3. Validate sequences.
    4. Store validated sequences & enable Run buttons if no blocking errors.
    """
    triggered = ctx.triggered_id

    # ── Resolve FASTA text ──
    fasta_text = ""

    if triggered == "useq-file-upload" or active_tab == "upload":
        if not upload_contents:
            raise PreventUpdate
        # Decode base64 file contents (dcc.Upload always returns base64)
        content_type, content_string = upload_contents.split(",", 1)
        decoded = base64.b64decode(content_string).decode("utf-8", errors="replace")
        fasta_text = decoded
    else:
        # Paste mode
        fasta_text = textarea_value or ""

    if not fasta_text.strip():
        alert = dbc.Alert(
            "Please paste or upload a FASTA sequence first.",
            color="secondary",
        )
        return alert, None, True, True

    # ── Parse ──
    sequences, parse_errors = _parse_fasta(fasta_text)

    # ── Validate ──
    validation_issues = _validate_sequences(sequences)

    # ── Detect virus ──
    if selected_virus_override and selected_virus_override != "auto":
        detected_virus = selected_virus_override
    else:
        detected_virus = _detect_virus(sequences) if sequences else None

    # ── Determine if blocking errors exist ──
    blocking = bool(
        parse_errors
        or any(
            any(k in issue for k in ("Too many", "illegal", "empty"))
            for issue in validation_issues
        )
        or len(sequences) == 0
    )

    # ── Build alert UI ──
    alert = _build_validation_alert(
        sequences, parse_errors, validation_issues, detected_virus
    )

    # ── Store & button state ──
    store_data = None
    btn_disabled = True

    if not blocking and sequences:
        store_data = {
            "sequences": sequences,
            "detected_virus": detected_virus,
            "validated_at": datetime.utcnow().isoformat(),
            "count": len(sequences),
        }
        btn_disabled = False

    return alert, store_data, btn_disabled, btn_disabled


# ── 5. Run Analysis buttons (Basic, Recomb, and Inline Recomb) ───────────
@callback(
    Output("useq-job-state",        "data",     allow_duplicate=True),
    Output("useq-job-id",           "data",     allow_duplicate=True),
    Output("useq-progress-section", "style",    allow_duplicate=True),
    Output("useq-progress-label",   "children", allow_duplicate=True),
    Output("useq-progress-bar",     "value",    allow_duplicate=True),
    Output("useq-poll-interval",    "disabled", allow_duplicate=True),
    Output("useq-results-section",  "style",    allow_duplicate=True),
    Input("useq-btn-run",           "n_clicks"),
    Input("useq-btn-recomb",        "n_clicks"),
    Input("useq-btn-run-recombination-inline", "n_clicks"),
    State("useq-validated-store",   "data"),
    State("useq-job-id",            "data"),
    prevent_initial_call=True,
)
def start_analysis(run_clicks, recomb_clicks, inline_clicks, validated_data, current_job_id):
    triggered = ctx.triggered_id
    if not triggered:
        raise PreventUpdate

    if triggered == "useq-btn-run":
        if not validated_data:
            raise PreventUpdate
        from pipeline_runner import dispatch_user_pipeline
        job_id = dispatch_user_pipeline(validated_data, run_recombination=False)
        return (
            "pending",
            job_id,
            {"display": "block"},
            "Submitting job…",
            5,
            False,
            {"display": "none"},
        )

    elif triggered == "useq-btn-recomb":
        if not validated_data:
            raise PreventUpdate
        from pipeline_runner import dispatch_user_pipeline
        job_id = dispatch_user_pipeline(validated_data, run_recombination=True)
        return (
            "pending",
            job_id,
            {"display": "block"},
            "Submitting job…",
            5,
            False,
            {"display": "none"},
        )

    elif triggered == "useq-btn-run-recombination-inline":
        if not current_job_id:
            raise PreventUpdate
        from pipeline_runner import dispatch_recombination_only
        success = dispatch_recombination_only(current_job_id)
        if not success:
            raise PreventUpdate
        return (
            "running",
            current_job_id,
            {"display": "block"},
            "Queuing recombination analysis…",
            80,
            False,
            dash.no_update,
        )

    raise PreventUpdate


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
    Output("useq-summary-table",         "children"),
    Output("useq-recombination-panel",   "children"),
    Output("useq-mutations-panel",       "children"),
    Output("useq-sequence-map-graph",    "figure"),
    Output("useq-phylo-tree-container",  "children"),
    Output("useq-recombination-inline-btn-container", "style"),
    Input("useq-results-store", "data"),
    prevent_initial_call=True,
)
def render_results(results):
    """
    Populates the four result panels from the pipeline output dict.

    Expected results schema
    -----------------------
    {
      "sequences": [
        {
          "id":          str,           # sequence ID
          "virus":       str,           # "HBV" | "HCV" | "HEV"
          "genotype":    str,           # e.g. "HBV-D"
          "is_recombinant": bool,
          "breakpoints": [[int, int]],  # genomic positions, may be []
          "nearest_ref": str,           # accession of nearest reference
          "epa_score":   float,
        },
        ...
      ],
      "newick": str,                    # Newick string of pruned subtree
      "sequence_map": {                 # per-sequence annotation data for Plotly
        "<seq_id>": {
          "orfs":        [[start, end, label], ...],
          "mutations":   [[pos, label, color], ...],
          "breakpoints": [[start, end], ...],
        }
      }
    }
    """
    import plotly.graph_objects as go

    if not results:
        empty = dbc.Alert("No results available.", color="secondary")
        return empty, empty, html.Div(), go.Figure(), html.Div(), {"display": "none"}

    seqs = results.get("sequences", [])
    newick = results.get("newick", "")
    seq_map = results.get("sequence_map", {})

    # ── 7a: Summary table ──
    if seqs:
        header = html.Thead(html.Tr([
            html.Th("Sequence ID"),
            html.Th("Virus"),
            html.Th("Genotype"),
            html.Th("Recombinant"),
            html.Th("Nearest Reference"),
            html.Th("EPA Score"),
        ]))
        rows = []
        for rec in seqs:
            status = rec.get("validation_status", "none")
            if status == "high_confidence":
                badge_color = "danger"
                badge_text = "Recombinant: High confidence"
            elif status == "needs_review":
                badge_color = "warning"
                badge_text = "Candidate recombinant: Needs review"
            else:
                badge_color = "success"
                badge_text = "No validated recombination detected"
            
            rows.append(html.Tr([
                html.Td(html.Code(rec.get("id", "—"))),
                html.Td(rec.get("virus", "—")),
                html.Td(dbc.Badge(rec.get("genotype", "—"),
                                  color="primary", pill=True)),
                html.Td(dbc.Badge(badge_text, color=badge_color, pill=True)),
                html.Td(html.Code(rec.get("nearest_ref", "—"))),
                html.Td(f"{rec.get('epa_score', 0):.3f}"),
            ]))
        summary_table = dbc.Table(
            [header, html.Tbody(rows)],
            bordered=True, hover=True, responsive=True, striped=True,
            size="sm",
        )
    else:
        summary_table = dbc.Alert("No sequence results returned.", color="warning")

    # ── 7b: Recombination panel ──
    recombination_run = results.get("recombination_run", False)
    btn_container_style = {"display": "none"}
    
    if not recombination_run:
        recombination_panel = dbc.Alert(
            [
                html.I(className="bi bi-info-circle-fill me-2"),
                "Recombination analysis (3Seq) was not run for this job."
            ],
            color="info", className="mb-0"
        )
        btn_container_style = {"display": "block"}
    else:
        recombinants = [r for r in seqs if r.get("validation_status", "none") in ("high_confidence", "needs_review")]
        if recombinants:
            rec_items = []
            for rec in recombinants:
                bps = rec.get("breakpoints", [])
                bp_text = (
                    ", ".join(f"{s}–{e} nt" for s, e in bps)
                    if bps else "Breakpoints not resolved"
                )
                status = rec.get("validation_status", "none")
                status_text = (
                    "High confidence" if status == "high_confidence"
                    else "Needs review"
                )
                status_color = "danger" if status == "high_confidence" else "warning"
                
                rec_items.append(
                    dbc.ListGroupItem([
                        html.Strong(rec["id"]),
                        html.Span(f" — {rec.get('genotype','?')}",
                                  className="text-muted ms-1"),
                        dbc.Badge(status_text, color=status_color, className="ms-2", pill=True),
                        html.Br(),
                        html.Small([
                            html.I(className="bi bi-scissors me-1"),
                            f"Breakpoints: {bp_text}"
                        ], className="text-muted"),
                    ])
                )
            recombination_panel = html.Div([
                dbc.Alert(
                    [html.I(className="bi bi-exclamation-triangle-fill me-2"),
                     f"{len(recombinants)} recombinant sequence(s) detected."],
                    color="warning", className="mb-2",
                ),
                dbc.ListGroup(rec_items, flush=True),
            ])
        else:
            recombination_panel = dbc.Alert(
                [html.I(className="bi bi-check-circle-fill me-2"),
                 "No validated recombination detected."],
                color="success",
            )

    # ── 7c: Mutations & Resistance panel ──
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

    # ── 7d: Sequence map (Plotly) ──
    fig = _build_sequence_map_figure(seqs, seq_map)

    # ── 7e: Phylogenetic tree container ──
    from phylo_plot import build_tree_figure
    if newick:
        phylo_fig = build_tree_figure(newick)
        phylo_div = dcc.Graph(figure=phylo_fig, config={"displayModeBar": False})
    else:
        phylo_div = html.Div(
            dbc.Alert("No tree data returned from pipeline.", color="secondary")
        )

    return summary_table, recombination_panel, mutations_panel, fig, phylo_div, btn_container_style


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

    n_seqs = len(seqs)
    row_height = 1.0
    fig = go.Figure()

    for i, rec in enumerate(seqs):
        y_center = i * 2        # vertical spacing
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

        # Mutation lollipops
        for pos, label, color in muts:
            fig.add_shape(
                type="line",
                x0=pos, x1=pos,
                y0=y_center + 0.35, y1=y_center + 0.7,
                line=dict(color=color, width=1.5),
            )
            fig.add_trace(go.Scatter(
                x=[pos], y=[y_center + 0.75],
                mode="markers+text",
                marker=dict(size=8, color=color),
                text=[label],
                textposition="top center",
                textfont=dict(size=8),
                hovertemplate=f"<b>{label}</b><br>Position: {pos}<extra></extra>",
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
        height=max(300, 120 + n_seqs * 120),
        xaxis=dict(
            title="Genomic position (nt)",
            range=[-200, max((rec.get("length", 3200) for rec in seqs), default=3200) + 200],
            gridcolor="rgba(0,0,0,0.06)",
        ),
        yaxis=dict(
            visible=False,
            range=[-1, n_seqs * 2 + 0.5],
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(t=20, b=60, l=160, r=20),
        showlegend=False,
    )
    return fig


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
