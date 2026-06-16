import os
import re
import pandas as pd
import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
from flask import Flask, redirect
from data_loader import load_and_preprocess_data
from cache_config import cache

# Create Flask app first
server = Flask(__name__)
server.secret_key = os.environ.get('SECRET_KEY', 'default-secret-key')

# Initialize flask-caching
cache.init_app(server, config={
    'CACHE_TYPE': 'FileSystemCache',
    'CACHE_DIR': 'results/cache',
    'CACHE_DEFAULT_TIMEOUT': 86400  # 24 hours
})

# Add redirect before creating Dash app - ONLY ONE REDIRECT
@server.route('/')
def redirect_to_dashboard():
    return redirect('/dashboard', code=302)

app = dash.Dash(
    __name__,
    server=server,  # Use our Flask app
    use_pages=True,
    external_stylesheets=[dbc.themes.BOOTSTRAP,
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css"],
    suppress_callback_exceptions=True,
)

def navbar():
    reg = {p["name"]: p["path"] for p in dash.page_registry.values()}
    order = ["Dashboard", "About", "Resources", "Contact"]
    items = [dbc.NavItem(dbc.NavLink(name, href=reg[name], active="exact"))
             for name in order if name in reg]
    return dbc.Navbar(
        dbc.Container([
            dbc.NavbarBrand("Hepatitis", className="fw-bold"),
            dbc.Nav(items, pills=True, navbar=True, className="me-auto"),
            html.Div([
                dcc.Input(
                    id="global-search-input", 
                    placeholder="Search Accession (e.g. KF779250.1)...", 
                    type="text", 
                    className="form-control me-2", 
                    style={"width": "260px", "borderRadius": "20px"}
                ),
                dbc.Button(
                    [html.I(className="fa fa-search me-1"), "Search"], 
                    id="global-search-btn", 
                    color="light", 
                    outline=True, 
                    className="my-2 my-sm-0 fw-semibold", 
                    style={"borderRadius": "20px"}
                ),
            ], className="d-flex align-items-center")
        ]),
        color="primary", dark=True, sticky="top", className="mb-4",
    )

app.layout = dbc.Container(
    [
        navbar(),
        dash.page_container,
        # Global Search Modal
        dbc.Modal(
            [
                dbc.ModalHeader(dbc.ModalTitle("Sequence Details", id="search-modal-title", className="fw-bold")),
                dbc.ModalBody(id="search-modal-body"),
                dbc.ModalFooter(
                    dbc.Button("Close", id="search-modal-close", className="ms-auto", color="secondary")
                ),
            ],
            id="search-modal",
            size="lg",
            is_open=False,
        ),
    ],
    fluid=True,
)

app.server.config["DATA_STORE"] = load_and_preprocess_data()

# Global search callback
@app.callback(
    Output("search-modal", "is_open"),
    Output("search-modal-title", "children"),
    Output("search-modal-body", "children"),
    Input("global-search-btn", "n_clicks"),
    Input("global-search-input", "n_submit"),
    Input("search-modal-close", "n_clicks"),
    State("global-search-input", "value"),
    State("search-modal", "is_open"),
    prevent_initial_call=True
)
def handle_global_search(search_clicks, search_submit, close_clicks, search_term, is_open):
    ctx = dash.callback_context
    if not ctx.triggered:
        return is_open, "", ""
    
    trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]
    
    if trigger_id == "search-modal-close":
        return False, "", ""
        
    if trigger_id in ["global-search-btn", "global-search-input"] and search_term:
        term = str(search_term).strip()
        if not term:
            return is_open, "", ""
            
        # Standardize search term
        from data_loader import normalize_accession_id
        normalized_term = normalize_accession_id(term)
        
        # Load from DATA_STORE
        store = app.server.config.get("DATA_STORE", {})
        if not store:
            return True, "Error", html.Div("Data store not initialized.")
            
        # Search across HBV, HCV, HEV
        match_row = None
        virus_type = None
        
        for v in ["hbv", "hcv", "hev"]:
            df = store.get(f"{v}_data", pd.DataFrame())
            if not df.empty and "ID" in df.columns:
                res = df[df["ID"].apply(normalize_accession_id) == normalized_term]
                if not res.empty:
                    match_row = res.iloc[0]
                    virus_type = v.upper()
                    break
        
        if match_row is None:
            return True, "Search Results", dbc.Alert(
                f"No sequence found with ID: '{term}' (Normalized: '{normalized_term}')",
                color="warning",
                className="d-flex align-items-center"
            )
            
        # Found a match!
        seq_id = match_row.get("ID", term)
        genotype = match_row.get("genotype", "Unknown")
        subgenotype = match_row.get("subgenotype", "N/A")
        country = match_row.get("Country_standard", "Unknown")
        year = match_row.get("Year", "Unknown")
        region = match_row.get("WHO_Regions", "Unknown")
        length = match_row.get("length", "Unknown")
        is_recomb = str(match_row.get("is_recombinant", "false")).lower() == "true"
        
        # Determine genotype color
        if virus_type == "HBV":
            color = "#08519c"
        elif virus_type == "HCV":
            color = "#e6550d"
        else:
            color = "#2ca25f"
            
        # Query recombination details
        recomb_details = None
        if is_recomb:
            recomb_df = store.get(f"{virus_type.lower()}_recombs", pd.DataFrame())
            if not recomb_df.empty and "ID" in recomb_df.columns:
                rec_res = recomb_df[recomb_df["ID"].apply(normalize_accession_id) == normalized_term]
                if not rec_res.empty:
                    recomb_row = rec_res.iloc[0]
                    parent_1 = recomb_row.get("parent_1", "Unknown")
                    parent_2 = recomb_row.get("parent_2", "Unknown")
                    bp_start = recomb_row.get("breakpoint_start", "Unknown")
                    bp_end = recomb_row.get("breakpoint_end", "Unknown")
                    p_val = recomb_row.get("p_value", "N/A")
                    methods = recomb_row.get("methods", "N/A")
                    recomb_class = recomb_row.get("recombination_class", "unknown").replace("_", " ").title()
                    
                    try:
                        p_val_float = float(p_val)
                        p_val_str = f"{p_val_float:.3e}"
                    except:
                        p_val_str = str(p_val)
                        
                    recomb_details = dbc.Card([
                        dbc.CardHeader([
                            html.I(className="fa fa-dna me-2 text-warning"),
                            html.Span("Recombination Profile", className="fw-bold text-warning")
                        ], className="bg-dark text-white border-bottom border-warning"),
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    html.Div([
                                        html.Small("Recombination Class", className="text-muted d-block"),
                                        dbc.Badge(recomb_class, color="warning", text_color="dark", className="px-3 py-2 fs-6 fw-bold")
                                    ], className="mb-3")
                                ], width=6),
                                dbc.Col([
                                    html.Div([
                                        html.Small("Breakpoints", className="text-muted d-block"),
                                        html.Span(f"{bp_start} – {bp_end} bp", className="fw-bold fs-5")
                                    ], className="mb-3")
                                ], width=6),
                            ]),
                            dbc.Row([
                                dbc.Col([
                                    html.Div([
                                        html.Small("Parent 1 Reference", className="text-muted d-block"),
                                        html.Span(str(parent_1).replace("ref_", ""), className="fw-bold")
                                    ], className="mb-2")
                                ], width=6),
                                dbc.Col([
                                    html.Div([
                                        html.Small("Parent 2 Reference", className="text-muted d-block"),
                                        html.Span(str(parent_2).replace("ref_", ""), className="fw-bold")
                                    ], className="mb-2")
                                ], width=6),
                            ], className="border-bottom pb-2 mb-2"),
                            dbc.Row([
                                dbc.Col([
                                    html.Small("Statistical Confidence (p-value)", className="text-muted d-block"),
                                    html.Span(p_val_str, className="font-monospace fw-bold text-danger")
                                ], width=6),
                                dbc.Col([
                                    html.Small("Detection Methods", className="text-muted d-block"),
                                    html.Span(methods, className="small")
                                ], width=6),
                            ])
                        ])
                    ], className="border-warning border-2 shadow-sm mb-4")
                    
        # Query mutations
        mutation_rows = []
        mut_df = store.get(f"{virus_type.lower()}_mut", pd.DataFrame())
        if not mut_df.empty and "ID" in mut_df.columns:
            mut_res = mut_df[mut_df["ID"].apply(normalize_accession_id) == normalized_term]
            if not mut_res.empty:
                mutation_rows = mut_res.to_dict("records")
                
        mutations_details = None
        if mutation_rows:
            mut_list = []
            for row in mutation_rows:
                mut_name = row.get("mutation", "Unknown")
                gene = row.get("gene", "N/A")
                drug = row.get("drug", "N/A")
                mut_type = row.get("type", "N/A").replace("_", " ").title()
                mut_list.append(
                    html.Tr([
                        html.Td(html.Span(mut_name, className="badge bg-danger fs-6 fw-bold")),
                        html.Td(gene, className="fw-semibold"),
                        html.Td(drug, className="text-wrap"),
                        html.Td(mut_type, className="small text-muted")
                    ])
                )
                
            mutations_details = dbc.Card([
                dbc.CardHeader([
                    html.I(className="fa fa-triangle-exclamation me-2 text-danger"),
                    html.Span("Drug Resistance Mutations Detected", className="fw-bold text-danger")
                ], className="bg-dark text-white border-bottom border-danger"),
                dbc.CardBody([
                    dbc.Table([
                        html.Thead([
                            html.Tr([
                                html.Th("Mutation"),
                                html.Th("Gene"),
                                html.Th("Drug(s) Affected"),
                                html.Th("Classification")
                            ])
                        ], className="table-dark"),
                        html.Tbody(mut_list)
                    ], bordered=True, hover=True, responsive=True, size="sm")
                ])
            ], className="border-danger border-2 shadow-sm mb-4")
        else:
            if virus_type in ["HBV", "HCV"]:
                mutations_details = dbc.Card([
                    dbc.CardBody([
                        html.I(className="fa fa-circle-check text-success me-2 fs-5"),
                        html.Span("No drug-resistance mutations detected.", className="fw-semibold text-success")
                    ], className="d-flex align-items-center bg-success-subtle border border-success rounded")
                ], className="shadow-sm mb-4")
            else:
                mutations_details = dbc.Card([
                    dbc.CardBody([
                        html.I(className="fa fa-circle-info text-secondary me-2 fs-5"),
                        html.Span("Drug resistance mutations profiling is not available for HEV.", className="text-secondary")
                    ], className="d-flex align-items-center bg-light border rounded")
                ], className="shadow-sm mb-4")

        # Format details layout
        body = html.Div([
            # General Metadata
            dbc.Card([
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            html.Div([
                                html.Small("Virus type", className="text-muted d-block"),
                                html.Span(virus_type, className="fw-bold fs-4", style={"color": color})
                            ], className="mb-3")
                        ], width=4),
                        dbc.Col([
                            html.Div([
                                html.Small("Genotype / Subgenotype", className="text-muted d-block"),
                                html.Span(f"{genotype} / {subgenotype}", className="fw-bold fs-4")
                            ], className="mb-3")
                        ], width=4),
                        dbc.Col([
                            html.Div([
                                html.Small("Genome Length", className="text-muted d-block"),
                                html.Span(f"{length} bp" if pd.notna(length) else "N/A", className="fw-bold fs-4")
                            ], className="mb-3")
                        ], width=4),
                    ], className="border-bottom pb-2 mb-3"),
                    dbc.Row([
                        dbc.Col([
                            html.Div([
                                html.Small("Country of Origin", className="text-muted d-block"),
                                html.Span(country, className="fw-semibold fs-5")
                            ])
                        ], width=4),
                        dbc.Col([
                            html.Div([
                                html.Small("Collection Year", className="text-muted d-block"),
                                html.Span(str(year) if pd.notna(year) else "N/A", className="fw-semibold fs-5")
                            ])
                        ], width=4),
                        dbc.Col([
                            html.Div([
                                html.Small("WHO Region", className="text-muted d-block"),
                                html.Span(region, className="fw-semibold fs-5")
                            ])
                        ], width=4),
                    ])
                ])
            ], className="shadow-sm border-start border-5 mb-4", style={"borderLeftColor": color}),
            
            # Recombination
            recomb_details if recomb_details else dbc.Card([
                dbc.CardBody([
                    html.I(className="fa fa-circle-check text-success me-2 fs-5"),
                    html.Span("Recombination status: No recombination events detected.", className="fw-semibold text-success")
                ], className="d-flex align-items-center bg-success-subtle border border-success rounded")
            ], className="shadow-sm mb-4"),
            
            # Mutations
            mutations_details
        ])
        
        return True, f"Sequence ID: {seq_id}", body
        
    return is_open, "", ""

if __name__ == "__main__":
    app.run(debug=True, port=8051)