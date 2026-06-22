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

def global_sidebar():
    reg = {p["name"]: p["path"] for p in dash.page_registry.values()}
    
    # Brand Header
    brand = html.Div(
        className="hep-sidebar-brand",
        children=[
            html.Div(html.I(className="fa-solid fa-dna"), className="hep-brand-icon"),
            html.Div([
                html.Div("HepTracker", className="hep-brand-title"),
                html.Div("GENOMICS DASHBOARD", className="hep-brand-subtitle"),
            ]),
        ],
    )
    
    # Dashboard Tools Section
    tools_section = html.Div(
        className="hep-sidebar-section",
        children=[
            html.Div("Dashboard tools:", className="hep-sidebar-label"),
            dbc.Button(
                [html.I(className="fa-solid fa-table-cells-large hep-nav-icon"), html.Span("Overview")],
                id="tab-overview",
                n_clicks=1,
                color="link",
                className="hep-nav-item hep-nav-active",
            ),
            dbc.Button(
                [html.I(className="fa-solid fa-dna hep-nav-icon"), html.Span("Mutations")],
                id="tab-mutations",
                n_clicks=0,
                color="link",
                className="hep-nav-item",
            ),
            dbc.Button(
                [html.I(className="fa-solid fa-wave-square hep-nav-icon"), html.Span("Epidemiology")],
                id="tab-epidemiology",
                n_clicks=0,
                color="link",
                className="hep-nav-item",
            ),
            dbc.Button(
                [html.I(className="fa-solid fa-flask hep-nav-icon"), html.Span("My Sequences")],
                id="tab-user-seq",
                n_clicks=0,
                color="link",
                className="hep-nav-item",
            ),
        ]
    )
    
    #Virus Selector Section
    virus_section = html.Div(
        className="hep-virus-panel",
        children=[
            html.Div("Virus selector:", className="hep-sidebar-label"),
            html.Div(
                [html.Span(className="hep-dot hep-dot-a"), html.Span("Hepatitis A")],
                className="hep-virus-item hep-virus-disabled",
            ),
            dbc.Button(
                [html.Span(className="hep-radio-dot"), html.Span("Hepatitis B")],
                id="btn-hbv",
                n_clicks=1,
                color="link",
                className="hep-virus-item hep-virus-active",
            ),
            dbc.Button(
                [html.Span(className="hep-dot hep-dot-c"), html.Span("Hepatitis C")],
                id="btn-hcv",
                n_clicks=0,
                color="link",
                className="hep-virus-item",
            ),
            html.Div(
                [html.Span(className="hep-dot hep-dot-d"), html.Span("Hepatitis D")],
                className="hep-virus-item hep-virus-disabled",
            ),
            dbc.Button(
                [html.Span(className="hep-dot hep-dot-e"), html.Span("Hepatitis E")],
                id="btn-hev",
                n_clicks=0,
                color="link",
                className="hep-virus-item",
            ),
        ]
    )
    
    # Pages Section
    pages_section = html.Div(
        className="hep-sidebar-section",
        style={"borderTop": "1px solid var(--border)", "paddingTop": "15px"},
        children=[
            html.Div("Pages:", className="hep-sidebar-label"),
            dcc.Link(
                [html.I(className="fa-solid fa-gauge hep-nav-icon"), html.Span("Dashboard")],
                id="link-dashboard",
                href=reg.get("Dashboard", "/dashboard"),
                className="hep-page-link hep-nav-item btn",
            ),
            dcc.Link(
                [html.I(className="fa-solid fa-circle-info hep-nav-icon"), html.Span("About")],
                id="link-about",
                href=reg.get("About", "/about"),
                className="hep-page-link hep-nav-item btn",
            ),
            dcc.Link(
                [html.I(className="fa-solid fa-book-open hep-nav-icon"), html.Span("Resources")],
                id="link-resources",
                href=reg.get("Resources", "/resources"),
                className="hep-page-link hep-nav-item btn",
            ),
            dcc.Link(
                [html.I(className="fa-solid fa-envelope hep-nav-icon"), html.Span("Contact")],
                id="link-contact",
                href=reg.get("Contact", "/contact"),
                className="hep-page-link hep-nav-item btn",
            ),
        ]
    )
    
    return html.Aside(
        id="sidebar-container",
        className="hep-sidebar",
        children=[
            brand,
            tools_section,
            virus_section,
            pages_section
        ]
    )

app.layout = html.Div(
    [
        dcc.Location(id="url", refresh=False),
        dcc.Store(id="selected-virus", data="HBV"),
        dcc.Store(id="dashboard-active-page-store", data=True),
        dcc.Store(id="active-tab-store", data="overview"),
        
        # Mobile top header (only visible on mobile/tablet)
        html.Div(
            className="mobile-header d-flex d-md-none align-items-center justify-content-between p-3",
            children=[
                dbc.Button(
                    html.I(className="fa-solid fa-bars"),
                    id="mobile-sidebar-toggle",
                    color="link",
                    className="text-white fs-3 p-0",
                ),
                html.Div(
                    "HepTracker",
                    className="fw-bold fs-4 text-white",
                    style={"fontFamily": "var(--serif)"}
                ),
                html.Div(style={"width": "30px"}), # Spacer to center
            ]
        ),
        
        # Sidebar Backdrop (for closing sidebar on tap on mobile)
        html.Div(id="sidebar-backdrop", className="sidebar-backdrop"),

        html.Div(
            className="app-shell-layout",
            children=[
                global_sidebar(),
                html.Div(
                    dash.page_container,
                    id="page-content-wrapper",
                    className="hep-main-content",
                ),
            ]
        ),

        dbc.Modal(
            [
                dbc.ModalHeader(
                    dbc.ModalTitle(
                        "Sequence Details",
                        id="search-modal-title",
                        className="fw-bold"
                    )
                ),
                dbc.ModalBody(id="search-modal-body"),
                dbc.ModalFooter(
                    dbc.Button(
                        "Close",
                        id="search-modal-close",
                        className="ms-auto",
                        color="secondary"
                    )
                ),
            ],
            id="search-modal",
            size="lg",
            is_open=False,
        ),
    ],
    className="app-shell",
)

#Global Navigation redirects
@app.callback(
    Output("url", "pathname"),
    Input("tab-overview", "n_clicks"),
    Input("tab-mutations", "n_clicks"),
    Input("tab-epidemiology", "n_clicks"),
    Input("tab-user-seq", "n_clicks"),
    Input("btn-hbv", "n_clicks"),
    Input("btn-hcv", "n_clicks"),
    Input("btn-hev", "n_clicks"),
    State("url", "pathname"),
    prevent_initial_call=True
)
def redirect_to_dashboard(*args):
    current_path = args[-1]
    if current_path not in ["/", "/dashboard"]:
        return "/dashboard"
    return dash.no_update

# Toggle mobile sidebar class
@app.callback(
    Output("sidebar-container", "className"),
    Input("mobile-sidebar-toggle", "n_clicks"),
    Input("sidebar-backdrop", "n_clicks"),
    Input("url", "pathname"),
    State("sidebar-container", "className"),
    prevent_initial_call=True
)
def toggle_sidebar_class(n_clicks_toggle, n_clicks_backdrop, pathname, current_className):
    ctx = dash.callback_context
    if not ctx.triggered:
        return "hep-sidebar"
    
    trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]
    
    if trigger_id in ["mobile-sidebar-toggle", "sidebar-backdrop"]:
        if "mobile-open" in current_className:
            return "hep-sidebar"
        else:
            return "hep-sidebar mobile-open"
            
    return "hep-sidebar"

# Highlight active Page links in sidebar
@app.callback(
    Output("link-dashboard", "className"),
    Output("link-about", "className"),
    Output("link-resources", "className"),
    Output("link-contact", "className"),
    Input("url", "pathname")
)
def update_sidebar_active_links(pathname):
    base_class = "hep-page-link hep-nav-item btn"
    active_class = "hep-page-link hep-nav-item btn hep-page-active"
    return (
        active_class if pathname in ["/", "/dashboard"] else base_class,
        active_class if pathname == "/about" else base_class,
        active_class if pathname == "/resources" else base_class,
        active_class if pathname == "/contact" else base_class,
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
        
        def get_base_id(accession_id):
            if not isinstance(accession_id, str):
                return accession_id
            return accession_id.split('.')[0].strip().upper()

        term_base = get_base_id(normalized_term)
        
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
                if res.empty:
                    res = df[df["ID"].apply(get_base_id) == term_base]
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
                if rec_res.empty:
                    rec_res = recomb_df[recomb_df["ID"].apply(get_base_id) == term_base]
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
            if mut_res.empty:
                mut_res = mut_df[mut_df["ID"].apply(get_base_id) == term_base]
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