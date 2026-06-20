"""
hep_theme.py
Shared dark-navy theme for the Hepatitis dashboard: brand colors, a registered
Plotly template ("hep_dark"), genotype/sequential color-scale builders, and
shared dash_table style constants. Import this once near the top of any page
module and call register_theme() before building figures.
"""

import plotly.graph_objects as go
import plotly.io as pio
import colorsys

# ---------------------------------------------------------------------------
# BRAND PALETTE
# ---------------------------------------------------------------------------
VIRUS_COLORS = {
    "HAV": "#F5A623",   # warm amber
    "HBV": "#E84057",   # vivid crimson
    "HCV": "#00D4AA",   # electric teal
    "HDV": "#A259FF",   # vivid violet
    "HEV": "#8BC34A",   # lime green
    "RECOMBINANT": "#A259FF",
    "DANGER": "#FF4D6D",
}

BG = "#0F1419"
PANEL = "#161D26"
PANEL_2 = "#1B232E"
BORDER = "rgba(255,255,255,0.10)"
TEXT = "#ECEFF2"
TEXT_DIM = "#8B98A5"
TEXT_FAINT = "#5C6773"

FONT_SANS = "IBM Plex Sans, Helvetica, Arial, sans-serif"
FONT_MONO = "IBM Plex Mono, monospace"

GENOTYPE_SHADE_OFFSETS = [0, -22, 24, -40, 42, -55, 14, -10, 30, -30]


def shade(hex_color, percent):
    """Lighten (positive percent) or darken (negative percent) a hex color."""
    num = int(hex_color.lstrip("#"), 16)
    r = (num >> 16) + round(2.55 * percent)
    g = ((num >> 8) & 0x00FF) + round(2.55 * percent)
    b = (num & 0x0000FF) + round(2.55 * percent)
    r = max(0, min(255, r))
    g = max(0, min(255, g))
    b = max(0, min(255, b))
    return "#{:02x}{:02x}{:02x}".format(r, g, b)


def sequential_scale(hex_color):
    """5-stop sequential ramp for choropleths/heatmaps on a dark basemap:
    low values fade toward the background, high values glow at full/peak
    saturation, rather than fading to white as a light-theme scale would."""
    return [
        [0.0, shade(hex_color, -75)],
        [0.25, shade(hex_color, -45)],
        [0.5, shade(hex_color, -15)],
        [0.75, hex_color],
        [1.0, shade(hex_color, 15)],
    ]


def hex_to_rgb(hex_str):
    """Convert hex string (e.g. '#E84057') to fractional RGB tuple (r, g, b) in [0, 1]."""
    hex_str = hex_str.lstrip('#')
    return tuple(int(hex_str[i:i+2], 16) / 255.0 for i in (0, 2, 4))


def rgb_to_hex(rgb):
    """Convert fractional RGB tuple to hex string."""
    return '#' + ''.join('{:02x}'.format(max(0, min(255, int(round(x * 255.0))))) for x in rgb)


def genotype_palette(prefix, base_hex, keys):
    """Builds {"<PREFIX>-<key>": color, ..., "Recombinant": color} where each
    genotype maps to a universal, highly distinct color palette. This ensures
    genotypes are easily differentiable and consistent across all viruses."""
    palette = {}
    
    universal_map = {
        "A": "#3B82F6", "B": "#F59E0B", "C": "#EF4444", "D": "#8B5CF6",
        "E": "#10B981", "F": "#EC4899", "G": "#2D9CDB", "H": "#00D4AA",
        "I": "#F5A623", "J": "#E84057",
        "1": "#3B82F6", "2": "#F59E0B", "3": "#EF4444", "4": "#8B5CF6",
        "5": "#10B981", "6": "#EC4899", "7": "#2D9CDB", "8": "#00D4AA"
    }
    
    fallback_colors = [
        "#3B82F6", "#F59E0B", "#EF4444", "#8B5CF6",
        "#10B981", "#EC4899", "#2D9CDB", "#00D4AA",
        "#F5A623", "#E84057"
    ]
    
    for i, key in enumerate(keys):
        suffix = key.split('-')[-1]
        base_char = suffix[0].upper() if suffix else ""
        
        if base_char in universal_map:
            color = universal_map[base_char]
        else:
            color = fallback_colors[i % len(fallback_colors)]
            
        palette[f"{prefix}-{key}"] = color
        
    palette["Recombinant"] = VIRUS_COLORS["RECOMBINANT"]
    return palette


# ---------------------------------------------------------------------------
# DASH_TABLE STYLE CONSTANTS
# ---------------------------------------------------------------------------
TABLE_CELL_STYLE = {
    "textAlign": "left",
    "padding": "10px",
    "overflow": "hidden",
    "textOverflow": "ellipsis",
    "maxWidth": 0,
    "backgroundColor": PANEL,
    "color": TEXT,
    "fontFamily": FONT_SANS,
    "fontSize": "13px",
    "border": f"1px solid {BORDER}",
}

TABLE_HEADER_STYLE = {
    "backgroundColor": PANEL_2,
    "color": TEXT,
    "fontWeight": "bold",
    "fontFamily": FONT_MONO,
    "fontSize": "12px",
    "textTransform": "uppercase",
    "letterSpacing": "0.04em",
    "border": f"1px solid {BORDER}",
}

TABLE_ODD_ROW_STYLE = {
    "if": {"row_index": "odd"},
    "backgroundColor": PANEL_2,
}


# ---------------------------------------------------------------------------
# PLOTLY TEMPLATE
# ---------------------------------------------------------------------------
def register_theme():
    """Registers the 'hep_dark' Plotly template and makes it the default."""
    template = go.layout.Template(
        layout=go.Layout(
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
            font=dict(family=FONT_SANS, color=TEXT_DIM, size=12),
            title=dict(font=dict(family=FONT_SANS, color=TEXT, size=16)),
            colorway=[
                VIRUS_COLORS["HBV"],
                VIRUS_COLORS["HCV"],
                VIRUS_COLORS["HEV"],
                VIRUS_COLORS["HAV"],
                VIRUS_COLORS["HDV"],
                "#4FAEFF",
                "#FF6B9D",
                "#FFD166",
            ],
            xaxis=dict(
                gridcolor="rgba(255,255,255,0.04)",
                linecolor="rgba(255,255,255,0.22)",
                zerolinecolor="rgba(255,255,255,0.15)",
                tickfont=dict(family=FONT_MONO, color=TEXT_DIM, size=11),
                title=dict(font=dict(family=FONT_SANS, color=TEXT_DIM)),
            ),
            yaxis=dict(
                gridcolor="rgba(255,255,255,0.04)",
                linecolor="rgba(255,255,255,0.22)",
                zerolinecolor="rgba(255,255,255,0.15)",
                tickfont=dict(family=FONT_MONO, color=TEXT_DIM, size=11),
                title=dict(font=dict(family=FONT_SANS, color=TEXT_DIM)),
            ),
            legend=dict(
                font=dict(family=FONT_SANS, color=TEXT_DIM, size=11),
                bgcolor="rgba(0,0,0,0)",
            ),
            hoverlabel=dict(
                bgcolor=PANEL_2,
                bordercolor=BORDER,
                font=dict(family=FONT_MONO, color=TEXT, size=12),
            ),
            geo=dict(
                bgcolor="rgba(0,0,0,0)",
                lakecolor="#07111A",
                landcolor="#101E2B",
                showland=True,
                showframe=False,
                showcoastlines=True,
                coastlinecolor="rgba(255,255,255,0.08)",
                countrycolor="rgba(255,255,255,0.06)",
                showocean=True,
                oceancolor="#07111A",
            ),
            coloraxis=dict(colorbar=dict(tickfont=dict(family=FONT_MONO, color=TEXT_DIM))),
        )
    )
    pio.templates["hep_dark"] = template
    pio.templates.default = "hep_dark"
