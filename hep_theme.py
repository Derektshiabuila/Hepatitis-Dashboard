"""
hep_theme.py
Shared light theme for the Hepatitis dashboard: brand colors, a registered
Plotly template ("hep_dark" registered for compatibility, tuned for light theme),
genotype/sequential color-scale builders, and shared dash_table style constants.
Import this once near the top of any page module and call register_theme() before building figures.
"""

import plotly.graph_objects as go
import plotly.io as pio
import colorsys

# ---------------------------------------------------------------------------
# BRAND PALETTE
# ---------------------------------------------------------------------------
VIRUS_COLORS = {
    "HAV": "#D97706",   # warm amber / orange
    "HBV": "#E11D48",   # vivid rose crimson
    "HCV": "#0D9488",   # deep teal
    "HDV": "#7C3AED",   # vivid violet
    "HEV": "#65A30D",   # lime green
    "RECOMBINANT": "#8B5CF6",
    "DANGER": "#EF4444",
}

BG = "#F5F6F8"
PANEL = "#FFFFFF"
PANEL_2 = "#F0F2F5"
BORDER = "rgba(15,23,42,0.10)"
TEXT = "#0F172A"
TEXT_DIM = "#5B6472"
TEXT_FAINT = "#8A94A6"

FONT_DISPLAY = '"Plus Jakarta Sans", sans-serif'
FONT_SANS = '"IBM Plex Sans", -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif'
FONT_MONO = '"IBM Plex Mono", Consolas, monospace'

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
    """5-stop sequential ramp for choropleths/heatmaps on a light basemap."""
    return [
        [0.0, "#F1F5F9"],
        [0.25, shade(hex_color, 65)],
        [0.5, shade(hex_color, 35)],
        [0.75, hex_color],
        [1.0, shade(hex_color, -20)],
    ]


def hex_to_rgb(hex_str):
    """Convert hex string (e.g. '#E84057') to fractional RGB tuple (r, g, b) in [0, 1]."""
    hex_str = hex_str.lstrip('#')
    return tuple(int(hex_str[i:i+2], 16) / 255.0 for i in (0, 2, 4))


def rgb_to_hex(rgb):
    """Convert fractional RGB tuple to hex string."""
    return '#' + ''.join('{:02x}'.format(max(0, min(255, int(round(x * 255.0))))) for x in rgb)


_GOLDEN_ANGLE = 0.6180339887  # fraction of a full turn


def _distinct_hue_colors(n, s=0.75, l=0.45, hue_start=0.58):
    """Generate n colors stepped around the color wheel by the golden angle,
    tuned for visibility on the light panel background."""
    colors = []
    for i in range(n):
        h = (hue_start + i * _GOLDEN_ANGLE) % 1.0
        r, g, b = colorsys.hls_to_rgb(h, l, s)
        colors.append(rgb_to_hex((r, g, b)))
    return colors


# Fixed, order-independent assignment so a given genotype letter/number always
# gets the same color across pages and re-renders.
_GENOTYPE_KEY_ORDER = list("ABCDEFGHIJKLMNOP") + [str(i) for i in range(1, 9)]
_GENOTYPE_COLOR_LOOKUP = dict(
    zip(_GENOTYPE_KEY_ORDER, _distinct_hue_colors(len(_GENOTYPE_KEY_ORDER)))
)


def genotype_palette(prefix, base_hex, keys):
    """Builds {"<PREFIX>-<key>": color, ..., "Recombinant": color} where each
    genotype maps to a universal, highly distinct color palette."""
    palette = {}
    fallback_colors = _distinct_hue_colors(max(len(keys), 8))

    for i, key in enumerate(keys):
        suffix = key.split('-')[-1]
        base_char = suffix[0].upper() if suffix else ""

        color = _GENOTYPE_COLOR_LOOKUP.get(base_char, fallback_colors[i % len(fallback_colors)])
        palette[f"{prefix}-{key}"] = color

    palette["Recombinant"] = VIRUS_COLORS["RECOMBINANT"]
    return palette


# ---------------------------------------------------------------------------
# DASH_TABLE STYLE CONSTANTS
# ---------------------------------------------------------------------------
TABLE_CELL_STYLE = {
    "textAlign": "left",
    "padding": "12px 14px",
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
    "fontWeight": "700",
    "fontFamily": FONT_SANS,
    "fontSize": "11px",
    "textTransform": "uppercase",
    "letterSpacing": "0.06em",
    "border": f"1px solid {BORDER}",
    "padding": "12px 14px",
}

TABLE_ODD_ROW_STYLE = {
    "if": {"row_index": "odd"},
    "backgroundColor": "#F8FAFC",
}


# ---------------------------------------------------------------------------
# PLOTLY TEMPLATE
# ---------------------------------------------------------------------------
def register_theme():
    """Registers the 'hep_dark' Plotly template (configured for light background) and makes it default."""
    template = go.layout.Template(
        layout=go.Layout(
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
            font=dict(family=FONT_SANS, color=TEXT_DIM, size=12),
            title=dict(font=dict(family=FONT_DISPLAY, color=TEXT, size=16, weight="bold")),
            colorway=[
                VIRUS_COLORS["HBV"],
                VIRUS_COLORS["HCV"],
                VIRUS_COLORS["HEV"],
                VIRUS_COLORS["HAV"],
                VIRUS_COLORS["HDV"],
                "#2563EB",
                "#DB2777",
                "#D97706",
            ],
            xaxis=dict(
                gridcolor="rgba(15,23,42,0.05)",
                linecolor="rgba(15,23,42,0.18)",
                zerolinecolor="rgba(15,23,42,0.12)",
                tickfont=dict(family=FONT_SANS, color=TEXT_DIM, size=11),
                title=dict(font=dict(family=FONT_SANS, color=TEXT_DIM, size=12, weight=600)),
            ),
            yaxis=dict(
                gridcolor="rgba(15,23,42,0.05)",
                linecolor="rgba(15,23,42,0.18)",
                zerolinecolor="rgba(15,23,42,0.12)",
                tickfont=dict(family=FONT_SANS, color=TEXT_DIM, size=11),
                title=dict(font=dict(family=FONT_SANS, color=TEXT_DIM, size=12, weight=600)),
            ),
            legend=dict(
                font=dict(family=FONT_SANS, color=TEXT_DIM, size=11),
                bgcolor="rgba(0,0,0,0)",
            ),
            hoverlabel=dict(
                bgcolor="#0F172A",
                bordercolor="rgba(255,255,255,0.1)",
                font=dict(family=FONT_SANS, color="#FFFFFF", size=12),
            ),
            geo=dict(
                bgcolor="rgba(0,0,0,0)",
                lakecolor="#E2E8F0",
                landcolor="#F1F5F9",
                showland=True,
                showframe=False,
                showcoastlines=True,
                coastlinecolor="rgba(15,23,42,0.12)",
                countrycolor="rgba(15,23,42,0.06)",
                showocean=True,
                oceancolor="#E2E8F0",
            ),
            coloraxis=dict(colorbar=dict(tickfont=dict(family=FONT_SANS, color=TEXT_DIM))),
        )
    )
    pio.templates["hep_dark"] = template
    pio.templates.default = "hep_dark"


# ---------------------------------------------------------------------------
# CENTRALIZED PLOTLY FIGURE TYPOGRAPHY & THEME HELPER
# ---------------------------------------------------------------------------
HEP_PLOT_FONT = dict(
    family="IBM Plex Sans, sans-serif",
    size=12,
    color="#475569",
)

HEP_TITLE_FONT = dict(
    family="Plus Jakarta Sans, sans-serif",
    size=16,
    color="#0F172A",
)

HEP_AXIS_TITLE_FONT = dict(
    family="IBM Plex Sans, sans-serif",
    size=12,
    color="#475569",
)

HEP_TICK_FONT = dict(
    family="IBM Plex Sans, sans-serif",
    size=11,
    color="#64748B",
)

HEP_LEGEND_FONT = dict(
    family="IBM Plex Sans, sans-serif",
    size=11,
    color="#475569",
)


def apply_heptracker_figure_style(fig):
    """Applies standardized HepTracker typography, grid, margin, and background styling to any Plotly figure."""
    if fig is None:
        return fig

    fig.update_layout(
        font=HEP_PLOT_FONT,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        legend=dict(
            font=HEP_LEGEND_FONT,
            title_font=HEP_LEGEND_FONT,
        ),
    )

    fig.update_xaxes(
        title_font=HEP_AXIS_TITLE_FONT,
        tickfont=HEP_TICK_FONT,
        gridcolor="#E2E8F0",
        linecolor="#CBD5E1",
    )

    fig.update_yaxes(
        title_font=HEP_AXIS_TITLE_FONT,
        tickfont=HEP_TICK_FONT,
        gridcolor="#E2E8F0",
        linecolor="#CBD5E1",
    )

    return fig

