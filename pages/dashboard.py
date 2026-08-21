import dash
from dash import html, register_page, dcc, Input, Output, State, dash_table, callback
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from datetime import datetime
import pycountry
from functools import lru_cache
from dash import callback_context
from concurrent.futures import ThreadPoolExecutor
import re
import json
from dash.exceptions import PreventUpdate
from plotly.subplots import make_subplots
import plotly.io as pio

from hep_theme import (
    VIRUS_COLORS, sequential_scale, genotype_palette,
    TABLE_HEADER_STYLE, TABLE_CELL_STYLE, TABLE_ODD_ROW_STYLE,
    register_theme, shade, _distinct_hue_colors, _GENOTYPE_COLOR_LOOKUP,
    apply_heptracker_figure_style,
)
register_theme()


register_page(__name__, path="/", name="Dashboard", order=0)
register_page(__name__, path="/dashboard", name="Dashboard")
# callback decorator imported directly from dash

# Import data loading functions
from data_loader import load_and_preprocess_data
from user_sequence_analysis import (
    USER_SEQ_STORES,
    user_seq_tab_button,
    user_seq_tab_content,
)

# ---------------------------------------------------------------------------
# CONSTANTS
# ---------------------------------------------------------------------------
MAX_SEQUENCES = 50
MIN_SEQUENCES = 1
VALID_NUCLEOTIDES = re.compile(r"^[ACGTNacgtnRYSWKMBDHVryswkmbdhv\-]+$")

VIRUS_KEYWORDS = {
    "HBV": ["hepadna", "hepatitis b", "Whbsag", "hbcag", "hbv"],
    "HCV": ["flaviviri", "hepatitis c", "hcv", "ns5b", "ns3"],
    "HEV": ["hepeviri", "hepatitis e", "hev", "orf2"],
}

# Colorscale from red to green
coverage_colorscale = [
    [0.0, "#1a9850"],   # Red (lowest coverage)
    [0.25, "#d9ef8b"],
    [0.5, "#fee08b"],
    [0.75, "#fc8d59"],
    [1.0, "#d73027"]    # Dark green (highest coverage)
]

def _stepped_colorscale(colors):
    """Build a Plotly colorscale with hard edges between N colors (no blending),
    so each bin renders as a solid block instead of a gradient — matches the
    reference legend style rather than Plotly's default smooth interpolation."""
    n = len(colors)
    scale = []
    for i, c in enumerate(colors):
        scale.append([i / n, c])
        scale.append([(i + 1) / n, c])
    return scale

SEQUENCE_COUNT_COLORS = ["#EAF2FB", "#C3DBF2", "#9AC2E8", "#6FA8DC", "#3D7FC4", "#1F5FA8", "#123F73", "#0A2647"]
SEQUENCE_COUNT_STEPPED_COLORSCALE = _stepped_colorscale(SEQUENCE_COUNT_COLORS)

# Initialize data store (will be loaded on first access)
data_store = None

def get_data_store():
    """Helper function to access the global data store"""
    global data_store
    if data_store is None:
        try:
            from data_loader import load_and_preprocess_data
            data_store = load_and_preprocess_data()
            print("Data store loaded successfully")
        except Exception as e:
            print(f"Error loading data: {e}")
            # Create empty data store structure to prevent further errors
            data_store = {
                'hbv_data': pd.DataFrame(),
                'hcv_data': pd.DataFrame(),
                'hev_data': pd.DataFrame(),
                'ihme_df': pd.DataFrame(),
                'population_df': pd.DataFrame(),
                'coord_lookup': {},
                'hbv_mut': pd.DataFrame(),
                'hcv_mut': pd.DataFrame(),
                'hev_mut': pd.DataFrame()
            }
    return data_store

# === CONFIG & CONSTANTS ======================================================
HBV_GENOTYPE_COLORS = genotype_palette(
    "HBV",
    VIRUS_COLORS["HBV"],
    list("ABCDEFGHIJ"),
)

HCV_GENOTYPE_COLORS = genotype_palette(
    "HCV",
    VIRUS_COLORS["HCV"],
    [str(i) for i in range(1, 9)],
)

HEV_GENOTYPE_COLORS = genotype_palette(
    "HEV",
    VIRUS_COLORS["HEV"],
    [str(i) for i in range(1, 9)],
)

BURDEN_MEASURE_FALLBACK = "Prevalence|Number"  # used if dropdown missing
        
# === HELPERS & FIGURE BUILDERS ==============================================
def calculate_epidemiology_summary(ihme_df, virus, latest_year=None, regions=None, countries=None):
    """
    Calculate key epidemiology metrics for summary cards
    Returns dict with all summary metrics
    """
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause = cause_lookup.get((virus or "HBV").upper())
    
    if ihme_df.empty:
        return {}
    
    # Get data for this cause
    cause_data = ihme_df[ihme_df["cause"] == cause].copy()
    
    if cause_data.empty:
        return {}
    
    # Determine latest year
    if latest_year is None:
        latest_year = int(cause_data["year"].max())
    
    summary = {
        "latest_year": latest_year,
        "virus": virus
    }
    
    # Filter for latest year
    latest_data = cause_data[cause_data["year"] == latest_year].copy()
    
    if latest_data.empty:
        return summary
    
    # Apply region/country filters
    if regions:
        latest_data = latest_data[latest_data["WHO_Regions"].isin(regions)]
    if countries:
        latest_data = latest_data[latest_data["Country_standard"].isin(countries)]
    
    # Calculate prevalence (Number) - DON'T filter by "All ages"
    prevalence_data = latest_data[
        (latest_data["measure"] == "Prevalence") &
        (latest_data["metric"] == "Number")
    ].copy()
    
    if not prevalence_data.empty:
        # Sum across all age groups and sexes
        total_prevalence = prevalence_data["val"].sum()
        summary["prevalence_total"] = total_prevalence
        
        # Calculate trend (compare with 5 years ago)
        prev_year = latest_year - 5
        if prev_year >= 1980:  # Ensure reasonable year
            prev_data = cause_data[
                (cause_data["year"] == prev_year) &
                (cause_data["measure"] == "Prevalence") &
                (cause_data["metric"] == "Number")
            ]
            if not prev_data.empty:
                prev_prevalence = prev_data["val"].sum()
                if prev_prevalence > 0:
                    trend = ((total_prevalence - prev_prevalence) / prev_prevalence) * 100
                    summary["prevalence_trend"] = trend
    
    # Calculate incidence (Number)
    incidence_data = latest_data[
        (latest_data["measure"] == "Incidence") &
        (latest_data["metric"] == "Number")
    ].copy()
    
    if not incidence_data.empty:
        total_incidence = incidence_data["val"].sum()
        summary["incidence_total"] = total_incidence
        
        # Incidence trend
        prev_year = latest_year - 5
        if prev_year >= 1980:
            prev_incidence_data = cause_data[
                (cause_data["year"] == prev_year) &
                (cause_data["measure"] == "Incidence") &
                (cause_data["metric"] == "Number")
            ]
            if not prev_incidence_data.empty:
                prev_incidence = prev_incidence_data["val"].sum()
                if prev_incidence > 0:
                    trend = ((total_incidence - prev_incidence) / prev_incidence) * 100
                    summary["incidence_trend"] = trend
    
    # Calculate deaths (Number)
    deaths_data = latest_data[
        (latest_data["measure"] == "Deaths") &
        (latest_data["metric"] == "Number")
    ].copy()
    
    if not deaths_data.empty:
        total_deaths = deaths_data["val"].sum()
        summary["deaths_total"] = total_deaths
        
        # Deaths trend
        prev_year = latest_year - 5
        if prev_year >= 1980:
            prev_deaths_data = cause_data[
                (cause_data["year"] == prev_year) &
                (cause_data["measure"] == "Deaths") &
                (cause_data["metric"] == "Number")
            ]
            if not prev_deaths_data.empty:
                prev_deaths = prev_deaths_data["val"].sum()
                if prev_deaths > 0:
                    trend = ((total_deaths - prev_deaths) / prev_deaths) * 100
                    summary["deaths_trend"] = trend
    
    # Calculate sex ratio (Male:Female) from prevalence data
    if not prevalence_data.empty:
        # Group by sex
        sex_totals = prevalence_data.groupby("sex")["val"].sum()
        
        if "Male" in sex_totals.index and "Female" in sex_totals.index:
            male_val = sex_totals["Male"]
            female_val = sex_totals["Female"]
            if female_val > 0:
                sex_ratio = male_val / female_val
                summary["sex_ratio"] = sex_ratio
    
    # Find top age group from prevalence data
    if not prevalence_data.empty:
        # Group by age
        age_totals = prevalence_data.groupby("age")["val"].sum()
        
        if not age_totals.empty:
            top_age = age_totals.idxmax()
            top_age_value = age_totals.max()
            
            total_prevalence = summary.get("prevalence_total", 0)
            if total_prevalence > 0:
                age_percentage = (top_age_value / total_prevalence) * 100
                summary["top_age_group"] = top_age
                summary["top_age_percentage"] = age_percentage
    
    # Find top region from prevalence data
    if not prevalence_data.empty:
        # Group by region
        region_totals = prevalence_data.groupby("WHO_Regions")["val"].sum()
        
        if not region_totals.empty:
            top_region = region_totals.idxmax()
            top_region_value = region_totals.max()
            
            total_prevalence = summary.get("prevalence_total", 0)
            if total_prevalence > 0:
                region_percentage = (top_region_value / total_prevalence) * 100
                summary["top_region"] = top_region
                summary["top_region_percentage"] = region_percentage
    
    # Calculate WHO 2030 progress (using incidence)
    if "incidence_total" in summary:
        # Get 2015 baseline
        baseline_2015 = cause_data[
            (cause_data["year"] == 2015) &
            (cause_data["measure"] == "Incidence") &
            (cause_data["metric"] == "Number")
        ]
        if not baseline_2015.empty:
            baseline = baseline_2015["val"].sum()
            current = summary["incidence_total"]
            if baseline > 0:
                reduction = ((baseline - current) / baseline) * 100
                target_reduction = 90  # WHO 2030 target
                progress = (reduction / target_reduction) * 100
                summary["who_progress"] = min(100, max(0, progress))
                summary["reduction_needed"] = max(0, target_reduction - reduction)
    
    return summary

def format_large_number(num):
    """Format large numbers with K, M, B suffixes"""
    if pd.isna(num) or num == "N/A" or not isinstance(num, (int, float, np.number)):
        return "N/A"
    
    try:
        num = float(num)
    except (ValueError, TypeError):
        return "N/A"
    
    if num >= 1_000_000_000:
        return f"{num/1_000_000_000:.1f}B"
    elif num >= 1_000_000:
        return f"{num/1_000_000:.1f}M"
    elif num >= 1_000:
        return f"{num/1_000:.1f}K"
    elif num >= 0:
        return f"{num:,.0f}"
    else:
        return f"{num:,.0f}"

def format_trend(trend_value):
    """Format trend with arrow and color"""
    if pd.isna(trend_value) or trend_value == "N/A":
        return html.Span("No trend data", className="text-muted")
    
    try:
        trend_value = float(trend_value)
    except (ValueError, TypeError):
        return html.Span("No trend data", className="text-muted")
    
    if trend_value > 0:
        return html.Span([
            html.I(className="bi bi-arrow-up text-danger me-1"),
            f"+{abs(trend_value):.1f}%",
        ], className="text-danger")
    elif trend_value < 0:
        return html.Span([
            html.I(className="bi bi-arrow-down text-success me-1"),
            f"-{abs(trend_value):.1f}%",
        ], className="text-success")
    else:
        return html.Span([
            html.I(className="bi bi-dash text-secondary me-1"),
            "0.0%",
        ], className="text-secondary")

def compute_gap_df(
    virus: str,
    filtered_seq_df: pd.DataFrame,
    ihme_df: pd.DataFrame,
    selected_years: list | tuple,
    who_regions: list | None,
    countries: list | None,
    ihme_metric_choice: str | None,
    sex,
    target_per_10k: float = 5.0
) -> pd.DataFrame:
    # Example placeholder logic
    if filtered_seq_df is None or filtered_seq_df.empty:
        return pd.DataFrame(columns=["Country_standard", "observed_sequences", "burden", "expected_sequences", "coverage_gap", "coverage_ratio", "WHO_Regions"])

    # --- count observed sequences ---
    obs = (
        filtered_seq_df.groupby("Country_standard")
        .size()
        .reset_index(name="observed_sequences")
    )

    # --- get IHME burden for same filters ---
    try:
        measure, metric = (ihme_metric_choice or BURDEN_MEASURE_FALLBACK).split("|")
    except ValueError:
        measure, metric = "Prevalence", "Number"

    cause_lookup = {"HBV": "Total burden related to hepatitis B", "HCV": "Total burden related to hepatitis C", "HEV": "Total burden related to hepatitis E"}
    cause = cause_lookup.get(virus.upper(), "")

    # Filter by cause, measure, metric
    base_burden = ihme_df[
        (ihme_df["cause"] == cause)
        & (ihme_df["measure"] == measure)
        & (ihme_df["metric"] == metric)
    ].copy()

    # Filter by sex
    if sex == "Both":
        base_burden = base_burden[base_burden["sex"].isin(["Male", "Female"])]
    else:
        base_burden = base_burden[base_burden["sex"] == sex]

    # Filter by years if specified
    if selected_years:
        y0, y1 = selected_years
        base_burden = base_burden[(base_burden["year"] >= y0) & (base_burden["year"] <= y1)]

    # Sum over all ages and sexes per country/year
    if not base_burden.empty:
        # Group by country and year to get the yearly totals
        yearly_burden = base_burden.groupby(["Country_standard", "year"], as_index=False)["val"].sum()
        # Take the latest year in the range for each country
        latest_years = yearly_burden.groupby("Country_standard")["year"].transform("max")
        latest_burden = yearly_burden[yearly_burden["year"] == latest_years]
        # Rename column
        burden = latest_burden.rename(columns={"val": "burden"})[["Country_standard", "burden"]]
    else:
        burden = pd.DataFrame(columns=["Country_standard", "burden"])

    # --- merge and compute expected ---
    df = pd.merge(obs, burden, on="Country_standard", how="outer").fillna(0)
    
    # FIXED: Handle zero burden cases properly
    df["expected_sequences"] = np.where(
        df["burden"] > 0,
        (df["burden"] / 10000.0) * target_per_10k,
        0
    )
    
    # FIXED: Calculate coverage ratio safely
    df["coverage_ratio"] = np.where(
        df["expected_sequences"] > 0,
        df["observed_sequences"] / df["expected_sequences"],
        np.where(df["observed_sequences"] > 0, np.inf, 0)  # Infinite if we have sequences but no expected
    )
    
    df["coverage_gap"] = df["expected_sequences"] - df["observed_sequences"]
    df.loc[df["coverage_gap"] < 0, "coverage_gap"] = 0  # clip negative
    
    who_map = (
        ihme_df[["Country_standard", "WHO_Regions"]]
        .dropna(subset=["Country_standard"])
        .drop_duplicates("Country_standard")
    )
    
    df = df.merge(who_map, on="Country_standard", how="left")
    
    # CRITICAL FIX: Ensure Country_standard column exists even if no data
    if "Country_standard" not in df.columns:
        df["Country_standard"] = None
    
    return df
    
def ihme_latest_by_country(ihme_df, virus, measure_metric, sex, regions=None, countries=None, years=None):
    try:
        ihme_measure, ihme_metric = (measure_metric or "Prevalence|Number").split("|")
    except Exception:
        ihme_measure, ihme_metric = "Prevalence", "Number"

    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause_filter = cause_lookup.get((virus or "HBV").upper(), "")
    
    # Check if cause exists
    if ihme_df is None or ihme_df.empty or "cause" not in ihme_df.columns or cause_filter not in ihme_df["cause"].values:
        return pd.DataFrame(columns=["Country_standard", "Metric_raw", "Metric", "year"])

    # Handle "Both" sexes by getting Male and Female separately
    if sex == "Both":
        # Get Male data
        male_data = ihme_df[
            (ihme_df["sex"] == "Male") &
            (ihme_df["cause"] == cause_filter) &
            (ihme_df["measure"] == ihme_measure) &
            (ihme_df["metric"] == ihme_metric)
        ].copy()
        
        # Get Female data
        female_data = ihme_df[
            (ihme_df["sex"] == "Female") &
            (ihme_df["cause"] == cause_filter) &
            (ihme_df["measure"] == ihme_measure) &
            (ihme_df["metric"] == ihme_metric)
        ].copy()
        
        # Combine
        base = pd.concat([male_data, female_data])
    else:
        # Use the sex as-is
        base = ihme_df[
            (ihme_df["sex"] == sex) &
            (ihme_df["cause"] == cause_filter) &
            (ihme_df["measure"] == ihme_measure) &
            (ihme_df["metric"] == ihme_metric)
        ].copy()

    if base.empty:
        return pd.DataFrame(columns=["Country_standard", "Metric_raw", "Metric", "year"])

    # Apply region filter
    if regions:
        base = base[base["WHO_Regions"].isin(regions)]
    
    # Apply country filter
    if countries:
        base = base[base["Country_standard"].isin(countries)]

    # Limit to selected range then pick latest year within that range
    if years and len(years) == 2:
        y0, y1 = int(years[0]), int(years[1])
        base = base[(base["year"] >= y0) & (base["year"] <= y1)]
    
    if base.empty:
        return pd.DataFrame(columns=["Country_standard", "Metric_raw", "Metric", "year"])

    # Get latest year
    latest_year = int(base["year"].max())
    latest = base[base["year"] == latest_year].copy()

    # Sum across age groups for each country (or mean for rates/percents)
    if not latest.empty:
        if ihme_metric in ["Percent", "Rate"]:
            latest = latest.groupby(["Country_standard", "year"], as_index=False)["val"].mean()
        else:
            latest = latest.groupby(["Country_standard", "year"], as_index=False)["val"].sum()
    
    latest["Metric_raw"] = pd.to_numeric(latest["val"], errors="coerce")
    latest["Metric"] = latest["Metric_raw"].apply(lambda x: np.log10(x) if (np.isfinite(x) and x > 0) else np.nan)

    return latest[["Country_standard", "Metric_raw", "Metric", "year"]]

_ACCESSION_RE = re.compile(r"([A-Z0-9]+\.\d+)")

def _to_df(obj) -> pd.DataFrame:
    """Coerce list/dict/DF/None into a DataFrame copy."""
    if obj is None:
        return pd.DataFrame()
    if isinstance(obj, pd.DataFrame):
        return obj.copy()
    if isinstance(obj, list):
        return pd.DataFrame(obj)
    if isinstance(obj, dict):
        # If it's a single record dict, wrap it
        if all(not isinstance(v, (list, tuple)) for v in obj.values()):
            return pd.DataFrame([obj])
        return pd.DataFrame(obj)
    return pd.DataFrame(obj)

def _rename_flex(df: pd.DataFrame) -> pd.DataFrame:
    """Case/alias-insensitive renames into canonical names used in plots."""
    if df.empty:
        return df
    col_map = {}
    for c in df.columns:
        lc = c.lower().strip()
        if lc in {"country_standard", "country_std", "country name", "country"}:
            col_map[c] = "Country_standard"
        elif lc in {"who_regions", "who region", "region", "who"}:
            col_map[c] = "WHO_Regions"
        elif lc in {"genotype", "genotypes", "geno"}:
            col_map[c] = "genotype"
        elif lc in {"year", "yr"}:
            col_map[c] = "Year"
        elif lc in {"date"}:
            col_map[c] = "Date"
        elif lc in {"taxa"}:
            col_map[c] = "Taxa"
        elif lc in {"id", "accession"}:
            col_map[c] = "ID"
        elif lc in {"population"}:
            col_map[c] = "Population"
    return df.rename(columns=col_map)

def _ensure_year(df: pd.DataFrame) -> pd.DataFrame:
    """Create numeric Year from Year or Date, drop rows without it."""
    if "Year" not in df.columns and "Date" in df.columns:
        df["Year"] = pd.to_numeric(df["Date"], errors="coerce")
    if "Year" not in df.columns:
        df["Year"] = pd.NA
    df["Year"] = pd.to_numeric(df["Year"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["Year"]).copy()
    df["Year"] = df["Year"].astype(int)
    return df

def _normalize_seq_df(df: pd.DataFrame,
                      required=("Country_standard", "Year"),
                      fill_region=True) -> pd.DataFrame:
    """Normalize an input (list/dict/df) to the columns our plots expect."""
    df = _to_df(df)
    df = _rename_flex(df)
    df = _ensure_year(df)
    if fill_region and "WHO_Regions" not in df.columns:
        df["WHO_Regions"] = "Unknown"
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(
            f"normalize: missing {missing}; present columns: {list(df.columns)}"
        )
    return df

def _extract_id(series: pd.Series) -> pd.Series:
    """Extract accession-like ID (e.g., AB123456.1) from Taxa strings."""
    return series.astype(str).str.extract(_ACCESSION_RE, expand=False)

def _enrich_mutation_df(mutation_df, sequence_df):
    """
    Enrich mutation dataframe with sequence-level metadata.
    Requires a stable 'ID' column for joining.
    """
    if mutation_df is None or mutation_df.empty or sequence_df is None or sequence_df.empty:
        return pd.DataFrame()

    # --- Defensive copy (CRITICAL to avoid UnboundLocalError) ---
    mut = mutation_df.copy()

    # --- Ensure we have an ID column ---
    if "ID" not in mut.columns:
        if "sample" in mut.columns:
            mut["ID"] = mut["sample"]
        elif "Taxa" in mut.columns:
            mut["ID"] = mut["Taxa"]
        else:
            return pd.DataFrame()

    # --- Normalize ID just in case ---
    mut["ID"] = mut["ID"].astype(str).str.strip()

    # --- Ensure sequence_df has ID ---
    if "ID" not in sequence_df.columns:
        return pd.DataFrame()

    seq = sequence_df.copy()
    seq["ID"] = seq["ID"].astype(str).str.strip()

    # --- Enrich mutations with sequence metadata ---
    enriched = mut.merge(
        seq.drop_duplicates(subset=["ID"]),
        on="ID",
        how="left",
        suffixes=("", "_seq")
    )

    return enriched



def merge_population_nearest(counts_df: pd.DataFrame,
                             pop_df: pd.DataFrame,
                             tol_years: int = 3) -> pd.DataFrame:
    if counts_df is None or len(counts_df) == 0:
        return counts_df.assign(Population=np.nan)

    # Copy and select req cols on population
    left = counts_df.copy()
    pop     = pop_df[["Country_standard", "Year", "Population"]].copy()

    # Coerce types consistently
    left["Country_standard"] = left["Country_standard"].astype(str)
    pop["Country_standard"]     = pop["Country_standard"].astype(str)

    # Coerce Year to numeric, drop NaNs, cast BOTH to EXACT SAME dtype (int64)
    left["Year"] = pd.to_numeric(left["Year"], errors="coerce")
    pop["Year"]     = pd.to_numeric(pop["Year"],  errors="coerce")
    left = left.dropna(subset=["Year"])
    pop     = pop.dropna(subset=["Year"])
    left["Year"] = left["Year"].astype("int64")
    pop["Year"]     = pop["Year"].astype("int64")

    # Ensure Population numeric
    pop["Population"] = pd.to_numeric(pop["Population"], errors="coerce")

    # Merge per country so keys are sorted within each group
    out = []
    for ctry, lgrp in left.groupby("Country_standard", sort=False):
        rgrp = pop[pop["Country_standard"] == ctry]
        if rgrp.empty:
            out.append(lgrp.assign(Population=np.nan))
            continue

        lgrp = lgrp.sort_values("Year", kind="mergesort")
        rgrp = rgrp.sort_values("Year", kind="mergesort")

        merged = pd.merge_asof(
            lgrp,
            rgrp,
            on="Year",                      # SAME dtype on both sides (int64)
            tolerance=int(tol_years),      # tolerance in years
            direction="nearest",
            allow_exact_matches=True,
        )
        out.append(merged)

    return pd.concat(out, ignore_index=True)

def merge_population_nearest_two_pass(counts_df: pd.DataFrame,
                                      pop_df: pd.DataFrame,
                                      tol_years_first: int = 3,
                                      tol_years_wide: int = 50) -> pd.DataFrame:

    first = merge_population_nearest(counts_df, pop_df, tol_years=tol_years_first)

    # If Population missing anywhere, try a wider tolerance and fill only those rows
    if "Population" in first.columns and first["Population"].isna().any():
        widened = merge_population_nearest(counts_df, pop_df, tol_years=tol_years_wide)
        first["Population"] = first["Population"].fillna(widened["Population"])

    return first

def _country_pie_heading(virus: str, years_text: str, top_n: int, has_filters: bool) -> str:
    scope = "current selection" if has_filters else "all data"
    return f"{virus} sequences by country · {years_text} — Top {top_n} ({scope})"
    
def _fmt_list(values, max_items=3, *, empty_label="All"):
    vals = [v for v in (values or []) if v]
    if not vals:
        return empty_label
    if len(vals) <= max_items:
        return ", ".join(vals)
    return f"{', '.join(vals[:max_items])} +{len(vals)-max_items} more"

def _mutations_heading(virus, years_text, filters_text, has_filters):
    scope = "current selection" if has_filters else "all data"
    v = (virus or "HBV").upper()
    prefix = "HBV resistance mutations" if v == "HBV" else "HCV mutations"
    facet  = "drug" if v == "HBV" else "gene"
    return f"{prefix} by {facet} — {filters_text} · {years_text} ({scope})"
    
def _label_declutter(fig, xs, ys, labels, top_n=14, min_gap_frac=0.04, y_key=None):
    """Label only the most notable points on a scatter, greedily skipping any
    candidate whose x sits too close to an already-placed label. Alternates
    label direction (up/down) to further reduce collisions. Avoids Plotly's
    lack of built-in text-collision handling, which otherwise stacks labels
    unreadably in dense clusters."""
    if not xs:
        return
    x_range = max(xs) - min(xs) or 1
    min_gap = x_range * min_gap_frac
    order = sorted(range(len(xs)), key=lambda i: -(y_key(i) if y_key else ys[i]))
    placed_x = []
    for rank, i in enumerate(order[:top_n * 3]):  # scan a bit beyond top_n to allow skips
        if len(placed_x) >= top_n:
            break
        if any(abs(xs[i] - px) < min_gap for px in placed_x):
            continue
        placed_x.append(xs[i])
        fig.add_annotation(
            x=xs[i], y=ys[i], text=labels[i], showarrow=True,
            arrowhead=0, arrowsize=0.6, arrowwidth=1, arrowcolor="rgba(100,116,139,0.6)",
            ax=0, ay=-22 if len(placed_x) % 2 == 0 else 22,
            font=dict(size=9, color="#334155"),
            bgcolor="rgba(255,255,255,0.85)",
        )

def calculate_global_mutation_maximum():
    data = get_data_store()
    
    max_percentage = 0
    
    # Check HBV mutations
    if not data['hbv_mut'].empty and not data['hbv_data'].empty:
        hbv_mutations = _enrich_mutation_df(data['hbv_mut'], data['hbv_data'])
        if not hbv_mutations.empty:
            hbv_counts = hbv_mutations.groupby("mutation")["ID"].nunique()
            hbv_total = len(data['hbv_data'])
            hbv_max = (hbv_counts.max() / hbv_total * 100) if hbv_total > 0 else 0
            max_percentage = max(max_percentage, hbv_max)
    
    # Check HCV mutations  
    if not data['hcv_mut'].empty and not data['hcv_data'].empty:
        hcv_mutations = _enrich_mutation_df(data['hcv_mut'], data['hcv_data'])
        if not hcv_mutations.empty:
            hcv_counts = hcv_mutations.groupby("mutation")["ID"].nunique()
            hcv_total = len(data['hcv_data'])
            hcv_max = (hcv_counts.max() / hcv_total * 100) if hcv_total > 0 else 0
            max_percentage = max(max_percentage, hcv_max)

    # Check HCV mutations  
    if not data['hev_mut'].empty and not data['hev_data'].empty:
        hev_mutations = _enrich_mutation_df(data['hev_mut'], data['hev_data'])
        if not hev_mutations.empty:
            hev_counts = hev_mutations.groupby("mutation")["ID"].nunique()
            hev_total = len(data['hev_data'])
            hev_max = (hev_counts.max() / hev_total * 100) if hev_total > 0 else 0
            max_percentage = max(max_percentage, hev_max)
    
    # Round up to nearest 10 and add some padding
    global_max = min(100, np.ceil(max_percentage / 10) * 10 + 10) if max_percentage > 0 else 50
    return global_max


def _add_single_box_legend(fig: go.Figure, colors: list[str], labels: list[str], x: float = 0.01, y: float = 0.02):
    """Add discrete single box swatches without legend headings to bottom-left corner of world map."""
    for color, label in zip(colors, labels):
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode="markers",
            marker=dict(
                symbol="square",
                size=12,
                color=color,
                line=dict(width=1, color="rgba(15, 23, 42, 0.4)")
            ),
            name=str(label),
            showlegend=True
        ))
    
    fig.update_layout(
        showlegend=True,
        xaxis=dict(visible=False, showgrid=False, zeroline=False),
        yaxis=dict(visible=False, showgrid=False, zeroline=False),
        legend=dict(
            title=dict(text=""),  # No legend title / heading as requested
            orientation="v",       # Vertical stacked single boxes
            x=x, xanchor="left",
            y=y, yanchor="bottom",
            bgcolor="rgba(255, 255, 255, 0.92)",
            bordercolor="rgba(15, 23, 42, 0.15)",
            borderwidth=1,
            font=dict(family="IBM Plex Sans, sans-serif", size=10, color="#1e293b"),
            itemsizing="constant",
            traceorder="normal"
        )
    )


def _bin_continuous_series(z_series, is_pct=False, max_bins=5):
    """Bin continuous numerical data into discrete steps for single box legend rendering."""
    s = pd.Series(z_series).dropna()
    s_pos = s[s > 0]
    if s_pos.empty:
        return pd.Series(index=z_series.index, dtype=float), [], []
    
    num_bins = min(max_bins, max(1, len(s_pos.unique())))
    if num_bins == 1:
        val = s_pos.iloc[0]
        lbl = f"{val:.1f}%" if is_pct else f"{val:,.0f}"
        return pd.Series(0, index=s.index), [SEQUENCE_COUNT_COLORS[0]], [lbl]
    
    ranks = s_pos.rank(method="first")
    bins = pd.qcut(ranks, q=num_bins, labels=False, duplicates="drop")
    
    color_indices = [int(i * (len(SEQUENCE_COUNT_COLORS) - 1) / max(1, num_bins - 1)) for i in range(num_bins)]
    step_colors = [SEQUENCE_COUNT_COLORS[idx] for idx in color_indices]
    
    def _fmt_short(num):
        if is_pct:
            return f"{num:.1f}%"
        if num >= 1_000_000:
            val = num / 1_000_000
            return f"{val:.1f}M" if val % 1 != 0 else f"{int(val)}M"
        if num >= 1_000:
            val = num / 1_000
            return f"{val:.1f}k" if val % 1 != 0 else f"{int(val)}k"
        return f"{int(num):,}"

    tick_labels = []
    for b in range(num_bins):
        sub = s_pos[bins == b]
        if not sub.empty:
            low, high = sub.min(), sub.max()
            low_str, high_str = _fmt_short(low), _fmt_short(high)
            lbl = low_str if low_str == high_str else f"{low_str} – {high_str}"
            tick_labels.append(lbl)
        else:
            tick_labels.append(f"Bin {b+1}")
            
    return bins, step_colors, tick_labels


def create_world_map(
    country_data: pd.DataFrame, 
    country_genotype_counts: pd.DataFrame, 
    coord_lookup: dict[str, dict[str, float]], 
    virus_type: str = "HBV", 
    display_mode: str = "raw",  # "raw", "PerMillion", or "ihme"
    map_title: str = "",  # Add this to know what metric we're showing
    height: int = 480
) -> go.Figure:
    
    if virus_type == "HBV":
        genotype_colors = coverage_colorscale
    elif virus_type == "HCV":
        genotype_colors = coverage_colorscale
    else:
        genotype_colors = coverage_colorscale
    
    df = country_data.copy()
    
    # Ensure we have Metric_raw column
    if "Metric_raw" not in df.columns:
        # Try to find alternative columns
        if "val" in df.columns:
            df["Metric_raw"] = df["val"]
            print(f"Map debug - Using 'val' column for Metric_raw")
        elif "value" in df.columns:
            df["Metric_raw"] = df["value"]
            print(f"Map debug - Using 'value' column for Metric_raw")
        else:
            df["Metric_raw"] = np.nan
            print(f"Map debug - No value column found")
    
    valid = df.dropna(subset=["Country_standard"]).copy()
    
    # Ensure Metric_raw is numeric
    valid["Metric_raw"] = pd.to_numeric(valid["Metric_raw"], errors="coerce")

    fig = go.Figure()
    if valid.empty:
        return _empty_world("No country data available for current filters")

    # EPIDEMIOLOGY MODE (IHME data)
    if display_mode == "ihme":
        print(f"Map debug - IN EPIDEMIOLOGY MODE")
        
        # Check if we have any data
        if valid["Metric_raw"].isna().all():
            print(f"Map debug - All Metric_raw values are NaN")
            fig = _empty_world("No valid IHME data (all values missing)")
            return fig
        
        # Filter out non-positive values for log scaling
        valid_nonzero = valid[valid["Metric_raw"] > 0].copy()
        
        if valid_nonzero.empty:
            # Check why - all zeros or all negative?
            zero_count = (valid["Metric_raw"] == 0).sum()
            negative_count = (valid["Metric_raw"] < 0).sum()
            nan_count = valid["Metric_raw"].isna().sum()
            
            print(f"Map debug - Zero values: {zero_count}, Negative: {negative_count}, NaN: {nan_count}")
            
            fig = _empty_world(f"No positive IHME data available. Zero values: {zero_count}")
            fig.add_annotation(
                text="IHME data may have zeros for some countries/years", 
                x=0.5, y=0.4, showarrow=False
            )
            return fig
        
        # Apply logarithmic transformation safely
        valid_nonzero["log_value"] = np.log10(valid_nonzero["Metric_raw"])
        
        # Determine appropriate min/max for color scale
        if valid_nonzero["log_value"].notna().any():
            vmin = float(valid_nonzero["log_value"].min())
            vmax = float(valid_nonzero["log_value"].max())
            
            # Add some padding
            vmin = max(0, vmin - 0.5)  # Don't go below 0
            vmax = vmax + 0.5
            
            print(f"Map debug - Log range: {vmin} to {vmax}")
        else:
            vmin, vmax = 3.0, 8.0  # Default range
            print(f"Map debug - Using default log range: {vmin} to {vmax}")
        
        # Use appropriate color scale based on virus type
        
        if virus_type == "HBV":
            colorscale = SEQUENCE_COUNT_STEPPED_COLORSCALE
        elif virus_type == "HCV":
            colorscale = SEQUENCE_COUNT_STEPPED_COLORSCALE
        else:
            colorscale = SEQUENCE_COUNT_STEPPED_COLORSCALE
        
        # Determine colorbar title from map_title
        if "Prevalence" in map_title:
            colorbar_title = "Prevalence (Log10)"
        elif "Incidence" in map_title:
            colorbar_title = "Incidence (Log10)"
        elif "Deaths" in map_title:
            colorbar_title = "Deaths (Log10)"
        else:
            colorbar_title = "Value (Log10)"
        
        # Bin continuous values for single-box legend
        bins, step_colors, tick_labels = _bin_continuous_series(
            valid_nonzero["Metric_raw"],
            is_pct=("rate" in map_title.lower() or "prevalence" in map_title.lower())
        )
        
        fig.add_trace(
            go.Choropleth(
                locations=valid_nonzero["Country_standard"],
                locationmode="country names",
                z=valid_nonzero["log_value"],
                zmin=vmin,
                zmax=vmax,
                colorscale=SEQUENCE_COUNT_STEPPED_COLORSCALE,
                showscale=False,
                marker_line_color="rgba(0,0,0,0.3)",
                marker_line_width=0.5,
                hovertemplate=(
                    "<b>%{location}</b><br>" +
                    "Value: %{customdata:,.0f}<br>" +
                    "Log10: %{z:.2f}<extra></extra>"
                ),
                customdata=valid_nonzero["Metric_raw"].values,
            )
        )
        _add_single_box_legend(fig, step_colors, tick_labels)
    
    # PER MILLION MODE
    elif display_mode == "PerMillion":
        bins_pm = [0, 0.01, 0.1, 1, 10, 100, 1000, 10000, float("inf")]
        labels_pm = ["0", "<0.1", "0.1–1", "1–10", "10–100", "100–1K", "1K–10K", "10K+"]
        valid["bin_pm"] = pd.cut(valid["Metric_raw"], bins=bins_pm, labels=labels_pm, include_lowest=True, right=False)
        bin_to_idx_pm = {lab: i for i, lab in enumerate(labels_pm)}
        valid["z_value_pm"] = valid["bin_pm"].map(bin_to_idx_pm)

        nonzero = valid[valid["Metric_raw"] > 0].copy()
        driving = nonzero if not nonzero.empty else valid
        z_numeric = driving["z_value_pm"].astype(float).to_numpy()
        locations = driving["Country_standard"]

        fig.add_trace(
            go.Choropleth(
                locations=locations,
                locationmode="country names",
                z=z_numeric,
                zmin=0,
                zmax=len(labels_pm) - 1,
                colorscale=SEQUENCE_COUNT_STEPPED_COLORSCALE,
                showscale=False,
                marker_line_color="rgba(0,0,0,0.3)",
                marker_line_width=0.5,
                hovertext=driving.apply(
                    lambda r: f"<b>{r['Country_standard']}</b><br>Per million: {float(r['Metric_raw']):,.2f}<br>Range: {r['bin_pm']}",
                    axis=1,
                ),
                hoverinfo="text",
            )
        )
        _add_single_box_legend(fig, SEQUENCE_COUNT_COLORS, labels_pm)
    
    # RAW COUNTS MODE (default)
    else:
        bins = [0, 1, 5, 20, 100, 500, 2000, 4000, float("inf")]
        labels = ["0", "1–4", "5–19", "20–99", "100–499", "500–1,999", "2,000–3,999", "4,000+"]
        valid["bin"] = pd.cut(valid["Metric_raw"], bins=bins, labels=labels, include_lowest=True, right=False)
        bin_to_idx = {lab: i for i, lab in enumerate(labels)}
        valid["z_value"] = valid["bin"].map(bin_to_idx)

        nonzero = valid[valid["Metric_raw"] > 0].copy()
        driving = nonzero if not nonzero.empty else valid
        z_numeric = driving["z_value"].astype(float).to_numpy()
        locations = driving["Country_standard"]

        colorscale = SEQUENCE_COUNT_STEPPED_COLORSCALE

        fig.add_trace(
            go.Choropleth(
                locations=locations,
                locationmode="country names",
                z=z_numeric,
                zmin=0,
                zmax=len(labels) - 1,
                colorscale=colorscale,
                showscale=False,
                marker_line_color="rgba(0,0,0,0.3)",
                marker_line_width=0.5,
                hovertext=driving.apply(
                    lambda r: f"<b>{r['Country_standard']}</b><br>Exact count: {float(r['Metric_raw']):.0f}<br>Range: {r['bin']}",
                    axis=1,
                ),
                hoverinfo="text",
            )
        )
        _add_single_box_legend(fig, SEQUENCE_COUNT_COLORS, labels)

    # genotype overlay markers (only show for sequence data, not IHME data)
    if display_mode != "ihme":
        lons, lats, texts, sizes = [], [], [], []
        for country in df["Country_standard"].dropna().unique():
            subset = country_genotype_counts[country_genotype_counts["Country_standard"] == country]
            subset = subset[subset["Count"] > 0]
            total = int(subset["Count"].sum()) if not subset.empty else 0
            coords = coord_lookup.get(country)
            if total <= 0 or not coords:
                continue
            lat = float(coords["latitude"])
            lon = float(coords["longitude"])
            genotype_text = "<br>".join(
                f"{row.genotype}: {row.Count} ({row.Count/total:.1%})" for _, row in subset.sort_values("Count", ascending=False).iterrows()
            )
            texts.append(f"<b>{country}</b><br>Total: {total}<br>{genotype_text}")
            lats.append(lat)
            lons.append(lon)
            sizes.append(10 + min(20, total ** 0.2))

        if lons:
            fig.add_trace(
                go.Scattergeo(
                    lon=lons,
                    lat=lats,
                    text=texts,
                    hoverinfo="text",
                    mode="markers",
                    marker=dict(size=sizes, color="lightgrey", opacity=0.7, line=dict(width=1.5, color="black")),
                    showlegend=False,
                )
            )

    # Apply to all modes after traces so fitbounds works
    fig.update_geos(
        projection_type="equirectangular",
        showcountries=True,
        countrycolor="rgba(0,0,0,0.2)",
        showsubunits=True,
        fitbounds=False,
        showframe=False,
        lataxis_range=[-65, 85],
        domain=dict(x=[0, 1], y=[0, 1]),
    )
    fig.update_layout(
        height=height,
        margin=dict(t=0, b=0, l=0, r=0),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(visible=False, showgrid=False, zeroline=False),
        yaxis=dict(visible=False, showgrid=False, zeroline=False),
    )
    return apply_heptracker_figure_style(fig)
    
def create_coverage_map(
    cov_df: pd.DataFrame,
    coord_lookup: dict[str, any],
    coords_df: pd.DataFrame | None = None,
    virus_type: str = "HBV",
    who_regions: list[str] | None = None,
    countries: list[str] | None = None,
    height: int = 480
) -> go.Figure:
    
    if cov_df is None or cov_df.empty:
        print("Warning: cov_df is empty")
        return apply_heptracker_figure_style(_empty_world("No data available"))
        return _empty_world("No coverage data available")
    
    valid = cov_df.copy()
    
    # AGGREGATE BY COUNTRY
    if "Country_standard" in valid.columns and "Coverage_ratio" in valid.columns:
        country_data = valid.groupby("Country_standard", as_index=False).agg({
            "Coverage_ratio": "mean",
            "Seq_count": "sum",
            "Est_infections_genotype": "sum"
        })
        
        country_data = country_data.rename(columns={
            "Coverage_ratio": "coverage_ratio",
            "Seq_count": "observed_sequences",
            "Est_infections_genotype": "expected_sequences"
        })
        
        country_data["coverage_gap"] = (country_data["expected_sequences"] - country_data["observed_sequences"]).clip(lower=0)
        
        valid = country_data.copy()
        
        print(f"📊 Coverage stats - Min: {valid['coverage_ratio'].min():.6%}, "
              f"Max: {valid['coverage_ratio'].max():.6%}, "
              f"Median: {valid['coverage_ratio'].median():.6%}")
    
    if "Country_standard" not in valid.columns:
        return _empty_world("No country information available in coverage data")
    
    valid["Country_standard"] = valid["Country_standard"].astype(str).str.strip()
    
    # Apply filters
    if who_regions and coords_df is not None and "WHO_Regions" in coords_df.columns:
        region_map = coords_df[["Country_standard", "WHO_Regions"]].drop_duplicates()
        valid = valid.merge(region_map, on="Country_standard", how="left")
        valid = valid[valid["WHO_Regions"].isin(who_regions)]
    
    if countries:
        valid = valid[valid["Country_standard"].isin(countries)]
    
    # Filter valid positive numbers and drop NaNs
    plot_data = valid.dropna(subset=["coverage_ratio"]).copy()
    plot_data = plot_data[plot_data["coverage_ratio"] > 0].copy()
    
    if plot_data.empty:
        return _empty_world("No coverage ratio data available for selected filters")
    
    # Calculate quantile bins using rank-based qcut for 100% fail-safe execution
    num_countries = len(plot_data)
    num_bins = min(8, max(1, num_countries))
    
    ranks = plot_data["coverage_ratio"].rank(method="first")
    plot_data["coverage_bin"] = pd.qcut(ranks, q=num_bins, labels=False, duplicates="drop")
    
    max_bin = int(plot_data["coverage_bin"].max()) if not plot_data["coverage_bin"].isna().all() else 0
    plot_data["z_value"] = plot_data["coverage_bin"].astype(float)
    
    # Build tick labels dynamically from bin min/max values
    tick_positions = list(range(max_bin + 1))
    tick_labels = []
    
    for b in tick_positions:
        bin_subset = plot_data[plot_data["coverage_bin"] == b]["coverage_ratio"]
        if not bin_subset.empty:
            low = bin_subset.min() * 100
            high = bin_subset.max() * 100
            if low == high:
                tick_labels.append(f"{low:.1f}%")
            else:
                tick_labels.append(f"{low:.1f}%–{high:.1f}%")
        else:
            tick_labels.append(f"Bin {b+1}")
    
    color_indices = [int(i * (len(SEQUENCE_COUNT_COLORS) - 1) / max(1, len(tick_labels) - 1)) for i in range(len(tick_labels))]
    step_colors = [SEQUENCE_COUNT_COLORS[idx] for idx in color_indices]

    fig = go.Figure()
    
    fig.add_trace(go.Choropleth(
        locations=plot_data["Country_standard"],
        locationmode="country names",
        z=plot_data["z_value"],
        zmin=0,
        zmax=len(tick_labels) - 1,
        colorscale=SEQUENCE_COUNT_STEPPED_COLORSCALE,
        showscale=False,
        marker_line_color="rgba(0,0,0,0.3)",
        marker_line_width=0.5,
        hovertemplate=(
            "<b>%{location}</b><br>" +
            "<br>" +
            "📊 Coverage Ratio: <b>%{customdata[0]:.4%}</b><br>" +
            "🔬 Observed Sequences: <b>%{customdata[1]:,.0f}</b><br>" +
            "📈 Expected Sequences: <b>%{customdata[2]:,.0f}</b><br>" +
            "⚠️ Sequencing Gap: <b>%{customdata[3]:,.0f}</b><br>" +
            "📊 Relative Rank: <b>%{customdata[4]}</b><br>" +
            "<extra></extra>"
        ),
        customdata=np.stack([
            plot_data["coverage_ratio"].values,
            plot_data["observed_sequences"].values,
            plot_data["expected_sequences"].values,
            plot_data["coverage_gap"].values,
            plot_data["coverage_bin"].apply(lambda x: f"Bin {int(x)+1}/{int(max_bin)+1}" if pd.notna(x) else "N/A").values
        ], axis=1)
    ))
    
    _add_single_box_legend(fig, step_colors, tick_labels)
    
    fig.update_geos(
        projection_type="equirectangular",
        showcountries=True,
        countrycolor="rgba(0,0,0,0.2)",
        showsubunits=True,
        fitbounds=False,
        showframe=False,
        lataxis_range=[-65, 85],
        domain=dict(x=[0, 1], y=[0, 1]),
    )
    fig.update_layout(
        height=height,
        margin=dict(t=0, b=0, l=0, r=0),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(visible=False, showgrid=False, zeroline=False),
        yaxis=dict(visible=False, showgrid=False, zeroline=False),
    )
    
    return apply_heptracker_figure_style(fig)
    
def get_clean_log_axis(ymin, ymax, is_per_million=False):
    """
    Generate clean range, tick values, and tick labels for Plotly log axes.

    Important:
    - yaxis.range must be log10 values.
    - yaxis.tickvals must be real/raw data values, not log10 values.
    """
    min_allowed = 0.01 if is_per_million else 1.0

    ymin = max(min_allowed, float(ymin))
    ymax = max(ymin * 1.1, float(ymax))

    lo_val = max(min_allowed, ymin / 1.5)
    hi_val = ymax * 1.5

    lo = np.log10(lo_val)
    hi = np.log10(hi_val)

    # Guarantee enough visual range
    if (hi - lo) < 1.5:
        mid = (lo + hi) / 2.0
        lo = mid - 0.75
        hi = mid + 0.75
        lo_val = 10.0 ** lo
        hi_val = 10.0 ** hi

    span_ratio = hi_val / lo_val

    start_dec = int(np.floor(np.log10(lo_val)))
    end_dec = int(np.ceil(np.log10(hi_val)))

    if span_ratio > 10000:
        multipliers = [1]
    elif span_ratio > 100:
        multipliers = [1, 3]
    else:
        multipliers = [1, 2, 5]

    tickvals = []
    ticktext = []

    for dec in range(start_dec, end_dec + 1):
        for mult in multipliers:
            val = (10.0 ** dec) * mult

            if lo_val <= val <= hi_val:
                # CRITICAL FIX:
                # Plotly log axis tickvals must be raw values, not np.log10(val)
                tickvals.append(val)

                if val >= 1e6:
                    ticktext.append(f"{val / 1e6:g}M")
                elif val >= 1e3:
                    ticktext.append(f"{val / 1e3:g}k")
                elif val >= 1:
                    ticktext.append(f"{val:g}")
                else:
                    ticktext.append(f"{val:g}")

    if not tickvals:
        tickvals = [ymin]
        ticktext = [f"{ymin:g}"]

    return {
        "range": [float(lo), float(hi)],
        "tickvals": tickvals,
        "ticktext": ticktext,
    }

def make_line_trend(
    filtered_df: pd.DataFrame,
    selected_virus: str,
) -> go.Figure:

    if selected_virus == "HBV":
        genotype_colors = HBV_GENOTYPE_COLORS
    elif selected_virus == "HCV":
        genotype_colors = HCV_GENOTYPE_COLORS
    else:
        genotype_colors = HEV_GENOTYPE_COLORS

    if filtered_df is None or filtered_df.empty:
        return _empty_plot("No sequence data available")

    line_data = (
        filtered_df.groupby(["Year", "genotype"])
        .size()
        .reset_index(name="Genome Sequences")
    )

    if line_data.empty:
        return _empty_plot("No sequence trend data available")

    smoothed_data = []

    for genotype in line_data["genotype"].dropna().unique():
        genotype_df = line_data[line_data["genotype"] == genotype].copy()
        genotype_df = genotype_df.sort_values("Year")

        genotype_df["Smoothed_Sequences"] = (
            genotype_df["Genome Sequences"]
            .rolling(window=3, min_periods=1, center=True)
            .mean()
        )

        smoothed_data.append(genotype_df)

    if not smoothed_data:
        return _empty_plot("No genotype trend data available")

    smoothed_df = pd.concat(smoothed_data, ignore_index=True)

    smoothed_df["Year"] = pd.to_numeric(smoothed_df["Year"], errors="coerce")
    smoothed_df["Smoothed_Sequences"] = pd.to_numeric(
        smoothed_df["Smoothed_Sequences"],
        errors="coerce"
    )

    smoothed_df = smoothed_df.dropna(subset=["Year", "Smoothed_Sequences"])
    smoothed_df = smoothed_df[smoothed_df["Smoothed_Sequences"] > 0]

    if smoothed_df.empty:
        return _empty_plot("No positive values available for log-scale trend")

    fig = px.line(
        smoothed_df,
        x="Year",
        y="Smoothed_Sequences",
        color="genotype",
        color_discrete_map=genotype_colors,
        markers=True,
        line_shape="linear",
    )

    fig.update_traces(
        mode="lines+markers",
        marker=dict(size=5),
        line=dict(width=2.5),
    )

    y_pos = smoothed_df["Smoothed_Sequences"].dropna()

    if len(y_pos) and y_pos.max() > 0:
        y_min = float(y_pos[y_pos > 0].min())
        y_max = float(y_pos.max())
        axis_config = get_clean_log_axis(y_min, y_max, is_per_million=False)
        lo, hi = axis_config["range"]
    else:
        axis_config = {
            "range": [0, 2],
            "tickvals": [1, 10, 100],
            "ticktext": ["1", "10", "100"],
        }
        lo, hi = 0, 2

    year_min = int(smoothed_df["Year"].min())
    year_max = int(smoothed_df["Year"].max())

    xaxis_dict = dict(
        title="Year",
        fixedrange=True,
        showticklabels=True,
        automargin=True,
        gridcolor="rgba(255,255,255,0.08)",
        linecolor="rgba(255,255,255,0.25)",
    )

    if year_min <= year_max:
        tickvals = list(range(year_min, year_max + 1, 4))
        xaxis_dict["tickmode"] = "array"
        xaxis_dict["tickvals"] = tickvals
        xaxis_dict["ticktext"] = [str(y) for y in tickvals]

    fig.update_layout(
        xaxis=xaxis_dict,
        yaxis=dict(
            title=dict(
                text="Number of Sequences (log scale)",
                standoff=24,
            ),
            type="log",
            autorange=False,
            range=axis_config["range"],
            tickmode="array",
            tickvals=axis_config["tickvals"],
            ticktext=axis_config["ticktext"],
            gridcolor="rgba(255,255,255,0.08)",
            linecolor="rgba(255,255,255,0.25)",
            zeroline=False,
            fixedrange=True,
            showticklabels=True,
            automargin=True,
        ),
        height=420,
        margin=dict(t=90, b=50, l=95, r=25),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="center",
            x=0.5,
        ),
        hovermode="closest",
        meta={"y_full_log_range": [lo, hi]},
    )

    return fig

def make_genotype_bar(
    filtered_df: pd.DataFrame,
    population_df: pd.DataFrame,
    selected_virus: str,
    display_mode: str,
) -> go.Figure:
    selected_virus = selected_virus or "HBV"
    
    if selected_virus.upper() == "HBV":
        genotype_colors = HBV_GENOTYPE_COLORS
    elif selected_virus.upper() == "HCV":
        genotype_colors = HCV_GENOTYPE_COLORS
    else:  # HEV
        genotype_colors = HEV_GENOTYPE_COLORS

    # 1) Aggregate counts at (Country, Year, genotype)
    if "Count" in filtered_df.columns:
        cyg = (
            filtered_df.groupby(["Country_standard", "Year", "genotype"], as_index=False)["Count"]
                       .sum()
        )
    elif "sequences" in filtered_df.columns:
        cyg = (
            filtered_df.groupby(["Country_standard", "Year", "genotype"], as_index=False)["sequences"]
                       .sum().rename(columns={"sequences": "Count"})
        )
    else:
        cyg = (
            filtered_df.groupby(["Country_standard", "Year", "genotype"], as_index=False)
                       .size().rename(columns={"size": "Count"})
        )

    # 2) Merge population (same helper you used before)
    cyg = merge_population_nearest_two_pass(
        cyg, population_df, tol_years_first=3, tol_years_wide=50
    )

    # 3) Aggregate to genotype totals and population denominators
    agg = (
        cyg.groupby("genotype", as_index=False)
           .agg(Total=("Count", "sum"), Pop=("Population", "sum"))
    )

    # Normalise genotype labels for consistent ordering and formatting
    def _norm(g):
        g = str(g).strip()
        if selected_virus.upper() == "HBV":
            if len(g) == 1 and g.isalpha():
                return f"HBV-{g.upper()}"
            elif g.startswith("HBV-"):
                return g
            elif g == "Recombinant":
                return "Recombinant"  # Keep as is for color mapping
            else:
                return f"HBV-{g}"
        elif selected_virus.upper() == "HCV":
            if g.isdigit() or (g.replace('.', '').isdigit() and g.count('.') <= 1):
                return f"HCV-{g}"
            elif g.startswith("HCV-"):
                return g
            elif g == "Recombinant":
                return "Recombinant"  # Keep as is for color mapping
            else:
                return f"HCV-{g}"
        else: # HEV
            if g.isdigit() or (g.replace('.', '').isdigit() and g.count('.') <= 1):
                return f"HEV-{g}"
            elif g.startswith("HEV-"):
                return g
            elif g == "Recombinant":
                return "Recombinant"
            else:
                return f"HCV-{g}"

    agg["genotype"] = agg["genotype"].apply(_norm)

    # 4) Compute per-million when requested
    agg["PerMillion"] = np.where(
        (agg["Pop"].notna()) & (agg["Pop"] > 0),
        (agg["Total"] / agg["Pop"]) * 1_000_000.0,
        np.nan
    )
    y_col = "PerMillion" if (display_mode or "raw") == "PerMillion" else "Total"
    y_title = "Sequences per Million" if y_col == "PerMillion" else "Number of Sequences"
    
    # Add this after calculating vals
    print(f"DEBUG - Virus: {selected_virus}, Max value: {max_val if 'max_val' in locals() else 'N/A'}, y_col: {y_col}")
    print(f"DEBUG - Agg values: {agg[y_col].tolist()}")

    # 5) Ensure stable order
    if selected_virus.upper() == "HBV":
        target_order = [f"HBV-{c}" for c in list("ABCDEFGHIJ")] + ["Recombinant"]
    elif selected_virus.upper() == "HCV":
        # HCV genotypes: 1-7 plus recombinant
        target_order = [f"HCV-{i}" for i in range(1, 9)] + ["Recombinant"]
    else:
        target_order = [f"HEV-{i}" for i in range(1, 9)] + ["Recombinant"]

    # Add missing categories as zero so the axis is complete
    base = pd.DataFrame({"genotype": target_order})
    agg = base.merge(agg, on="genotype", how="left").fillna({y_col: 0, "Total": 0})

    # Remove any genotypes that don't exist in our data (all zeros)
    agg = agg[agg[y_col] > 0] if len(agg[agg[y_col] > 0]) > 0 else agg

    # 6) Prepare log range (avoid log(0))
    vals = pd.to_numeric(agg[y_col], errors="coerce").dropna()
    if len(vals) == 0 or vals.sum() == 0:
        # Empty data case
        ymin, ymax = 0.1, 10
    else:
        # CRITICAL FIX: Force ymin to 1.0 (or 0.1 for PerMillion) to prevent the Plotly log-bar SVG clipping/cut-off bug
        ymin = 0.1 if y_col == "PerMillion" else 1.0
        ymax = float(vals.max()) if len(vals) else 10
    
    low_pad = max(0.1, ymin / 1.2)
    high_pad = max(1.0, ymax * 1.3)
    log_range = [np.log10(low_pad), np.log10(high_pad)]

    # 7) Create custom color mapping for the actual genotypes present
    present_genotypes = agg["genotype"].tolist()
    custom_color_map = {}
    
    for genotype in present_genotypes:
        # Look for exact match first
        if genotype in genotype_colors:
            custom_color_map[genotype] = genotype_colors[genotype]
        # For recombinant - check if it's just "Recombinant" and use the color from dictionary
        elif genotype == "Recombinant" and "Recombinant" in genotype_colors:
            custom_color_map[genotype] = genotype_colors["Recombinant"]
        # For HBV genotypes without prefix
        elif selected_virus.upper() == "HBV" and len(genotype) == 1 and genotype.isalpha():
            hbv_key = f"HBV-{genotype}"
            if hbv_key in genotype_colors:
                custom_color_map[genotype] = genotype_colors[hbv_key]
        # For HCV genotypes without prefix
        elif selected_virus.upper() == "HCV" and genotype.isdigit():
            hcv_key = f"HCV-{genotype}"
            if hcv_key in genotype_colors:
                custom_color_map[genotype] = genotype_colors[hcv_key]
        # for HEV genotypes without prefix
        elif selected_virus.upper() == "HEV" and genotype.isdigit():
            hev_key = f"HEV-{genotype}"
            if hev_key in genotype_colors:
                custom_color_map[genotype] = genotype_colors[hev_key]
        else:
            # Fallback color if not found
            custom_color_map[genotype] = "#CCCCCC"

    # 8) Plot with custom color mapping
    bar_fig = px.bar(
        agg,
        x="genotype",
        y=y_col,
        color="genotype",
        category_orders={"genotype": [g for g in target_order if g in agg["genotype"].values]},
        color_discrete_map=custom_color_map,  # Use the custom mapping
    )

    # Hover text
    unit = " per million" if y_col == "PerMillion" else ""
    bar_fig.update_traces(
        marker_cornerradius=6,
        hovertemplate="<b>%{x}</b><br>%{y:,}" + unit + " sequences<extra></extra>"
    )

    # Layout & legend
    axis_config = get_clean_log_axis(ymin, ymax, is_per_million=(y_col == "PerMillion"))
    bar_fig.update_layout(
        #title=f"Total {selected_virus.upper()} Sequences by genotype",
        xaxis_title="genotype",
        xaxis=dict(fixedrange=True),
        yaxis=dict(
            title=y_title + " (log scale)",
            type="log",
            autorange=False,
            range=axis_config["range"],
            tickmode="array",
            tickvals=axis_config["tickvals"],
            ticktext=axis_config["ticktext"],
            gridcolor="rgba(0,0,0,0.08)",
            zeroline=False,
            linecolor="rgba(0,0,0,0.25)",
            fixedrange=True,
            constrain="domain",
        ),
        bargap=0.25,
        
        
        height=420,
        margin=dict(t=70, b=40, l=40, r=20),
        legend=dict(orientation="h", y=1.02, x=0.5, xanchor="center", yanchor="bottom"),
    )

    return bar_fig
    
def make_country_pie(df, selected_regions=None, selected_countries=None, selected_years=None,
                     virus="HBV", top_n=10):
    df = _normalize_seq_df(df, required=("Country_standard", "Year"))

    # Apply region/country filters first
    if selected_regions and "WHO_Regions" in df.columns:
        df = df[df["WHO_Regions"].isin(selected_regions)]
    if selected_countries:
        df = df[df["Country_standard"].isin(selected_countries)]

    # Compute the full year span *after* region/country filters
    year_series = pd.to_numeric(df["Year"], errors="coerce")
    data_ymin = int(year_series.min()) if year_series.notna().any() else None
    data_ymax = int(year_series.max()) if year_series.notna().any() else None

    # Apply year filter (if any) and build the label
    if selected_years and len(selected_years) == 2 and all(v is not None for v in selected_years):
        y0, y1 = int(selected_years[0]), int(selected_years[1])
        years_text = f"{y0}–{y1}"
        df = df[(year_series >= y0) & (year_series <= y1)]
        # A years filter counts only if it narrows the full span
        years_filter_active = (data_ymin is not None and data_ymax is not None) and (y0 > data_ymin or y1 < data_ymax)
    else:
        years_text = "All years"
        years_filter_active = False

    # Determine if *any* filters are active
    has_filters = bool(
        (selected_regions and len(selected_regions) > 0) or
        (selected_countries and len(selected_countries) > 0) or
        years_filter_active
    )

    if df.empty:
        fig = px.pie(pd.DataFrame({"Country": ["No data"], "Count": [1]}),
                     names="Country", values="Count", hole=0.6)
        fig.update_traces(textinfo="none", hoverinfo="skip", showlegend=False)
        fig.update_layout(margin=dict(l=10, r=10, t=40, b=40))
        title = _country_pie_heading(virus, years_text, top_n, has_filters)
        return fig, title

    # --- counts & figure (unchanged from your working version) ---
    vc = (df["Country_standard"].fillna("Unknown")
          .value_counts(dropna=False).rename_axis("Country").reset_index(name="Count")
          .sort_values("Count", ascending=False))
    top = vc.head(top_n).copy()
    other = vc.iloc[top_n:]
    other_countries, other_count = int(other.shape[0]), int(other["Count"].sum())
    total = int(vc["Count"].sum())
    top["SharePct"] = (top["Count"] / total) * 100

    fig = px.pie(top, names="Country", values="Count", hole=0.55,
                 category_orders={"Country": list(top["Country"])})
    fig.update_traces(
        text=[f"{p:.1f}%" for p in top["SharePct"]],
        texttemplate="%{text}", textinfo="text", textposition="inside", sort=False,
        marker=dict(line=dict(color="#fff", width=1)),
        hovertemplate="<b>%{label}</b><br>Sequences: %{value:,}<br>"
                      "Share of total: %{customdata:.1f}%<extra></extra>",
        customdata=top["SharePct"], showlegend=True,
    )
    fig.update_layout(
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="center", x=0.5, bgcolor="rgba(0,0,0,0)",
                    font=dict(size=11), itemwidth=30, title_text=""),
        margin=dict(l=20, r=20, t=90, b=80), autosize=True
    )
    footnote = (f"+ {other_countries} more countries ({other_count:,} sequences not shown)"
                if other_count > 0 and other_countries > 0
                else "All countries shown (complete data)")
    fig.add_annotation(text=footnote, x=0.5, y=-0.12, xref="paper", yref="paper",
                       showarrow=False, font=dict(size=14), align="center")

    title = _country_pie_heading(virus, years_text, min(top_n, len(vc)), has_filters)
    return fig, title

##### HEV just has genotype data from glue
def make_mutation_bar(
    mutation_df: pd.DataFrame,
    total_sequences: int,
    selected_virus: str,
    selected_filter: list | None = None,
    years_range: list | tuple | None = None,
    data_span: list | tuple | None = None,
    other_filters_active: bool = False,
) -> tuple[go.Figure, str]:
    import numpy as np
    v = (selected_virus or "HBV").upper()
    is_hbv = (v == "HBV")
    color = "#3182bd" if is_hbv else "#d94801"
    col_name = "drug" if is_hbv else "gene"

    # years text + active
    if years_range and len(years_range) == 2 and None not in years_range:
        y0, y1 = int(years_range[0]), int(years_range[1])
        years_text = f"{y0}–{y1}"
        if data_span and len(data_span) == 2 and None not in data_span:
            ymin, ymax = int(data_span[0]), int(data_span[1])
            years_active = (y0 > ymin) or (y1 < ymax)
        else:
            years_active = True
    else:
        years_text = "All years"
        years_active = False

    # filter by drug/gene
    filters_text = "All drugs" if is_hbv else "All genes"
    if selected_filter:
        mutation_df = mutation_df[mutation_df[col_name].notna()]
        selected_filter_lower = [str(s).strip().lower() for s in selected_filter]
        mutation_df = mutation_df[
            mutation_df[col_name].astype(str).str.strip().str.lower().isin(selected_filter_lower)
        ]
        filters_text = _fmt_list(
            selected_filter, max_items=3,
            empty_label=("All drugs" if is_hbv else "All genes")
        )

    has_filters = bool(other_filters_active or years_active or (selected_filter and len(selected_filter) > 0))

    fig = go.Figure()
    if mutation_df.empty or not total_sequences:
        fig.update_layout(
            xaxis={"visible": False}, yaxis={"visible": False},
            annotations=[{"text": "No mutations found for current selection",
                          "xref": "paper", "yref": "paper", "x": 0.5, "y": 0.5,
                          "showarrow": False, "font": {"size": 16}}],
            height=450,  
            margin=dict(t=40, b=0, l=0, r=0)
        )
        title = _mutations_heading(v, years_text, filters_text, has_filters)
        return fig, title

    # unique sequences per mutation
    mutation_counts = (
        mutation_df.groupby("mutation")["ID"].nunique().reset_index()
        .rename(columns={"mutation": "Mutation", "ID": "Unique_Sequences"})
    )
    mutation_counts["Proportion"] = (mutation_counts["Unique_Sequences"] / total_sequences * 100.0).round(2)
    
    # FILTER OUT MUTATIONS WITH 0% - KEY FIX
    mutation_counts = mutation_counts[mutation_counts["Proportion"] > 0.1]
    
    # If no mutations left after filtering, return empty plot
    if mutation_counts.empty:
        fig.update_layout(
            xaxis={"visible": False}, yaxis={"visible": False},
            annotations=[{"text": "No mutations with >0% frequency found",
                          "xref": "paper", "yref": "paper", "x": 0.5, "y": 0.5,
                          "showarrow": False, "font": {"size": 16}}],
            height=450,  
            margin=dict(t=40, b=0, l=0, r=0)
        )
        title = _mutations_heading(v, years_text, filters_text, has_filters)
        return fig, title
    
    mutation_counts = mutation_counts.sort_values("Proportion", ascending=False).head(20)
    
    # Calculate global y-axis maximum
    GLOBAL_Y_MAX = calculate_global_mutation_maximum()
    
    # Get current maximum from the data
    current_max = mutation_counts["Proportion"].max() if not mutation_counts.empty else 0
    
    # Determine final y_max - use whichever is larger: global max or current max + padding
    if current_max > GLOBAL_Y_MAX:
        # If current data exceeds global max, round up to nearest 10 above current max
        y_max = min(100, np.ceil(current_max / 10) * 10 + 10)
    else:
        # Use the global maximum for consistent scaling
        y_max = GLOBAL_Y_MAX

    fig = px.bar(
        mutation_counts,
        x="Mutation", y="Proportion",
        labels={"Proportion": "Sequences with Mutation (%)", "Mutation": "Mutation"},
        color_discrete_sequence=[color],
    )
    fig.update_traces(
        marker_cornerradius=6,
        hovertemplate="<b>%{x}</b><br>Percentage: %{y:.1f}%<br>"
                      "Count: %{customdata[0]}<br>"
                      "Total sequences: %{customdata[1]}<extra></extra>",
        customdata=np.stack([
            mutation_counts["Unique_Sequences"].astype(int).to_numpy(),
            np.full(len(mutation_counts), int(total_sequences))
        ], axis=1),
        texttemplate="%{y:.1f}%", textposition="outside", cliponaxis=False
    )
    
    # Update layout with consistent y-axis range
    fig.update_layout(
        height=450, 
        margin=dict(t=16, b=0, l=0, r=0),
         
        
        xaxis=dict(
            showgrid=False, 
            linecolor="rgba(0,0,0,0.25)", 
            tickangle=-45, 
            automargin=True
        ),
        yaxis=dict(
            title="Sequences with Mutation (%)", 
            gridcolor="rgba(0,0,0,0.08)",
            zeroline=False, 
            range=[-2, y_max],  # Consistent upper limit
            fixedrange=True,
            # Consistent tick marks for both viruses
            tickmode='linear',
            tick0=0,
            dtick=10
        )
    )

    title = _mutations_heading(v, years_text, filters_text, has_filters)
    return fig, title
    
def make_coverage_bar(
    gap_df: pd.DataFrame,
    selected_virus: str,
    order: str = "lowest",
    top_n: int = 20
) -> go.Figure:
    virus = (selected_virus or "HBV").upper()
    
    if gap_df.empty:
        return go.Figure().update_layout(
            height=450,
            annotations=[dict(text="No data available for current filters",
                              x=0.5, y=0.5, showarrow=False)],
            
            paper_bgcolor="white"
        )

    # --- normalize expected columns & dtypes ---
    # allow for different casings or missing columns
    colmap = {c.lower(): c for c in gap_df.columns}
    def col(name):    # case-insensitive getter
        return colmap.get(name.lower(), name)

    for c in ["observed_sequences", "expected_sequences", "coverage_gap"]:
        if col(c) in gap_df.columns:
            gap_df[col(c)] = pd.to_numeric(gap_df[col(c)], errors="coerce")

    # compute Coverage_Ratio if missing (case-insensitively)
    if "coverage_ratio" in colmap:
        gap_df["Coverage_Ratio"] = pd.to_numeric(gap_df[col("coverage_ratio")], errors="coerce")
    else:
        # create it from observed/expected
        obs = gap_df[col("observed_sequences")] if col("observed_sequences") in gap_df.columns else np.nan
        exp = gap_df[col("expected_sequences")] if col("expected_sequences") in gap_df.columns else np.nan
        gap_df["Coverage_Ratio"] = np.where(pd.to_numeric(exp, errors="coerce") > 0,
                                            pd.to_numeric(obs, errors="coerce") / pd.to_numeric(exp, errors="coerce"),
                                            np.nan)

    # basic clean
    if col("coverage_gap") in gap_df.columns:
        gap_df = gap_df.dropna(subset=[col("coverage_gap")])
    gap_df["Coverage_Ratio"] = pd.to_numeric(gap_df["Coverage_Ratio"], errors="coerce")

    # --- select / order top N ---
    if order == "lowest":
        # "under-sequenced" = biggest additional genomes needed
        gap_df = gap_df.sort_values(col("coverage_gap"), ascending=False).head(top_n)
        gap_df = gap_df.sort_values(col("coverage_gap"), ascending=True)  # for horizontal bar (small→large)
        x_col = col("coverage_gap")
        x_label = "Estimated additional genomes needed"
        title = f"Top {len(gap_df)} under-sequenced countries ({virus})"
    else:
        # "best covered" = highest coverage ratio
        gap_df = gap_df.sort_values("Coverage_Ratio", ascending=False).head(top_n)
        gap_df = gap_df.sort_values("Coverage_Ratio", ascending=True)      # small→large left→right
        x_col = "Coverage_Ratio"
        x_label = "Coverage ratio (Observed / Expected)"
        title = f"Top {len(gap_df)} best-covered countries ({virus})"

    if virus == "HBV":
        color = "#3182bd"
    elif virus == "HCV":
        color = "#d94801"
    else:
        color = "#d94801"

    fig = px.bar(
        gap_df,
        x=x_col,
        y="Country_standard",
        orientation="h",
        color_discrete_sequence=[color],
        labels={x_col: x_label, "Country_standard": ""},
        title=title,
    )
    fig.update_layout(
        height=450,
        margin=dict(t=60, r=10, b=30, l=10),
        xaxis=dict(title=x_label, gridcolor="rgba(0,0,0,0.08)", linecolor="rgba(0,0,0,0.25)", range=[0, None], fixedrange=True),
        yaxis=dict(title="", fixedrange=True),
        
        
    )
    return fig
    
def _empty_world(message: str) -> go.Figure:
    """Create an empty world map with a message."""
    fig = go.Figure()
    fig.update_layout(
        geo=dict(
            showframe=False,
            showcoastlines=True,
            projection_type='equirectangular'
        ),
        annotations=[dict(
            text=message,
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=16)
        )],
        height=500
    )
    return fig


#Time Series with Projections
def get_forecast_log_axis(values):
    """
    Clean custom ticks for Plotly log axes.

    Plotly rule:
    - yaxis.range uses log10 values
    - yaxis.tickvals must use real/raw values
    """
    values = pd.to_numeric(pd.Series(values), errors="coerce")
    values = values.replace([np.inf, -np.inf], np.nan).dropna()
    values = values[values > 0]

    if values.empty:
        return {
            "range": [6, 8],
            "tickvals": [1_000_000, 10_000_000, 100_000_000],
            "ticktext": ["1M", "10M", "100M"],
        }

    ymin = float(values.min())
    ymax = float(values.max())

    lo_val = max(1, ymin / 1.4)
    hi_val = max(lo_val * 1.2, ymax * 1.4)

    lo = np.floor(np.log10(lo_val))
    hi = np.ceil(np.log10(hi_val))

    tickvals = []
    ticktext = []

    for dec in range(int(lo), int(hi) + 1):
        for mult in [1, 2, 5]:
            val = mult * (10 ** dec)

            if lo_val <= val <= hi_val:
                tickvals.append(val)

                if val >= 1_000_000_000:
                    ticktext.append(f"{val / 1_000_000_000:g}B")
                elif val >= 1_000_000:
                    ticktext.append(f"{val / 1_000_000:g}M")
                elif val >= 1_000:
                    ticktext.append(f"{val / 1_000:g}k")
                else:
                    ticktext.append(f"{val:g}")

    if not tickvals:
        tickvals = [ymin]
        ticktext = [f"{ymin:g}"]

    return {
        "range": [float(np.log10(lo_val)), float(np.log10(hi_val))],
        "tickvals": tickvals,
        "ticktext": ticktext,
    }


def create_forecast_chart(
    ihme_df,
    selected_virus,
    sex,
    selected_regions=None,
    selected_countries=None
):
    """Show historical burden trends with simple linear projections and clean log-axis labels."""

    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }

    selected_virus = (selected_virus or "HBV").upper()
    cause_filter = cause_lookup.get(selected_virus)

    if ihme_df is None or ihme_df.empty:
        return _empty_plot("No IHME burden data loaded")

    if cause_filter not in ihme_df["cause"].values:
        return _empty_plot(f"No {selected_virus} burden data available for forecasting")

    cause_data = ihme_df[ihme_df["cause"] == cause_filter].copy()
    available_metrics = cause_data["metric"].dropna().unique()

    metric_to_use = None
    for metric in ["Number", "Rate", "Percent"]:
        if metric in available_metrics:
            metric_to_use = metric
            break

    if not metric_to_use:
        return _empty_plot(f"No suitable metric found. Available: {list(available_metrics)}")

    if sex == "Both":
        male_data = ihme_df[
            (ihme_df["cause"] == cause_filter) &
            (ihme_df["metric"] == metric_to_use) &
            (ihme_df["sex"] == "Male")
        ].copy()

        female_data = ihme_df[
            (ihme_df["cause"] == cause_filter) &
            (ihme_df["metric"] == metric_to_use) &
            (ihme_df["sex"] == "Female")
        ].copy()

        base = pd.concat([male_data, female_data], ignore_index=True)
    else:
        base = ihme_df[
            (ihme_df["sex"] == sex) &
            (ihme_df["cause"] == cause_filter) &
            (ihme_df["metric"] == metric_to_use)
        ].copy()

    if selected_regions:
        base = base[base["WHO_Regions"].isin(selected_regions)]

    if selected_countries:
        base = base[base["Country_standard"].isin(selected_countries)]

    if base.empty:
        return _empty_plot(f"No data available for {selected_virus} with current filters")

    base["year"] = pd.to_numeric(base["year"], errors="coerce")
    base["val"] = pd.to_numeric(base["val"], errors="coerce")
    base = base.dropna(subset=["year", "val"])
    base = base[base["val"] > 0]

    if base.empty:
        return _empty_plot("No positive values available for log-scale forecasting")

    latest_data_year = int(base["year"].max())

    fig = go.Figure()

    measure_colors = {
        "Prevalence": "#4FAEFF",
        "Incidence": "#FFD166",
        "Deaths": "#FF4D6D",
    }

    all_y_values = []
    data_found = False

    for measure in ["Prevalence", "Incidence", "Deaths"]:
        measure_data = base[base["measure"] == measure].copy()

        if measure_data.empty:
            continue

        yearly_data = (
            measure_data
            .groupby("year", as_index=False)["val"]
            .sum()
            .sort_values("year")
        )

        yearly_data = yearly_data[yearly_data["val"] > 0]

        if yearly_data.empty:
            continue

        years = yearly_data["year"].values.astype(float)
        values = yearly_data["val"].values.astype(float)

        all_y_values.extend(values[values > 0].tolist())

        if len(yearly_data) < 2:
            fig.add_trace(
                go.Scatter(
                    x=years,
                    y=values,
                    name=f"{measure} ({len(yearly_data)} point)",
                    mode="markers",
                    marker=dict(size=8, color=measure_colors.get(measure)),
                    hovertemplate=f"{measure}: %{{y:,.0f}}<extra></extra>",
                )
            )
            data_found = True
            continue

        try:
            n = len(years)
            sum_x = np.sum(years)
            sum_y = np.sum(values)
            sum_xy = np.sum(years * values)
            sum_x2 = np.sum(years * years)

            denominator = n * sum_x2 - sum_x * sum_x

            fig.add_trace(
                go.Scatter(
                    x=years,
                    y=values,
                    name=f"{measure} Historical",
                    line=dict(color=measure_colors.get(measure), width=3),
                    mode="lines+markers",
                    marker=dict(size=6),
                    hovertemplate=f"{measure}: %{{y:,.0f}}<extra></extra>",
                )
            )

            data_found = True

            if denominator != 0 and latest_data_year < 2030:
                m = (n * sum_xy - sum_x * sum_y) / denominator
                b = (sum_y - m * sum_x) / n

                future_years = np.arange(latest_data_year + 1, 2031)
                forecast_pred = m * future_years + b

                # Log axes cannot display zero/negative values.
                forecast_pred = np.where(forecast_pred > 0, forecast_pred, np.nan)

                positive_forecast = forecast_pred[np.isfinite(forecast_pred) & (forecast_pred > 0)]
                all_y_values.extend(positive_forecast.tolist())

                if len(positive_forecast) > 0:
                    fig.add_trace(
                        go.Scatter(
                            x=future_years,
                            y=forecast_pred,
                            name=f"{measure} Forecast",
                            line=dict(
                                color=measure_colors.get(measure),
                                width=2,
                                dash="dash",
                            ),
                            mode="lines+markers",
                            marker=dict(size=5),
                            hovertemplate=f"{measure} forecast: %{{y:,.0f}}<extra></extra>",
                        )
                    )

        except Exception:
            fig.add_trace(
                go.Scatter(
                    x=years,
                    y=values,
                    name=f"{measure}",
                    line=dict(color=measure_colors.get(measure), width=2),
                    mode="markers",
                    marker=dict(size=6),
                    hovertemplate=f"{measure}: %{{y:,.0f}}<extra></extra>",
                )
            )
            data_found = True

    if not data_found:
        available_measures = base["measure"].dropna().unique()
        return _empty_plot(
            f"No forecast data available for {selected_virus}. "
            f"Available measures: {list(available_measures)}. "
            f"Using metric: {metric_to_use}"
        )

    if latest_data_year:
        fig.add_vline(
            x=latest_data_year,
            line_dash="dot",
            line_color="#FF4D6D",
            line_width=1.5,
            annotation_text="Data limit",
            annotation_position="top left",
        )

    fig.add_vline(
        x=2030,
        line_dash="dot",
        line_color="#8BC34A",
        line_width=1.5,
        annotation_text="WHO 2030",
        annotation_position="top right",
    )

    axis_config = get_forecast_log_axis(all_y_values)

    fig.update_layout(
        xaxis=dict(
            title="Year",
            fixedrange=True,
            gridcolor="rgba(255,255,255,0.08)",
            linecolor="rgba(255,255,255,0.25)",
            showgrid=True,
            zeroline=False,
        ),
        yaxis=dict(
            title=dict(
                text=f"Value ({metric_to_use})",
                standoff=22,
            ),
            type="log",
            autorange=False,
            range=axis_config["range"],
            tickmode="array",
            tickvals=axis_config["tickvals"],
            ticktext=axis_config["ticktext"],
            fixedrange=True,
            gridcolor="rgba(255,255,255,0.08)",
            linecolor="rgba(255,255,255,0.25)",
            showgrid=True,
            zeroline=False,
            automargin=True,
        ),
        height=400,
        hovermode="x unified",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="center",
            x=0.5,
        ),
        margin=dict(t=60, b=45, l=90, r=25),
    )

    return fig

#Mutation Timeline
def create_mutation_timeline(mutation_df, sequence_df, selected_virus, top_mutations=10):
    """Show emergence and spread of key mutations over time"""
    
    # Enrich mutation data with temporal information
    enriched_mutations = _enrich_mutation_df(mutation_df, sequence_df)
    
    if enriched_mutations.empty:
        return _empty_plot("No mutation data available")
    
    # Get top mutations by frequency
    top_muts = (enriched_mutations.groupby("mutation")["ID"]
                .nunique()
                .nlargest(top_mutations)
                .index.tolist())
    
    # Aggregate by year and mutation
    timeline_data = (enriched_mutations[enriched_mutations["mutation"].isin(top_muts)]
                    .groupby(["Year", "mutation"])
                    .size()
                    .reset_index(name="count"))
    
    # Calculate cumulative prevalence
    yearly_totals = sequence_df.groupby("Year").size().reset_index(name="total_sequences")
    timeline_data = timeline_data.merge(yearly_totals, on="Year", how="left")
    timeline_data["prevalence_pct"] = (timeline_data["count"] / timeline_data["total_sequences"]) * 100
    
    timeline_data = timeline_data.sort_values(["mutation", "Year"])

    fig = px.line(
        timeline_data,
        x="Year",
        y="prevalence_pct",
        color="mutation",
        markers=True,
        hover_data={"count": True, "prevalence_pct": ":.2f"},
    )

    # Scale marker size by sample count (per mutation, so each trace's own
    # min/max sets its point sizes) instead of a single flat gray line
    # underneath differently-colored dots.
    for trace in fig.data:
        mutation_counts = timeline_data.loc[
            timeline_data["mutation"] == trace.name, "count"
        ]
        if mutation_counts.empty:
            continue
        c_min, c_max = mutation_counts.min(), mutation_counts.max()
        if c_max > c_min:
            sizes = 6 + (mutation_counts - c_min) / (c_max - c_min) * 14
        else:
            sizes = [10] * len(mutation_counts)
        trace.update(
            line=dict(width=2),
            marker=dict(size=list(sizes), line=dict(width=0)),
            opacity=0.9,
        )

    fig.update_layout(
        yaxis_title="Prevalence (%)",
        yaxis=dict(range=[0, None], fixedrange=True, rangemode="nonnegative"),
        xaxis_title="Year",
        xaxis=dict(fixedrange=True),
        height=420,
        hovermode="closest",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
    )

    return fig

#Transmission Cluster Map
def create_transmission_clusters(sequence_df, selected_virus, genetic_distance_threshold=0.05):    
    if sequence_df.empty:
        return _empty_plot("No sequence data available for cluster analysis")
    
    # Use temporal and geographic patterns since we don't have genetic distance data
    cluster_data = sequence_df.copy()
    
    # Group by country and count sequences per year
    country_year_counts = (cluster_data.groupby(["Country_standard", "Year"])
                          .size()
                          .reset_index(name="sequence_count"))
    
    # Identify countries with increasing sequence counts (potential outbreaks)
    clusters = []
    for country in country_year_counts["Country_standard"].unique():
        country_data = country_year_counts[country_year_counts["Country_standard"] == country].sort_values("Year")
        
        if len(country_data) >= 2:
            # Calculate year-over-year growth
            country_data["growth"] = country_data["sequence_count"].pct_change()
            
            # Flag as cluster if significant growth detected
            significant_growth = country_data[country_data["growth"] > 0.5]  # 50% growth threshold
            
            for _, row in significant_growth.iterrows():
                clusters.append({
                    "Country_standard": country,
                    "Year": row["Year"],
                    "sequences": row["sequence_count"],
                    "growth_pct": row["growth"] * 100,
                    "cluster_type": "Emerging" if row["growth"] > 1.0 else "Growing"
                })
    
    if not clusters:
        # Fallback: show countries with highest recent sequencing activity
        recent_year = cluster_data["Year"].max()
        recent_data = cluster_data[cluster_data["Year"] == recent_year]
        if not recent_data.empty:
            country_counts = recent_data["Country_standard"].value_counts().head(10)
            for country, count in country_counts.items():
                clusters.append({
                    "Country_standard": country,
                    "Year": recent_year,
                    "sequences": count,
                    "growth_pct": 0,
                    "cluster_type": "Active Sequencing"
                })
    
    clusters_df = pd.DataFrame(clusters)
    
    if clusters_df.empty:
        return _empty_plot("No transmission patterns detected with current data")
    
    # Create cluster map with different colors for cluster types
    fig = px.scatter_geo(clusters_df,
                        locations="Country_standard",
                        locationmode="country names",
                        size="sequences",
                        color="cluster_type",
                        hover_name="Country_standard",
                        hover_data={
                            "Year": True, 
                            "sequences": True,
                            "growth_pct": ":.1f",
                            "cluster_type": True
                        },
                        color_discrete_map={
                            "Emerging": "#ff4444",
                            "Growing": "#ffaa00", 
                            "Active Sequencing": "#44ff44"
                        },
                        title=f"{selected_virus} Transmission Patterns and Sequencing Activity")
    
    fig.update_geos(
        projection_type="natural earth",
        showcountries=True,
        countrycolor="rgba(0,0,0,0.2)",
        fitbounds="locations"
    )
    
    fig.update_layout(
        height=500,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="center",
            x=0.5
        )
    )
    
    return fig

#Replace Pie Chart with Treemap:
def create_country_treemap(df, selected_virus, selected_regions=None):    
    if df.empty:
        return _empty_plot("No data available")
    
    # Aggregate data
    if selected_regions and "WHO_Regions" in df.columns:
        df = df[df["WHO_Regions"].isin(selected_regions)]
    
    region_country_counts = (df.groupby(["WHO_Regions", "Country_standard"])
                            .size()
                            .reset_index(name="count"))
    
    fig = px.treemap(region_country_counts,
                     path=["WHO_Regions", "Country_standard"],
                     values="count",
                     color="count",
                     color_continuous_scale='Blues',
                     title=f"{selected_virus} Sequences by Region and Country")
    
    fig.update_layout(height=400)
    fig.update_traces(
        hovertemplate='<b>%{label}</b><br>Sequences: %{value:,}<br>Parent: %{parent}'
    )
    
    return fig

#Country barchart
def create_country_stacked_bar(df, selected_virus, selected_regions=None, selected_countries=None, top_n=10):    
    if df.empty:
        return _empty_plot("No data available")
    
    # Apply region/country filters
    if selected_regions and "WHO_Regions" in df.columns:
        df = df[df["WHO_Regions"].isin(selected_regions)]
    if selected_countries:
        df = df[df["Country_standard"].isin(selected_countries)]
    
    # Get top countries by total sequences
    country_totals = df.groupby('Country_standard').size().sort_values(ascending=False)
    top_countries = country_totals.head(top_n).index
    
    # Filter to top countries
    top_df = df[df['Country_standard'].isin(top_countries)].copy()
    
    # Aggregate by country and genotype
    stacked_data = (top_df.groupby(['Country_standard', 'genotype'])
                    .size()
                    .reset_index(name='count'))
    
    if stacked_data.empty:
        return _empty_plot("No data available for stacked bar chart")
    
    # Get genotype colors based on virus type
    if selected_virus.upper() == "HBV":
        genotype_colors = HBV_GENOTYPE_COLORS
    elif selected_virus.upper() == "HCV":
        genotype_colors = HCV_GENOTYPE_COLORS
    else:
        genotype_colors = HEV_GENOTYPE_COLORS
    
    # Create stacked bar chart
    fig = px.bar(
        stacked_data,
        x='Country_standard',
        y='count',
        color='genotype',
        labels={'count': 'Number of Sequences', 'Country_standard': 'Country'},
        color_discrete_map=genotype_colors
    )
    fig.update_traces(
        marker_cornerradius=6
    )
    
    # Update layout for better readability
    fig.update_layout(
        height=500,
        xaxis_title="Country",
        xaxis=dict(fixedrange=True, categoryorder='total descending'),
        yaxis_title="Number of Sequences",
        yaxis=dict(range=[0, None], fixedrange=True, rangemode="nonnegative"),
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02,
            title="genotype"
        ),
        margin=dict(r=150),  # Add margin for legend
    )
    
    # Update hover template
    fig.update_traces(
        hovertemplate="<b>%{x}</b><br>genotype: %{fullData.name}<br>Sequences: %{y:,}<extra></extra>"
    )
    
    return fig

#Priority Setting Tool
def create_priority_calculator(gap_df, ihme_df, selected_virus, weights=None):
    priority_data = gap_df.copy() if not gap_df.empty else pd.DataFrame()
    
    if priority_data.empty:
        return _empty_plot("No data available for priority calculation"), priority_data

    
    if "Seq_count" in priority_data.columns:
        priority_data["observed_sequences"] = priority_data["Seq_count"]
    
    if "Coverage_ratio" in priority_data.columns:
        priority_data["coverage_gap"] = 1 - priority_data["Coverage_ratio"]
    
    default_weights = {
        "burden": 0.4,
        "coverage_gap": 0.3, 
        "population": 0.2,
        "neighbor_sequencing": 0.1
    }
    weights = weights or default_weights
    
    # Calculate priority scores
    priority_data = gap_df.copy()
    
    # Normalize metrics (handle missing columns safely)
    for metric in ["burden", "coverage_gap", "observed_sequences"]:
        if metric in priority_data.columns:
            col_min = priority_data[metric].min()
            col_max = priority_data[metric].max()
            if col_max > col_min:  # Avoid division by zero
                priority_data[f"{metric}_norm"] = (priority_data[metric] - col_min) / (col_max - col_min)
            else:
                priority_data[f"{metric}_norm"] = 0.5  # Default value if all values are same
    
    # Proper population normalization
    population_col = priority_data.get("Population", pd.Series([1] * len(priority_data)))
    if hasattr(population_col, 'max'):    # Check if it's a Series with max method
        pop_max = population_col.max()
        population_norm = population_col / max(pop_max, 1)    # Avoid division by zero
    else:
        population_norm = 0     # Fallback if Population is not available
    
    # Calculate composite score
    priority_data["priority_score"] = (
        weights["burden"] * priority_data.get("burden_norm", 0) +
        weights["coverage_gap"] * priority_data.get("coverage_gap_norm", 0) +
        weights["population"] * population_norm
    )
    
    # Rank countries
    priority_data = priority_data.sort_values("priority_score", ascending=False)
    priority_data["rank"] = range(1, len(priority_data) + 1)
    
    # Create interactive table
    fig = go.Figure(data=[go.Table(
        header=dict(values=["Rank", "Country", "Priority Score", "Burden", "Coverage Gap", "Sequences"],
                    fill_color='paleturquoise',
                    align='left'),
        cells=dict(values=[priority_data["rank"], 
                          priority_data["Country_standard"],
                          priority_data["priority_score"].round(3),
                          priority_data.get("burden", 0).round(0),
                          priority_data.get("coverage_gap", 0).round(0),
                          priority_data.get("observed_sequences", 0)],
                   align='left'))
    ])
    
    fig.update_layout(
        height=400
    )
    
    return fig, priority_data

#---burden-lines-helper-----

MASTER_TICKS = np.array([3e5, 6e5, 1e6, 2e6, 5e6, 1e7, 2e7, 5e7,
                         1e8, 2e8, 5e8], dtype=float)
MASTER_TEXT     = ["300k","600k","1M","2M","5M","10M","20M","50M",
                "100M","200M","500M"]

def virus_log_axis(selected_virus: str, values: np.ndarray) -> dict:
    """Return {range, tickvals, ticktext} for a single log axis tuned per virus."""
    v = (selected_virus or "HBV").upper()
    # sensible presets reflecting your data ranges
    presets = {
        "HBV": dict(floor=3e5, cap=6e8),   # 70M–400M + deaths ~500–700k
        "HCV": dict(floor=2e5, cap=3e8),   # 9M–200M + deaths ~200–600k
    }
    p = presets.get(v, dict(floor=2e5, cap=9e8))

    vals = np.asarray(values, float)
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if vals.size == 0:
        lo, hi = p["floor"], p["cap"]
    else:
        lo = max(p["floor"], np.nanmin(vals) * 0.9)        # small padding
        hi = min(p["cap"],     np.nanmax(vals) * 1.1)

    # snap to decades so the axis looks clean
    lo = 10 ** np.floor(np.log10(lo))
    hi = 10 ** np.ceil( np.log10(hi))

    mask = (MASTER_TICKS >= lo) & (MASTER_TICKS <= hi)
    tickvals = np.log10(MASTER_TICKS[mask]).tolist()
    ticktext = [t for t,m in zip(MASTER_TEXT, mask) if m]

    return dict(
        range=[float(np.log10(lo)), float(np.log10(hi))],
        tickvals=tickvals,
        ticktext=ticktext
    )

def make_burden_lines(
    ihme_df: pd.DataFrame,
    selected_virus: str,
    selected_continents: list = None,
    selected_countries: list = None,
    selected_years: list = None
) -> tuple[go.Figure, str]:
    
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause_filter = cause_lookup.get((selected_virus or "HBV").upper())
    
    base = ihme_df[
        (ihme_df["sex"] == sex) &
        (ihme_df["age"] == "All ages") &
        (ihme_df["cause"] == cause_filter)
    ].copy()
    
    # Apply year filtering if provided
    if selected_years and len(selected_years) == 2:
        y0, y1 = map(int, selected_years)
        base = base[(base["year"] >= y0) & (base["year"] <= y1)]
    
    if selected_continents:
        base = base[base["WHO_Regions"].isin(selected_continents)]
    if selected_countries:
        base = base[base["Country_standard"].isin(selected_countries)]

    def series(measure, metric="Number"):
        s = base[(base["measure"] == measure) & (base["metric"] == metric)].copy()
        if s.empty:
            return pd.DataFrame(columns=["year", f"{measure}_{metric}"])
        s["val"] = pd.to_numeric(s["val"], errors="coerce")
        s = s.dropna(subset=["val"])
        return (s.groupby("year", as_index=False)["val"]
                  .sum().rename(columns={"val": f"{measure}_{metric}"}))

    d_prev = series("Prevalence", "Number")
    d_inc  = series("Incidence",  "Number")
    d_dea  = series("Deaths",      "Number")

    from functools import reduce
    dfs = [d for d in [d_prev, d_inc, d_dea] if not d.empty]
    if dfs:
        burden_df = reduce(lambda L,R: pd.merge(L,R,on="year",how="outer"), dfs).sort_values("year")
    else:
        burden_df = pd.DataFrame(columns=["year","Prevalence_Number","Incidence_Number","Deaths_Number"])

    # Ensure numeric and mask non-positives for log plot
    for col in ["Prevalence_Number","Incidence_Number","Deaths_Number"]:
        if col in burden_df:
            burden_df[col] = pd.to_numeric(burden_df[col], errors="coerce")

    def ymask(col):
        return burden_df[col].mask(~(burden_df[col] > 0), None)

    # Build the figure
    burden_fig = go.Figure()

    # Prevalence
    if "Prevalence_Number" in burden_df:
        burden_fig.add_trace(go.Scatter(
            x=burden_df["year"], y=ymask("Prevalence_Number"),
            mode="lines+markers", name="Prevalence",
            hovertemplate="Year: %{x}<br>Prevalence: %{y:,.0f}<extra></extra>",
            line=dict(width=3)
        ))

    # Incidence
    if "Incidence_Number" in burden_df:
        burden_fig.add_trace(go.Scatter(
            x=burden_df["year"], y=ymask("Incidence_Number"),
            mode="lines+markers", name="Incidence",
            hovertemplate="Year: %{x}<br>Incidence: %{y:,.0f}<extra></extra>",
            line=dict(width=3)
        ))

    # Deaths
    if "Deaths_Number" in burden_df:
        burden_fig.add_trace(go.Scatter(
            x=burden_df["year"], y=ymask("Deaths_Number"),
            mode="lines+markers", name="Deaths",
            hovertemplate="Year: %{x}<br>Deaths: %{y:,.0f}<extra></extra>",
            line=dict(width=3)
        ))

    # Get values for axis scaling
    vals_for_axis = []
    for col in ("Prevalence_Number","Incidence_Number","Deaths_Number"):
        if col in burden_df.columns:
            vals_for_axis.append(burden_df[col].to_numpy())
    vals_for_axis = np.concatenate(vals_for_axis) if vals_for_axis else np.array([])
    
    # Use the helper function for consistent axis scaling
    axis = virus_log_axis(selected_virus, vals_for_axis)
    
    # Layout with single y-axis
    burden_fig.update_layout(
        height=400,
        margin=dict(t=50, b=20, l=60, r=20),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        xaxis=dict(
            title="Year",
            dtick=3, tickmode="linear",
            gridcolor="rgba(0,0,0,0.06)", showgrid=True,
            fixedrange=True
        ),
        yaxis=dict(
            title="Number of Cases (log scale)",
            type="log",
            autorange=False,
            range=axis["range"],
            tickmode="array",
            tickvals=axis["tickvals"],
            ticktext=axis["ticktext"],
            gridcolor="rgba(0,0,0,0.05)",
            zeroline=False,
            fixedrange=True,
            constrain="domain",
        ),
        
        paper_bgcolor="white"
    )
    
    burden_title = f"Global burdens in {selected_virus}: Prevalence, Incidence & Deaths"
    
    return burden_fig, burden_title
    
def _empty_plot(message):
    """Create an empty plot with a message"""
    fig = go.Figure()
    fig.update_layout(
        xaxis={"visible": False},
        yaxis={"visible": False},
        annotations=[{
            "text": message,
            "xref": "paper", "yref": "paper",
            "x": 0.5, "y": 0.5, "showarrow": False,
            "font": {"size": 14, "color": "#64748B", "family": "IBM Plex Sans, sans-serif"}
        }],
        height=300
    )
    return apply_heptracker_figure_style(fig)

def _df_to_json(df: pd.DataFrame) -> str:
    if df is None or df.empty:
        return json.dumps({"columns": [], "data": []})
    return json.dumps({"columns": list(df.columns), "data": df.to_dict("records")})

def _df_from_json(payload: str) -> pd.DataFrame:
    try:
        obj = json.loads(payload or "{}")
        return pd.DataFrame(obj.get("data", []), columns=obj.get("columns", []))
    except Exception:
        return pd.DataFrame()
        
# Mutations map helper
def create_mutation_map(mutation_df, coord_lookup, virus_type, mutation_type=None):
    if mutation_df.empty:
        return _empty_world("No mutation data available")
    
    # Group by country
    country_mutations = mutation_df.groupby("Country_standard").size().reset_index(name="mutation_count")
    
    # Merge with coordinates
    fig = go.Figure()
    
    if not country_mutations.empty:
        # Create choropleth
        fig.add_trace(go.Choropleth(
            locations=country_mutations["Country_standard"],
            locationmode="country names",
            z=country_mutations["mutation_count"],
            colorscale="Reds",
            marker_line_color="rgba(0,0,0,0.3)",
            marker_line_width=0.5,
            hovertemplate="<b>%{location}</b><br>Mutations: %{z}<extra></extra>"
        ))
    
    fig.update_geos(
        projection_type="natural earth",
        showcountries=True,
        countrycolor="rgba(0,0,0,0.2)",
        showsubunits=True,
        fitbounds="locations"
    )
    
    fig.update_layout(
        height=500,
        title=f"{virus_type} Mutation Distribution Map" + (f" - {mutation_type}" if mutation_type else ""),
        margin=dict(t=50, b=30, l=10, r=10)
    )
    
    return fig

def create_drug_resistance_profile(mutation_df, virus_type):
    if mutation_df.empty or "drug" not in mutation_df.columns:
        return _empty_plot("No drug resistance data available")
    
    # Filter for antiviral resistance
    resistance_mutations = mutation_df[mutation_df["type"] == "antiviral_resistance"]
    
    if resistance_mutations.empty:
        return _empty_plot("No antiviral resistance mutations")
    
    # Group by drug and count unique samples
    drug_resistance = (resistance_mutations.groupby("drug")["ID"]
                      .nunique()
                      .reset_index()
                      .rename(columns={"ID": "sample_count"}))
    
    drug_resistance = drug_resistance.sort_values("sample_count", ascending=True)
    
    # Create horizontal bar chart
    fig = px.bar(
        drug_resistance,
        y="drug",
        x="sample_count",
        orientation="h",
        color="sample_count",
        color_continuous_scale="Viridis",
        title=f"{virus_type} Drug Resistance Profile"
    )
    
    fig.update_traces(
        hovertemplate="<b>%{y}</b><br>Samples with resistance: %{x}<extra></extra>"
    )
    
    fig.update_layout(
        height=400,
        xaxis_title="Number of Samples with Resistance",
        yaxis_title="Drug",
        coloraxis_showscale=False,
        
        paper_bgcolor="white"
    )
    
    return fig

def build_filter_bar_with_logo():
    """Sticky top navbar (Row 1: site brand + page navigation) with Data Control Bar (Row 2: filters + search + export) below."""
    header_bar = html.Div(
        className="hep-top-navbar sticky-top mb-3 shadow-sm",
        id="hep-header-bar",
        children=[
            dbc.Container(
                [
                    dbc.Row(
                        [
                            # Left: Brand Logo & Title
                            dbc.Col(
                                [
                                    html.A(
                                        href="/dashboard",
                                        className="d-flex align-items-center text-decoration-none py-1",
                                        children=[
                                            html.Div(
                                                html.Img(
                                                    src="/assets/Hepatitis_genome2.jpg",
                                                    alt="HepTracker Logo",
                                                    className="w-100 h-100 rounded-3",
                                                    style={"objectFit": "cover"}
                                                ),
                                                className="brand-icon-bg p-1 rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0 shadow-sm border border-secondary border-opacity-25 bg-dark",
                                                style={"width": "62px", "height": "62px", "overflow": "hidden"},
                                            ),
                                            html.Div([
                                                html.Span("HepTracker", className="fw-extrabold lh-1 d-block mb-1", style={"fontSize": "1.6rem", "letterSpacing": "-0.4px", "color": "#38BDF8"}),
                                                html.Span("Hepatitis Genomic Surveillance Dashboard", className="fw-bold d-block", style={"fontSize": "0.78rem", "letterSpacing": "1.1px", "color": "#A8B4C7"}),
                                            ]),
                                        ],
                                    )
                                ],
                                xs=12, md=7,
                                className="d-flex align-items-center mb-2 mb-md-0",
                            ),
                            # Right: Navigation Links
                            dbc.Col(
                                [
                                    html.Div(
                                        className="d-flex align-items-center justify-content-md-end gap-2",
                                        children=[
                                            dcc.Link(
                                                [html.I(className="fa-solid fa-gauge me-1"), "Dashboard"],
                                                href="/dashboard",
                                                className="hep-nav-header-link active",
                                            ),
                                            dcc.Link(
                                                [html.I(className="fa-solid fa-book-open me-1"), "Resources"],
                                                href="/resources",
                                                className="hep-nav-header-link",
                                            ),
                                            dcc.Link(
                                                [html.I(className="fa-solid fa-circle-info me-1"), "About"],
                                                href="/about",
                                                className="hep-nav-header-link",
                                            ),
                                        ],
                                    )
                                ],
                                xs=12, md=5,
                            ),
                        ],
                        className="align-items-center py-2 px-3",
                    )
                ],
                fluid=True,
            )
        ],
    )

    return [header_bar]


def build_section_selector():
    """Section selector tabs aligned with the main dashboard content width."""
    return html.Div(
        className="hep-section-selector-container",
        children=[
            html.Div(
                dbc.ButtonGroup(
                    [
                        dbc.Button(
                            [html.I(className="fa-solid fa-map-location-dot me-2"), "Surveillance Map"],
                            id="tab-overview",
                            n_clicks=1,
                            className="hep-section-tab-btn active",
                        ),
                        dbc.Button(
                            [html.I(className="fa-solid fa-chart-line me-2"), "Epidemiology"],
                            id="tab-epidemiology",
                            n_clicks=0,
                            className="hep-section-tab-btn",
                        ),
                        dbc.Button(
                            [html.I(className="fa-solid fa-chart-pie me-2"), "Genotypes"],
                            id="tab-genotypes",
                            n_clicks=0,
                            className="hep-section-tab-btn",
                        ),
                        dbc.Button(
                            [html.I(className="fa-solid fa-dna me-2"), "Mutations"],
                            id="tab-mutations",
                            n_clicks=0,
                            className="hep-section-tab-btn",
                        ),
                        dbc.Button(
                            [html.I(className="fa-solid fa-vial me-2"), "My sequences"],
                            id="tab-user-seq",
                            n_clicks=0,
                            className="hep-section-tab-btn",
                        ),
                    ],
                    className="w-100 justify-content-center flex-wrap gap-2",
                ),
                className="mx-auto",
                style={"maxWidth": "1400px"},
            )
        ],
    )


def build_overview_summary_cards():
    """Horizontal summary cards shown above filters/plots (Optional)."""
    return html.Div()


def build_epi_summary_cards():
    """Horizontal summary cards for Epidemiology tab."""
    return html.Div()


def make_priority_table(priority_data):
    """Render priority ranking as a Dash DataTable instead of a Plotly go.Table."""
    if priority_data is None or priority_data.empty:
        return html.Div(
            "No priority data available for current filters",
            className="hep-empty-table text-muted p-4 text-center",
        )

    display_df = priority_data.copy()

    columns_needed = [
        "rank",
        "Country_standard",
        "priority_score",
        "burden",
        "coverage_gap",
        "observed_sequences",
    ]

    available_cols = [c for c in columns_needed if c in display_df.columns]
    display_df = display_df[available_cols].copy()

    rename_map = {
        "rank": "Rank",
        "Country_standard": "Country",
        "priority_score": "Priority Score",
        "burden": "Burden",
        "coverage_gap": "Coverage Gap",
        "observed_sequences": "Sequences",
    }
    display_df = display_df.rename(columns=rename_map)

    if "Priority Score" in display_df.columns:
        display_df["Priority Score"] = pd.to_numeric(display_df["Priority Score"], errors="coerce").round(3)

    for col in ["Burden", "Coverage Gap", "Sequences"]:
        if col in display_df.columns:
            display_df[col] = pd.to_numeric(display_df[col], errors="coerce").apply(
                lambda x: f"{x:,.0f}" if pd.notna(x) else "N/A"
            )

    return dash_table.DataTable(
        data=display_df.to_dict("records"),
        columns=[{"name": col, "id": col} for col in display_df.columns],
        page_size=8,
        sort_action="native",
        style_table={
            "overflowX": "auto",
            "backgroundColor": "transparent",
        },
        style_header=TABLE_HEADER_STYLE,
        style_cell=TABLE_CELL_STYLE,
        style_data_conditional=[
            TABLE_ODD_ROW_STYLE,
            {
                "if": {"column_id": "Rank"},
                "color": VIRUS_COLORS["HCV"],
                "fontWeight": "800",
            },
            {
                "if": {"column_id": "Priority Score"},
                "color": VIRUS_COLORS["HEV"],
                "fontWeight": "700",
            },
        ],
    )
# === APP SETUP ==============================================================
external_stylesheets = [
    "https://cdn.jsdelivr.net/npm/bootstrap-icons@1.10.5/font/bootstrap-icons.css",
    dbc.themes.BOOTSTRAP
]

# === APP LAYOUT ==============================================================
def create_dashboard_layout():
    stores = html.Div(
        [
            dcc.Store(id="filtered-store"),
            dcc.Store(id="gap-store"),
            dcc.Store(id="computed-metrics-store"),
            dcc.Store(id="ihme-latest-store"),
            dcc.Store(id="priority-data-store"),
            dcc.Download(id="priority-download"),
            dcc.Download(id="download-mutations"),
            dcc.Download(id="download-mutation-report"),
            html.Div(id="download-trigger", style={"display": "none"}),
            dcc.Download(id="download-data"),
            *USER_SEQ_STORES(),
        ],
        style={"display": "none"},
    )

    main_children = [
        *build_filter_bar_with_logo(),

        # Section Selector Tabs at the top of main dashboard layout
        build_section_selector(),

        # === TAB 1: SURVEILLANCE MAP (OVERVIEW CONTENT) ===
        html.Div(
            id="overview-content",
            className="pb-5",
            style={"paddingBottom": "80px"},
            children=[
                # LEVEL 1: WHAT DOES THIS DATASET CONTAIN? — Summary KPI Cards (Full Width)
                dbc.Row(
                    [
                        # KPI 1: Total Sequences
                        dbc.Col(
                            dbc.Card(
                                dbc.CardBody([
                                    html.Div(
                                        [
                                            html.Div(
                                                html.I(className="fa-solid fa-dna text-primary fs-4"),
                                                className="p-2 bg-primary-subtle text-primary rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                                                style={"width": "44px", "height": "44px"}
                                            ),
                                            html.Div([
                                                html.Div("TOTAL SEQUENCES", className="text-muted small fw-bold tracking-wider mb-1", style={"fontSize": "0.68rem"}),
                                                html.H3(id="indicator-total", className="fw-extrabold text-dark mb-0 fs-3"),
                                            ])
                                        ],
                                        className="d-flex align-items-center"
                                    )
                                ]),
                                className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
                            ),
                            xs=12, sm=6, lg=3, className="mb-3 mb-lg-0"
                        ),
                        # KPI 2: Countries Represented
                        dbc.Col(
                            dbc.Card(
                                dbc.CardBody([
                                    html.Div(
                                        [
                                            html.Div(
                                                html.I(className="fa-solid fa-earth-americas text-success fs-4"),
                                                className="p-2 bg-success-subtle text-success rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                                                style={"width": "44px", "height": "44px"}
                                            ),
                                            html.Div([
                                                html.Div("COUNTRIES", className="text-muted small fw-bold tracking-wider mb-1", style={"fontSize": "0.68rem"}),
                                                html.H3(id="indicator-countries", className="fw-extrabold text-dark mb-0 fs-3"),
                                            ])
                                        ],
                                        className="d-flex align-items-center"
                                    )
                                ]),
                                className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
                            ),
                            xs=12, sm=6, lg=3, className="mb-3 mb-lg-0"
                        ),
                        # KPI 3: Unique Genotypes
                        dbc.Col(
                            dbc.Card(
                                dbc.CardBody([
                                    html.Div(
                                        [
                                            html.Div(
                                                html.I(className="fa-solid fa-microscope text-info fs-4"),
                                                className="p-2 bg-info-subtle text-info rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                                                style={"width": "44px", "height": "44px"}
                                            ),
                                            html.Div([
                                                html.Div("GENOTYPES", className="text-muted small fw-bold tracking-wider mb-1", style={"fontSize": "0.68rem"}),
                                                html.H3(id="indicator-genotypes", className="fw-extrabold text-dark mb-0 fs-3"),
                                            ])
                                        ],
                                        className="d-flex align-items-center"
                                    )
                                ]),
                                className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
                            ),
                            xs=12, sm=6, lg=3, className="mb-3 mb-lg-0"
                        ),
                        # KPI 4: Date Range
                        dbc.Col(
                            dbc.Card(
                                dbc.CardBody([
                                    html.Div(
                                        [
                                            html.Div(
                                                html.I(className="fa-solid fa-calendar-days text-warning fs-4"),
                                                className="p-2 bg-warning-subtle text-warning rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                                                style={"width": "44px", "height": "44px"}
                                            ),
                                            html.Div([
                                                html.Div("DATE RANGE", className="text-muted small fw-bold tracking-wider mb-1", style={"fontSize": "0.68rem"}),
                                                html.H3(id="indicator-years", children="1963–2025", className="fw-extrabold text-dark mb-0 fs-3"),
                                            ])
                                        ],
                                        className="d-flex align-items-center"
                                    )
                                ]),
                                className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
                            ),
                            xs=12, sm=6, lg=3, className="mb-3 mb-lg-0"
                        ),
                    ],
                    className="mb-4 g-3 mx-auto",
                    style={"maxWidth": "1400px"}
                ),

                # LEVEL 2: WHERE ARE THE SEQUENCES AND HOW CAN I INTERROGATE THEM? — Compact Filters (Left) + Map (Right)
                dbc.Row(
                    [
                        # LEFT: Compact Vertical Filter Card Aligned with Map
                        dbc.Col(
                            [
                                dbc.Card(
                                    [
                                        dbc.CardBody([
                                            html.Div([
                                                html.I(className="fa-solid fa-sliders text-primary me-2"),
                                                html.Span("Filters", className="fw-extrabold small tracking-wider text-dark"),
                                            ], className="d-flex align-items-center mb-3 pb-2 border-bottom"),

                                            # Virus Dropdown
                                            html.Div([
                                                html.Label("VIRUS", className="hep-filter-label fw-bold small text-muted mb-1"),
                                                dcc.Dropdown(
                                                    id="virus-dropdown",
                                                    options=[
                                                        {"label": "Hepatitis B (HBV)", "value": "HBV"},
                                                        {"label": "Hepatitis C (HCV)", "value": "HCV"},
                                                        {"label": "Hepatitis E (HEV)", "value": "HEV"},
                                                        {"label": "Hepatitis A (HAV)", "value": "HAV", "disabled": True},
                                                        {"label": "Hepatitis D (HDV)", "value": "HDV", "disabled": True},
                                                    ],
                                                    value="HBV",
                                                    clearable=False,
                                                    className="hep-dropdown mb-3",
                                                ),
                                            ]),

                                            # Region Dropdown
                                            html.Div([
                                                html.Label("REGION", className="hep-filter-label fw-bold small text-muted mb-1"),
                                                dcc.Dropdown(
                                                    id="continent-dropdown",
                                                    multi=True,
                                                    placeholder="All regions",
                                                    className="hep-dropdown mb-3",
                                                ),
                                            ]),

                                            # Country Dropdown
                                            html.Div([
                                                html.Label("COUNTRY", className="hep-filter-label fw-bold small text-muted mb-1"),
                                                dcc.Dropdown(
                                                    id="country-dropdown",
                                                    multi=True,
                                                    placeholder="All countries",
                                                    className="hep-dropdown mb-3",
                                                ),
                                            ]),

                                            # Genotype Dropdown
                                            html.Div([
                                                html.Label("GENOTYPE", className="hep-filter-label fw-bold small text-muted mb-1"),
                                                dcc.Dropdown(
                                                    id="genotype-dropdown",
                                                    multi=True,
                                                    placeholder="All genotypes",
                                                    className="hep-dropdown mb-3",
                                                ),
                                            ]),

                                            # Year Range Slider & Inputs
                                            html.Div([
                                                html.Div([
                                                    html.Label("YEAR RANGE", className="hep-filter-label fw-bold small text-muted mb-1"),
                                                    html.Div([
                                                        dcc.Input(
                                                            id="year-start-input",
                                                            type="number",
                                                            min=1963,
                                                            max=2025,
                                                            value=1963,
                                                            className="form-control form-control-sm text-center px-1 hep-year-input",
                                                            style={"width": "58px", "fontSize": "0.78rem"}
                                                        ),
                                                        html.Span("–", className="text-muted small px-1 fw-bold"),
                                                        dcc.Input(
                                                            id="year-end-input",
                                                            type="number",
                                                            min=1963,
                                                            max=2025,
                                                            value=2025,
                                                            className="form-control form-control-sm text-center px-1 hep-year-input",
                                                            style={"width": "58px", "fontSize": "0.78rem"}
                                                        ),
                                                    ], className="d-flex align-items-center me-1"),
                                                ], className="d-flex justify-content-between align-items-center mb-1"),
                                                dcc.RangeSlider(
                                                    id="year-range-slider",
                                                    min=1963,
                                                    max=2025,
                                                    step=1,
                                                    value=[1963, 2025],
                                                    marks={1963: "1963", 1990: "1990", 2010: "2010", 2025: "2025"},
                                                    className="hep-range-slider mb-3",
                                                ),
                                            ]),

                                            # Search Bar
                                            html.Div([
                                                html.Label("SEARCH", className="hep-filter-label fw-bold small text-muted mb-1"),
                                                dbc.InputGroup([
                                                    dbc.Input(
                                                        id="global-search-input",
                                                        type="text",
                                                        placeholder="GenBank ID...",
                                                        className="form-control-sm",
                                                    ),
                                                    dbc.Button(
                                                        html.I(className="fa-solid fa-magnifying-glass"),
                                                        id="global-search-btn",
                                                        color="primary",
                                                        size="sm",
                                                    ),
                                                ], className="mb-3"),
                                            ]),

                                            # Action Buttons: Export & Reset
                                            html.Div([
                                                dbc.Button(
                                                    [html.I(className="fa-solid fa-download me-1"), "Export"],
                                                    id="btn-download-data",
                                                    color="outline-secondary",
                                                    size="sm",
                                                    className="w-50 me-1 fw-bold",
                                                ),
                                                dbc.Button(
                                                    [html.I(className="fa-solid fa-rotate-left me-1"), "Reset"],
                                                    id="btn-reset-filters",
                                                    color="outline-danger",
                                                    size="sm",
                                                    className="w-50 ms-1 fw-bold",
                                                ),
                                            ], className="d-flex justify-content-between pt-2 border-top mt-auto"),
                                        ], className="p-3 d-flex flex-column h-100 justify-content-between")
                                    ],
                                    className="h-100 shadow-sm border-0 bg-white rounded-3",
                                    id="common-filters"
                                )
                            ],
                            xs=12,
                            lg=4,
                            xl=3,
                            className="mb-4 mb-lg-0 d-flex flex-column",
                        ),
                        # RIGHT: Global Distribution Map
                        dbc.Col(
                            [
                                dbc.Card(
                                    [
                                        dbc.CardBody(
                                            [
                                                # Compact map header
                                                dbc.Row(
                                                    [
                                                        # LEFT: title + update date
                                                        dbc.Col(
                                                            [
                                                                html.H4(
                                                                    id="dashboard-page-title",
                                                                    className="mb-0 fw-bold fs-6 text-dark"
                                                                ),
                                                                html.Small(
                                                                    id="dashboard-page-subtitle",
                                                                    className="text-secondary",
                                                                    style={"fontSize": "0.76rem"}
                                                                ),
                        
                                                                # Keep callback target but don't display
                                                                html.Div(
                                                                    id="map-title-sub",
                                                                    style={"display": "none"}
                                                                ),
                                                                html.Div(
                                                                    id="map-title-main",
                                                                    style={"display": "none"}
                                                                ),
                                                            ],
                                                            width=True,
                                                        ),
                        
                                                        # RIGHT: controls
                                                        dbc.Col(
                                                            [
                                                                html.Div(
                                                                    [
                                                                        dcc.RadioItems(
                                                                            id="map-mode",
                                                                            options=[
                                                                                {
                                                                                    "label": "Sequences",
                                                                                    "value": "sequences"
                                                                                },
                                                                                {
                                                                                    "label": "Coverage",
                                                                                    "value": "coverage"
                                                                                },
                                                                                {
                                                                                    "label": "Epidemiology",
                                                                                    "value": "epidemiology"
                                                                                },
                                                                            ],
                                                                            value="sequences",
                                                                            className="hep-btn-group",
                                                                            inputClassName="btn-check",
                                                                            labelClassName="btn btn-outline-secondary btn-sm",
                                                                        ),
                        
                                                                        html.Div(
                                                                            dcc.RadioItems(
                                                                                id="ihme-metric-type",
                                                                                options=[
                                                                                    {
                                                                                        "label": "Deaths",
                                                                                        "value": "Deaths|Number"
                                                                                    },
                                                                                    {
                                                                                        "label": "Incidence",
                                                                                        "value": "Incidence|Number"
                                                                                    },
                                                                                    {
                                                                                        "label": "Prevalence",
                                                                                        "value": "Prevalence|Number"
                                                                                    },
                                                                                ],
                                                                                value="Prevalence|Number",
                                                                                className="hep-btn-group",
                                                                                inputClassName="btn-check",
                                                                                labelClassName="btn btn-outline-secondary btn-sm",
                                                                            ),
                                                                            id="epidemiology-controls",
                                                                            style={"display": "none"},
                                                                        ),
                        
                                                                        dcc.RadioItems(
                                                                            id="display-mode",
                                                                            options=[
                                                                                {
                                                                                    "label": "Raw",
                                                                                    "value": "raw"
                                                                                },
                                                                                {
                                                                                    "label": "Per M",
                                                                                    "value": "PerMillion"
                                                                                },
                                                                            ],
                                                                            value="raw",
                                                                            className="hep-btn-group",
                                                                            inputClassName="btn-check",
                                                                            labelClassName="btn btn-outline-secondary btn-sm",
                                                                        ),
                                                                    ],
                                                                    className=(
                                                                        "d-flex align-items-center "
                                                                        "justify-content-end gap-2 flex-wrap"
                                                                    ),
                                                                )
                                                            ],
                                                            width="auto",
                                                        ),
                                                    ],
                                                    className="align-items-center g-2 mb-1",
                                                ),
                        
                                                # Map
                                                dcc.Loading(
                                                    dcc.Graph(
                                                        id="genotype-map",
                                                        config={
                                                            "displayModeBar": True,
                                                            "displaylogo": False
                                                        },
                                                        style={"height": "480px"},
                                                    ),
                                                    type="circle",
                                                ),
                                            ],
                                            className="px-3 pt-2 pb-2",
                                        )
                                    ],
                                    className="shadow-sm border-0 bg-white rounded-3 h-100",
                                )
                            ],
                            xs=12,
                            lg=8,
                            xl=9,
                            className="d-flex flex-column",
                        ),
                    ],
                    className="mb-4 g-3 mx-auto align-items-stretch",
                    style={"maxWidth": "1400px"}
                ),

                # ROW 4: Sequences by Genotype (Left) + Mutation Types Pie (Right)
                dbc.Row(
                    [
                        # LEFT: Sequences by Genotype Bar Chart
                        dbc.Col(
                            [
                                dbc.Card(
                                    [
                                        dbc.CardBody(
                                            [
                                                html.H5(
                                                    "Sequences by Genotype (Regional Distribution)",
                                                    className="hep-card-title mb-3",
                                                ),
                                                dcc.Loading(
                                                    dcc.Graph(
                                                        id="region-genotype-bar-chart",
                                                        config={"displayModeBar": False},
                                                        style={"height": "350px"},
                                                        className="w-100",
                                                    ),
                                                    type="circle",
                                                ),
                                            ],
                                            className="p-3 d-flex flex-column h-100",
                                        )
                                    ],
                                    className="h-100 shadow-sm border-0 bg-white rounded-3",
                                )
                            ],
                            xs=12,
                            lg=8,
                            className="mb-4 mb-lg-0",
                        ),
                        # RIGHT: Mutation Types Pie Chart
                        dbc.Col(
                            [
                                dbc.Card(
                                    [
                                        dbc.CardBody(
                                            [
                                                html.H5(
                                                    "Mutation Types & Variant Profile",
                                                    className="hep-card-title mb-3",
                                                ),
                                                dcc.Loading(
                                                    dcc.Graph(
                                                        id="homepage-mutation-pie-graph",
                                                        config={"displayModeBar": False},
                                                        style={"height": "350px"},
                                                        className="w-100",
                                                    ),
                                                    type="circle",
                                                ),
                                            ],
                                            className="p-3 d-flex flex-column h-100",
                                        )
                                    ],
                                    className="h-100 shadow-sm border-0 bg-white rounded-3",
                                )
                            ],
                            xs=12,
                            lg=4,
                        ),
                    ],
                    className="mb-4 g-3 mx-auto align-items-stretch",
                    style={"maxWidth": "1400px"},
                ),

                # ROW 5: Burden Forecast (Left) + Sequencing Priority Ranking (Right)
                dbc.Row(
                    [
                        # LEFT: Burden Forecast
                        dbc.Col(
                            [
                                dbc.Card(
                                    [
                                        dbc.CardBody(
                                            [
                                                html.H5(
                                                    "Burden Forecast with Projections",
                                                    className="hep-card-title mb-3",
                                                ),
                                                dcc.Loading(
                                                    dcc.Graph(
                                                        id="homepage-forecast-summary-graph",
                                                        config={"displayModeBar": True, "displaylogo": False},
                                                        style={"height": "350px"},
                                                        className="w-100",
                                                    ),
                                                    type="circle",
                                                ),
                                            ],
                                            className="p-3 d-flex flex-column h-100",
                                        )
                                    ],
                                    className="h-100 shadow-sm border-0 bg-white rounded-3",
                                )
                            ],
                            xs=12,
                            lg=6,
                            className="mb-4 mb-lg-0",
                        ),
                        # RIGHT: Sequencing Priority Ranking Table
                        dbc.Col(
                            [
                                dbc.Card(
                                    [
                                        dbc.CardBody(
                                            [
                                                html.Div(
                                                    [
                                                        html.H5(
                                                            "Sequencing Priority Ranking",
                                                            className="hep-card-title mb-0",
                                                        ),
                                                        dbc.Button(
                                                            "Full Ranking →",
                                                            id="btn-goto-priority-tab",
                                                            color="link",
                                                            className="p-0 text-decoration-none small fw-bold text-primary",
                                                        ),
                                                    ],
                                                    className="d-flex justify-content-between align-items-center mb-3",
                                                ),
                                                dcc.Loading(
                                                    html.Div(
                                                        id="homepage-priority-table-summary",
                                                        style={"minHeight": "300px"},
                                                    ),
                                                    type="circle",
                                                ),
                                            ],
                                            className="p-3 d-flex flex-column h-100",
                                        )
                                    ],
                                    className="h-100 shadow-sm border-0 bg-white rounded-3",
                                )
                            ],
                            xs=12,
                            lg=6,
                        ),
                    ],
                    className="mb-4 g-3 mx-auto align-items-stretch",
                    style={"maxWidth": "1400px"},
                ),
            ],
        ),
        # === TAB 2: EPIDEMIOLOGY CONTENT ===
        html.Div(id="epidemiology-content", className="pb-5", style={"display": "none", "paddingBottom": "80px"}, children=[
            
            # =========================================================================
            # ROW 1: SUMMARY KPIs (4 Compact Cards: 25% | 25% | 25% | 25%)
            # =========================================================================
            dbc.Row([
                # Card 1: Prevalence
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody([
                            html.Div([
                                html.Div(
                                    html.I(className="fa-solid fa-chart-line text-primary fs-4"),
                                    className="p-2 bg-primary-subtle text-primary rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                                    style={"width": "44px", "height": "44px"}
                                ),
                                html.Div([
                                    html.Div("PREVALENCE", className="hep-filter-label mb-1"),
                                    html.Div(id="epi-prevalence-total", className="hep-kpi-value mb-0"),
                                    html.Div(id="epi-prevalence-trend", className="hep-meta-text mt-1")
                                ])
                            ], className="d-flex align-items-center")
                        ], className="p-3"),
                        className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
                    ),
                    xs=12, sm=6, lg=3, className="mb-3 mb-lg-0"
                ),

                # Card 2: Incidence Rate
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody([
                            html.Div([
                                html.Div(
                                    html.I(className="fa-solid fa-virus text-danger fs-4"),
                                    className="p-2 bg-danger-subtle text-danger rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                                    style={"width": "44px", "height": "44px"}
                                ),
                                html.Div([
                                    html.Div("INCIDENCE RATE", className="hep-filter-label mb-1"),
                                    html.Div(id="epi-incidence-total", className="hep-kpi-value mb-0"),
                                    html.Div(id="epi-incidence-trend", className="hep-meta-text mt-1")
                                ])
                            ], className="d-flex align-items-center")
                        ], className="p-3"),
                        className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
                    ),
                    xs=12, sm=6, lg=3, className="mb-3 mb-lg-0"
                ),

                # Card 3: Mortality / Deaths
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody([
                            html.Div([
                                html.Div(
                                    html.I(className="fa-solid fa-heart-pulse text-warning fs-4"),
                                    className="p-2 bg-warning-subtle text-warning rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                                    style={"width": "44px", "height": "44px"}
                                ),
                                html.Div([
                                    html.Div("MORTALITY / DEATHS", className="hep-filter-label mb-1"),
                                    html.Div(id="epi-deaths-total", className="hep-kpi-value mb-0"),
                                    html.Div(id="epi-deaths-trend", className="hep-meta-text mt-1")
                                ])
                            ], className="d-flex align-items-center")
                        ], className="p-3"),
                        className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
                    ),
                    xs=12, sm=6, lg=3, className="mb-3 mb-lg-0"
                ),

                # Card 4: Genomic Surveillance Coverage
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody([
                            html.Div([
                                html.Div(
                                    html.I(className="fa-solid fa-shield-halved text-success fs-4"),
                                    className="p-2 bg-success-subtle text-success rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                                    style={"width": "44px", "height": "44px"}
                                ),
                                html.Div([
                                    html.Div("GENOMIC COVERAGE", className="hep-filter-label mb-1"),
                                    html.Div(id="epi-coverage-percent", className="hep-kpi-value mb-0"),
                                    html.Div(id="epi-coverage-status", className="hep-meta-text mt-1 fw-semibold")
                                ])
                            ], className="d-flex align-items-center")
                        ], className="p-3"),
                        className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
                    ),
                    xs=12, sm=6, lg=3, className="mb-3 mb-lg-0"
                ),
            ], className="g-3 mb-4"),

            # Hidden targets for callback outputs to maintain strict callback integrity
            html.Div([
                html.Div(id="epi-sex-ratio"),
                html.Div(id="epi-top-age-group"),
                html.Div(id="epi-age-percentage"),
                html.Div(id="epi-top-region"),
                html.Div(id="epi-region-percentage"),
            ], style={"display": "none"}),


            # =========================================================================
            # ROW 2: EPIDEMIOLOGY OVER TIME & CARE CASCADE (Side-by-Side 50% | 50%)
            # =========================================================================
            dbc.Row([
                # Left 50% (6 cols): Disease Burden Over Time with Filters
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("Disease Burden Over Time", className="hep-card-title mb-1"),
                                html.Small("Historical trajectory (1990–2025) with 2015 baseline comparison", className="hep-meta-text d-block mb-3"),
                                
                                # Interactive Filter Controls Bar
                                html.Div([
                                    html.Div([
                                        html.Label("Metric:", className="hep-filter-label me-1 mb-0"),
                                        dcc.Dropdown(
                                            id="epi-trend-metric",
                                            options=[
                                                {"label": "Prevalence", "value": "Prevalence"},
                                                {"label": "Incidence", "value": "Incidence"},
                                                {"label": "Deaths / Mortality", "value": "Deaths"}
                                            ],
                                            value="Prevalence",
                                            clearable=False,
                                            searchable=False,
                                            className="hep-control-dropdown hep-control-sm me-2",
                                        ),
                                    ], className="d-flex align-items-center me-2 mb-2 mb-sm-0"),
                                    
                                    html.Div([
                                        html.Label("Measure:", className="hep-filter-label me-1 mb-0"),
                                        dcc.Dropdown(
                                            id="epi-trend-measure",
                                            options=[
                                                {"label": "Rate", "value": "Rate"},
                                                {"label": "Number", "value": "Number"}
                                            ],
                                            value="Rate",
                                            clearable=False,
                                            searchable=False,
                                            className="hep-control-dropdown hep-control-xs me-2",
                                        ),
                                    ], className="d-flex align-items-center me-2 mb-2 mb-sm-0"),
                                    
                                    
                                ], className="d-flex align-items-center flex-wrap gap-1 mb-3 bg-light p-2 rounded-2"),
                            ]),
                            dcc.Loading(dcc.Graph(id="epi-disease-burden-trend-graph", style={"height": "340px"}), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], lg=6, className="mb-4 mb-lg-0"),

                # Right 50% (6 cols): Diagnosis & Treatment Cascade (2022)
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("Diagnosis & Treatment Cascade (2022)", className="hep-card-title mb-1"),
                                html.Small("Global care cascade: living with infection, diagnosed, and treated", className="hep-meta-text d-block mb-3"),
                            ]),
                            dcc.Loading(dcc.Graph(id="epi-cascade-bars-graph", style={"height": "395px"}), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], lg=6),
            ], className="g-3 mb-4 align-items-stretch"),


            # =========================================================================
            # ROW 3: GEOGRAPHIC COMPARISON (Side-by-Side 50% | 50%)
            # =========================================================================
            dbc.Row([
                # Left 50%: Disease Burden & Care Response Map
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    html.H5("Disease Burden & Care Response", className="hep-card-title mb-1"),
                                ], width=True),
                                dbc.Col([
                                    html.Div([
                                        dcc.RadioItems(
                                            id="epi-map-category",
                                            options=[
                                                {"label": "Disease Burden", "value": "burden"},
                                                {"label": "Program Response", "value": "response"},
                                            ],
                                            value="burden",
                                            className="hep-btn-group me-2",
                                            inputClassName="btn-check",
                                            labelClassName="btn btn-outline-secondary btn-sm",
                                        ),
                                        dcc.Dropdown(
                                            id="epi-metric-dropdown",
                                            options=[
                                                {"label": "People living with infection", "value": "livingwith_num"},
                                                {"label": "New infections", "value": "new_infections_num"},
                                                {"label": "Deaths", "value": "deaths_num"},
                                                {"label": "Prevalence %", "value": "prevalence_pct"},
                                                {"label": "Incidence rate", "value": "incidence_rate"},
                                                {"label": "Mortality rate", "value": "mortality_rate"},
                                            ],
                                            value="livingwith_num",
                                            clearable=False,
                                            searchable=False,
                                            className="hep-control-dropdown hep-control-md d-inline-block",
                                        )
                                    ], className="d-flex align-items-center flex-wrap gap-2 mb-3")
                                ], width="auto", className="d-flex align-items-center")
                            ], className="g-2 align-items-center mb-2"),
                            dcc.Loading(dcc.Graph(id="gho-burden-map", style={"height": "330px"}), type="circle")
                        ], className="p-3")
                    ], id="epi-controls-row", className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], lg=6, className="mb-4 mb-lg-0"),

                # Right 50%: Genomic Surveillance Coverage Map
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("Genomic Surveillance Coverage Map", className="hep-card-title mb-1"),
                                html.Small("Sequences available relative to estimated disease burden", className="hep-meta-text mb-3 d-block"),
                            ]),
                            dcc.Loading(dcc.Graph(id="epi-coverage-map-graph", style={"height": "330px"}), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], lg=6),
            ], className="g-3 mb-4 align-items-stretch"),


            # =========================================================================
            # ROW 4: SURVEILLANCE GAP — DECISION-SUPPORT SECTION (58% | 42%)
            # =========================================================================
            dbc.Row([
                # Left 58% (7 cols): Burden vs Sequencing Scatter Plot
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    html.H5("Burden vs. Sequencing Coverage (Surveillance Gap)", className="hep-card-title mb-1"),
                                    html.Small("Countries in bottom-right quadrant (High Burden + Low Coverage) represent critical surveillance priorities", className="hep-meta-text d-block mb-0"),
                                ], width=True),
                                dbc.Col([
                                    html.Div(
                                        id="hbv-epi-trends-controls",
                                        children=[
                                            dcc.Dropdown(
                                                id="hbv-y-metric",
                                                options=[
                                                    {"label": "Prevalence %", "value": "hbv_prevalence_pct"},
                                                    {"label": "New Infections", "value": "hbv_new_infections_num"},
                                                    {"label": "Living with HBV", "value": "hbv_livingwith_num"},
                                                ],
                                                value="hbv_prevalence_pct",
                                                clearable=False,
                                                searchable=False,
                                                className="hep-control-dropdown hep-control-sm",
                                            )
                                        ],
                                        style={"display": "flex", "alignItems": "center"}
                                    ),
                                    html.Div(
                                        id="hcv-epi-trends-controls",
                                        children=[
                                            dcc.Dropdown(
                                                id="hcv-x-metric",
                                                options=[
                                                    {"label": "Diagnosis Rate %", "value": "hcv_diagnosis_rate_pct"},
                                                    {"label": "Treatment Rate %", "value": "hcv_treatment_rate_diagnosed_pct"},
                                                ],
                                                value="hcv_diagnosis_rate_pct",
                                                clearable=False,
                                                searchable=False,
                                                className="hep-control-dropdown hep-control-sm me-2",
                                            ),
                                            dcc.Dropdown(
                                                id="hcv-y-metric",
                                                options=[
                                                    {"label": "Prevalence %", "value": "hcv_prevalence_pct"},
                                                    {"label": "New Infections", "value": "hcv_new_infections_num"},
                                                    {"label": "Living with HCV", "value": "hcv_livingwith_num"},
                                                ],
                                                value="hcv_prevalence_pct",
                                                clearable=False,
                                                searchable=False,
                                                className="hep-control-dropdown hep-control-sm",
                                            ),
                                        ],
                                        style={"display": "none", "alignItems": "center"}
                                    ),
                                ], width="auto", className="d-flex align-items-center")
                            ], className="g-2 align-items-center mb-3"),
                            
                            # Hidden toggling containers to preserve existing HBV/HCV callbacks
                            html.Div(id="hbv-epi-trends-row", children=[
                                dcc.Loading(dcc.Graph(id="hbv-epi-trend-plot", style={"height": "360px", "width": "100%"}), type="circle")
                            ]),

                            html.Div(id="hcv-epi-trends-row", children=[
                                dcc.Loading(dcc.Graph(id="hcv-epi-trend-plot", style={"height": "360px", "width": "100%"}), type="circle")
                            ]),

                            # Keep burden-coverage-scatter in DOM for fallback callbacks
                            html.Div(dcc.Graph(id="burden-coverage-scatter", style={"display": "none"}))
                        ], className="p-3 style-overflow-hidden", style={"overflow": "hidden"})
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100 style-overflow-hidden", style={"overflow": "hidden"})
                ], lg=7, className="mb-4 mb-lg-0"),

                # Right 42% (5 cols): Priority / Gap Ranking Table
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("Surveillance Priorities", className="hep-card-title mb-0"),
                                dbc.Button(
                                    "Full Ranking →",
                                    id="epi-btn-goto-priority-tab",
                                    color="link",
                                    className="p-0 text-decoration-none small fw-bold text-primary",
                                ),
                            ], className="d-flex justify-content-between align-items-center mb-3"),
                            
                            html.Div(id="hbv-epi-priority-table"),
                            html.Div(id="hcv-epi-priority-table"),
                        ], className="p-3 d-flex flex-column h-100")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], lg=5),
            ], className="g-3 mb-4 align-items-stretch"),


            # =========================================================================
            # ROW 5: WHO 2030 ELIMINATION PROGRESS (Full Width 100%)
            # =========================================================================
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.Div([
                                    html.H5("WHO 2030 Elimination Progress", className="hep-card-title mb-1"),
                                    html.Small("Historical trajectory (2015–2025), forecast to 2030, and comparison against WHO 90% reduction target", className="hep-meta-text")
                                ]),
                                html.Div([
                                    dbc.Badge([
                                        html.I(className="fa-solid fa-circle-exclamation me-1"),
                                        html.Span(id="epi-2030-progress", children="● OFF TRACK")
                                    ], color="danger", className="px-3 py-2 fs-6 rounded-pill shadow-sm")
                                ], className="mt-2 mt-sm-0")
                            ], className="d-flex justify-content-between align-items-center mb-3 flex-wrap"),
                            dcc.Loading(dcc.Graph(id="epi-forecast-summary-graph", style={"height": "350px"}), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], width=12)
            ], className="mb-4"),


            # =========================================================================
            # ROW 6: WHO GHO COUNTRY ELIMINATION COMPARISON TABLE (Full Width 100%)
            # =========================================================================
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.Div([
                                    html.H5("Country Elimination Progress & GHO Profile", className="hep-card-title d-inline mb-0"),
                                    html.Small(" Detailed epidemiological metrics, diagnosis rates, and treatment coverage by country", className="hep-meta-text d-block mt-1")
                                ]),
                                html.Div([
                                    dbc.Button("Download Data", id="btn-download-epi", color="secondary", size="sm", className="fw-bold px-3"),
                                ], className="d-flex align-items-center mt-2 mt-sm-0"),
                                dcc.Download(id="download-epi-data")
                            ], className="mb-3 d-flex justify-content-between align-items-center flex-wrap"),
                            dcc.Loading(html.Div(id="gho-data-table"), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], width=12)
            ], className="mb-4")
        ]),

        # === TAB 3: GENOTYPES CONTENT ===
        html.Div(id="genotypes-content", className="pb-5", style={"display": "none", "paddingBottom": "80px"}, children=[
            
            # =========================================================================
            # ROW 1: SUMMARY KPIs (4 Compact Cards: 25% | 25% | 25% | 25%)
            # =========================================================================
            dbc.Row([
                # Card 1: Genotypes Observed
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody([
                            html.Div([
                                html.Div(
                                    html.I(className="fa-solid fa-dna text-primary fs-4"),
                                    className="p-2 bg-primary-subtle text-primary rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                                    style={"width": "44px", "height": "44px"}
                                ),
                                html.Div([
                                    html.Div("GENOTYPES OBSERVED", className="hep-filter-label mb-1"),
                                    html.Div(id="geno-kpi-observed", className="hep-kpi-value mb-0", children="--"),
                                    html.Div("Distinct genotypes in dataset", className="hep-meta-text mt-1")
                                ])
                            ], className="d-flex align-items-center")
                        ], className="p-3"),
                        className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
                    ),
                    xs=12, sm=6, lg=3, className="mb-3 mb-lg-0"
                ),
                # Card 2: Dominant Genotype
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody([
                            html.Div([
                                html.Div(
                                    html.I(className="fa-solid fa-crown text-warning fs-4"),
                                    className="p-2 bg-warning-subtle text-warning rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                                    style={"width": "44px", "height": "44px"}
                                ),
                                html.Div([
                                    html.Div("DOMINANT GENOTYPE", className="hep-filter-label mb-1"),
                                    html.Div(id="geno-kpi-dominant", className="hep-kpi-value mb-0", children="--"),
                                    html.Div("Most represented among sequences", className="hep-meta-text mt-1")
                                ])
                            ], className="d-flex align-items-center")
                        ], className="p-3"),
                        className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
                    ),
                    xs=12, sm=6, lg=3, className="mb-3 mb-lg-0"
                ),
                # Card 3: Most Diverse Region
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody([
                            html.Div([
                                html.Div(
                                    html.I(className="fa-solid fa-chart-pie text-info fs-4"),
                                    className="p-2 bg-info-subtle text-info rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                                    style={"width": "44px", "height": "44px"}
                                ),
                                html.Div([
                                    html.Div("MOST DIVERSE REGION", className="hep-filter-label mb-1"),
                                    html.Div(id="geno-kpi-diversity", className="hep-kpi-value mb-0", children="--"),
                                    html.Div("Highest Shannon diversity index", className="hep-meta-text mt-1")
                                ])
                            ], className="d-flex align-items-center")
                        ], className="p-3"),
                        className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
                    ),
                    xs=12, sm=6, lg=3, className="mb-3 mb-lg-0"
                ),
                # Card 4: Recent Genotype Shift
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody([
                            html.Div([
                                html.Div(
                                    html.I(className="fa-solid fa-arrow-trend-up text-success fs-4"),
                                    className="p-2 bg-success-subtle text-success rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                                    style={"width": "44px", "height": "44px"}
                                ),
                                html.Div([
                                    html.Div("RECENT GENOTYPE SHIFT", className="hep-filter-label mb-1"),
                                    html.Div(id="geno-kpi-shift", className="hep-kpi-value mb-0", children="--"),
                                    html.Div("Largest recent increase (>=2020 vs <2020)", className="hep-meta-text mt-1")
                                ])
                            ], className="d-flex align-items-center")
                        ], className="p-3"),
                        className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
                    ),
                    xs=12, sm=6, lg=3, className="mb-3 mb-lg-0"
                ),
            ], className="g-3 mb-4"),


            # =========================================================================
            # ROW 2: GLOBAL GENOTYPE DISTRIBUTION MAP (100% Full Width)
            # =========================================================================
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.Div([
                                    html.H5("Global Genotype Distribution", className="hep-card-title mb-1"),
                                ]),
                                html.Div([
                                    dcc.RadioItems(
                                        id="geno-map-mode",
                                        options=[
                                            {"label": "Dominant Genotype", "value": "dominant"},
                                            {"label": "Genotype Frequency", "value": "frequency"},
                                            {"label": "Diversity Index", "value": "diversity"},
                                        ],
                                        value="dominant",
                                        className="hep-btn-group me-2",
                                        inputClassName="btn-check",
                                        labelClassName="btn btn-outline-secondary btn-sm",
                                    ),
                                    dcc.Dropdown(
                                        id="geno-map-genotype-filter",
                                        options=[{"label": "All Genotypes", "value": "ALL"}],
                                        value="ALL",
                                        clearable=False,
                                        searchable=False,
                                        className="hep-control-dropdown hep-control-xs d-inline-block",
                                    )
                                ], className="d-flex align-items-center flex-wrap gap-2 mt-2 mt-sm-0")
                            ], className="d-flex justify-content-between align-items-center mb-3 flex-wrap"),
                            dcc.Loading(dcc.Graph(id="geno-distribution-map-graph", style={"height": "480px"}), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], width=12)
            ], className="mb-4"),


            # =========================================================================
            # ROW 3: REGIONAL COMPOSITION + DIVERSITY RANKING (67% | 33%)
            # =========================================================================
            dbc.Row([
                # Left 67% (8 cols): Genotype × Region Heatmap
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    html.H5("Genotype × WHO Region Heatmap", className="hep-card-title mb-1"),
                                    html.Small("Distribution of genotypes across global WHO regions", className="hep-meta-text d-block mb-0")
                                ], width=True),
                                dbc.Col(
                                    dcc.RadioItems(
                                        id="geno-heatmap-mode",
                                        options=[
                                            {"label": "Within-region %", "value": "pct"},
                                            {"label": "Sequence count", "value": "count"},
                                        ],
                                        value="pct",
                                        className="hep-btn-group",
                                        inputClassName="btn-check",
                                        labelClassName="btn btn-outline-secondary btn-sm",
                                    ),
                                    width="auto"
                                )
                            ], className="g-2 align-items-center mb-3"),
                            dcc.Loading(dcc.Graph(id="geno-region-heatmap-graph", style={"height": "360px"}), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], lg=8, className="mb-4 mb-lg-0"),

                # Right 33% (4 cols): Genotype Diversity Ranking
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    html.H5("Genotype Diversity Ranking", className="hep-card-title mb-1"),
                                    html.Small("Shannon Diversity Index (H') by region or country", className="hep-meta-text d-block mb-0")
                                ], width=True),
                                dbc.Col(
                                    dcc.RadioItems(
                                        id="geno-diversity-scope",
                                        options=[
                                            {"label": "WHO Regions", "value": "region"},
                                            {"label": "Countries", "value": "country"},
                                        ],
                                        value="region",
                                        className="hep-btn-group",
                                        inputClassName="btn-check",
                                        labelClassName="btn btn-outline-secondary btn-sm",
                                    ),
                                    width="auto"
                                )
                            ], className="g-2 align-items-center mb-3"),
                            dcc.Loading(dcc.Graph(id="geno-diversity-ranking-graph", style={"height": "320px"}), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], lg=4),
            ], className="g-3 mb-4 align-items-stretch"),


            # =========================================================================
            # ROW 4: GENOTYPE COMPOSITION THROUGH TIME (100% Full Width)
            # =========================================================================
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5(id="line-title-main", className="hep-card-title mb-1"),
                                html.Small("Temporal trajectory of genotype proportions (1963–2026)", className="hep-meta-text")
                            ], className="mb-3"),
                            dcc.Loading(dcc.Graph(id="line-chart", style={"height": "360px"}), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], width=12),           
            ], className="mb-4"),


            # =========================================================================
            # ROW 5: EMERGING / CHANGING GENOTYPE SIGNALS (67% | 33%)
            # =========================================================================
            dbc.Row([
                # Left 67% (8 cols): Recent vs Historical Frequency Scatter
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("Recent vs. Historical Genotype Shift", className="hep-card-title mb-1"),
                                html.Small("Comparing genotype frequencies: Historical (pre-2020) vs Recent (>=2020). Points above diagonal indicate expanding genotypes.", className="hep-meta-text")
                            ], className="mb-3"),
                            dcc.Loading(dcc.Graph(id="geno-scatter-signals-graph", style={"height": "350px"}), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], lg=8, className="mb-4 mb-lg-0"),

                # Right 33% (4 cols): Recent Genotype Signals Table
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("Recent Genotype Signals", className="hep-card-title mb-1"),
                                html.Small("Genotypes with highest relative shift in recent surveillance", className="hep-meta-text")
                            ], className="mb-3"),
                            dcc.Loading(html.Div(id="geno-signals-table", style={"minHeight": "300px"}), type="circle")
                        ], className="p-3 d-flex flex-column h-100")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], lg=4),
            ], className="g-3 mb-4 align-items-stretch"),


            # =========================================================================
            # ROW 6: GENOTYPE DETAIL & COUNTRY MATRIX (100% Full Width)
            # =========================================================================
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5(id="bar-title-main", className="hep-card-title mb-3"),
                            dcc.Loading(
                                dcc.Graph(id="genotype-bar-chart", style={"height": "340px"}),
                                type="circle"
                            )
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], width=6),
                
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col(
                                    html.H5("Sequences by Country and Genotype", className="hep-card-title mb-0"),
                                    width=True
                                ),
                                dbc.Col(
                                    html.Div([
                                        html.Label("Top N:", className="hep-filter-label me-2 mb-0"),
                                        dcc.Dropdown(
                                            id="top-countries-count",
                                            options=[{"label": str(i), "value": i} for i in [10, 15, 20, 25]],
                                            value=10,
                                            clearable=False,
                                            searchable=False,
                                            className="hep-control-dropdown hep-control-xs d-inline-block",
                                        )
                                    ], className="d-flex align-items-center"),
                                    width="auto"
                                )
                            ], className="g-2 align-items-center mb-3"),
                            dcc.Loading(
                                dcc.Graph(id="country-barchart", style={"height": "340px"}),
                                type="circle"
                            )
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], width=6)
            ], className="g-3 mb-4"),
        ]),

        # === TAB 4: MUTATIONS CONTENT ===
        html.Div(id="mutations-content", className="pb-5", style={"display": "none", "paddingBottom": "80px"}, children=[
            
            # =========================================================================
            # ROW 1: SUMMARY KPIs (4 Compact Cards: 25% | 25% | 25% | 25%)
            # =========================================================================
            dbc.Row([
                # Card 1: Unique Mutations
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody([
                            html.Div([
                                html.Div(
                                    html.I(className="fa-solid fa-bug text-danger fs-4"),
                                    className="p-2 bg-danger-subtle text-danger rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                                    style={"width": "44px", "height": "44px"}
                                ),
                                html.Div([
                                    html.Div("MUTATIONS TRACKED", className="hep-filter-label mb-1"),
                                    html.Div(id="mut-kpi-unique", className="hep-kpi-value mb-0", children="--"),
                                    html.Div("Total distinct variants identified", className="hep-meta-text mt-1")
                                ])
                            ], className="d-flex align-items-center")
                        ], className="p-3"),
                        className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
                    ),
                    xs=12, sm=6, lg=3, className="mb-3 mb-lg-0"
                ),
                # Card 2: High-Frequency Mutations
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody([
                            html.Div([
                                html.Div(
                                    html.I(className="fa-solid fa-fire text-warning fs-4"),
                                    className="p-2 bg-warning-subtle text-warning rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                                    style={"width": "44px", "height": "44px"}
                                ),
                                html.Div([
                                    html.Div("HIGH-FREQ MUTATIONS", className="hep-filter-label mb-1"),
                                    html.Div(id="mut-kpi-high-freq", className="hep-kpi-value mb-0", children="--"),
                                    html.Div(">5% sequence frequency share", className="hep-meta-text mt-1")
                                ])
                            ], className="d-flex align-items-center")
                        ], className="p-3"),
                        className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
                    ),
                    xs=12, sm=6, lg=3, className="mb-3 mb-lg-0"
                ),
                # Card 3: Countries Affected
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody([
                            html.Div([
                                html.Div(
                                    html.I(className="fa-solid fa-globe text-primary fs-4"),
                                    className="p-2 bg-primary-subtle text-primary rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                                    style={"width": "44px", "height": "44px"}
                                ),
                                html.Div([
                                    html.Div("COUNTRIES AFFECTED", className="hep-filter-label mb-1"),
                                    html.Div(id="mut-kpi-countries", className="hep-kpi-value mb-0", children="--"),
                                    html.Div("Countries with detected mutations", className="hep-meta-text mt-1")
                                ])
                            ], className="d-flex align-items-center")
                        ], className="p-3"),
                        className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
                    ),
                    xs=12, sm=6, lg=3, className="mb-3 mb-lg-0"
                ),
                # Card 4: Emerging Mutation Signals
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody([
                            html.Div([
                                html.Div(
                                    html.I(className="fa-solid fa-bolt text-success fs-4"),
                                    className="p-2 bg-success-subtle text-success rounded-3 me-3 d-flex align-items-center justify-content-center flex-shrink-0",
                                    style={"width": "44px", "height": "44px"}
                                ),
                                html.Div([
                                    html.Div("EMERGING SIGNALS", className="hep-filter-label mb-1"),
                                    html.Div(id="mut-kpi-signals", className="hep-kpi-value mb-0", children="--"),
                                    html.Div("Recent increasing trend (>=2020 vs <2020)", className="hep-meta-text mt-1")
                                ])
                            ], className="d-flex align-items-center")
                        ], className="p-3"),
                        className="hep-kpi-card shadow-sm border-0 bg-white rounded-3 h-100"
                    ),
                    xs=12, sm=6, lg=3, className="mb-3 mb-lg-0"
                ),
            ], className="g-3 mb-4"),


            # =========================================================================
            # ROW 2: MUTATION EXPLORER TABLE & TOP FILTERS (100% Full Width)
            # =========================================================================
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("Mutation Explorer", className="hep-card-title mb-1"),
                                html.Small("Filter, search, and interrogate functional mutation variants across genes and drugs", className="hep-meta-text mb-3 d-block")
                            ]),

                            # Filters row
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Mutation Type:", className="hep-filter-label me-2 mb-0"),
                                    dcc.Dropdown(
                                        id="mutation-type-filter",
                                        options=[
                                            {"label": "All Types", "value": "all"},
                                            {"label": "Antiviral Resistance", "value": "antiviral_resistance"},
                                            {"label": "Vaccine Escape (HBV)", "value": "vaccine_escape"},
                                            {"label": "Substitutions of Interest", "value": "substitution_of_interest"},
                                            {"label": "No Resistance", "value": "no_resistance"}
                                        ],
                                        value="all", clearable=False, searchable=False, className="hep-control-dropdown hep-control-sm",
                                    )
                                ], width="auto", className="d-flex align-items-center me-2"),
                                dbc.Col([
                                    html.Label("Gene/Drug:", className="hep-filter-label me-2 mb-0"),
                                    dcc.Dropdown(id="mutation-category-filter", options=[],
                                                  placeholder="Select category...", className="hep-control-dropdown hep-control-sm")
                                ], width="auto", className="d-flex align-items-center me-2"),
                                dbc.Col([
                                    html.Label("Show Top:", className="hep-filter-label me-2 mb-0"),
                                    dcc.Dropdown(
                                        id="mutation-top-n",
                                        options=[{"label": str(i), "value": i} for i in [10, 20, 30, 50]],
                                        value=20, clearable=False, searchable=False, className="hep-control-dropdown hep-control-xs",
                                    )
                                ], width="auto", className="d-flex align-items-center me-2"),
                                dbc.Col([
                                    html.Label("Group By:", className="hep-filter-label me-2 mb-0"),
                                    dcc.RadioItems(
                                        id="mutation-group-by",
                                        options=[{"label": " Type", "value": "type"}, {"label": " Drug", "value": "drug"}, {"label": " Gene", "value": "gene"}],
                                        value="type", inline=True, className="hep-btn-group"
                                    )
                                ], width="auto", className="d-flex align-items-center me-2"),
                                dbc.Col([
                                    dbc.Button("Download CSV", id="btn-download-mutations", color="secondary",
                                               size="sm", className="float-end fw-bold")
                                ], width="auto", className="ms-auto")
                            ], className="g-2 mb-3 align-items-center"),

                            dcc.Loading(html.Div(id="mutation-details-table"), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], width=12)
            ], className="mb-4"),


            # =========================================================================
            # ROW 3: MUTATION DISTRIBUTION & MAP (67% | 33% Side-by-Side)
            # =========================================================================
            dbc.Row([
                # Left 67% (8 cols): Mutation Prevalence Over Time
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    html.H5("Mutation Prevalence Over Time", className="hep-card-title mb-1"),
                                    html.Small("Temporal trajectory of mutation frequency among available sequences", className="hep-meta-text d-block mb-0")
                                ], width=True),
                                dbc.Col(
                                    html.Div([
                                        html.Label("Top:", className="hep-filter-label me-2 mb-0"),
                                        dcc.Dropdown(
                                            id="top-mutations-count",
                                            options=[{"label": str(i), "value": i} for i in [5, 10, 15, 20]],
                                            value=10, clearable=False, searchable=False, className="hep-control-dropdown hep-control-xs",
                                        )
                                    ], className="d-flex align-items-center"),
                                    width="auto"
                                )
                            ], className="g-2 align-items-center mb-3"),
                            dcc.Loading(dcc.Graph(id="mutation-timeline", style={"height": "360px"}), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], lg=8, className="mb-4 mb-lg-0"),

                # Right 33% (4 cols): Geographic Distribution Map of Selected Mutation
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("Geographic Distribution", className="hep-card-title mb-1"),
                                html.Small("Country-level frequency share of mutations", className="hep-meta-text mb-3 d-block")
                            ]),
                            dcc.Loading(dcc.Graph(id="mut-geo-map-graph", style={"height": "360px"}), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], lg=4),
            ], className="g-3 mb-4 align-items-stretch"),


            # =========================================================================
            # ROW 4: PROTEIN POSITION / HOTSPOTS LOLLIPOP PLOT (100% Full Width)
            # =========================================================================
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    html.H5("Protein / Genome Position Hotspots (Lollipop Plot)", className="hep-card-title mb-1"),
                                    html.Small("Amino-acid position vs mutation frequency (≥5%) across protein targets", className="hep-meta-text d-block mb-0")
                                ], width=True),
                                dbc.Col(
                                    dcc.RadioItems(
                                        id="mut-protein-selector",
                                        options=[
                                            {"label": "All Targets", "value": "ALL"},
                                            {"label": "HBsAg", "value": "HBsAg"},
                                            {"label": "Polymerase", "value": "Polymerase"},
                                            {"label": "Core", "value": "Core"},
                                            {"label": "X", "value": "X"},
                                        ],
                                        value="ALL",
                                        className="hep-btn-group",
                                        inputClassName="btn-check",
                                        labelClassName="btn btn-outline-secondary btn-sm",
                                    ),
                                    width="auto"
                                )
                            ], className="g-2 align-items-center mb-3"),
                            dcc.Loading(dcc.Graph(id="mut-lollipop-graph", style={"height": "360px"}), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], width=12)
            ], className="mb-4"),


            # =========================================================================
            # ROW 5: FUNCTIONAL SIGNIFICANCE (67% | 33%)
            # =========================================================================
            dbc.Row([
                # Left 67% (8 cols): Functional Category Summary Bar
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("Functional Significance Summary", className="hep-card-title mb-1"),
                                html.Small("Distribution of mutations by clinical & biological impact (Drug Resistance, Vaccine Escape, Substitutions)", className="hep-meta-text mb-3 d-block")
                            ]),
                            dcc.Loading(dcc.Graph(id="mutation-frequency-chart", style={"height": "350px"}), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], lg=8, className="mb-4 mb-lg-0"),

                # Right 33% (4 cols): High-Priority Mutation Table
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("High-Priority Mutations", className="hep-card-title mb-1"),
                                html.Small("Curated antiviral resistance & vaccine escape variants", className="hep-meta-text mb-3 d-block")
                            ]),
                            dcc.Loading(html.Div(id="mut-priority-table", style={"minHeight": "300px"}), type="circle")
                        ], className="p-3 d-flex flex-column h-100")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], lg=4),
            ], className="g-3 mb-4 align-items-stretch"),


            # =========================================================================
            # ROW 6: GENOTYPE ASSOCIATION + CO-OCCURRENCE (67% | 33%)
            # =========================================================================
            dbc.Row([
                # Left 67% (8 cols): Mutation x Genotype Heatmap
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("Mutation × Genotype Association Heatmap", className="hep-card-title mb-1"),
                                html.Small("Prevalence of specific mutations across viral genotypes", className="hep-meta-text mb-3 d-block")
                            ]),
                            dcc.Loading(dcc.Graph(id="mut-genotype-heatmap-graph", style={"height": "350px"}), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], lg=8, className="mb-4 mb-lg-0"),

                # Right 33% (4 cols): Common Co-occurrences Table
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("Co-Occurring Mutations", className="hep-card-title mb-1"),
                                html.Small("Frequent mutation pairs co-detected in single genomes", className="hep-meta-text mb-3 d-block")
                            ]),
                            dcc.Loading(html.Div(id="mut-cooccurrence-table", style={"minHeight": "300px"}), type="circle")
                        ], className="p-3 d-flex flex-column h-100")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], lg=4),
            ], className="g-3 mb-4 align-items-stretch"),


            # =========================================================================
            # ROW 7: EMERGING MUTATION SIGNALS (100% Full Width)
            # =========================================================================
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("Emerging Mutation Signals", className="hep-card-title mb-1"),
                                html.Small("Historical (pre-2020) vs Recent (>=2020) frequency comparison. Points above diagonal indicate expanding resistance/escape variants.", className="hep-meta-text mb-3 d-block")
                            ]),
                            dcc.Loading(dcc.Graph(id="mut-emerging-signals-graph", style={"height": "360px"}), type="circle")
                        ], className="p-3")
                    ], className="shadow-sm border-0 bg-white rounded-3 h-100")
                ], width=12)
            ], className="mb-4")
        ]),
        
        # === TAB 5: USER SEQUENCE SUBMISSION ===
        html.Div(
            user_seq_tab_content(),
            className="mx-auto px-2",
            style={
                "maxWidth": "1400px",
                "width": "100%",
            },
        ),



        # Footer
        html.Footer([
            dbc.Container([
                dbc.Row([
                    dbc.Col([
                        html.Div([
                            html.Img(src="/assets/ceri_logo.png", 
                                    className="footer-logo mx-2",
                                    style={"height": "45px", "objectFit": "contain"}),
                            html.Img(src="/assets/CRICK.png", 
                                    className="footer-logo mx-2",
                                    style={"height": "45px", "objectFit": "contain"}),
                            html.Img(src="/assets/AHRI_logo.png", 
                                    className="footer-logo mx-2",
                                    style={"height": "45px", "objectFit": "contain"})
                        ], className="d-flex justify-content-center align-items-center flex-wrap")
                    ], width=12, className="mb-3")
                ]),
                
                dbc.Row([
                    dbc.Col([
                        html.Div([
                            html.P([
                                f"© {datetime.now().year} Hepatitis Virus Sequence Dashboard. ",
                                html.Span("All rights reserved.", className="text-muted")
                            ], className="mb-1"),
                            html.P("Developed by Derek Tshiabuila, Vagner Fonseca, Eduan Wilkinson, Tulio de Oliveira", 
                                  className="mb-1 text-muted"),
                            html.A("GitHub Repository", 
                                  href="https://github.com/CERI-KRISP/Hepatitis-Dashboard.git", 
                                  target="_blank", 
                                  className="text-decoration-none text-primary",
                                  style={"fontSize": "0.9rem"})
                        ], className="text-center")
                    ], width=12)
                ])
            ], fluid=True)
        ], className="bg-light py-4 mt-5 border-top", 
           style={"marginTop": "2rem !important"}),
    ]
    return html.Div(
        id="hep-dashboard-shell-container",
        className="hep-dashboard-shell px-4 py-3",
        children=[
            stores,
            html.Main(
                style={"width": "100%", "flex": "1"},
                children=main_children
            )
        ]
    )

    
# === STATE STORES ==============================================================
# — Year bounds & dropdown option lists (fast) —
layout = create_dashboard_layout()  # contains the year-slider

# Callbacks from user_sequence_analysis register automatically at import time.

@callback(Output("dl-toast", "is_open"), Input("btn-download-data", "n_clicks"), prevent_initial_call=True)
def _show_toast(n): return True



# — Filtered sequence dataframe store —
@callback(
    Output("filtered-store", "data"),
    Input("selected-virus", "data"),
    Input("year-range-slider", "value"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
    Input("genotype-dropdown", "value"),
)
def compute_filtered_store(virus, year_range, regions, countries, genotypes):
    data = get_data_store()  # UPDATED
    if data['hbv_data'].empty and data['hcv_data'].empty:
        return _df_to_json(pd.DataFrame())

    if virus == "HBV":
        base = data['hbv_data']
    elif virus == "HCV":
        base = data['hcv_data']
    elif virus == "HEV":
        base = data['hev_data']
    else:
        base = data['hbv_data']
        
    if base.empty:
        return _df_to_json(pd.DataFrame())

    # default bounds
    ymin, ymax = int(base["Year"].min()), int(base["Year"].max())
    y0 = int(year_range[0]) if year_range and len(year_range) == 2 else ymin
    y1 = int(year_range[1]) if year_range and len(year_range) == 2 else ymax
    
    if y0 > y1:
        y0, y1 = y1, y0

    df = base[(base["Year"] >= y0) & (base["Year"] <= y1)].copy()
    if regions:
        df = df[df["WHO_Regions"].isin(regions)]
    if countries:
        df = df[df["Country_standard"].isin(countries)]
    if genotypes:
        df = df[df["genotype"].isin(genotypes)]

    light = df[["ID", "Country_standard", "WHO_Regions", "Year", "genotype"]].copy()
    return _df_to_json(light)


# — Burden-adjusted coverage store —
@callback(
    Output("gap-store", "data"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
    Input("ihme-metric-type", "value"),
)
def update_gap_store(filtered_json, virus, ihme_metric_choice):
    return compute_gap_from_filtered(
        filtered_json=filtered_json,
        virus=virus,
        ihme_metric_choice=ihme_metric_choice,
        sex="Both",
    )
def compute_gap_from_filtered(
    filtered_json,
    virus,
    ihme_metric_choice,
    sex,   # ✅ ADD THIS
):
    seq_df = _df_from_json(filtered_json)
    data = get_data_store()

    if seq_df.empty:
        return _df_to_json(
            pd.DataFrame(columns=["Country_standard", "coverage_gap"])
        )

    gap = compute_gap_df(
        virus=(virus or "HBV"),
        filtered_seq_df=seq_df,
        ihme_df=data["ihme_df"],
        selected_years=[seq_df["Year"].min(), seq_df["Year"].max()],
        who_regions=None,
        countries=None,
        ihme_metric_choice=ihme_metric_choice or BURDEN_MEASURE_FALLBACK,
        sex=sex,                 # ✅ NOW VALID
        target_per_10k=5.0,
    )

    # tidy & drop Unknown
    if "Country_standard" in gap.columns:
        gap = gap[gap["Country_standard"] != "Unknown"]

    return _df_to_json(gap)


# — Latest IHME per country store —
@callback(
    Output("ihme-latest-store", "data"),
    Input("selected-virus", "data"),
    Input("ihme-metric-type", "value"),
    Input("year-range-slider", "value"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
)
def compute_ihme_latest_store(virus, metric, year_range, regions, countries):
    sex = "Both"
    data = get_data_store()
    
    # default bounds
    if virus == "HBV":
        base = data['hbv_data']
    elif virus == "HCV":
        base = data['hcv_data']
    elif virus == "HEV":
        base = data['hev_data']
    else:
        base = data['hbv_data']
        
    if base.empty:
        ymin, ymax = 1963, 2024
    else:
        ymin, ymax = int(base["Year"].min()), int(base["Year"].max())
        
    y0 = int(year_range[0]) if year_range and len(year_range) == 2 else ymin
    y1 = int(year_range[1]) if year_range and len(year_range) == 2 else ymax
    
    if y0 > y1:
        y0, y1 = y1, y0
        
    years = [y0, y1]
    
    df = ihme_latest_by_country(
        ihme_df=data["ihme_df"],
        virus=(virus or "HBV"),
        measure_metric=(metric or BURDEN_MEASURE_FALLBACK),
        sex=sex,  # This is critical - make sure it's not "Both" if your data doesn't support it
        regions=regions, 
        countries=countries, 
        years=years
    )
    
    return _df_to_json(df)
    

# — Indicator values —
@callback(
    Output("indicator-total", "children"),
    Output("indicator-countries", "children"),
    Output("indicator-genotypes", "children"),
    Output("indicator-years", "children"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
)
def update_indicators(filtered_json, selected_virus):
    df = _df_from_json(filtered_json)
    if df is None or df.empty:
        return "0", "0", "0", "N/A"

    total_genomes = len(df)
    unique_countries = df.get("Country_standard", pd.Series(dtype="object")).nunique()

    g = df.get("genotype", pd.Series(dtype="object")).astype("string").str.strip()
    recomb_mask = g.str.contains(r"recomb", case=False, na=False)
    base_genotype_count = g[~recomb_mask].dropna().nunique()

    date_col = None
    for c in ["Year", "year", "Release_Date", "Collection_Date", "Date"]:
        if c in df.columns:
            date_col = c
            break

    year_range_str = "1963–2025"
    if date_col:
        try:
            s = df[date_col].dropna()
            numeric_years = pd.to_numeric(s, errors="coerce").dropna()
            numeric_years = numeric_years[(numeric_years >= 1900) & (numeric_years <= 2030)]
            if not numeric_years.empty:
                year_range_str = f"{int(numeric_years.min())}–{int(numeric_years.max())}"
            else:
                dates = pd.to_datetime(s, errors="coerce").dt.year.dropna()
                dates = dates[(dates >= 1900) & (dates <= 2030)]
                if not dates.empty:
                    year_range_str = f"{int(dates.min())}–{int(dates.max())}"
        except Exception:
            pass

    return (
        f"{total_genomes:,}",
        f"{unique_countries:,}",
        f"{base_genotype_count:,}",
        year_range_str,
    )


@callback(
    Output("dashboard-page-title", "children"),
    Output("dashboard-page-subtitle", "children"),
    Input("selected-virus", "data"),
    Input("active-tab-store", "data"),         # ← REPLACES the 4 style Inputs
    Input("filtered-store", "data"),
)
def update_dashboard_heading(virus, active_tab, filtered_json):
    virus = virus or "HBV"
    virus_label = {
        "HBV": "Hepatitis B",
        "HCV": "Hepatitis C",
        "HEV": "Hepatitis E",
    }.get(virus, virus)

    section_map = {
        "mutations":    "Mutations",
        "epidemiology": "Epidemiology",
        "user-seq":     "My Sequences",
    }
    section = section_map.get(active_tab, "Overview")

    df = _df_from_json(filtered_json)
    sequence_count = len(df) if df is not None and not df.empty else 0

    data = get_data_store()
    if virus == "HBV":
        mut_df = data.get("hbv_mut", pd.DataFrame())
    elif virus == "HCV":
        mut_df = data.get("hcv_mut", pd.DataFrame())
    else:
        mut_df = data.get("hev_mut", pd.DataFrame())
    mutation_count = mut_df["mutation"].nunique() if not mut_df.empty and "mutation" in mut_df.columns else 0

    title = f"{virus_label} {section}"
    if section == "Overview":
        last_updated = "2026-08-02"
        if df is not None and not df.empty:
            for dcol in ["Release_Date", "Collection_Date", "Date", "date", "created_at"]:
                if dcol in df.columns:
                    valid_dates = pd.to_datetime(df[dcol], errors="coerce").dropna()
                    if not valid_dates.empty:
                        last_updated = valid_dates.max().strftime("%Y-%m-%d")
                        break
        subtitle = f"Last updated: {last_updated}"
    elif section == "Mutations":
        subtitle = f"{mutation_count:,} detected mutation markers · {sequence_count:,} filtered sequences"
    elif section == "Epidemiology":
        subtitle = "Burden, prevalence, incidence, deaths, and sequencing coverage"
    else:
        subtitle = "Upload or analyze your own sequences"

    return title, subtitle


@callback(
    Output("priority-insight-countries", "children"),
    Output("priority-insight-reason", "children"),
    Input("gap-store", "data"),
    Input("selected-virus", "data"),
)
def update_priority_insight(gap_json, virus):
    gap_df = _df_from_json(gap_json)
    
    if gap_df.empty or "coverage_gap" not in gap_df.columns or "Country_standard" not in gap_df.columns:
        return "No priority data available", "Apply filters to see priority countries"
    
    # Get top 3 countries with largest coverage gaps
    top_priority = gap_df.nlargest(3, "coverage_gap")[["Country_standard", "coverage_gap"]]
    
    countries_list = []
    for _, row in top_priority.iterrows():
        countries_list.append(html.Span([
            html.Strong(row["Country_standard"]),
            html.Small(f" ({int(row['coverage_gap']):,} seq needed)", className="ms-2 text-muted")
        ], className="d-block"))
    
    # Reason for priority
    total_gap = gap_df["coverage_gap"].sum()
    reason = f"Total sequencing gap: {int(total_gap):,} genomes needed"
    
    return countries_list, reason
    
@callback(
    Output("tab-epidemiology", "n_clicks", allow_duplicate=True),
    Output("tab-mutations", "n_clicks", allow_duplicate=True),
    Input("btn-quick-forecast", "n_clicks"),
    Input("btn-quick-priority", "n_clicks"),
    Input("btn-quick-timeline", "n_clicks"),
    prevent_initial_call=True
)
def navigate_from_quick_buttons(forecast_clicks, priority_clicks, timeline_clicks):
    ctx = callback_context
    if not ctx.triggered:
        return dash.no_update, dash.no_update
    
    trigger = ctx.triggered[0]
    val = trigger.get('value')
    if val is None or val == 0:
        return dash.no_update, dash.no_update
        
    button_id = trigger['prop_id'].split('.')[0]
    
    if button_id == "btn-quick-forecast":
        return 1, dash.no_update  # Navigate to epidemiology tab
    elif button_id == "btn-quick-timeline":
        return dash.no_update, 1  # Navigate to mutations tab
    elif button_id == "btn-quick-priority":
        # Could create a new tab or use existing
        return 1, dash.no_update  # Navigate to epidemiology tab for now
    
    return dash.no_update, dash.no_update

# === FIGURE CALLBACKS ==============================================================
@callback(
    [Output("mutation-section-title", "children"),
     Output("mutation-section-content", "children"),
     Output("mutation-summary-row", "style")],
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
    Input("year-range-slider", "value"),
)
def update_mutation_section(filtered_json, selected_virus, regions, countries, year_range):
    """Updates the mutation section based on selected virus and filters"""
    
    data = get_data_store()
    selected_virus = selected_virus or "HBV"
    
    # default bounds
    if selected_virus == "HBV":
        base = data['hbv_data']
    elif selected_virus == "HCV":
        base = data['hcv_data']
    else:
        base = data['hev_data']
        
    if base.empty:
        ymin, ymax = 1963, 2024
    else:
        ymin, ymax = int(base["Year"].min()), int(base["Year"].max())
        
    y0 = int(year_range[0]) if year_range and len(year_range) == 2 else ymin
    y1 = int(year_range[1]) if year_range and len(year_range) == 2 else ymax
    
    if y0 > y1:
        y0, y1 = y1, y0
        
    years = [y0, y1]
    filtered_df = _df_from_json(filtered_json)
    selected_virus = selected_virus or "HBV"
    
    # For HEV: Show a message instead of mutation data
    if selected_virus == "HEV":
        title = "Molecular Analysis (HEV)"
        content = get_hev_no_mutation_content()
        style = {"display": "block"}
        return title, content, style
    
    # For HBV and HCV: Get mutation data
    if selected_virus == "HBV":
        mutation_data = data["hbv_mut"]
        seq_data = data["hbv_data"]
        mutation_type_title = "HBV Vaccine Escape"
        virus_name = "HBV"
    else:  # HCV
        mutation_data = data["hcv_mut"]
        seq_data = data["hcv_data"]
        mutation_type_title = "HCV Key Substitutions"
        virus_name = "HCV"
    
    # Default values for empty data
    if filtered_df.empty or mutation_data.empty:
        title = f"Mutation Analysis Summary ({virus_name})"
        content = dbc.Alert(
            "No mutation data available for current filters",
            color="warning",
            className="text-center"
        )
        style = {"display": "block"}
        return title, content, style
    
    # Calculate basic stats
    total_samples = len(filtered_df)
    
    # Filter mutations to match current filtered sequences
    if "ID" in mutation_data.columns and "ID" in filtered_df.columns:
        filtered_mutations = mutation_data[mutation_data["ID"].isin(filtered_df["ID"].unique())]
    else:
        filtered_mutations = mutation_data
    
    if filtered_mutations.empty:
        title = f"Mutation Analysis Summary ({virus_name})"
        content = dbc.Alert(
            f"No mutations found in {total_samples:,} filtered samples",
            color="info",
            className="text-center"
        )
        style = {"display": "block"}
        return title, content, style
    
    # Calculate metrics
    samples_with_mutations = filtered_mutations["ID"].nunique() if "ID" in filtered_mutations.columns else 0
    mutation_coverage_pct = (samples_with_mutations / total_samples * 100) if total_samples > 0 else 0
    
    # Count by type
    type_counts = {}
    type_percentages = {}
    
    if "type" in filtered_mutations.columns:
        type_counts = filtered_mutations.groupby("type")["ID"].nunique().to_dict()
        for mut_type, count in type_counts.items():
            type_percentages[mut_type] = (count / total_samples * 100) if total_samples > 0 else 0
    
    # Antiviral resistance
    antiviral_resistance = type_counts.get("antiviral_resistance", 0)
    antiviral_percentage = (antiviral_resistance / total_samples * 100) if total_samples > 0 else 0
    
    # Vaccine escape or key substitutions
    if selected_virus == "HBV":
        vaccine_escape = type_counts.get("vaccine_escape", 0)
        escape_percentage = (vaccine_escape / total_samples * 100) if total_samples > 0 else 0
    else:
        key_substitutions = type_counts.get("substitution_of_interest", 0)
        escape_percentage = (key_substitutions / total_samples * 100) if total_samples > 0 else 0
        vaccine_escape = key_substitutions
    
    # Find most common mutation
    if "mutation" in filtered_mutations.columns:
        mutation_counts = filtered_mutations.groupby("mutation")["ID"].nunique()
        if not mutation_counts.empty:
            top_mutation = mutation_counts.idxmax()
            top_mutation_count = mutation_counts.max()
            top_mutation_percentage = (top_mutation_count / total_samples * 100) if total_samples > 0 else 0
        else:
            top_mutation = "None"
            top_mutation_percentage = 0
    else:
        top_mutation = "N/A"
        top_mutation_percentage = 0
    
    # Create pie chart
    pie_fig = create_mutation_type_pie_chart(type_counts, selected_virus)
    
    # Generate insights
    insights = generate_mutation_insights(
        total_samples=total_samples,
        samples_with_mutations=samples_with_mutations,
        type_counts=type_counts,
        type_percentages=type_percentages,
        antiviral_resistance=antiviral_resistance,
        antiviral_percentage=antiviral_percentage,
        virus=selected_virus
    )
    
    # Data range info
    if not filtered_df.empty and "Year" in filtered_df.columns:
        min_year = int(filtered_df["Year"].min())
        max_year = int(filtered_df["Year"].max())
        data_range = f"Data: {min_year}–{max_year}"
        if years and len(years) == 2:
            data_range += f" (filtered from {years[0]}–{years[1]})"
    else:
        data_range = "Year data not available"
    
    # Format text
    mutation_coverage_text = f"of {total_samples:,} samples ({mutation_coverage_pct:.1f}%)"
    resistance_percentage_text = f"({antiviral_percentage:.1f}% of samples)" if antiviral_resistance > 0 else ""
    escape_percentage_text = f"({escape_percentage:.1f}% of samples)" if vaccine_escape > 0 else ""
    top_mutation_percentage_text = f"{top_mutation_percentage:.1f}% of samples" if top_mutation != "None" and top_mutation_percentage > 0 else ""
    
    # Create the content
    title = f"Mutation Analysis Summary ({virus_name})"
    content = html.Div([
        # First row: Key mutation metrics
        dbc.Row([
            # Card 1: Total Samples with Mutations
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            html.I(className="bi bi-virus text-primary", 
                                   style={"fontSize": "2rem", "marginBottom": "10px"}),
                        ], className="text-center mb-2"),
                        html.H4(f"{samples_with_mutations:,}", className="card-title text-center mb-1"),
                        html.H6("Samples with Mutations", className="card-subtitle text-center text-muted mb-2"),
                        html.Small(mutation_coverage_text, className="text-center d-block text-success")
                    ])
                ], className="text-center h-100 border-start border-5 border-primary shadow-sm")
            ], width=3),
            
            # Card 2: Antiviral Resistance
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            html.I(className="bi bi-shield-exclamation text-danger", 
                                   style={"fontSize": "2rem", "marginBottom": "10px"}),
                        ], className="text-center mb-2"),
                        html.H4(f"{antiviral_resistance:,}", className="card-title text-center mb-1"),
                        html.H6("Antiviral Resistance", className="card-subtitle text-center text-muted mb-2"),
                        html.Small(resistance_percentage_text, className="text-center d-block")
                    ])
                ], className="text-center h-100 border-start border-5 border-danger shadow-sm")
            ], width=3),
            
            # Card 3: Vaccine Escape (HBV) / Key Substitutions (HCV)
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            html.I(className="bi bi-prescription2 text-warning", 
                                   style={"fontSize": "2rem", "marginBottom": "10px"}),
                        ], className="text-center mb-2"),
                        html.H4(f"{vaccine_escape:,}", className="card-title text-center mb-1"),
                        html.H6(mutation_type_title, className="card-subtitle text-center text-muted mb-2"),
                        html.Small(escape_percentage_text, className="text-center d-block")
                    ])
                ], className="text-center h-100 border-start border-5 border-warning shadow-sm")
            ], width=3),
            
            # Card 4: Most Common Mutation
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            html.I(className="bi bi-bar-chart-fill text-success", 
                                   style={"fontSize": "2rem", "marginBottom": "10px"}),
                        ], className="text-center mb-2"),
                        html.H4(top_mutation[:20] + ("..." if len(top_mutation) > 20 else ""), 
                               className="card-title text-center mb-1"),
                        html.H6("Most Common Mutation", className="card-subtitle text-center text-muted mb-2"),
                        html.Small(top_mutation_percentage_text, className="text-center d-block")
                    ])
                ], className="text-center h-100 border-start border-5 border-success shadow-sm")
            ], width=3),
        ], className="g-3 mb-4"),
        
        # Second row: Mutation distribution and quick insights
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H6("Mutation Types Distribution", className="card-subtitle mb-3"),
                        dcc.Graph(figure=pie_fig, config={'displayModeBar': False})
                    ])
                ], className="h-100 shadow-sm")
            ], width=6),
            
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H6("Quick Insights", className="card-subtitle mb-3"),
                        html.Ul([
                            html.Li(insights[0], className="mb-2"),
                            html.Li(insights[1], className="mb-2"),
                            html.Li(insights[2], className="mb-2"),
                        ], className="list-unstyled"),
                        html.Div([
                            html.Small("Based on current filters", className="text-muted"),
                            html.Br(),
                            html.Small(data_range, className="text-muted")
                        ], className="mt-3")
                    ])
                ], className="h-100 shadow-sm")
            ], width=6),
        ], className="g-3 mb-3"),
        
        # Third row: Action buttons
        dbc.Row([
            dbc.Col([
                dbc.ButtonGroup([
                    dbc.Button(
                        [html.I(className="bi bi-arrow-right me-2"), "View Detailed Mutation Analysis"],
                        id="hbv-btn-go-to-mutations",
                        color="primary",
                        size="lg"
                    ),
                    dbc.Button(
                        [html.I(className="bi bi-download me-2"), "Download Mutation Report"],
                        id="hbv-btn-download-mutation-report",
                        color="secondary",
                        size="lg"
                    ),
                ], className="w-100")
            ], width=12)
        ])
    ])
    
    style = {"display": "block"}
    return title, content, style

def get_hbv_mutation_content():
    """Returns the original HBV mutation summary content"""
    return html.Div([
        # First row: Key mutation metrics
        dbc.Row([
            # Card 1: Total Samples with Mutations
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            html.I(className="bi bi-virus text-primary", 
                                   style={"fontSize": "2rem", "marginBottom": "10px"}),
                        ], className="text-center mb-2"),
                        html.H4(id="mutation-samples-count", className="card-title text-center mb-1"),
                        html.H6("Samples with Mutations", className="card-subtitle text-center text-muted mb-2"),
                        html.Small(id="mutation-coverage-percent", className="text-center d-block text-success")
                    ])
                ], className="text-center h-100 border-start border-5 border-primary shadow-sm")
            ], width=3),
            
            # Card 2: Antiviral Resistance
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            html.I(className="bi bi-shield-exclamation text-danger", 
                                   style={"fontSize": "2rem", "marginBottom": "10px"}),
                        ], className="text-center mb-2"),
                        html.H4(id="antiviral-resistance-count", className="card-title text-center mb-1"),
                        html.H6("Antiviral Resistance", className="card-subtitle text-center text-muted mb-2"),
                        html.Small(id="resistance-percentage", className="text-center d-block")
                    ])
                ], className="text-center h-100 border-start border-5 border-danger shadow-sm")
            ], width=3),
            
            # Card 3: Vaccine Escape (for HBV)
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            html.I(className="bi bi-prescription2 text-warning", 
                                   style={"fontSize": "2rem", "marginBottom": "10px"}),
                        ], className="text-center mb-2"),
                        html.H4(id="vaccine-escape-count", className="card-title text-center mb-1"),
                        html.H6("Vaccine Escape", className="card-subtitle text-center text-muted mb-2"),
                        html.Small(id="escape-percentage", className="text-center d-block")
                    ])
                ], className="text-center h-100 border-start border-5 border-warning shadow-sm")
            ], width=3),
            
            # Card 4: Most Common Mutation
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            html.I(className="bi bi-bar-chart-fill text-success", 
                                   style={"fontSize": "2rem", "marginBottom": "10px"}),
                        ], className="text-center mb-2"),
                        html.H4(id="top-mutation-name", className="card-title text-center mb-1"),
                        html.H6("Most Common Mutation", className="card-subtitle text-center text-muted mb-2"),
                        html.Small(id="top-mutation-percentage", className="text-center d-block")
                    ])
                ], className="text-center h-100 border-start border-5 border-success shadow-sm")
            ], width=3),
        ], className="g-3 mb-4"),
        
        # Second row: Mutation distribution and quick insights
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H6("Mutation Types Distribution", className="card-subtitle mb-3"),
                        dcc.Loading(
                            dcc.Graph(id="mutation-type-pie", config={'displayModeBar': False}),
                            type="circle"
                        )
                    ])
                ], className="h-100 shadow-sm")
            ], width=6),
            
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H6("Quick Insights", className="card-subtitle mb-3"),
                        html.Ul([
                            html.Li(id="insight-1", className="mb-2"),
                            html.Li(id="insight-2", className="mb-2"),
                            html.Li(id="insight-3", className="mb-2"),
                        ], className="list-unstyled"),
                        html.Div([
                            html.Small("Based on current filters", className="text-muted"),
                            html.Br(),
                            html.Small(id="mutation-data-range", className="text-muted")
                        ], className="mt-3")
                    ])
                ], className="h-100 shadow-sm")
            ], width=6),
        ], className="g-3 mb-3"),
        
        # Third row: Action buttons
        dbc.Row([
            dbc.Col([
                dbc.ButtonGroup([
                    dbc.Button(
                        [html.I(className="bi bi-arrow-right me-2"), "View Detailed Mutation Analysis"],
                        id="btn-go-to-mutations",
                        color="primary",
                        size="lg"
                    ),
                    dbc.Button(
                        [html.I(className="bi bi-download me-2"), "Download Mutation Report"],
                        id="btn-download-mutation-report",
                        color="secondary",
                        size="lg"
                    ),
                ], className="w-100")
            ], width=12)
        ])
    ])

def get_hcv_mutation_content():
    """Returns HCV mutation summary content (similar to HBV but with HCV-specific titles)"""
    return html.Div([
        # Same structure as HBV but with HCV-specific IDs/titles
        dbc.Row([
            # Card 1: Total Samples with Mutations
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            html.I(className="bi bi-virus text-primary", 
                                   style={"fontSize": "2rem", "marginBottom": "10px"}),
                        ], className="text-center mb-2"),
                        html.H4(id="hcv-mutation-samples-count", className="card-title text-center mb-1"),
                        html.H6("Samples with Mutations", className="card-subtitle text-center text-muted mb-2"),
                        html.Small(id="hcv-mutation-coverage-percent", className="text-center d-block text-success")
                    ])
                ], className="text-center h-100 border-start border-5 border-primary shadow-sm")
            ], width=3),
            
            # Card 2: Antiviral Resistance
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            html.I(className="bi bi-shield-exclamation text-danger", 
                                   style={"fontSize": "2rem", "marginBottom": "10px"}),
                        ], className="text-center mb-2"),
                        html.H4(id="hcv-antiviral-resistance-count", className="card-title text-center mb-1"),
                        html.H6("Antiviral Resistance", className="card-subtitle text-center text-muted mb-2"),
                        html.Small(id="hcv-resistance-percentage", className="text-center d-block")
                    ])
                ], className="text-center h-100 border-start border-5 border-danger shadow-sm")
            ], width=3),
            
            # Card 3: Key Substitutions (for HCV - different title)
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            html.I(className="bi bi-prescription2 text-warning", 
                                   style={"fontSize": "2rem", "marginBottom": "10px"}),
                        ], className="text-center mb-2"),
                        html.H4(id="hcv-key-substitutions-count", className="card-title text-center mb-1"),
                        html.H6("Key Substitutions", className="card-subtitle text-center text-muted mb-2"),  # Different title
                        html.Small(id="hcv-substitutions-percentage", className="text-center d-block")
                    ])
                ], className="text-center h-100 border-start border-5 border-warning shadow-sm")
            ], width=3),
            
            # Card 4: Most Common Mutation
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            html.I(className="bi bi-bar-chart-fill text-success", 
                                   style={"fontSize": "2rem", "marginBottom": "10px"}),
                        ], className="text-center mb-2"),
                        html.H4(id="hcv-top-mutation-name", className="card-title text-center mb-1"),
                        html.H6("Most Common Mutation", className="card-subtitle text-center text-muted mb-2"),
                        html.Small(id="hcv-top-mutation-percentage", className="text-center d-block")
                    ])
                ], className="text-center h-100 border-start border-5 border-success shadow-sm")
            ], width=3),
        ], className="g-3 mb-4"),
        
        # Second row: Mutation distribution and quick insights
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H6("Mutation Types Distribution", className="card-subtitle mb-3"),
                        dcc.Loading(
                            dcc.Graph(id="hcv-mutation-type-pie", config={'displayModeBar': False}),
                            type="circle"
                        )
                    ])
                ], className="h-100 shadow-sm")
            ], width=6),
            
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H6("Quick Insights", className="card-subtitle mb-3"),
                        html.Ul([
                            html.Li(id="hcv-insight-1", className="mb-2"),
                            html.Li(id="hcv-insight-2", className="mb-2"),
                            html.Li(id="hcv-insight-3", className="mb-2"),
                        ], className="list-unstyled"),
                        html.Div([
                            html.Small("Based on current filters", className="text-muted"),
                            html.Br(),
                            html.Small(id="hcv-mutation-data-range", className="text-muted")
                        ], className="mt-3")
                    ])
                ], className="h-100 shadow-sm")
            ], width=6),
        ], className="g-3 mb-3"),
        
        # Third row: Action buttons
        dbc.Row([
            dbc.Col([
                dbc.ButtonGroup([
                    dbc.Button(
                        [html.I(className="bi bi-arrow-right me-2"), "View Detailed Mutation Analysis"],
                        id="hcv-btn-go-to-mutations",
                        color="primary",
                        size="lg"
                    ),
                    dbc.Button(
                        [html.I(className="bi bi-download me-2"), "Download Mutation Report"],
                        id="hcv-btn-download-mutation-report",
                        color="secondary",
                        size="lg"
                    ),
                ], className="w-100")
            ], width=12)
        ])
    ])

def get_hev_no_mutation_content():
    """Returns a message when mutation data is not available for HEV"""
    return dbc.Alert(
        [
            html.H4("Mutation Analysis Not Available", className="alert-heading mb-3"),
            html.P([
                "Detailed mutation analysis is not available for HEV in this dashboard. ",
                html.B("HEV typically has limited therapeutic resistance mutations "),
                "compared to HBV and HCV due to differences in treatment approaches."
            ], className="mb-3"),
            html.Hr(),
            html.P("What you can explore instead:", className="fw-bold mb-2"),
            html.Ul([
                html.Li("Epidemiology and burden data for HEV"),
                html.Li("Genotype distribution and geographical patterns"),
                html.Li("Sequencing coverage and gap analysis"),
                html.Li("Outbreak detection patterns"),
            ], className="mb-3"),
            html.Div([
                html.Small("HEV is primarily a zoonotic virus with genotypes 1-4 having different epidemiological patterns.", 
                          className="text-muted")
            ])
        ],
        color="info",
        className="text-center py-4"
    )

# - Line trend -
@callback(
    Output("line-chart", "figure"),
    Output("line-title-main", "children"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
)
def render_line(filtered_json, virus):
    df = _df_from_json(filtered_json)
    if df.empty:
        empty = go.Figure()
        empty.update_layout(title="No data", xaxis={"visible": True}, yaxis={"visible": True})
        return empty, "No Data"
    selected_virus = virus or "HBV"
    return make_line_trend(df, selected_virus), f"{selected_virus.upper()} Whole Genomes Per Year"


# - genotype bar -
@callback(
    Output("genotype-bar-chart", "figure"),
    Output("bar-title-main", "children"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
    Input("display-mode", "value"),
)
def render_genotype_bar(filtered_json, virus, display_mode):
    df = _df_from_json(filtered_json)
    if df.empty:
        empty = go.Figure()
        empty.update_layout(title="No data", xaxis={"visible": False}, yaxis={"visible": False})
        return empty, "No Data"
    
    data = get_data_store()
    fig = make_genotype_bar(
        filtered_df=df,
        population_df=data["population_df"],
        selected_virus=(virus or "HBV").upper(),
        display_mode=(display_mode or "raw"),
    )
    return fig, f"Total {(virus or 'HBV').upper()} Sequences by genotype"

# - Map callback with three distinct modes -
@callback(
    Output("genotype-map", "figure"),
    Output("map-title-main", "children"),
    Output("map-title-sub", "children"),
    Output("epidemiology-controls", "style"),  # Show/hide IHME metric selector
    Input("filtered-store", "data"),
    Input("gap-store", "data"), 
    Input("ihme-latest-store", "data"),
    Input("selected-virus", "data"),
    Input("display-mode", "value"),         # raw / PerMillion (for sequences mode only)
    Input("map-mode", "value"),             # sequences / coverage / epidemiology
    Input("ihme-metric-type", "value"),     # For epidemiology mode
    State("continent-dropdown", "value"),
    State("country-dropdown", "value"),
)
def render_map(filtered_json, gap_json, ihme_json, virus, display_mode, map_mode, 
               ihme_metric, regions, countries):
    selected_virus = (virus or "HBV")
    filtered = _df_from_json(filtered_json)
    gap = _df_from_json(gap_json)
    ihme = _df_from_json(ihme_json)
    data = get_data_store()
    
    # Show/hide epidemiology controls based on map mode - DEFINE THIS AT THE START
    epi_controls_style = {"display": "block"} if map_mode == "epidemiology" else {"display": "none"}

    # MODE 1: Coverage map
    if map_mode == "coverage":
        # Get the coverage dataframe based on virus type
        if selected_virus == "HBV":
            cov = data.get("cov_hbv", pd.DataFrame())
        elif selected_virus == "HCV":
            cov = data.get("cov_hcv", pd.DataFrame())
        else:  # HEV
            cov = data.get("cov_hev", pd.DataFrame())
        
        # FIXED: Debug prints using the correct variable name 'cov' not 'cov_df'
        print(f"\n{'='*50}")
        print(f"COVERAGE MAP DEBUG - Starting")
        print(f"COVERAGE MAP DEBUG - cov type: {type(cov)}")
        if cov is not None and not cov.empty:
            print(f"COVERAGE MAP DEBUG - cov shape: {cov.shape}")
            print(f"COVERAGE MAP DEBUG - cov columns: {cov.columns.tolist()}")
            print(f"COVERAGE MAP DEBUG - cov first 2 rows:\n{cov.head(2)}")
        else:
            print(f"COVERAGE MAP DEBUG - cov is empty or None!")
        print(f"{'='*50}\n")
    
        fig = create_coverage_map(
            cov,
            data["coord_lookup"],
            data.get("coords"),
            virus_type=selected_virus,
            who_regions=regions,
            countries=countries
        )
    
        title = f"{selected_virus} Burden-adjusted Sequencing Coverage"
        subtitle = "Coverage = sequences / estimated infections"
    
        return fig, title, subtitle, epi_controls_style

    # MODE 2: Epidemiology map (IHME data)
    elif map_mode == "epidemiology":
        print(f"EPIDEMIOLOGY MAP DEBUG - Starting")
        print(f"EPIDEMIOLOGY MAP DEBUG - ihme shape: {ihme.shape if hasattr(ihme, 'shape') else 'No shape'}")
        print(f"EPIDEMIOLOGY MAP DEBUG - ihme columns: {ihme.columns.tolist() if not ihme.empty else 'Empty'}")
        print(f"EPIDEMIOLOGY MAP DEBUG - ihme sample:")
        print(ihme.head() if not ihme.empty else 'Empty dataframe')
        
        if ihme.empty:
            print(f"EPIDEMIOLOGY MAP DEBUG - ihme is empty!")
            fig = _empty_world("No IHME epidemiology data available for current filters")
            title = f"{selected_virus} Epidemiology Map"
            subtitle = f"({ihme_metric}) - No data"
            return fig, title, subtitle, epi_controls_style
        
        # Parse metric for title
        try:
            measure, metric_type = (ihme_metric or "Prevalence|Number").split("|")
            metric_display = f"{measure} ({metric_type})"
        except:
            measure, metric_type = "Prevalence", "Number"
            metric_display = "Prevalence (Number)"
        
        print(f"EPIDEMIOLOGY MAP DEBUG - Creating map with measure: {measure}")
        
        # Create the epidemiology map
        fig = create_world_map(
            country_data=ihme,
            country_genotype_counts=pd.DataFrame(),  # No genotype markers for epidemiology
            coord_lookup=data["coord_lookup"],
            virus_type=selected_virus,
            display_mode="ihme",
            map_title=measure  # Pass the measure for colorbar title
        )
        
        title = f"{selected_virus} Epidemiology Map"
        subtitle = f"{metric_display}"
        return fig, title, subtitle, epi_controls_style

    # MODE 3: Sequences map (default)
    else:
        # Prepare data for the sequence map
        if display_mode == "PerMillion":
            # Calculate per million values
            if not filtered.empty and "Population" in data and not data["population_df"].empty:
                pop_data = data["population_df"]
                country_data = filtered.groupby("Country_standard").size().reset_index(name="count")
                country_data = country_data.merge(
                    pop_data[["Country_standard", "Population"]].drop_duplicates(),
                    on="Country_standard",
                    how="left"
                )
                country_data["Metric_raw"] = (country_data["count"] / country_data["Population"]) * 1_000_000
            else:
                country_data = filtered.groupby("Country_standard").size().reset_index(name="Metric_raw")
        else:
            # Raw counts
            country_data = filtered.groupby("Country_standard").size().reset_index(name="Metric_raw")
        
        # Prepare genotype counts for sequence map
        country_genotype_counts = filtered.groupby(["Country_standard", "genotype"]).size().reset_index(name="Count")
        
        # Create the sequence map
        fig = create_world_map(
            country_data,
            country_genotype_counts,
            data["coord_lookup"],
            virus_type=selected_virus,
            display_mode=display_mode or "raw"
        )
        
        # Set titles based on display mode
        if display_mode == "PerMillion":
            title = f"{selected_virus} whole-genome sequence map"
            subtitle = "Sequence count per million by country"
        else:
            title = f"{selected_virus} whole-genome sequence map" 
            subtitle = "Sequence count by country"
        
        return fig, title, subtitle, epi_controls_style
    
# - Mutation barplot -
@callback(
    Output("mutation-barplot", "figure"),
    Output("mutation-title", "children"),
    Input("filtered-store", "data"),
    Input("mutation-filter-dropdown", "value"),
    Input("selected-virus", "data"),
)
def render_mutation_bar(filtered_json, selected_filter, virus):
    seq_df = _df_from_json(filtered_json)
    selected_virus = (virus or "HBV").upper()
    data = get_data_store()  # UPDATED

    # Total sequences (denominator)
    total_sequences = len(seq_df)

    # Choose mutation table + full source sequences
    if selected_virus == "HBV":
        mut = data["hbv_mut"].copy()
        seq_source = data["hbv_data"]
        facet_col = "drug"  # Use 'drug' for HBV
    else:
        mut = data["hcv_mut"].copy()
        seq_source = data["hcv_data"]
        facet_col = "gene"  # Use 'gene' for HCV
        
    # For HCV, also allow 'drug' as alternative facet column
    if selected_virus == "HCV" and "drug" in mut.columns and facet_col not in mut.columns:
        facet_col = "drug"

    # Enrich mutation rows with Country/Region/Year/genotype
    enriched = _enrich_mutation_df(mut, seq_source)

    # Base span (ignore year filter): match Region/Country/genotype only
    if not seq_df.empty:
        keys_base = ["Country_standard", "WHO_Regions", "genotype"]
        base_for_span = enriched.merge(seq_df[keys_base].drop_duplicates(), on=keys_base, how="inner")
    else:
        base_for_span = enriched

    # Actual overlap for current selection (includes Year)
    if not seq_df.empty:
        keys_full = ["Country_standard", "WHO_Regions", "Year", "genotype"]
        overlap = enriched.merge(seq_df[keys_full].drop_duplicates(), on=keys_full, how="inner")
    else:
        overlap = enriched

    # Apply filter by drug/gene
    if selected_filter:
        vals = [str(s).strip().lower() for s in selected_filter]
        if facet_col in enriched.columns:
            enriched = enriched[enriched[facet_col].astype(str).str.strip().str.lower().isin(vals)]

    # Compute years_range for the *current* selection
    if not seq_df.empty and "Year" in seq_df.columns and seq_df["Year"].notna().any():
        y0, y1 = int(np.nanmin(seq_df["Year"])), int(np.nanmax(seq_df["Year"]))
        years_range = [y0, y1]
    else:
        years_range = None

    # Compute data_span from base_for_span (for "does years actually narrow the span?")
    if not base_for_span.empty and base_for_span["Year"].notna().any():
        ymin, ymax = int(np.nanmin(base_for_span["Year"])), int(np.nanmax(base_for_span["Year"]))
        data_span = [ymin, ymax]
    else:
        data_span = None

    # Other filters active? (is current selection a strict subset of all seqs)
    try:
        all_ids = set(seq_source["ID"]) if "ID" in seq_source.columns else set()
        sel_ids = set(seq_df["ID"]) if "ID" in seq_df.columns else set()
        other_active = bool(all_ids and len(sel_ids) < len(all_ids))
    except Exception:
        other_active = False

    # Build the figure + title in one place
    fig, title = make_mutation_bar(
        mutation_df=overlap,
        total_sequences=total_sequences,
        selected_virus=selected_virus,
        selected_filter=selected_filter,
        years_range=years_range,
        data_span=data_span,
        other_filters_active=other_active,
    )

    return fig, title

# === TAB NAVIGATION CALLBACKS ===
# === TAB NAVIGATION CALLBACKS ===
@callback(
    Output("tab-overview", "className"),
    Output("tab-epidemiology", "className"),
    Output("tab-genotypes", "className"),
    Output("tab-mutations", "className"),
    Output("tab-user-seq", "className"),
    Output("tab-mutations", "disabled"),
    Output("active-tab-store", "data"),
    
    # Inputs
    Input("tab-overview", "n_clicks"),
    Input("tab-epidemiology", "n_clicks"),
    Input("tab-genotypes", "n_clicks"),
    Input("tab-mutations", "n_clicks"),
    Input("tab-user-seq", "n_clicks"),
    Input("btn-goto-priority-tab", "n_clicks"),
    Input("selected-virus", "data"),
    State("url", "pathname"),
    prevent_initial_call=True
)
def update_tab_highlights(overview_clicks, epidemiology_clicks, genotypes_clicks,
                          mutations_clicks, user_seq_clicks, priority_goto_clicks, selected_virus, pathname):
    ctx = callback_context
    mutations_disabled = selected_virus == "HEV" if selected_virus else False

    if pathname not in ["/", "/dashboard"]:
        return ("hep-section-tab-btn", "hep-section-tab-btn", "hep-section-tab-btn",
                "hep-section-tab-btn" + (" disabled" if mutations_disabled else ""),
                "hep-section-tab-btn", mutations_disabled, "overview")

    trigger_id = ctx.triggered_id
    active_tab = "overview"
    if trigger_id == "tab-overview":
        active_tab = "overview"
    elif trigger_id in ["tab-epidemiology", "btn-goto-priority-tab"]:
        active_tab = "epidemiology"
    elif trigger_id == "tab-genotypes":
        active_tab = "genotypes"
    elif trigger_id == "tab-mutations" and not mutations_disabled:
        active_tab = "mutations"
    elif trigger_id == "tab-user-seq":
        active_tab = "user-seq"
    else:
        active_tab = "overview"

    classes = []
    for t in ["overview", "epidemiology", "genotypes", "mutations", "user-seq"]:
        if t == "mutations" and mutations_disabled:
            classes.append("hep-section-tab-btn disabled")
        elif t == active_tab:
            classes.append("hep-section-tab-btn active")
        else:
            classes.append("hep-section-tab-btn")

    return (*classes, mutations_disabled, active_tab)


# Client-side tab switching callback
dash.clientside_callback(
    """
    function(activeTab, pathname, overviewId) {
        if (pathname !== "/" && pathname !== "/dashboard") {
            return [
                {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "none"}
            ];
        }
        
        const showBlock = {"display": "block"};
        const hide = {"display": "none"};
        const active = activeTab || "overview";
        
        return [
            active === "overview" ? showBlock : hide,
            active === "epidemiology" ? showBlock : hide,
            active === "genotypes" ? showBlock : hide,
            active === "mutations" ? showBlock : hide,
            active === "user-seq" ? showBlock : hide
        ];
    }
    """,
    Output("overview-content", "style"),
    Output("epidemiology-content", "style"),
    Output("genotypes-content", "style"),
    Output("mutations-content", "style"),
    Output("user-seq-content", "style"),
    Input("active-tab-store", "data"),
    Input("url", "pathname"),
    Input("overview-content", "id"),
    prevent_initial_call=True
)




def create_mutation_type_pie_chart(type_counts, virus):
    """Create a pie chart showing mutation type distribution"""
    if not type_counts:
        return _empty_pie_chart("No mutation types")
    
    # Prepare data for pie chart
    labels = []
    values = []
    colors = []
    
    # Color mapping for different mutation types
    color_map = {
        "antiviral_resistance": "#e41a1c",  # Red
        "vaccine_escape": "#377eb8",  # Blue
        "substitution_of_interest": "#4daf4a",  # Green
        "no_resistance": "#999999",  # Gray
    }
    
    for mut_type, count in type_counts.items():
        if count > 0:
            labels.append(mut_type.replace("_", " ").title())
            values.append(count)
            colors.append(color_map.get(mut_type, "#ff7f00"))
    
    if not values:
        return _empty_pie_chart("No mutations found")
    
    # Create pie chart
    fig = go.Figure(data=[go.Pie(
        labels=labels,
        values=values,
        hole=0.4,
        marker=dict(colors=colors),
        textinfo='label+percent',
        textposition='inside',
        hovertemplate='<b>%{label}</b><br>Count: %{value}<br>Percentage: %{percent}<extra></extra>'
    )])
    
    fig.update_layout(
        height=200,
        margin=dict(t=10, b=10, l=10, r=10),
        showlegend=False,
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
    )
    
    return fig

def _empty_pie_chart(message):
    fig = go.Figure()
    fig.update_layout(
        xaxis={"visible": False},
        yaxis={"visible": False},
        annotations=[{
            "text": message,
            "xref": "paper", "yref": "paper",
            "x": 0.5, "y": 0.5, "showarrow": False,
            "font": {"size": 12, "color": "gray"}
        }],
        height=200,
        margin=dict(t=10, b=10, l=10, r=10),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
    )
    return fig

def generate_mutation_insights(total_samples, samples_with_mutations, type_counts, 
                              type_percentages, antiviral_resistance, antiviral_percentage, virus):
    insights = []
    
    # Insight 1: Mutation prevalence
    if total_samples > 0:
        mutation_prevalence = (samples_with_mutations / total_samples) * 100
        if mutation_prevalence > 50:
            insights.append(html.Span([
                html.I(className="bi bi-exclamation-triangle-fill text-warning me-1"),
                f"High mutation prevalence ({mutation_prevalence:.1f}%) detected"
            ]))
        elif mutation_prevalence > 20:
            insights.append(html.Span([
                html.I(className="bi bi-info-circle-fill text-info me-1"),
                f"Moderate mutation prevalence ({mutation_prevalence:.1f}%)"
            ]))
        else:
            insights.append(html.Span([
                html.I(className="bi bi-check-circle-fill text-success me-1"),
                f"Low mutation prevalence ({mutation_prevalence:.1f}%)"
            ]))
    else:
        insights.append(html.Span("No samples available for analysis", className="text-muted"))
    
    # Insight 2: Antiviral resistance
    if antiviral_resistance > 0:
        if antiviral_percentage > 10:
            insights.append(html.Span([
                html.I(className="bi bi-shield-exclamation text-danger me-1"),
                f"Significant antiviral resistance ({antiviral_percentage:.1f}%)"
            ]))
        elif antiviral_percentage > 5:
            insights.append(html.Span([
                html.I(className="bi bi-shield text-warning me-1"),
                f"Moderate antiviral resistance ({antiviral_percentage:.1f}%)"
            ]))
        else:
            insights.append(html.Span([
                html.I(className="bi bi-shield-check text-success me-1"),
                f"Low antiviral resistance ({antiviral_percentage:.1f}%)"
            ]))
    else:
        insights.append(html.Span([
            html.I(className="bi bi-shield-check text-success me-1"),
            "No antiviral resistance detected"
        ]))
    
    # Insight 3: Mutation type distribution
    if type_counts:
        dominant_type = max(type_counts.items(), key=lambda x: x[1])[0] if type_counts else None
        if dominant_type:
            dominant_percentage = type_percentages.get(dominant_type, 0)
            type_display = dominant_type.replace("_", " ").title()
            
            if dominant_type == "antiviral_resistance" and dominant_percentage > 5:
                insights.append(html.Span([
                    html.I(className="bi bi-activity text-danger me-1"),
                    f"Dominant: {type_display} ({dominant_percentage:.1f}%)"
                ]))
            elif dominant_type == "vaccine_escape" and virus == "HBV":
                insights.append(html.Span([
                    html.I(className="bi bi-exclamation-diamond text-warning me-1"),
                    f"Vaccine escape mutations present"
                ]))
            else:
                insights.append(html.Span([
                    html.I(className="bi bi-clipboard-data text-info me-1"),
                    f"Most common: {type_display}"
                ]))
        else:
            insights.append(html.Span("No dominant mutation type", className="text-muted"))
    else:
        insights.append(html.Span("No mutation type data available", className="text-muted"))
    
    # Ensure we always return exactly 3 insights
    while len(insights) < 3:
        insights.append(html.Span("No insight available", className="text-muted"))
    
    return insights[:3]


@callback(
    Output("mutation-category-filter", "options"),
    Input("selected-virus", "data"),
    Input("mutation-type-filter", "value"),
)
def update_mutation_category_options(virus, mutation_type):
    data = get_data_store()
    selected_virus = virus or "HBV"
    
    # Get mutation data
    mutation_data = data["hbv_mut"] if selected_virus == "HBV" else data["hcv_mut"]
    
    if mutation_data.empty:
        return []
    
    # Filter by mutation type if specified
    if mutation_type != "all" and "type" in mutation_data.columns:
        mutation_data = mutation_data[mutation_data["type"] == mutation_type]
    
    # For HBV: use drug column; for HCV: use gene column
    if selected_virus == "HBV" and "drug" in mutation_data.columns:
        categories = sorted(mutation_data["drug"].dropna().unique())
        label_prefix = "Drug: "
    elif selected_virus == "HCV" and "gene" in mutation_data.columns:
        categories = sorted(mutation_data["gene"].dropna().unique())
        label_prefix = "Gene: "
    else:
        return []
    
    # Create options
    options = [{"label": "All Categories", "value": "all"}]
    options.extend([
        {"label": f"{label_prefix}{cat}", "value": cat}
        for cat in categories
    ])
    
    return options


@callback(
    Output("mutation-frequency-chart", "figure"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
    Input("mutation-type-filter", "value"),
    Input("mutation-category-filter", "value"),
    Input("mutation-top-n", "value"),
)
def update_mutation_frequency_chart(filtered_json, virus, mutation_type, category, top_n):
    data = get_data_store()
    selected_virus = virus or "HBV"
    
    # Get filtered data
    filtered_df = _df_from_json(filtered_json)
    
    # Get mutation data
    mutation_data = data["hbv_mut"] if selected_virus == "HBV" else data["hcv_mut"]
    
    if mutation_data.empty or filtered_df.empty:
        return _empty_plot("No mutation data available")
    
    # Filter mutations to match current filtered sequences
    if "ID" in mutation_data.columns and "ID" in filtered_df.columns:
        filtered_mutations = mutation_data[mutation_data["ID"].isin(filtered_df["ID"].unique())]
    else:
        filtered_mutations = mutation_data
    
    # Apply additional filters
    if mutation_type != "all" and "type" in filtered_mutations.columns:
        filtered_mutations = filtered_mutations[filtered_mutations["type"] == mutation_type]
    
    if category != "all" and category:
        if selected_virus == "HBV" and "drug" in filtered_mutations.columns:
            filtered_mutations = filtered_mutations[filtered_mutations["drug"] == category]
        elif selected_virus == "HCV" and "gene" in filtered_mutations.columns:
            filtered_mutations = filtered_mutations[filtered_mutations["gene"] == category]
    
    if filtered_mutations.empty:
        return _empty_plot("No mutations match the selected filters")
    
    # Count mutations
    mutation_counts = (filtered_mutations.groupby("mutation")["ID"]
                      .nunique()
                      .reset_index()
                      .rename(columns={"ID": "count"}))
    
    # Calculate percentage
    total_samples = filtered_df["ID"].nunique() if "ID" in filtered_df.columns else len(filtered_df)
    mutation_counts["percentage"] = (mutation_counts["count"] / total_samples * 100).round(2)
    
    # Sort and take top N
    mutation_counts = mutation_counts.sort_values("percentage", ascending=False).head(top_n)
    
    # Create bar chart
    fig = px.bar(
        mutation_counts,
        x="mutation",
        y="percentage",
        labels={"percentage": "Samples with Mutation (%)", "mutation": "Mutation"},
        color="percentage",
        color_continuous_scale="Viridis",
        title=f"{selected_virus} Top {len(mutation_counts)} Mutations"
    )
    
    fig.update_traces(
        marker_cornerradius=6,
        hovertemplate="<b>%{x}</b><br>Samples: %{y:.1f}%<br>Count: %{customdata}<extra></extra>",
        customdata=mutation_counts["count"]
    )
    
    fig.update_layout(
        height=400,
        xaxis_tickangle=-45,
        coloraxis_showscale=False,
        
        #paper_bgcolor="white"
    )
    
    return fig



@callback(
    Output("mutation-details-table", "children"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
    Input("mutation-type-filter", "value"),
    Input("mutation-category-filter", "value"),
)
def update_mutation_details_table(filtered_json, virus, mutation_type, category):
    data = get_data_store()
    selected_virus = virus or "HBV"
    
    # Get filtered data
    filtered_df = _df_from_json(filtered_json)
    
    # Get mutation data
    mutation_data = data["hbv_mut"] if selected_virus == "HBV" else data["hcv_mut"]
    
    if mutation_data.empty or filtered_df.empty:
        return html.P("No mutation data available", className="text-muted")
    
    # Filter mutations to match current filtered sequences
    if "ID" in mutation_data.columns and "ID" in filtered_df.columns:
        filtered_mutations = mutation_data[mutation_data["ID"].isin(filtered_df["ID"].unique())]
    else:
        filtered_mutations = mutation_data
    
    # Apply additional filters
    if mutation_type != "all" and "type" in filtered_mutations.columns:
        filtered_mutations = filtered_mutations[filtered_mutations["type"] == mutation_type]
    
    if category != "all" and category:
        if selected_virus == "HBV" and "drug" in filtered_mutations.columns:
            filtered_mutations = filtered_mutations[filtered_mutations["drug"] == category]
        elif selected_virus == "HCV" and "gene" in filtered_mutations.columns:
            filtered_mutations = filtered_mutations[filtered_mutations["gene"] == category]
    
    if filtered_mutations.empty:
        return html.P("No mutations match the selected filters", className="text-muted")
    
    # Select columns to display
    columns_to_show = ["ID", "mutation", "type", "detected"]
    if selected_virus == "HBV" and "drug" in filtered_mutations.columns:
        columns_to_show.append("drug")
    elif selected_virus == "HCV" and "gene" in filtered_mutations.columns:
        columns_to_show.append("gene")
    
    # Keep only existing columns
    columns_to_show = [col for col in columns_to_show if col in filtered_mutations.columns]
    
    # Create DataTable with elevated design system styling
    table = dash.dash_table.DataTable(
        data=filtered_mutations[columns_to_show].to_dict('records'),
        columns=[{"name": col.replace("_", " ").title(), "id": col} for col in columns_to_show],
        page_size=10,
        style_table={'overflowX': 'auto', 'borderRadius': '10px'},
        style_cell=TABLE_CELL_STYLE,
        style_header=TABLE_HEADER_STYLE,
        style_data_conditional=[
            TABLE_ODD_ROW_STYLE,
            {
                "if": {
                    "column_id": "type",
                    "filter_query": '{type} eq "antiviral_resistance"'
                },
                "backgroundColor": "#FFE4E6",
                "color": "#E11D48",
                "fontWeight": "700",
            },
            {
                "if": {
                    "column_id": "type",
                    "filter_query": '{type} eq "vaccine_escape"'
                },
                "backgroundColor": "#FEF3C7",
                "color": "#D97706",
                "fontWeight": "700",
            },
            {
                "if": {
                    "column_id": "type",
                    "filter_query": '{type} eq "substitution_of_interest"'
                },
                "backgroundColor": "#DBEAFE",
                "color": "#2563EB",
                "fontWeight": "700",
            },
        ],
        filter_action="native",
        sort_action="native",
        export_format="csv"
    )
    
    return table

@callback(
    Output("download-mutation-report", "data"),
    Input("btn-download-mutation-report", "n_clicks"),
    State("filtered-store", "data"),
    State("selected-virus", "data"),
    State("continent-dropdown", "value"),
    State("country-dropdown", "value"),
    State("year-range-slider", "value"),
    prevent_initial_call=True
)
def download_mutation_report(n_clicks, filtered_json, virus, regions, countries, year_range):
    if not n_clicks:
        return dash.no_update
    
    data = get_data_store()
    selected_virus = virus or "HBV"
    
    # default bounds
    if selected_virus == "HBV":
        base = data['hbv_data']
    elif selected_virus == "HCV":
        base = data['hcv_data']
    else:
        base = data['hev_data']
        
    if base.empty:
        ymin, ymax = 1963, 2024
    else:
        ymin, ymax = int(base["Year"].min()), int(base["Year"].max())
        
    y0 = int(year_range[0]) if year_range and len(year_range) == 2 else ymin
    y1 = int(year_range[1]) if year_range and len(year_range) == 2 else ymax
    
    if y0 > y1:
        y0, y1 = y1, y0
        
    years = [y0, y1]
    
    # Get filtered data
    filtered_df = _df_from_json(filtered_json)
    
    # Get mutation data
    if selected_virus == "HBV":
        mutation_data = data["hbv_mut"]
    else:
        mutation_data = data["hcv_mut"]
    
    # Filter mutations
    if "ID" in mutation_data.columns and "ID" in filtered_df.columns:
        filtered_mutations = mutation_data[mutation_data["ID"].isin(filtered_df["ID"].unique())]
    else:
        filtered_mutations = mutation_data
    
    # Create comprehensive report
    report_lines = []
    
    # Header
    report_lines.append(f"{selected_virus} Mutation Analysis Report")
    report_lines.append("=" * 50)
    report_lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report_lines.append(f"Virus: {selected_virus}")
    report_lines.append(f"Regions: {regions if regions else 'All'}")
    report_lines.append(f"Countries: {countries if countries else 'All'}")
    report_lines.append(f"Year Range: {f'{years[0]} - {years[1]}' if years else 'All'}")
    report_lines.append("")
    
    # Summary statistics
    total_samples = len(filtered_df)
    samples_with_mutations = filtered_mutations["ID"].nunique() if "ID" in filtered_mutations.columns else 0
    
    report_lines.append("SUMMARY STATISTICS")
    report_lines.append("-" * 30)
    report_lines.append(f"Total Samples: {total_samples:,}")
    report_lines.append(f"Samples with Mutations: {samples_with_mutations:,}")
    
    if total_samples > 0:
        mutation_coverage = (samples_with_mutations / total_samples) * 100
        report_lines.append(f"Mutation Coverage: {mutation_coverage:.1f}%")
    
    # Mutation type breakdown
    if "type" in filtered_mutations.columns:
        report_lines.append("")
        report_lines.append("MUTATION TYPE DISTRIBUTION")
        report_lines.append("-" * 30)
        
        type_counts = filtered_mutations.groupby("type")["ID"].nunique()
        for mut_type, count in type_counts.items():
            percentage = (count / total_samples * 100) if total_samples > 0 else 0
            report_lines.append(f"{mut_type.replace('_', ' ').title()}: {count:,} ({percentage:.1f}%)")
    
    # Top mutations
    if "mutation" in filtered_mutations.columns:
        report_lines.append("")
        report_lines.append("TOP 10 MUTATIONS")
        report_lines.append("-" * 30)
        
        mutation_counts = filtered_mutations.groupby("mutation")["ID"].nunique()
        top_mutations = mutation_counts.nlargest(10)
        
        for mutation, count in top_mutations.items():
            percentage = (count / total_samples * 100) if total_samples > 0 else 0
            report_lines.append(f"{mutation}: {count:,} ({percentage:.1f}%)")
    
    # Create filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{selected_virus}_mutation_report_{timestamp}.txt"
    
    # Return text file
    return dcc.send_string("\n".join(report_lines), filename)

@callback(
    Output("download-mutations", "data"),
    Input("btn-download-mutations", "n_clicks"),
    State("filtered-store", "data"),
    State("selected-virus", "data"),
    State("mutation-type-filter", "value"),
    State("mutation-category-filter", "value"),
    prevent_initial_call=True
)
def download_mutation_data(n_clicks, filtered_json, virus, mutation_type, category):
    if not n_clicks:
        return dash.no_update
    
    data = get_data_store()
    selected_virus = virus or "HBV"
    
    # Get filtered data
    filtered_df = _df_from_json(filtered_json)
    
    # Get mutation data
    mutation_data = data["hbv_mut"] if selected_virus == "HBV" else data["hcv_mut"]
    
    if mutation_data.empty or filtered_df.empty:
        return dcc.send_string("No mutation data available", "no_mutation_data.txt")
    
    # Filter mutations
    if "ID" in mutation_data.columns and "ID" in filtered_df.columns:
        filtered_mutations = mutation_data[mutation_data["ID"].isin(filtered_df["ID"].unique())]
    else:
        filtered_mutations = mutation_data
    
    # Apply additional filters
    if mutation_type != "all" and "type" in filtered_mutations.columns:
        filtered_mutations = filtered_mutations[filtered_mutations["type"] == mutation_type]
    
    if category != "all" and category:
        if selected_virus == "HBV" and "drug" in filtered_mutations.columns:
            filtered_mutations = filtered_mutations[filtered_mutations["drug"] == category]
        elif selected_virus == "HCV" and "gene" in filtered_mutations.columns:
            filtered_mutations = filtered_mutations[filtered_mutations["gene"] == category]
    
    if filtered_mutations.empty:
        return dcc.send_string("No mutations match the selected filters", "no_mutations_filtered.txt")
    
    # Create filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{selected_virus}_mutations_{timestamp}.csv"
    
    # Return CSV
    return dcc.send_data_frame(filtered_mutations.to_csv, filename, index=False)

# Burden vs. Sequencing Correlation Callback
@callback(
    Output("burden-coverage-scatter", "figure"),
    Input("selected-virus", "data"),
    Input("ihme-metric-type", "value"),
    Input("filtered-store", "data"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
)
def update_burden_coverage_scatter(virus, metric_choice, filtered_json, regions, countries):
    data = get_data_store()
    ihme_df = data["ihme_df"]
    
    seq_df = _df_from_json(filtered_json)
    if seq_df.empty or ihme_df.empty:
        return _empty_plot("No sequence data or GBD data matching criteria")
        
    # Get active years and genotypes from seq_df
    y0, y1 = seq_df["Year"].min(), seq_df["Year"].max()
    active_genotypes = seq_df["genotype"].unique()
    
    # Get original sequence database for this virus to calculate global/regional counts
    if virus == "HBV":
        base_seq = data['hbv_data']
    elif virus == "HCV":
        base_seq = data['hcv_data']
    elif virus == "HEV":
        base_seq = data['hev_data']
    else:
        base_seq = data['hbv_data']
        
    if base_seq.empty:
        return _empty_plot("No sequence data available")
        
    # Filter sequences ignoring country constraint
    filtered_seq = base_seq[
        (base_seq["Year"] >= y0) &
        (base_seq["Year"] <= y1)
    ]
    if regions:
        filtered_seq = filtered_seq[filtered_seq["WHO_Regions"].isin(regions)]
    if len(active_genotypes) > 0:
        filtered_seq = filtered_seq[filtered_seq["genotype"].isin(active_genotypes)]
        
    # Get sequence counts by country
    seq_counts = filtered_seq.groupby("Country_standard").size().reset_index(name="sequence_count")
    
    # Parse GBD measure & metric directly from metric_choice (e.g. Prevalence|Number)
    try:
        measure, metric_type = metric_choice.split("|")
    except:
        measure, metric_type = "Prevalence", "Number"
        
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause = cause_lookup.get((virus or "HBV").upper())
    
    # Filter GBD data (Both sexes, sum of age groups)
    male_data = ihme_df[
        (ihme_df["sex"] == "Male") &
        (ihme_df["cause"] == cause) &
        (ihme_df["measure"] == measure) &
        (ihme_df["metric"] == metric_type)
    ].copy()
    
    female_data = ihme_df[
        (ihme_df["sex"] == "Female") &
        (ihme_df["cause"] == cause) &
        (ihme_df["measure"] == measure) &
        (ihme_df["metric"] == metric_type)
    ].copy()
    
    burden_data = pd.concat([male_data, female_data], ignore_index=True)
    
    if burden_data.empty:
        return _empty_plot("No GBD burden data found")
        
    # Apply region filter only (ignore country dropdown for background dots)
    if regions:
        burden_data = burden_data[burden_data["WHO_Regions"].isin(regions)]
        
    if burden_data.empty:
        return _empty_plot("No GBD burden data matching criteria")
        
    # Pick the latest year in GBD data
    latest_year = int(burden_data["year"].max())
    burden_latest = burden_data[burden_data["year"] == latest_year]
    
    # Aggregate burden by country (mean for percents/rates, sum for counts)
    if metric_type in ["Percent", "Rate"]:
        country_burden = burden_latest.groupby("Country_standard")["val"].mean().reset_index(name="burden_val")
    else:
        country_burden = burden_latest.groupby("Country_standard")["val"].sum().reset_index(name="burden_val")
        
    # Merge GBD burden and sequence counts
    merged = pd.merge(country_burden, seq_counts, on="Country_standard", how="inner")
    
    if merged.empty:
        return _empty_plot("No overlapping data between GBD burden and sequencing counts")
        
    virus_color = VIRUS_COLORS.get((virus or "HBV").upper(), "#E84057")
    metric_label = f"{measure} ({metric_type})"
    
    # Separate into selected (highlighted) and unselected (background) countries
    if countries:
        selected_merged = merged[merged["Country_standard"].isin(countries)]
        unselected_merged = merged[~merged["Country_standard"].isin(countries)]
    else:
        selected_merged = merged
        unselected_merged = pd.DataFrame(columns=merged.columns)
        
    fig = go.Figure()
    
    # 1. Background / Unselected countries
    if not unselected_merged.empty:
        fig.add_trace(go.Scatter(
            x=unselected_merged["burden_val"],
            y=unselected_merged["sequence_count"],
            mode="markers",
            name="Other Countries",
            marker=dict(
                size=10,
                color="rgba(255, 255, 255, 0.15)",
                line=dict(width=1, color="rgba(255, 255, 255, 0.25)")
            ),
            hovertext=unselected_merged["Country_standard"],
            hovertemplate="<b>%{hovertext}</b><br>" + metric_label + ": %{x:,.2f}<br>Sequences: %{y:,}<extra></extra>"
        ))
        
    # 2. Highlighted / Selected countries
    if not selected_merged.empty:
        fig.add_trace(go.Scatter(
            x=selected_merged["burden_val"],
            y=selected_merged["sequence_count"],
            mode="markers",
            name="Selected Country" if countries else "Countries",
            marker=dict(
                size=12,
                color=virus_color,
                opacity=0.9,
                line=dict(width=1.5, color="white")
            ),
            hovertext=selected_merged["Country_standard"],
            hovertemplate="<b>%{hovertext}</b><br>" + metric_label + ": %{x:,.2f}<br>Sequences: %{y:,}<extra></extra>"
        ))
        
    # 3. Add global trend line (calculated on ALL countries matching active filters)
    if len(merged) > 1:
        try:
            x = np.log10(merged["burden_val"] + 1e-9)
            y = np.log10(merged["sequence_count"] + 1e-9)
            coefficients = np.polyfit(x, y, 1)
            polynomial = np.poly1d(coefficients)
            x_line = np.linspace(x.min(), x.max(), 100)
            y_line = polynomial(x_line)
            
            fig.add_trace(go.Scatter(
                x=10**x_line,
                y=10**y_line,
                mode='lines',
                name='Global Trend Line' if not regions else 'Regional Trend Line',
                line=dict(color='#FF4D6D', dash='dash', width=2),
                hovertemplate='Trend Line<extra></extra>'
            ))
        except Exception:
            pass
            
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="rgba(255,255,255,0.7)"),
        xaxis=dict(
            title=metric_label,
            type="log",
            gridcolor="rgba(255,255,255,0.08)",
            linecolor="rgba(255,255,255,0.25)",
            showgrid=True,
            zeroline=False,
        ),
        yaxis=dict(
            title="Number of Sequences",
            type="log",
            gridcolor="rgba(255,255,255,0.08)",
            linecolor="rgba(255,255,255,0.25)",
            showgrid=True,
            zeroline=False,
        ),
        title=dict(
            text=f"Burden vs. Sequences Correlation ({virus or 'HBV'}, {latest_year})",
            font=dict(color="white", size=14)
        ),
        height=400,
        margin=dict(l=40, r=20, t=40, b=40),
        hovermode="closest",
        showlegend=True
    )
    
    return fig
    
@callback(
    Output("mutation-timeline", "figure"),
    Input("selected-virus", "data"),
    Input("filtered-store", "data"),
    Input("top-mutations-count", "value"),
)
def update_mutation_timeline(virus, filtered_json, top_n):
    selected_virus = virus or "HBV"
    data = get_data_store()

    # --- Select correct mutation + sequence tables ---
    if selected_virus.upper() == "HBV":
        mutation_data = data["hbv_mut"]
        sequence_df = data["hbv_data"]   # ✅ MUST be per-sequence
    else:
        mutation_data = data["hcv_mut"]
        sequence_df = data["hcv_data"]   # ✅ MUST be per-sequence

    return create_mutation_timeline(
        mutation_data,
        sequence_df,
        selected_virus,
        top_n
    )

# Country barchart
@callback(
    Output("country-barchart", "figure"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
    Input("top-countries-count", "value"),  # ADD THIS
)
def update_country_stacked_bar_callback(filtered_json, virus, regions, countries, top_n):
    df = _df_from_json(filtered_json)
    return create_country_stacked_bar(df, virus or "HBV", regions, countries, top_n=top_n or 15)

# Update priority ranking to be more responsive
@callback(
    Output("priority-data-store", "data"),
    Input("gap-store", "data"),
    Input("selected-virus", "data"),
)
def update_priority_ranking_responsive(gap_json, virus):
    gap_df = _df_from_json(gap_json)
    data = get_data_store()
    
    # Use default weights instead of user inputs
    weights = {
        "burden": 0.4,
        "coverage_gap": 0.3,
        "population": 0.2,
        "neighbor_sequencing": 0.1
    }
    
    fig, priority_df = create_priority_calculator(
        gap_df, data["ihme_df"], virus or "HBV", weights
    )
    
    return _df_to_json(priority_df)
    
# === EPIDEMIOLOGY OVER TIME & CARE CASCADE FIGURES ===

def make_disease_burden_trend_plot(ihme_df, virus="HBV", metric="Prevalence", measure="Rate", geo="GLOBAL"):
    if ihme_df is None or ihme_df.empty:
        return _empty_plot("No IHME burden trend data available")

    cause_map = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause = cause_map.get(virus, cause_map["HBV"])
    df = ihme_df[
        (ihme_df["cause"] == cause) &
        (ihme_df["measure"] == metric) &
        (ihme_df["metric"] == measure)
    ].copy()

    if geo != "GLOBAL" and "WHO_Regions" in df.columns:
        df = df[df["WHO_Regions"] == geo]

    if df.empty:
        return _empty_plot(f"No burden data found for {virus} ({metric}, {measure}, {geo})")

    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    df["val"] = pd.to_numeric(df["val"], errors="coerce")
    df = df.dropna(subset=["year", "val"])

    if measure == "Rate":
        trend_df = df.groupby("year")["val"].mean().reset_index()
    else:
        trend_df = df.groupby("year")["val"].sum().reset_index()

    trend_df = trend_df.sort_values("year")

    if trend_df.empty:
        return _empty_plot("No annual trend points available")

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=trend_df["year"],
        y=trend_df["val"],
        mode="lines+markers",
        name=f"{metric} ({measure})",
        line=dict(color="#0284c7", width=2.5),
        marker=dict(size=6, color="#0284c7"),
        hovertemplate="Year %{x}: <b>%{y:,.1f}</b><extra></extra>"
    ))

    # Add 2015 baseline horizontal dashed line
    val_2015 = trend_df[trend_df["year"] == 2015]["val"].values
    if len(val_2015) > 0:
        y_2015 = val_2015[0]
        max_x = int(trend_df["year"].max())
        fig.add_shape(
            type="line",
            x0=1990, x1=max_x,
            y0=y_2015, y1=y_2015,
            line=dict(color="#f59e0b", width=1.5, dash="dash"),
        )
        fig.add_annotation(
            x=max_x, y=y_2015,
            text="2015 baseline",
            showarrow=False,
            font=dict(color="#d97706", size=11, weight="bold"),
            xanchor="right", yanchor="bottom"
        )

    fig.update_layout(
        xaxis=dict(
            title="Year",
            gridcolor="rgba(0,0,0,0.06)",
            range=[1989, int(trend_df["year"].max()) + 1],
            showline=True,
            linecolor="rgba(0,0,0,0.2)"
        ),
        yaxis=dict(
            title=f"{metric} ({measure})",
            gridcolor="rgba(0,0,0,0.06)",
            showline=True,
            linecolor="rgba(0,0,0,0.2)"
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(t=20, b=40, l=60, r=30),
        showlegend=False
    )
    return fig


def make_cascade_bars_plot(who_gho_df, virus="HBV", regions=None, countries=None):
    if who_gho_df is None or who_gho_df.empty:
        return _empty_plot("No WHO GHO care cascade data available")

    v = (virus or "HBV").lower()
    df22 = who_gho_df[who_gho_df["year"] == 2022].copy()

    if regions:
        df22 = df22[df22["WHO_Regions"].isin(regions)]
    if countries:
        df22 = df22[df22["Country_standard"].isin(countries)]

    living_col = f"{v}_livingwith_num"
    diag_col = f"{v}_diagnosed_num"
    treat_col = f"{v}_treatment_num" if v == "hbv" else f"{v}_treatment_cumulative_num"

    living = df22[living_col].sum() if living_col in df22.columns else 0
    diag = df22[diag_col].sum() if diag_col in df22.columns else 0
    treat = df22[treat_col].sum() if treat_col in df22.columns else 0

    living_pct = 100.0 if living > 0 else 0.0
    diag_pct = (diag / living * 100.0) if living > 0 else 0.0
    treat_pct = (treat / living * 100.0) if living > 0 else 0.0

    def fmt_val(num):
        if num >= 1_000_000:
            return f"{num / 1_000_000:.1f}M"
        if num >= 1_000:
            return f"{num / 1_000:.1f}k"
        return f"{int(num):,}"

    categories = ["Treated", "Diagnosed", "Living with infection"]
    pct_vals = [treat_pct, diag_pct, living_pct]
    text_labels = [
        f"  {fmt_val(treat)}  ({treat_pct:.1f}%)",
        f"  {fmt_val(diag)}  ({diag_pct:.1f}%)",
        f"  {fmt_val(living)}  ({living_pct:.1f}%)",
    ]

    fig = go.Figure(go.Bar(
        x=pct_vals,
        y=categories,
        orientation="h",
        marker=dict(
            color=["#10b981", "#f59e0b", "#0284c7"],
            line=dict(width=0)
        ),
        text=text_labels,
        textposition="auto",
        textfont=dict(color="white", size=12, weight="bold"),
        hovertemplate="<b>%{y}</b><br>Percent of total: %{x:.1f}%<extra></extra>"
    ))

    fig.update_layout(
        xaxis=dict(
            title="Percentage of Living with Infection (%)",
            range=[0, 118],
            gridcolor="rgba(0,0,0,0.06)",
            showline=True,
            linecolor="rgba(0,0,0,0.2)"
        ),
        yaxis=dict(
            title="",
            tickfont=dict(size=12, color="#1e293b", weight="bold"),
            showline=True,
            linecolor="rgba(0,0,0,0.2)"
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(t=20, b=40, l=145, r=25),
        showlegend=False
    )
    return fig


@callback(
    Output("epi-disease-burden-trend-graph", "figure"),
    Input("selected-virus", "data"),
    Input("epi-trend-metric", "value"),
    Input("epi-trend-measure", "value"),
    #Input("epi-trend-geo", "value"),
)
def update_epi_disease_burden_trend_graph(virus, metric, measure):
    data = get_data_store()
    ihme_df = data.get("ihme_df", pd.DataFrame())
    return make_disease_burden_trend_plot(
        ihme_df=ihme_df,
        virus=virus or "HBV",
        metric=metric or "Prevalence",
        measure=measure or "Rate",
        #geo=geo or "GLOBAL"
    )


@callback(
    Output("epi-cascade-bars-graph", "figure"),
    Input("selected-virus", "data"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
)
def update_epi_cascade_bars_graph(virus, regions, countries):
    data = get_data_store()
    who_gho_df = data.get("who_gho_df", pd.DataFrame())
    return make_cascade_bars_plot(
        who_gho_df=who_gho_df,
        virus=virus or "HBV",
        regions=regions,
        countries=countries
    )


def format_big_number(val):
    if pd.isna(val) or val is None or val == 0:
        return "N/A"
    if val >= 1_000_000:
        return f"{val / 1_000_000:.1f}M"
    if val >= 1_000:
        return f"{val / 1_000:.1f}k"
    return f"{int(val):,}"


@callback(
    Output("epi-metric-dropdown", "options"),
    Output("epi-metric-dropdown", "value"),
    Input("epi-map-category", "value"),
    Input("selected-virus", "data"),
    State("epi-metric-dropdown", "value"),
)
def update_epi_metric_dropdown_options(category, virus, current_val):
    category = category or "burden"
    virus = virus or "HBV"
    
    if category == "burden":
        options = [
            {"label": "People living with infection", "value": "livingwith_num"},
            {"label": "New infections", "value": "new_infections_num"},
            {"label": "Deaths", "value": "deaths_num"},
            {"label": "Prevalence %", "value": "prevalence_pct"},
            {"label": "Incidence rate", "value": "incidence_rate"},
            {"label": "Mortality rate", "value": "mortality_rate"},
        ]
    else: # "response"
        options = [
            {"label": "Diagnosis rate %", "value": "diagnosis_rate_pct"},
            {"label": "Treatment rate %", "value": "treatment_rate_diagnosed_pct"},
        ]
        if virus == "HBV":
            options.append({"label": "HBV vaccine coverage % (HepB3)", "value": "vaccine_hepb3_coverage_pct"})
        elif virus == "HCV":
            options.append({"label": "DAA treatment / cure coverage %", "value": "daa_treatment_coverage_pct"})
        
    valid_vals = [opt["value"] for opt in options]
    default_val = options[0]["value"] if options else "livingwith_num"
    val = current_val if current_val in valid_vals else default_val
    
    return options, val


@callback(
    Output("epi-card-livingwith", "children"),
    Output("epi-card-newinfections", "children"),
    Output("epi-card-deaths", "children"),
    Output("epi-card-diag-treat", "children"),
    Input("selected-virus", "data"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
)
def update_epi_summary_cards(virus, regions, countries):
    if virus not in ["HBV", "HCV"]:
        return "N/A", "N/A", "N/A", "N/A"
        
    data = get_data_store()
    who_gho_df = data.get("who_gho_df", pd.DataFrame())
    
    if who_gho_df.empty:
        return "N/A", "N/A", "N/A", "N/A"
        
    df = who_gho_df[who_gho_df["year"] == 2022]
    
    if regions:
        df = df[df["WHO_Regions"].isin(regions)]
    if countries:
        df = df[df["Country_standard"].isin(countries)]
        
    if df.empty:
        return "N/A", "N/A", "N/A", "N/A"
        
    prefix = virus.lower() + "_"
    
    livingwith_col = f"{prefix}livingwith_num"
    livingwith_val = df[livingwith_col].sum() if livingwith_col in df.columns else 0
    
    newinfections_col = f"{prefix}new_infections_num"
    newinfections_val = df[newinfections_col].sum() if newinfections_col in df.columns else 0
    
    deaths_col = f"{prefix}deaths_num"
    deaths_val = df[deaths_col].sum() if deaths_col in df.columns else 0
    
    diagnosed_col = f"{prefix}diagnosed_num"
    treatment_col = f"{prefix}treatment_num" if virus == "HBV" else f"{prefix}treatment_cumulative_num"
    
    diag_rate_col = f"{prefix}diagnosis_rate_pct"
    treat_rate_col = f"{prefix}treatment_rate_diagnosed_pct"
    
    sum_livingwith = df[livingwith_col].sum() if livingwith_col in df.columns else 0
    sum_diagnosed = df[diagnosed_col].sum() if diagnosed_col in df.columns else 0
    sum_treated = df[treatment_col].sum() if treatment_col in df.columns else 0
    
    if sum_livingwith > 0 and sum_diagnosed > 0:
        diag_rate = (sum_diagnosed / sum_livingwith) * 100.0
    else:
        diag_rate = df[diag_rate_col].mean() if diag_rate_col in df.columns else np.nan
        
    if sum_diagnosed > 0 and sum_treated > 0:
        treat_rate = (sum_treated / sum_diagnosed) * 100.0
    else:
        treat_rate = df[treat_rate_col].mean() if treat_rate_col in df.columns else np.nan
        
    diag_str = f"{diag_rate:.1f}%" if not pd.isna(diag_rate) else "N/A"
    treat_str = f"{treat_rate:.1f}%" if not pd.isna(treat_rate) else "N/A"
    diag_treat_str = f"{diag_str} / {treat_str}"
    
    return (
        format_big_number(livingwith_val),
        format_big_number(newinfections_val),
        format_big_number(deaths_val),
        diag_treat_str
    )


@callback(
    Output("gho-burden-map", "figure"),
    Input("selected-virus", "data"),
    Input("epi-metric-dropdown", "value"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
)
def update_gho_burden_map(virus, metric, regions, countries):
    if virus not in ["HBV", "HCV"]:
        return _empty_plot("WHO GHO country profile data is only available for Hepatitis B and Hepatitis C.")
        
    if not metric:
        return _empty_plot("Select a metric to view the map")
        
    data = get_data_store()
    ihme_df = data.get("ihme_df", pd.DataFrame())
    who_gho_df = data.get("who_gho_df", pd.DataFrame())
    
    metric_labels = {
        "livingwith_num": "People living with infection",
        "new_infections_num": "New infections",
        "deaths_num": "Deaths",
        "prevalence_pct": "Prevalence %",
        "incidence_rate": "Incidence rate (per 100k)",
        "mortality_rate": "Mortality rate (per 100k)",
        "diagnosis_rate_pct": "Diagnosis rate %",
        "treatment_rate_diagnosed_pct": "Treatment rate %",
        "vaccine_hepb3_coverage_pct": "HBV vaccine coverage % (HepB3)",
        "daa_treatment_coverage_pct": "DAA treatment / cure coverage %"
    }
    metric_label = metric_labels.get(metric, metric)
    
    # Handle IHME rates (incidence_rate, mortality_rate)
    if metric in ["incidence_rate", "mortality_rate"]:
        if ihme_df.empty:
            return _empty_plot("No IHME data available")
            
        cause_map = {
            "HBV": "Total burden related to hepatitis B",
            "HCV": "Total burden related to hepatitis C",
            "HEV": "Total burden related to hepatitis E"
        }
        m_name = "Incidence" if metric == "incidence_rate" else "Deaths"
        df_ihme = ihme_df[
            (ihme_df["cause"] == cause_map.get(virus, cause_map["HBV"])) &
            (ihme_df["measure"] == m_name) &
            (ihme_df["metric"] == "Rate") &
            (ihme_df["year"] == 2021)
        ].copy()
        
        if regions and "WHO_Regions" in df_ihme.columns:
            df_ihme = df_ihme[df_ihme["WHO_Regions"].isin(regions)]
        if countries and "Country_standard" in df_ihme.columns:
            df_ihme = df_ihme[df_ihme["Country_standard"].isin(countries)]
            
        if df_ihme.empty:
            return _empty_plot("No data available for current IHME filters")
            
        df_ihme["val_mean"] = df_ihme.groupby("Country_standard")["val"].transform("mean")
        df_plot = df_ihme.drop_duplicates(subset=["Country_standard"]).copy()
        
        locations = df_plot["Country_standard"]
        z_vals = df_plot["val_mean"]
        color_scale = "YlOrRd"
        
    else: # WHO GHO metrics
        if who_gho_df.empty:
            return _empty_plot("No WHO GHO data available")
            
        prefix = virus.lower() + "_"
        col_name = prefix + metric
        
        if metric == "daa_treatment_coverage_pct":
            col_name = prefix + "treatment_rate_diagnosed_pct"
            
        if col_name not in who_gho_df.columns:
            return _empty_plot(f"Metric '{metric}' not found for virus {virus}")
            
        map_year = 2022
        if metric == "vaccine_hepb3_coverage_pct":
            non_null_years = who_gho_df[who_gho_df[col_name].notna()]["year"]
            map_year = int(non_null_years.max()) if not non_null_years.empty else 2022
            
        df = who_gho_df[who_gho_df["year"] == map_year].copy()
        
        if regions:
            df = df[df["WHO_Regions"].isin(regions)]
        if countries:
            df = df[df["Country_standard"].isin(countries)]
            
        if df.empty or df[col_name].isna().all():
            return _empty_plot(f"No data available for year {map_year} with current filters")
            
        locations = df["Country_standard"]
        z_vals = df[col_name]
        
    is_pct = ("pct" in metric or "rate" in metric)
    bins, step_colors, tick_labels = _bin_continuous_series(z_vals, is_pct=is_pct)
        
    fig = go.Figure(data=go.Choropleth(
        locations=locations,
        locationmode="country names",
        z=bins,
        colorscale=_stepped_colorscale(step_colors) if step_colors else "YlOrRd",
        showscale=False,
        marker_line_color="rgba(0,0,0,0.3)",
        marker_line_width=0.5,
        hovertemplate="<b>%{location}</b><br>" + metric_label + ": <b>%{customdata:,.2f}</b><extra></extra>",
        customdata=z_vals,
    ))
    
    if step_colors and tick_labels:
        _add_single_box_legend(fig, step_colors, tick_labels)
    
    fig.update_geos(
        projection_type="equirectangular",
        showcountries=True,
        countrycolor="rgba(0,0,0,0.2)",
        showsubunits=True,
        fitbounds=False,
        showframe=False,
        lataxis_range=[-65, 85],
        domain=dict(x=[0, 1], y=[0, 1]),
    )
    fig.update_layout(
        height=330,
        margin=dict(t=0, b=0, l=0, r=0),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(visible=False, showgrid=False, zeroline=False),
        yaxis=dict(visible=False, showgrid=False, zeroline=False),
    )
    return fig


@callback(
    Output("gho-cascade-chart", "figure"),
    Input("selected-virus", "data"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
)
def update_gho_cascade_chart(virus, regions, countries):
    if virus not in ["HBV", "HCV"]:
        return _empty_plot("WHO GHO country profile data is only available for Hepatitis B and Hepatitis C.")
        
    data = get_data_store()
    who_gho_df = data.get("who_gho_df", pd.DataFrame())
    
    if who_gho_df.empty:
        return _empty_plot("No WHO GHO data available")
        
    df = who_gho_df[who_gho_df["year"] == 2022]
    
    if regions:
        df = df[df["WHO_Regions"].isin(regions)]
    if countries:
        df = df[df["Country_standard"].isin(countries)]
        
    if df.empty:
        return _empty_plot("No data available for current filters")
        
    prefix = virus.lower() + "_"
    
    livingwith_col = f"{prefix}livingwith_num"
    diagnosed_col = f"{prefix}diagnosed_num"
    treatment_col = f"{prefix}treatment_num" if virus == "HBV" else f"{prefix}treatment_cumulative_num"
    
    living_with = df[livingwith_col].sum() if livingwith_col in df.columns else 0
    diagnosed = df[diagnosed_col].sum() if diagnosed_col in df.columns else 0
    treated = df[treatment_col].sum() if treatment_col in df.columns else 0
    
    if pd.isna(living_with) or living_with == 0:
        return _empty_plot("No cascade numbers available for current filters")
        
    stages = ["Living with infection", "Diagnosed", "Treated"]
    values = [living_with, diagnosed, treated]
    values = [val if not pd.isna(val) else 0 for val in values]
    
    virus_color = VIRUS_COLORS.get(virus, "#E84057")
    colors = [virus_color, shade(virus_color, -20), shade(virus_color, -40)]
    
    fig = go.Figure(go.Funnel(
        y=stages,
        x=values,
        textinfo="value+percent initial",
        marker={"color": colors, "line": {"width": [1, 1, 1], "color": ["#161D26", "#161D26", "#161D26"]}},
        connector={"line": {"color": "rgba(255,255,255,0.05)", "width": 1}}
    ))
    
    fig.update_layout(
        title=dict(
            text=f"Hepatitis {virus[-1]} Care Cascade (2022)",
            font=dict(size=14, color="#ECEFF2"),
            x=0.5, y=0.95,
            xanchor="center"
        ),
        margin=dict(l=40, r=40, t=50, b=20),
        height=350,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )
    
    return fig


def make_epi_priority_table(virus, regions, countries):
    data = get_data_store()
    
    # 1. Get sequence data and compute gap
    if virus == "HBV":
        seq_df = data["hbv_data"]
    else:
        seq_df = data["hcv_data"]
        
    y0, y1 = int(seq_df["Year"].min()) if not seq_df.empty else 2000, int(seq_df["Year"].max()) if not seq_df.empty else 2021
    
    # Filter sequences by regions and countries for observed count
    filtered_seq = seq_df.copy()
    if regions:
        filtered_seq = filtered_seq[filtered_seq["WHO_Regions"].isin(regions)]
    if countries:
        filtered_seq = filtered_seq[filtered_seq["Country_standard"].isin(countries)]
        
    gap_df = compute_gap_df(
        virus=virus,
        filtered_seq_df=filtered_seq,
        ihme_df=data["ihme_df"],
        selected_years=[y0, y1],
        who_regions=None,
        countries=None,
        ihme_metric_choice="Prevalence|Number",
        sex="Both"
    )
    
    # 2. Run priority calculator
    weights = {
        "burden": 0.4,
        "coverage_gap": 0.3,
        "population": 0.2,
        "neighbor_sequencing": 0.1
    }
    _, priority_df = create_priority_calculator(gap_df, data["ihme_df"], virus, weights)
    
    if priority_df is None or priority_df.empty:
        return html.Div(
            "No priority data available for current filters",
            className="hep-empty-table",
        )
        
    # 3. Merge with WHO GHO 2022 data for diagnosis and treatment rates
    who_gho_df = data.get("who_gho_df", pd.DataFrame())
    if not who_gho_df.empty:
        who_2022 = who_gho_df[who_gho_df["year"] == 2022].copy()
        prefix = virus.lower() + "_"
        diag_col = f"{prefix}diagnosis_rate_pct"
        treat_col = f"{prefix}treatment_rate_diagnosed_pct"
        
        # Select columns to merge
        who_cols = ["Country_standard"]
        if diag_col in who_2022.columns:
            who_cols.append(diag_col)
        if treat_col in who_2022.columns:
            who_cols.append(treat_col)
            
        merged = pd.merge(
            priority_df,
            who_2022[who_cols],
            on="Country_standard",
            how="left"
        )
    else:
        merged = priority_df.copy()
        prefix = virus.lower() + "_"
        diag_col = f"{prefix}diagnosis_rate_pct"
        treat_col = f"{prefix}treatment_rate_diagnosed_pct"
        merged[diag_col] = np.nan
        merged[treat_col] = np.nan

    # 4. Filter priority table by user continent/country selections
    if regions:
        merged = merged[merged["WHO_Regions"].isin(regions)]
    if countries:
        merged = merged[merged["Country_standard"].isin(countries)]
        
    # Sort and rank
    merged = merged.sort_values("priority_score", ascending=False)
    if not merged.empty:
        merged["rank"] = range(1, len(merged) + 1)
    else:
        merged["rank"] = []
    
    # Select and rename columns
    display_df = merged.copy()
    columns_needed = [
        "rank",
        "Country_standard",
        "burden",
        diag_col,
        treat_col,
    ]
    
    available_cols = [c for c in columns_needed if c in display_df.columns]
    display_df = display_df[available_cols].copy()
    
    rename_map = {
        "rank": "Rank",
        "Country_standard": "Country",
        "burden": "Living with infection",
        diag_col: "Diagnosed %",
        treat_col: "Treatment %",
    }
    display_df = display_df.rename(columns=rename_map)
    
    # Format columns
    for col in ["Living with infection"]:
        if col in display_df.columns:
            display_df[col] = pd.to_numeric(display_df[col], errors="coerce").apply(
                lambda x: f"{x:,.0f}" if pd.notna(x) else "N/A"
            )
            
    for col in ["Diagnosed %", "Treatment %"]:
        if col in display_df.columns:
            display_df[col] = pd.to_numeric(display_df[col], errors="coerce").apply(
                lambda x: f"{x:.1f}%" if pd.notna(x) else "N/A"
            )
            
    # Render table
    from hep_theme import TABLE_HEADER_STYLE, TABLE_CELL_STYLE, TABLE_ODD_ROW_STYLE
    
    compact_cell_style = TABLE_CELL_STYLE.copy()
    compact_cell_style.update({"padding": "8px 8px", "fontSize": "12px", "maxWidth": None})
    
    compact_header_style = TABLE_HEADER_STYLE.copy()
    compact_header_style.update({"padding": "8px 8px", "fontSize": "11px"})

    return dash_table.DataTable(
        data=display_df.to_dict("records"),
        columns=[{"name": col, "id": col} for col in display_df.columns],
        page_size=6,
        sort_action="native",
        style_table={
            "overflowX": "hidden",
            "backgroundColor": "transparent",
            "width": "100%",
            "minWidth": "100%",
        },
        style_header=compact_header_style,
        style_cell=compact_cell_style,
        style_cell_conditional=[
            {
                "if": {"column_id": "Rank"},
                "width": "40px",
                "textAlign": "center",
            },
            {
                "if": {"column_id": "Country"},
                "width": "35%",
                "whiteSpace": "normal",
                "fontWeight": "600",
            },
            {
                "if": {"column_id": "Living with infection"},
                "width": "27%",
            },
            {
                "if": {"column_id": "Diagnosed %"},
                "width": "19%",
            },
            {
                "if": {"column_id": "Treatment %"},
                "width": "19%",
            },
        ],
        style_data_conditional=[TABLE_ODD_ROW_STYLE],
    )


def make_hbv_epi_scatter_plot(regions, countries, y_metric):
    data = get_data_store()
    who_gho_df = data.get("who_gho_df", pd.DataFrame())
    
    fig = go.Figure()
    
    if not who_gho_df.empty:
        df_2022 = who_gho_df[who_gho_df["year"] == 2022].copy()
        
        if regions:
            df_2022 = df_2022[df_2022["WHO_Regions"].isin(regions)]
        if countries:
            df_2022 = df_2022[df_2022["Country_standard"].isin(countries)]
            
        x_col = "hbv_vaccine_hepb3_coverage_pct"
        y_col = y_metric
        
        df_plot = df_2022.dropna(subset=[x_col, y_col]).copy()
        
        if not df_plot.empty:
            y_label_map = {
                "hbv_prevalence_pct": "Prevalence %",
                "hbv_new_infections_num": "New Infections",
                "hbv_livingwith_num": "Living with HBV"
            }
            y_label = y_label_map.get(y_metric, "Disease Burden")
            
            fig.add_trace(
                go.Scatter(
                    x=df_plot[x_col],
                    y=df_plot[y_col],
                    mode="markers",
                    text=df_plot["Country_standard"],
                    marker=dict(
                        color="#2EC4B6", 
                        size=10,
                        line=dict(width=1, color="white"),
                        opacity=0.8
                    ),
                    name="Countries",
                    hovertemplate="<b>%{text}</b><br>HepB3 Vaccine Coverage: %{x:.1f}%<br>" + y_label + ": %{y:,.2f}<extra></extra>"
                )
            )
            
            if len(df_plot) > 1:
                import numpy as np
                try:
                    x = df_plot[x_col].values
                    y = df_plot[y_col].values
                    idx = np.isfinite(x) & np.isfinite(y)
                    if np.sum(idx) > 1:
                        coef = np.polyfit(x[idx], y[idx], 1)
                        poly1d_fn = np.poly1d(coef)
                        x_range = np.linspace(x.min(), x.max(), 100)
                        y_fit = poly1d_fn(x_range)
                        fig.add_trace(
                            go.Scatter(
                                x=x_range,
                                y=y_fit,
                                mode="lines",
                                name="Trend Line",
                                line=dict(color="#FF9F1C", width=2, dash="dash"),
                                hoverinfo="skip"
                            )
                        )
                except Exception as e:
                    pass
                    
            fig.update_layout(
                xaxis_title="HepB3 Vaccine Coverage (%)",
                yaxis_title=y_label,
            )
            
    fig.update_layout(
        title=None,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="IBM Plex Sans, sans-serif", color="#475569", size=11),
        xaxis=dict(
            gridcolor="rgba(15, 23, 42, 0.08)",
            linecolor="rgba(15, 23, 42, 0.2)",
            showgrid=True,
            zeroline=False,
            title=dict(font=dict(color="#0f172a", size=11, weight="bold")),
            tickfont=dict(color="#475569", size=10),
        ),
        yaxis=dict(
            gridcolor="rgba(15, 23, 42, 0.08)",
            linecolor="rgba(15, 23, 42, 0.2)",
            showgrid=True,
            zeroline=False,
            title=dict(font=dict(color="#0f172a", size=11, weight="bold")),
            tickfont=dict(color="#475569", size=10),
        ),
        height=360,
        margin=dict(l=55, r=25, t=15, b=45),
        hovermode="closest",
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="top",
            y=0.98,
            xanchor="right",
            x=0.98,
            bgcolor="rgba(255, 255, 255, 0.85)",
            bordercolor="rgba(15, 23, 42, 0.12)",
            borderwidth=1,
            font=dict(color="#1e293b", size=10)
        )
    )
    return fig


def make_hcv_epi_scatter_plot(regions, countries, x_metric, y_metric):
    data = get_data_store()
    who_gho_df = data.get("who_gho_df", pd.DataFrame())
    
    fig = go.Figure()
    
    if not who_gho_df.empty:
        df_2022 = who_gho_df[who_gho_df["year"] == 2022].copy()
        
        if regions:
            df_2022 = df_2022[df_2022["WHO_Regions"].isin(regions)]
        if countries:
            df_2022 = df_2022[df_2022["Country_standard"].isin(countries)]
            
        x_col = x_metric
        y_col = y_metric
        
        df_plot = df_2022.dropna(subset=[x_col, y_col]).copy()
        
        if not df_plot.empty:
            x_label_map = {
                "hcv_diagnosis_rate_pct": "Diagnosis Rate %",
                "hcv_treatment_rate_diagnosed_pct": "Treatment Rate %"
            }
            y_label_map = {
                "hcv_prevalence_pct": "Prevalence %",
                "hcv_new_infections_num": "New Infections",
                "hcv_livingwith_num": "Living with HCV"
            }
            x_label = x_label_map.get(x_metric, "Program Response")
            y_label = y_label_map.get(y_metric, "Disease Burden")
            
            fig.add_trace(
                go.Scatter(
                    x=df_plot[x_col],
                    y=df_plot[y_col],
                    mode="markers",
                    text=df_plot["Country_standard"],
                    marker=dict(
                        color="#4FAEFF", 
                        size=10,
                        line=dict(width=1, color="white"),
                        opacity=0.8
                    ),
                    name="Countries",
                    hovertemplate="<b>%{text}</b><br>" + x_label + ": %{x:.1f}%<br>" + y_label + ": %{y:,.2f}<extra></extra>"
                )
            )
            
            if len(df_plot) > 1:
                import numpy as np
                try:
                    x = df_plot[x_col].values
                    y = df_plot[y_col].values
                    idx = np.isfinite(x) & np.isfinite(y)
                    if np.sum(idx) > 1:
                        coef = np.polyfit(x[idx], y[idx], 1)
                        poly1d_fn = np.poly1d(coef)
                        x_range = np.linspace(x.min(), x.max(), 100)
                        y_fit = poly1d_fn(x_range)
                        fig.add_trace(
                            go.Scatter(
                                x=x_range,
                                y=y_fit,
                                mode="lines",
                                name="Trend Line",
                                line=dict(color="#FFD166", width=2, dash="dash"),
                                hoverinfo="skip"
                            )
                        )
                except Exception as e:
                    pass
                    
            fig.update_layout(
                xaxis_title=x_label,
                yaxis_title=y_label,
            )
            
    fig.update_layout(
        title=None,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="IBM Plex Sans, sans-serif", color="#475569", size=11),
        xaxis=dict(
            gridcolor="rgba(15, 23, 42, 0.08)",
            linecolor="rgba(15, 23, 42, 0.2)",
            showgrid=True,
            zeroline=False,
            title=dict(font=dict(color="#0f172a", size=11, weight="bold")),
            tickfont=dict(color="#475569", size=10),
        ),
        yaxis=dict(
            gridcolor="rgba(15, 23, 42, 0.08)",
            linecolor="rgba(15, 23, 42, 0.2)",
            showgrid=True,
            zeroline=False,
            title=dict(font=dict(color="#0f172a", size=11, weight="bold")),
            tickfont=dict(color="#475569", size=10),
        ),
        height=360,
        margin=dict(l=55, r=25, t=15, b=45),
        hovermode="closest",
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="top",
            y=0.98,
            xanchor="right",
            x=0.98,
            bgcolor="rgba(255, 255, 255, 0.85)",
            bordercolor="rgba(15, 23, 42, 0.12)",
            borderwidth=1,
            font=dict(color="#1e293b", size=10)
        )
    )
    return fig


@callback(
    Output("hbv-epi-trend-plot", "figure"),
    Output("hbv-epi-priority-table", "children"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
    Input("hbv-y-metric", "value"),
)
def update_hbv_epi_row(regions, countries, y_metric):
    fig = make_hbv_epi_scatter_plot(regions, countries, y_metric)
    table = make_epi_priority_table("HBV", regions, countries)
    return fig, table


@callback(
    Output("hcv-epi-trend-plot", "figure"),
    Output("hcv-epi-priority-table", "children"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
    Input("hcv-x-metric", "value"),
    Input("hcv-y-metric", "value"),
)
def update_hcv_epi_row(regions, countries, x_metric, y_metric):
    fig = make_hcv_epi_scatter_plot(regions, countries, x_metric, y_metric)
    table = make_epi_priority_table("HCV", regions, countries)
    return fig, table


@callback(
    Output("hbv-epi-trends-controls", "style"),
    Output("hcv-epi-trends-controls", "style"),
    Output("hbv-epi-trends-row", "style"),
    Output("hcv-epi-trends-row", "style"),
    Output("hbv-epi-priority-table", "style"),
    Output("hcv-epi-priority-table", "style"),
    Input("selected-virus", "data"),
)
def toggle_epi_trends_rows(virus):
    if virus == "HBV":
        return (
            {"display": "flex", "alignItems": "center"},
            {"display": "none"},
            {"display": "block"},
            {"display": "none"},
            {"display": "block"},
            {"display": "none"}
        )
    elif virus == "HCV":
        return (
            {"display": "none"},
            {"display": "flex", "alignItems": "center"},
            {"display": "none"},
            {"display": "block"},
            {"display": "none"},
            {"display": "block"}
        )
    else:
        return (
            {"display": "none"},
            {"display": "none"},
            {"display": "none"},
            {"display": "none"},
            {"display": "none"},
            {"display": "none"}
        )





@callback(
    Output("epi-coverage-map-graph", "figure"),
    Input("selected-virus", "data"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
    Input("filtered-store", "data"),
)
def update_epi_coverage_map(virus, regions, countries, filtered_json):
    data = get_data_store()
    v = virus or "HBV"
    df = _df_from_json(filtered_json)
    
    if df.empty or "Country_standard" not in df.columns:
        df = data.get("hbv_data" if v == "HBV" else "hcv_data", pd.DataFrame())
        
    if regions and "WHO_Regions" in df.columns:
        df = df[df["WHO_Regions"].isin(regions)]
    if countries and "Country_standard" in df.columns:
        df = df[df["Country_standard"].isin(countries)]
        
    if df.empty or "Country_standard" not in df.columns:
        return _empty_world("No sequence data available for selected filters")
        
    ctry_counts = df.groupby("Country_standard").size().reset_index(name="Metric_raw")
    geno_counts = df.groupby(["Country_standard", "genotype"]).size().reset_index(name="Count") if "genotype" in df.columns else pd.DataFrame()
    
    return create_world_map(
        country_data=ctry_counts,
        country_genotype_counts=geno_counts,
        coord_lookup=data["coord_lookup"],
        virus_type=v,
        display_mode="raw",
        map_title="Sequence Count Map",
        height=330
    )


@callback(
    Output("gho-data-table", "children"),
    Input("selected-virus", "data"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
)
def update_gho_data_table(virus, regions, countries):
    if virus not in ["HBV", "HCV"]:
        return html.P("Select Hepatitis B or C to view WHO GHO detailed country profiles.", className="text-muted text-center py-4")
        
    data = get_data_store()
    who_gho_df = data.get("who_gho_df", pd.DataFrame())
    
    if who_gho_df.empty:
        return html.P("No WHO GHO data available.", className="text-muted text-center py-4")
        
    df = who_gho_df[who_gho_df["year"] == 2022].copy()
    
    if regions:
        df = df[df["WHO_Regions"].isin(regions)]
    if countries:
        df = df[df["Country_standard"].isin(countries)]
        
    if df.empty:
        return html.P("No data available matching filters.", className="text-muted text-center py-4")
        
    prefix = virus.lower() + "_"
    
    cols_map = {
        "Country_standard": "Country",
        "WHO_Regions": "WHO Region",
        f"{prefix}livingwith_num": "Living with infection",
        f"{prefix}new_infections_num": "New infections",
        f"{prefix}deaths_num": "Deaths",
        f"{prefix}prevalence_pct": "Prevalence %",
        f"{prefix}diagnosis_rate_pct": "Diagnosis rate %",
        f"{prefix}treatment_rate_diagnosed_pct": "Treatment rate %"
    }
    
    if virus == "HBV":
        cols_map[f"{prefix}vaccine_hepb3_coverage_pct"] = "HBV vaccine coverage % (HepB3)"
        
    available_cols = [c for c in cols_map.keys() if c in df.columns]
    table_df = df[available_cols].copy()
    
    for col in table_df.columns:
        if col in ["Country_standard", "WHO_Regions"]:
            continue
        if "num" in col or "livingwith" in col or "infections" in col or "deaths" in col:
            table_df[col] = table_df[col].apply(lambda x: f"{int(x):,}" if not pd.isna(x) and not np.isinf(x) else "N/A")
        elif "pct" in col or "rate" in col:
            table_df[col] = table_df[col].apply(lambda x: f"{x:.1f}%" if not pd.isna(x) and not np.isinf(x) else "N/A")
            
    table_df = table_df.rename(columns=cols_map)
    
    from hep_theme import TABLE_HEADER_STYLE, TABLE_CELL_STYLE, TABLE_ODD_ROW_STYLE
    
    table = dash.dash_table.DataTable(
        data=table_df.to_dict('records'),
        columns=[{"name": col, "id": col} for col in table_df.columns],
        page_size=10,
        style_table={'overflowX': 'auto'},
        style_cell=TABLE_CELL_STYLE,
        style_header=TABLE_HEADER_STYLE,
        style_data_conditional=[TABLE_ODD_ROW_STYLE],
        filter_action="native",
        sort_action="native",
    )
    
    return table


@callback(
    Output("priority-download", "data"),
    Input("priority-download-btn", "n_clicks"),
    State("priority-data-store", "data"),
    State("selected-virus", "data"),
    prevent_initial_call=True,
)
def download_priority_table(n_clicks, priority_json, virus):
    if not n_clicks:
        return dash.no_update
    
    priority_df = _df_from_json(priority_json)
    if priority_df.empty:
        return dcc.send_string("No priority data available", "priority_data_empty.txt")
    
    # Create comprehensive CSV with priority data
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{virus}_sequencing_priority_ranking_{timestamp}.csv"
    
    return dcc.send_data_frame(priority_df.to_csv, filename, index=False)
    
@callback(
    Output("year-range-slider", "min"),
    Output("year-range-slider", "max"),
    Output("year-range-slider", "value", allow_duplicate=True),
    Output("year-range-slider", "marks"),
    Output('continent-dropdown', 'options'),
    Output('country-dropdown', 'options'),
    Output('genotype-dropdown', 'options'),
    Input("selected-virus", "data"),
    prevent_initial_call='initial_duplicate'
)
def init_controls(virus):
    data = get_data_store()
    
    if virus == "HBV":
        base = data['hbv_data']
    elif virus == "HCV":
        base = data['hcv_data']
    elif virus == "HEV":
        base = data['hev_data']
    else:
        base = data['hbv_data']
    
    if base.empty:
        return 1963, 2024, [1963, 2024], {}, [], [], []
    
    y0, y1 = int(base["Year"].min()), int(base["Year"].max())
    
    # Generate year marks dynamically
    marks = {y0: str(y0), y1: str(y1)}
    start_decade = ((y0 // 10) + 1) * 10
    for y in range(start_decade, y1, 10):
        if y - y0 >= 3 and y1 - y >= 3:
            marks[y] = str(y)
            
    # sort keys
    marks = {k: marks[k] for k in sorted(marks.keys())}
    
    cont_opts = [{"label": r, "value": r} for r in sorted(base["WHO_Regions"].dropna().unique()) if r!="Unknown"]
    country_opts = [{"label": c, "value": c} for c in sorted(base["Country_standard"].dropna().unique()) if c!="Unknown"]
    geno_opts = [{"label": g, "value": g} for g in sorted(base["genotype"].dropna().unique())]
    return y0, y1, [y0, y1], marks, cont_opts, country_opts, geno_opts


@callback(
    Output("continent-dropdown", "style"),
    Output("continent-dropdown", "className"),
    Input("continent-dropdown", "value")
)
def update_continent_dropdown_style(selected_values):
    if not selected_values:
        return {"--selected-regions-summary": '""'}, "hep-dropdown"
    
    count = len(selected_values)
    summary_text = f"'{count} region selected'" if count == 1 else f"'{count} regions selected'"
    return {"--selected-regions-summary": summary_text}, "hep-dropdown has-selections"



@callback(
    Output("epi-prevalence-total", "children"),
    Output("epi-prevalence-trend", "children"),
    Output("epi-incidence-total", "children"),
    Output("epi-incidence-trend", "children"),
    Output("epi-deaths-total", "children"),
    Output("epi-deaths-trend", "children"),
    Output("epi-coverage-percent", "children"),
    Output("epi-coverage-status", "children"),
    Output("epi-sex-ratio", "children"),
    Output("epi-top-age-group", "children"),
    Output("epi-age-percentage", "children"),
    Output("epi-top-region", "children"),
    Output("epi-region-percentage", "children"),
    Output("epi-2030-progress", "children"),
    Input("selected-virus", "data"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
    Input("filtered-store", "data"),  # For sequence coverage calculation
    Input("gap-store", "data"),  # For coverage gap data
)
def update_epidemiology_summary(virus, regions, countries, filtered_json, gap_json):
    data = get_data_store()
    ihme_df = data["ihme_df"]
    
    # Default values
    defaults = ["N/A"] * 14
    defaults[1] = html.Span("No trend data", className="text-muted")  # prevalence trend
    defaults[3] = html.Span("No trend data", className="text-muted")  # incidence trend
    defaults[5] = html.Span("No trend data", className="text-muted")  # deaths trend
    defaults[7] = "No coverage data"  # coverage status
    defaults[10] = ""  # age percentage (empty string)
    defaults[12] = ""  # region percentage (empty string)
    
    # Calculate epidemiology summary
    summary = calculate_epidemiology_summary(
        ihme_df=ihme_df,
        virus=virus or "HBV",
        regions=regions,
        countries=countries
    )
    
    # Format prevalence
    prevalence_total = summary.get("prevalence_total")
    prevalence_total_display = format_large_number(prevalence_total) if prevalence_total is not None else "N/A"
    prevalence_trend = summary.get("prevalence_trend")
    prevalence_trend_display = format_trend(prevalence_trend)
    
    # Format incidence
    incidence_total = summary.get("incidence_total")
    incidence_total_display = format_large_number(incidence_total) if incidence_total is not None else "N/A"
    incidence_trend = summary.get("incidence_trend")
    incidence_trend_display = format_trend(incidence_trend)
    
    # Format deaths
    deaths_total = summary.get("deaths_total")
    deaths_total_display = format_large_number(deaths_total) if deaths_total is not None else "N/A"
    deaths_trend = summary.get("deaths_trend")
    deaths_trend_display = format_trend(deaths_trend)
    
    # Calculate sequencing coverage percentage
    filtered_df = _df_from_json(filtered_json)
    gap_df = _df_from_json(gap_json)
    
    coverage_percent = "N/A"
    coverage_status = "No coverage data"
    
    if not gap_df.empty and "observed_sequences" in gap_df.columns and "expected_sequences" in gap_df.columns:
        total_observed = gap_df["observed_sequences"].sum()
        total_expected = gap_df["expected_sequences"].sum()
        if total_expected > 0:
            coverage_pct = (total_observed / total_expected) * 100
            coverage_percent = f"{coverage_pct:.1f}%"
            
            if coverage_pct >= 100:
                coverage_status = "Adequate coverage"
            elif coverage_pct >= 50:
                coverage_status = "Moderate coverage"
            else:
                coverage_status = "Low coverage"
    
    # Format sex ratio
    sex_ratio = summary.get("sex_ratio")
    if sex_ratio is not None and isinstance(sex_ratio, (int, float)):
        sex_ratio_display = f"{sex_ratio:.1f}:1"
    else:
        sex_ratio_display = "N/A"
    
    # Format top age group
    top_age = summary.get("top_age_group", "N/A")
    age_percentage = summary.get("top_age_percentage")
    if age_percentage is not None and isinstance(age_percentage, (int, float)):
        age_percentage_display = f"({age_percentage:.1f}% of total)"
    else:
        age_percentage_display = ""
    
    # Format top region
    top_region = summary.get("top_region", "N/A")
    region_percentage = summary.get("top_region_percentage")
    if region_percentage is not None and isinstance(region_percentage, (int, float)):
        region_percentage_display = f"({region_percentage:.1f}% of total)"
    else:
        region_percentage_display = ""
    
    # Format WHO 2030 progress
    who_progress = summary.get("who_progress")
    reduction_needed = summary.get("reduction_needed")
    
    who_progress_display = "N/A"
    if who_progress is not None and isinstance(who_progress, (int, float)):
        who_progress_display = f"{who_progress:.0f}%"
        if reduction_needed is not None and isinstance(reduction_needed, (int, float)):
            who_progress_display += f" ({reduction_needed:.0f}% to go)"
    
    return [
        prevalence_total_display,  # epi-prevalence-total
        prevalence_trend_display,  # epi-prevalence-trend
        incidence_total_display,   # epi-incidence-total
        incidence_trend_display,   # epi-incidence-trend
        deaths_total_display,      # epi-deaths-total
        deaths_trend_display,      # epi-deaths-trend
        coverage_percent,          # epi-coverage-percent
        coverage_status,           # epi-coverage-status
        sex_ratio_display,         # epi-sex-ratio
        top_age,                   # epi-top-age-group
        age_percentage_display,    # epi-age-percentage
        top_region,                # epi-top-region
        region_percentage_display, # epi-region-percentage
        who_progress_display,      # epi-2030-progress
    ]



# === ACTIONS AND DOWNLOADS ==============================================================
# - Main data download with Taxa -
def _pick_mode(series: pd.Series):
    m = series.mode(dropna=True)
    if not m.empty:
        return m.iloc[0]
    s = series.dropna()
    return s.iloc[0] if not s.empty else None

def _build_keys(df: pd.DataFrame, is_main: bool) -> pd.DataFrame:
    # Drop duplicated column NAMES if any (keeps first occurrence)
    df = df.loc[:, ~df.columns.duplicated()].copy()

    # Country key (prefer Country_standard, fallback to Country)
    if "Country_standard" in df.columns and "Country" in df.columns:
        df["Country_key"] = df["Country_standard"].fillna(df["Country"])
    elif "Country_standard" in df.columns:
        df["Country_key"] = df["Country_standard"]
    elif "Country" in df.columns:
        df["Country_key"] = df["Country"]
    else:
        df["Country_key"] = pd.NA

    # Year key (accept Year or Date; coerce to numeric)
    year_num = None
    if "Year" in df.columns:
        year_num = pd.to_numeric(df["Year"], errors="coerce")
    if "Date" in df.columns:
        date_num = pd.to_numeric(df["Date"], errors="coerce")
        year_num = year_num.fillna(date_num) if year_num is not None else date_num
    df["Year_key"] = year_num if year_num is not None else pd.NA

    # genotype key
    df["Genotype_key"] = df["genotype"] if "genotype" in df.columns else pd.NA

    # Normalize types for robust grouping/merging
    for c in ["Country_key", "Genotype_key"]:
        df[c] = df[c].astype("string")

    # Ensure expected columns exist on main DF (for downstream CSV)
    if is_main:
        for c in ["Country_standard", "WHO_Regions", "Year", "genotype"]:
            if c not in df.columns:
                df[c] = pd.NA

    return df

# ----------------------------------------------------------------------
# Callback
# ----------------------------------------------------------------------
@callback(
    Output("download-data", "data"),
    Input("btn-download-data", "n_clicks"),
    State("filtered-store", "data"),
    State("selected-virus", "data"),
    State("year-range-slider", "value"),
    State("continent-dropdown", "value"),
    State("country-dropdown", "value"),
    State("genotype-dropdown", "value"),
    prevent_initial_call=True,
)
def download_main_data_with_taxa(n_clicks, filtered_json, virus, year_range, regions, countries, genotypes):
    # Only act on actual clicks
    if not n_clicks:
        raise PreventUpdate

    # Parse filtered data safely
    try:
        filtered_df = _df_from_json(filtered_json)
    except Exception:
        filtered_df = pd.DataFrame()

    selected_virus = (virus or "HBV").upper()
    data = get_data_store()  # UPDATED

    # If empty → still return a small CSV so the click always downloads something
    if filtered_df.empty:
        empty_df = pd.DataFrame(columns=["Country_standard", "WHO_Regions", "Year", "genotype", "Taxa"])
        return dcc.send_data_frame(
            empty_df.to_csv,
            f"{selected_virus}_no_data_available.csv",
            index=False
        )

    # Prepare main DF with join keys
    download_df = _build_keys(filtered_df.copy(), is_main=True)

    # Get summary data by virus
    summary_key = "hbv_summary_raw" if selected_virus == "HBV" else "hcv_summary_raw"
    summary_data = data.get(summary_key, pd.DataFrame())  # UPDATED

    if not summary_data.empty:
        summary_clean = _build_keys(summary_data.copy(), is_main=False)

        # Try progressively less-specific keys
        merge_strategies = [
            (["Genotype_key", "Country_key", "Year_key"], "exact_match"),
            (["Genotype_key", "Country_key"], "genotype_country"),
            (["Genotype_key"], "genotype_only"),
            (["Country_key"], "country_only"),
        ]

        merged_successfully = False
        base_len = len(download_df)

        for merge_cols, _strategy in merge_strategies:
            # Collapse to ONE Taxa per key to avoid many-to-many merges
            taxa_map = (
                summary_clean
                .dropna(subset=["Taxa"])
                .groupby(merge_cols, dropna=False)["Taxa"]
                .agg(_pick_mode)
                .reset_index()
            )

            if taxa_map.empty or taxa_map["Taxa"].isna().all():
                continue

            try:
                merged = download_df.merge(taxa_map, on=merge_cols, how="left")
            except Exception:
                continue

            # Guard against accidental row multiplication
            if len(merged) > base_len * 2:
                # Skip this strategy if it explodes rows
                continue

            download_df = merged
            merged_successfully = True
            break

        if not merged_successfully:
            download_df["Taxa"] = "No Taxa match found"
    else:
        download_df["Taxa"] = "Summary data not available"

    # Put Taxa first if present
    if "Taxa" in download_df.columns:
        ordered = ["Taxa"] + [c for c in download_df.columns if c != "Taxa"]
        download_df = download_df.reindex(columns=ordered)

    # Construct filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{selected_virus}_sequence_data_with_taxa_{timestamp}.csv"

    # Return CSV
    return dcc.send_data_frame(download_df.to_csv, filename, index=False)
    
def merge_taxa_information(main_df, summary_df, virus_type):
    """
    Merge Taxa information from summary data into main dataframe
    """
    if summary_df.empty or main_df.empty:
        main_df["Taxa"] = "No summary data available" if summary_df.empty else "No matches found"
        return main_df
    
    summary_clean = summary_df.rename(columns={
        "Country": "Country_standard",
        "Date": "Year"
    }).copy()
    
    summary_clean["Year"] = pd.to_numeric(summary_clean["Year"], errors="coerce")
    
    # Try multiple merge strategies
    strategies = [
        (["genotype", "Country_standard", "Year"], "exact_match"),
        (["genotype", "Country_standard"], "genotype_country"), 
        (["genotype"], "genotype_only"),
        (["Country_standard"], "country_only"),
    ]
    
    for merge_cols, strategy in strategies:
        if all(col in main_df.columns and col in summary_clean.columns for col in merge_cols):
            try:
                merged = main_df.merge(
                    summary_clean[merge_cols + ["Taxa"]].drop_duplicates(),
                    on=merge_cols,
                    how="left"
                )
                if merged["Taxa"].notna().any():
                    print(f"Taxa merge successful with {strategy}")
                    return merged
            except Exception as e:
                print(f"Taxa merge failed with {strategy}: {e}")
                continue
    
    main_df["Taxa"] = "No match found"
    return main_df
    
@callback(
    Output("download-data", "data", allow_duplicate=True),
    Input("btn-download", "n_clicks"),
    State("filtered-store", "data"),
    State("selected-virus", "data"),
    prevent_initial_call=True,
)
def download_main_data_simple(n_clicks, filtered_json, virus):
    if not n_clicks or n_clicks == 0:
        return dash.no_update
    
    filtered_df = _df_from_json(filtered_json)
    selected_virus = virus or "HBV"
    summary_key = "hbv_summary_raw" if selected_virus == "HBV" else "hcv_summary_raw"
    summary_data = data_store.get(summary_key, pd.DataFrame())
    
    # Merge Taxa information
    download_df = merge_taxa_information(filtered_df, summary_data, selected_virus)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{selected_virus}_sequence_data_with_taxa_{timestamp}.csv"
    
    return dcc.send_data_frame(download_df.to_csv, filename, index=False)
    


# === DYNAMIC ACCENT SWITCHING ===
@callback(
    Output("hep-dashboard-shell-container", "style"),
    Input("selected-virus", "data")
)
def update_dynamic_accent(virus):
    accent_map = {
        "HBV": "#E84057",
        "HCV": "#00D4AA",
        "HEV": "#8BC34A",
        "HAV": "#F5A623",
        "HDV": "#A259FF",
    }
    glow_map = {
        "HBV": "rgba(232, 64, 87, 0.25)",
        "HCV": "rgba(0, 212, 170, 0.25)",
        "HEV": "rgba(139, 195, 74, 0.25)",
        "HAV": "rgba(245, 166, 35, 0.25)",
        "HDV": "rgba(162, 89, 255, 0.25)",
    }
    selected = virus or "HBV"
    return {
        "--current-accent": accent_map.get(selected, "#E84057"),
        "--current-glow": glow_map.get(selected, "rgba(232, 64, 87, 0.25)"),
    }


# - Mutation section toggle -
@callback(
    Output("hbv-mutation-section", "style"),
    Output("hcv-mutation-section", "style"),
    Input("selected-virus", "data")
)
def toggle_mutation_sections(selected_virus):
    selected_virus = selected_virus or "HBV"
    if selected_virus == "HBV":
        return {"display": "block"}, {"display": "none"}
    elif selected_virus == "HCV":
        return {"display": "none"}, {"display": "block"}
    else:
        return {"display": "none"}, {"display": "none"}
        
# - Mutation filter options -
# DROP-IN REPLACEMENT
@callback(
    Output("mutation-filter-dropdown", "options"),
    Input("selected-virus", "data"),
    prevent_initial_call=False,
)
def update_mutation_filter_options(virus):
    data = get_data_store() or {}  # UPDATED
    selected = (virus or "HBV").upper()

    # Pick the right mutations table
    mut = data.get("hbv_mut") if selected == "HBV" else data.get("hcv_mut")
    if mut is None or len(mut) == 0:
        return []  # nothing to show yet

    # Column names differ: HBV uses 'drug', HCV uses 'gene'
    col = "drug" if selected == "HBV" else "gene"
    if col not in mut.columns:
        # Be defensive: try to discover a plausible column
        for candidate in ["drug", "gene", "Drug", "Gene"]:
            if candidate in mut.columns:
                col = candidate
                break
        else:
            return []

    opts = sorted(pd.Series(mut[col]).dropna().astype(str).unique())
    return [{"label": v, "value": v} for v in opts]


# === VIRUS DROPDOWN & RESET FILTER CALLBACKS ===
@callback(
    Output("selected-virus", "data", allow_duplicate=True),
    Input("virus-dropdown", "value"),
    prevent_initial_call=True
)
def sync_virus_dropdown_to_store(virus_val):
    return virus_val or "HBV"


@callback(
    Output("virus-dropdown", "value", allow_duplicate=True),
    Input("selected-virus", "data"),
    prevent_initial_call=True
)
def sync_store_to_virus_dropdown(virus_data):
    return virus_data or "HBV"


@callback(
    Output("virus-dropdown", "value", allow_duplicate=True),
    Output("year-range-slider", "value", allow_duplicate=True),
    Output("year-start-input", "value", allow_duplicate=True),
    Output("year-end-input", "value", allow_duplicate=True),
    Output("continent-dropdown", "value", allow_duplicate=True),
    Output("country-dropdown", "value", allow_duplicate=True),
    Output("genotype-dropdown", "value", allow_duplicate=True),
    Input("btn-reset-filters", "n_clicks"),
    prevent_initial_call=True
)
def reset_all_filters(n_clicks):
    if n_clicks:
        return "HBV", [1963, 2024], 1963, 2024, None, None, None
    return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update


# === YEAR SLIDER <-> MANUAL YEAR INPUTS SYNC ===
@callback(
    Output("year-start-input", "value"),
    Output("year-end-input", "value"),
    Input("year-range-slider", "value"),
    prevent_initial_call=False,
)
def sync_slider_to_year_inputs(range_val):
    if not range_val or len(range_val) != 2:
        return dash.no_update, dash.no_update
    return range_val[0], range_val[1]


@callback(
    Output("year-range-slider", "value", allow_duplicate=True),
    Input("year-start-input", "value"),
    Input("year-end-input", "value"),
    State("year-range-slider", "min"),
    State("year-range-slider", "max"),
    State("year-range-slider", "value"),
    prevent_initial_call=True,
)
def sync_year_inputs_to_slider(start_val, end_val, slider_min, slider_max, current_range):
    if start_val is None or end_val is None:
        return dash.no_update
    slider_min = slider_min if slider_min is not None else 1963
    slider_max = slider_max if slider_max is not None else 2024
    
    try:
        s = max(slider_min, min(int(start_val), slider_max))
        e = min(slider_max, max(int(end_val), slider_min))
        if s > e:
            s, e = e, s
        new_range = [s, e]
        if current_range and current_range == new_range:
            return dash.no_update
        return new_range
    except (ValueError, TypeError):
        return dash.no_update


# === REGION GENOTYPE HORIZONTAL BAR CHART CALLBACK ===
@callback(
    Output("region-genotype-bar-chart", "figure"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data")
)
def update_region_genotype_bar_chart(filtered_json, selected_virus):
    if not filtered_json:
        fig = go.Figure()
        fig.update_layout(
            annotations=[dict(text="No region data available", showarrow=False, font=dict(size=13, color="#64748B"))],
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)"
        )
        return fig

    try:
        df = _df_from_json(filtered_json)
    except Exception:
        df = pd.DataFrame()

    if df.empty or "WHO_Regions" not in df.columns or "genotype" not in df.columns:
        fig = go.Figure()
        fig.update_layout(
            annotations=[dict(text="No region data", showarrow=False, font=dict(size=13, color="#64748B"))],
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)"
        )
        return fig

    df_clean = df[df["WHO_Regions"].notna() & (df["WHO_Regions"] != "Unknown")].copy()
    if df_clean.empty:
        df_clean = df[df["WHO_Regions"].notna()].copy()

    if df_clean.empty:
        fig = go.Figure()
        fig.update_layout(
            annotations=[dict(text="No region data available", showarrow=False, font=dict(size=13, color="#64748B"))],
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)"
        )
        return fig

    region_abbrev = {
        "Eastern Mediterranean": "Mediterranean",
        "Western Pacific": "Western Pacific",
        "South-East Asia": "South-East Asia",
        "European": "Europe",
        "African": "Africa",
        "Americas": "Americas"
    }
    df_clean["WHO_Regions_Short"] = df_clean["WHO_Regions"].map(lambda x: region_abbrev.get(str(x), str(x)))
    grouped = df_clean.groupby(["WHO_Regions_Short", "genotype"]).size().reset_index(name="count")
    region_totals = df_clean.groupby("WHO_Regions_Short").size().sort_values(ascending=True)

    virus_str = (selected_virus or "HBV").upper()
    if virus_str == "HCV":
        color_map = HCV_GENOTYPE_COLORS
    elif virus_str == "HEV":
        color_map = HEV_GENOTYPE_COLORS
    else:
        color_map = HBV_GENOTYPE_COLORS

    fig = px.bar(
        grouped,
        y="WHO_Regions_Short",
        x="count",
        color="genotype",
        orientation="h",
        category_orders={"WHO_Regions_Short": region_totals.index.tolist()},
        labels={"WHO_Regions_Short": "Region", "count": "Sequences", "genotype": "Genotype"},
        color_discrete_map=color_map
    )

    fig.update_layout(
        barmode="stack",
        margin=dict(l=85, r=20, t=55, b=75),
        height=420,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        legend=dict(
            orientation="h",
            yanchor="top",
            y=10,
            xanchor="left",
            x=0,
            title=dict(text=""),
            font=dict(size=9, color="#475569")
        ),
        xaxis=dict(
            title=dict(text="Sequences", font=dict(size=11, color="#475569")),
            showgrid=True,
            gridcolor="#E2E8F0",
            tickfont=dict(color="#475569", size=11),
            tickangle=0
        ),
        yaxis=dict(
            title="",
            showgrid=False,
            tickfont=dict(color="#0F172A", size=11)
        )
    )

    return fig


# === HOMEPAGE EXECUTIVE SUMMARY CALLBACKS ===

@callback(
    Output("homepage-mutation-pie-graph", "figure"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
)
def update_homepage_mutation_pie_callback(filtered_json, selected_virus):
    df = _df_from_json(filtered_json)
    selected_virus = selected_virus or "HBV"
    if df is None or df.empty:
        return _empty_pie_chart("No sequence data")

    data = get_data_store()
    if selected_virus == "HBV":
        mut_df = data.get("hbv_mut", pd.DataFrame())
    elif selected_virus == "HCV":
        mut_df = data.get("hcv_mut", pd.DataFrame())
    elif selected_virus == "HEV":
        mut_df = data.get("hev_mut", pd.DataFrame())
    else:
        mut_df = pd.DataFrame()

    if mut_df.empty or "ID" not in mut_df.columns or "type" not in mut_df.columns:
        return _empty_pie_chart("No mutation markers")

    filtered_ids = set(df["ID"].dropna().unique()) if "ID" in df.columns else set()
    filtered_mut = mut_df[mut_df["ID"].isin(filtered_ids)] if filtered_ids else mut_df

    if filtered_mut.empty:
        return _empty_pie_chart("No mutations found")

    type_counts = filtered_mut.groupby("type")["ID"].nunique().to_dict()
    fig = create_mutation_type_pie_chart(type_counts, selected_virus)
    fig.update_layout(height=320, showlegend=True, margin=dict(t=20, b=20, l=10, r=10))
    return fig


@callback(
    Output("homepage-forecast-summary-graph", "figure"),
    Output("epi-forecast-summary-graph", "figure"),
    Input("selected-virus", "data"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
)
def update_homepage_forecast_callback(virus, regions, countries):
    data = get_data_store()
    ihme_df = data.get("ihme_df", pd.DataFrame())
    if ihme_df.empty:
        empty_fig = _empty_plot("No historical GBD data available")
        return empty_fig, empty_fig

    fig = create_forecast_chart(
        ihme_df=ihme_df,
        selected_virus=virus or "HBV",
        sex="Both",
        selected_regions=regions,
        selected_countries=countries
    )
    fig.update_layout(
        height=320,
        margin=dict(t=20, b=20, l=20, r=20),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="#1e293b"),
    )
    return fig, fig


@callback(
    Output("homepage-priority-table-summary", "children"),
    Input("priority-data-store", "data")
)
def update_homepage_priority_table_callback(priority_json):
    priority_df = _df_from_json(priority_json)
    if priority_df is None or priority_df.empty:
        return html.Div("No priority data available", className="text-muted p-3 text-center")

    display_df = priority_df.copy()
    columns_needed = ["rank", "Country_standard", "priority_score", "burden", "coverage_gap", "observed_sequences"]
    available_cols = [c for c in columns_needed if c in display_df.columns]
    display_df = display_df[available_cols].head(5).copy()

    rename_map = {
        "rank": "Rank",
        "Country_standard": "Country",
        "priority_score": "Priority Score",
        "burden": "Burden",
        "coverage_gap": "Coverage Gap",
        "observed_sequences": "Sequences",
    }
    display_df = display_df.rename(columns=rename_map)

    if "Priority Score" in display_df.columns:
        display_df["Priority Score"] = pd.to_numeric(display_df["Priority Score"], errors="coerce").round(2)
    for col in ["Burden", "Coverage Gap", "Sequences"]:
        if col in display_df.columns:
            display_df[col] = pd.to_numeric(display_df[col], errors="coerce").apply(
                lambda x: f"{x:,.0f}" if pd.notna(x) else "N/A"
            )

    return dash_table.DataTable(
        data=display_df.to_dict("records"),
        columns=[{"name": col, "id": col} for col in display_df.columns],
        style_table={"overflowX": "auto"},
        style_cell={
            "fontFamily": "IBM Plex Sans, sans-serif",
            "fontSize": "12px",
            "padding": "8px 10px",
            "textAlign": "left",
        },
        style_header={
            "backgroundColor": "#f8fafc",
            "fontWeight": "700",
            "color": "#334155",
            "borderBottom": "2px solid #e2e8f0",
        },
    )


# =============================================================================
# GENOTYPES TAB CALLBACKS
# =============================================================================
@callback(
    Output("geno-kpi-observed", "children"),
    Output("geno-kpi-dominant", "children"),
    Output("geno-kpi-diversity", "children"),
    Output("geno-kpi-shift", "children"),
    Output("geno-map-genotype-filter", "options"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
)
def update_genotype_kpi_summary(filtered_json, virus):
    df = _df_from_json(filtered_json)
    selected_virus = virus or "HBV"
    
    if df.empty or "genotype" not in df.columns:
        opts = [{"label": "All Genotypes", "value": "ALL"}]
        return "--", "--", "--", "--", opts
    
    df_valid = df[df["genotype"].notna() & (df["genotype"] != "Unknown")].copy()
    if df_valid.empty:
        opts = [{"label": "All Genotypes", "value": "ALL"}]
        return "--", "--", "--", "--", opts
        
    # 1. Observed Count
    observed_gts = sorted(df_valid["genotype"].astype(str).unique())
    observed_count = len(observed_gts)
    opts = [{"label": "All Genotypes", "value": "ALL"}] + [{"label": f"Genotype {g}", "value": g} for g in observed_gts]
    
    # 2. Dominant Genotype
    gt_counts = df_valid["genotype"].value_counts()
    if not gt_counts.empty:
        dom_gt = gt_counts.index[0]
        dom_pct = (gt_counts.iloc[0] / len(df_valid)) * 100
        dom_display = f"{dom_gt} ({dom_pct:.1f}%)"
    else:
        dom_display = "N/A"
        
    # 3. Most Diverse Region (Shannon Index)
    def calc_shannon(series):
        counts = series.value_counts()
        props = counts / counts.sum()
        return -np.sum(props * np.log(props + 1e-12))
        
    if "WHO_Regions" in df_valid.columns:
        region_div = df_valid.groupby("WHO_Regions")["genotype"].apply(calc_shannon)
        if not region_div.empty:
            top_reg = region_div.idxmax()
            top_shannon = region_div.max()
            diversity_display = f"{top_reg} ({top_shannon:.2f})"
        else:
            diversity_display = "N/A"
    else:
        diversity_display = "N/A"
        
    # 4. Recent Shift (>= 2020 vs < 2020)
    if "Year" in df_valid.columns:
        df_valid["Year_num"] = pd.to_numeric(df_valid["Year"], errors="coerce")
        hist = df_valid[df_valid["Year_num"] < 2020]["genotype"].value_counts(normalize=True)
        rec = df_valid[df_valid["Year_num"] >= 2020]["genotype"].value_counts(normalize=True)
        
        diffs = (rec - hist).dropna()
        if not diffs.empty:
            max_shift_gt = diffs.idxmax()
            shift_val = diffs.max() * 100
            shift_display = f"{max_shift_gt} ↑ {shift_val:+.1f}%"
        else:
            shift_display = "Stable"
    else:
        shift_display = "N/A"
        
    return f"{observed_count}", dom_display, diversity_display, shift_display, opts


@callback(
    Output("geno-distribution-map-graph", "figure"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
    Input("geno-map-mode", "value"),
    Input("geno-map-genotype-filter", "value"),
)
def update_geno_distribution_map(filtered_json, virus, map_mode, selected_gt):
    df = _df_from_json(filtered_json)
    if df.empty or "Country_standard" not in df.columns or "genotype" not in df.columns:
        return _empty_world("No genotype geographic data available")
        
    df_valid = df[df["genotype"].notna() & (df["genotype"] != "Unknown") & df["Country_standard"].notna()].copy()
    if df_valid.empty:
        return _empty_world("No genotype geographic data available")
        
    fig = go.Figure()
    
    if map_mode == "dominant":
        # Group by country, find most frequent genotype
        ctry_dom = df_valid.groupby("Country_standard")["genotype"].agg(lambda s: s.value_counts().index[0]).reset_index()
        unique_gts = sorted(ctry_dom["genotype"].unique())
        
        # Color scale index
        gt_to_idx = {g: i for i, g in enumerate(unique_gts)}
        ctry_dom["z_value"] = ctry_dom["genotype"].map(gt_to_idx)
        
        colors = _distinct_hue_colors(max(len(unique_gts), 1))
        custom_scale = []
        n_colors = len(unique_gts)
        for i, col in enumerate(colors):
            low = i / max(1, n_colors)
            high = (i + 1) / max(1, n_colors)
            custom_scale.append([low, col])
            custom_scale.append([high, col])
            
        fig.add_trace(go.Choropleth(
            locations=ctry_dom["Country_standard"],
            locationmode="country names",
            z=ctry_dom["z_value"],
            colorscale=custom_scale if len(unique_gts) > 1 else [[0, colors[0]], [1, colors[0]]],
            showscale=False,
            hovertemplate="<b>%{location}</b><br>Dominant Genotype: <b>%{text}</b><extra></extra>",
            text=ctry_dom["genotype"],
            marker_line_color="rgba(0,0,0,0.2)",
            marker_line_width=0.5,
        ))
        
    elif map_mode == "frequency":
        if selected_gt and str(selected_gt).upper() != "ALL":
            # Frequency of selected genotype per country
            target_gt = selected_gt
            counts = df_valid.groupby("Country_standard").agg(
                total=("genotype", "count"),
                target=("genotype", lambda s: (s == target_gt).sum())
            ).reset_index()
            counts["freq_pct"] = (counts["target"] / counts["total"]) * 100
            c_title = f"Genotype {target_gt} Share (%)"
            hover_lbl = f"Genotype {target_gt} Share"
        else:
            # Frequency (%) of the dominant genotype per country when "ALL" is selected
            def dom_freq(s):
                vc = s.value_counts()
                return (vc.iloc[0] / vc.sum()) * 100 if not vc.empty else 0
            
            counts = df_valid.groupby("Country_standard")["genotype"].agg(dom_freq).reset_index(name="freq_pct")
            c_title = "Dominant Genotype Share (%)"
            hover_lbl = "Dominant Genotype Share"

        fig.add_trace(go.Choropleth(
            locations=counts["Country_standard"],
            locationmode="country names",
            z=counts["freq_pct"],
            colorscale="Teal",
            colorbar=dict(
                title=dict(text=c_title, side="top", font=dict(size=10, color="#334155")),
                orientation="h",
                x=0.03,
                y=0.06,
                len=0.32,
                thickness=10,
                tickfont=dict(size=8, color="#475569")
            ),
            hovertemplate="<b>%{location}</b><br>" + hover_lbl + ": <b>%{z:.1f}%</b><extra></extra>",
            marker_line_color="rgba(0,0,0,0.2)",
            marker_line_width=0.5,
        ))
        
    else: # Diversity mode (Shannon index per country)
        def calc_shannon_ctry(s):
            counts = s.value_counts()
            props = counts / counts.sum()
            return -np.sum(props * np.log(props + 1e-12))
            
        ctry_div = df_valid.groupby("Country_standard")["genotype"].apply(calc_shannon_ctry).reset_index(name="shannon")
        
        fig.add_trace(go.Choropleth(
            locations=ctry_div["Country_standard"],
            locationmode="country names",
            z=ctry_div["shannon"],
            colorscale="Viridis",
            colorbar=dict(
                title=dict(text="Shannon Index (H')", side="top", font=dict(size=10, color="#334155")),
                orientation="h",
                x=0.03,
                y=0.06,
                len=0.32,
                thickness=10,
                tickfont=dict(size=8, color="#475569")
            ),
            hovertemplate="<b>%{location}</b><br>Shannon Diversity Index: <b>%{z:.2f}</b><extra></extra>",
            marker_line_color="rgba(0,0,0,0.2)",
            marker_line_width=0.5,
        ))
        
    fig.update_geos(
        projection_type="equirectangular",
        showcountries=True,
        countrycolor="rgba(0,0,0,0.2)",
        showsubunits=True,
        fitbounds=False,
        lataxis_range=[-65, 85],
        domain=dict(x=[0, 1], y=[0, 1]),
    )
    fig.update_layout(
        height=460, margin=dict(t=30, b=0, l=0, r=0), paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)"
    )
    return fig


@callback(
    Output("geno-region-heatmap-graph", "figure"),
    Output("geno-diversity-ranking-graph", "figure"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
    Input("geno-heatmap-mode", "value"),
    Input("geno-diversity-scope", "value"),
)
def update_geno_heatmap_and_diversity(filtered_json, virus, heat_mode, div_scope):
    df = _df_from_json(filtered_json)
    if df.empty or "genotype" not in df.columns:
        return _empty_plot("No genotype data available"), _empty_plot("No diversity data available")
        
    df_valid = df[df["genotype"].notna() & (df["genotype"] != "Unknown")].copy()
    if df_valid.empty:
        return _empty_plot("No genotype data available"), _empty_plot("No diversity data available")
        
    # 1. HEATMAP (Genotype x WHO Region)
    if "WHO_Regions" in df_valid.columns:
        pivot = pd.crosstab(df_valid["genotype"], df_valid["WHO_Regions"])
        if heat_mode == "pct":
            pivot_pct = pivot.div(pivot.sum(axis=0), axis=1) * 100
            z_data = pivot_pct.values
            hover_template = "Genotype: <b>%{y}</b><br>Region: <b>%{x}</b><br>Share: <b>%{z:.1f}%</b><extra></extra>"
            colorbar_title = "Share (%)"
        else:
            z_data = pivot.values
            hover_template = "Genotype: <b>%{y}</b><br>Region: <b>%{x}</b><br>Sequences: <b>%{z:,}</b><extra></extra>"
            colorbar_title = "Sequences"
            
        fig_heat = go.Figure(go.Heatmap(
            z=z_data,
            x=pivot.columns.tolist(),
            y=pivot.index.tolist(),
            colorscale="YlGnBu",
            colorbar=dict(title=dict(text=colorbar_title, side="top")),
            hovertemplate=hover_template
        ))
        fig_heat.update_layout(
            height=360, margin=dict(t=20, b=40, l=60, r=20),
            paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
            xaxis=dict(title="WHO Region"), yaxis=dict(title="Genotype")
        )
    else:
        fig_heat = _empty_plot("WHO Region information missing")
        
    # 2. DIVERSITY RANKING
    def calc_shannon_val(s):
        counts = s.value_counts()
        props = counts / counts.sum()
        return -np.sum(props * np.log(props + 1e-12))
        
    col_scope = "WHO_Regions" if div_scope == "region" and "WHO_Regions" in df_valid.columns else "Country_standard"
    if col_scope in df_valid.columns:
        div_scores = df_valid.groupby(col_scope)["genotype"].apply(calc_shannon_val).reset_index(name="shannon")
        div_scores = div_scores.sort_values("shannon", ascending=True).tail(10)
        
        fig_div = go.Figure(go.Bar(
            x=div_scores["shannon"],
            y=div_scores[col_scope],
            orientation="h",
            marker=dict(color="#3B82F6", cornerradius=4),
            hovertemplate="<b>%{y}</b><br>Shannon Diversity Index: <b>%{x:.2f}</b><extra></extra>"
        ))
        fig_div.update_layout(
            height=320, margin=dict(t=20, b=30, l=120, r=20),
            paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
            xaxis=dict(title="Shannon Index (H')"), yaxis=dict(title="")
        )
    else:
        fig_div = _empty_plot("No diversity hierarchy available")
        
    return fig_heat, fig_div


@callback(
    Output("geno-scatter-signals-graph", "figure"),
    Output("geno-signals-table", "children"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
)
def update_geno_signals(filtered_json, virus):
    df = _df_from_json(filtered_json)
    if df.empty or "genotype" not in df.columns or "Year" not in df.columns:
        return _empty_plot("No temporal genotype signals available"), html.Div("No signals data", className="text-muted p-3")
        
    df_valid = df[df["genotype"].notna() & (df["genotype"] != "Unknown")].copy()
    df_valid["Year_num"] = pd.to_numeric(df_valid["Year"], errors="coerce")
    df_valid = df_valid.dropna(subset=["Year_num"])
    
    if df_valid.empty:
        return _empty_plot("No temporal genotype signals available"), html.Div("No signals data", className="text-muted p-3")
        
    hist_df = df_valid[df_valid["Year_num"] < 2020]
    rec_df = df_valid[df_valid["Year_num"] >= 2020]
    
    hist_props = (hist_df["genotype"].value_counts() / len(hist_df) * 100) if not hist_df.empty else pd.Series(dtype=float)
    rec_props = (rec_df["genotype"].value_counts() / len(rec_df) * 100) if not rec_df.empty else pd.Series(dtype=float)
    
    all_gts = sorted(list(set(hist_props.index).union(set(rec_props.index))))
    
    signals_data = []
    for g in all_gts:
        h_val = float(hist_props.get(g, 0.0))
        r_val = float(rec_props.get(g, 0.0))
        diff = r_val - h_val
        rec_count = int((rec_df["genotype"] == g).sum())
        
        status = "Increasing" if diff > 1.0 else ("Decreasing" if diff < -1.0 else "Stable")
        signals_data.append({
            "Genotype": g,
            "Historical %": round(h_val, 1),
            "Recent %": round(r_val, 1),
            "Shift %": round(diff, 1),
            "Recent Count": rec_count,
            "Status": status
        })
        
    sig_df = pd.DataFrame(signals_data)
    
    # 1. SCATTER PLOT
    fig_scatter = go.Figure()
    max_val = max(sig_df["Historical %"].max(), sig_df["Recent %"].max(), 10.0) if not sig_df.empty else 100.0
    
    # Diagonal reference y=x
    fig_scatter.add_trace(go.Scatter(
        x=[0, max_val], y=[0, max_val],
        mode="lines", name="No Change (y=x)",
        line=dict(color="rgba(148,163,184,0.5)", dash="dash")
    ))
    
    if not sig_df.empty:
        fig_scatter.add_trace(go.Scatter(
            x=sig_df["Historical %"],
            y=sig_df["Recent %"],
            mode="markers+text",
            text=sig_df["Genotype"],
            textposition="top center",
            marker=dict(
                size=12,
                color=sig_df["Shift %"],
                colorscale="Tropic",
                showscale=True,
                colorbar=dict(title=dict(text="Shift %", side="top"))
            ),
            hovertemplate="<b>Genotype %{text}</b><br>Historical Share: %{x:.1f}%<br>Recent Share: %{y:.1f}%<br>Shift: %{marker.color:+.1f}%<extra></extra>"
        ))
        
    fig_scatter.update_layout(
        height=350, margin=dict(t=30, b=40, l=50, r=20),
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(title="Historical Frequency % (pre-2020)"),
        yaxis=dict(title="Recent Frequency % (>=2020)")
    )
    
    # 2. SIGNALS TABLE
    table_df = sig_df.sort_values("Shift %", ascending=False).head(6)
    
    table_comp = dash_table.DataTable(
        data=table_df.to_dict("records"),
        columns=[{"name": col, "id": col} for col in table_df.columns],
        style_table={"overflowX": "auto"},
        style_cell={"fontFamily": "IBM Plex Sans, sans-serif", "fontSize": "11px", "padding": "6px 8px", "textAlign": "left"},
        style_header={"backgroundColor": "#f8fafc", "fontWeight": "700", "fontSize": "11px"},
        style_data_conditional=[
            {"if": {"filter_query": '{Status} = "Increasing"', "column_id": "Status"}, "color": "#10B981", "fontWeight": "bold"},
            {"if": {"filter_query": '{Status} = "Decreasing"', "column_id": "Status"}, "color": "#EF4444", "fontWeight": "bold"},
        ]
    )
    
    return fig_scatter, table_comp


# =============================================================================
# MUTATIONS TAB CALLBACKS
# =============================================================================
from itertools import combinations

def _ensure_column(mut_df, seq_df, col):
    """Ensures `col` exists in mut_df by bringing it from seq_df on 'ID' if needed, avoiding duplicate columns."""
    if mut_df.empty:
        return mut_df
    if col in mut_df.columns:
        return mut_df
    if "ID" in mut_df.columns and "ID" in seq_df.columns and col in seq_df.columns:
        right = seq_df[["ID", col]].drop_duplicates()
        return mut_df.merge(right, on="ID", how="inner")
    return mut_df


@callback(
    Output("mut-kpi-unique", "children"),
    Output("mut-kpi-high-freq", "children"),
    Output("mut-kpi-countries", "children"),
    Output("mut-kpi-signals", "children"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
)
def update_mutation_kpi_summary(filtered_json, virus):
    data = get_data_store()
    selected_virus = virus or "HBV"
    df = _df_from_json(filtered_json)
    
    mut_data = data.get("hbv_mut") if selected_virus == "HBV" else data.get("hcv_mut")
    if mut_data is None or mut_data.empty or df.empty:
        return "--", "--", "--", "--"
        
    if "ID" in mut_data.columns and "ID" in df.columns:
        filtered_mut = mut_data[mut_data["ID"].isin(df["ID"].unique())]
    else:
        filtered_mut = mut_data
        
    if filtered_mut.empty or "mutation" not in filtered_mut.columns:
        return "0", "0", "0", "0"
        
    # 1. Unique Mutations
    unique_muts = filtered_mut["mutation"].dropna().unique()
    unique_count = len(unique_muts)
    
    # 2. High Frequency (>5% share)
    total_seqs = df["ID"].nunique() if "ID" in df.columns else max(len(df), 1)
    counts = filtered_mut.groupby("mutation")["ID"].nunique() if "ID" in filtered_mut.columns else filtered_mut["mutation"].value_counts()
    freq_pcts = (counts / total_seqs) * 100
    high_freq_count = (freq_pcts > 5.0).sum()
    
    # 3. Countries Affected
    if "ID" in filtered_mut.columns and "Country_standard" in df.columns and "ID" in df.columns:
        mut_ids = set(filtered_mut["ID"].dropna().unique())
        affected_ctries = df[df["ID"].isin(mut_ids)]["Country_standard"].dropna().nunique()
    elif "Country_standard" in df.columns:
        affected_ctries = df["Country_standard"].dropna().nunique()
    else:
        affected_ctries = 0
        
    # 4. Emerging Signals (>=2020 vs <2020)
    merged = _ensure_column(filtered_mut, df, "Year")
    if not merged.empty and "Year" in merged.columns and "ID" in merged.columns:
        merged["Year_num"] = pd.to_numeric(merged["Year"], errors="coerce")
        
        hist_ids = df[df["Year"] < 2020]["ID"].unique() if ("Year" in df.columns and "ID" in df.columns) else []
        rec_ids = df[df["Year"] >= 2020]["ID"].unique() if ("Year" in df.columns and "ID" in df.columns) else []
        
        hist_n = max(len(hist_ids), 1)
        rec_n = max(len(rec_ids), 1)
        
        hist_counts = merged[merged["ID"].isin(hist_ids)].groupby("mutation")["ID"].nunique() / hist_n * 100
        rec_counts = merged[merged["ID"].isin(rec_ids)].groupby("mutation")["ID"].nunique() / rec_n * 100
        
        diff = (rec_counts - hist_counts).dropna()
        emerging_count = (diff >= 1.0).sum()
    else:
        emerging_count = 0
        
    return f"{unique_count:,}", f"{high_freq_count}", f"{affected_ctries}", f"{emerging_count}"


@callback(
    Output("mut-geo-map-graph", "figure"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
)
def update_mut_geo_map(filtered_json, virus):
    data = get_data_store()
    selected_virus = virus or "HBV"
    df = _df_from_json(filtered_json)
    mut_data = data.get("hbv_mut") if selected_virus == "HBV" else data.get("hcv_mut")
    
    if mut_data is None or mut_data.empty or df.empty:
        return _empty_world("No geographic mutation data available")
        
    if "ID" in mut_data.columns and "ID" in df.columns:
        filtered_mut = mut_data[mut_data["ID"].isin(df["ID"].unique())]
    else:
        filtered_mut = mut_data
        
    merged = _ensure_column(filtered_mut, df, "Country_standard")
    if merged.empty or "Country_standard" not in merged.columns:
        return _empty_world("No geographic mutation data available")
        
    ctry_muts = merged.groupby("Country_standard").agg(
        total_mutations=("mutation", "count"),
        unique_mutations=("mutation", "nunique")
    ).reset_index()
    
    fig = go.Figure(go.Choropleth(
        locations=ctry_muts["Country_standard"],
        locationmode="country names",
        z=ctry_muts["total_mutations"],
        colorscale="Purples",
        colorbar=dict(title=dict(text="Mutations", side="top", font=dict(size=10)), orientation="h", x=0.02, y=0.02, len=0.34, thickness=8, tickfont=dict(size=8)),
        hovertemplate="<b>%{location}</b><br>Detected Mutations: <b>%{z:,}</b><br>Unique Variants: <b>%{text:,}</b><extra></extra>",
        text=ctry_muts["unique_mutations"],
        marker_line_color="rgba(0,0,0,0.2)",
        marker_line_width=0.5,
    ))
    
    fig.update_geos(
        projection_type="equirectangular",
        showcountries=True,
        countrycolor="rgba(0,0,0,0.2)",
        showsubunits=True,
        fitbounds=False,
        lataxis_range=[-65, 85],
        domain=dict(x=[0, 1], y=[0, 1]),
    )
    fig.update_layout(
        height=230, margin=dict(t=10, b=0, l=0, r=0), paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)"
    )
    return fig


@callback(
    Output("mut-lollipop-graph", "figure"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
    Input("mut-protein-selector", "value"),
)
def update_mut_lollipop(filtered_json, virus, protein_filter):
    data = get_data_store()
    selected_virus = virus or "HBV"
    df = _df_from_json(filtered_json)
    mut_data = data.get("hbv_mut") if selected_virus == "HBV" else data.get("hcv_mut")
    
    if mut_data is None or mut_data.empty or df.empty or "mutation" not in mut_data.columns:
        return _empty_plot("No protein position data available")
        
    if "ID" in mut_data.columns and "ID" in df.columns:
        filtered_mut = mut_data[mut_data["ID"].isin(df["ID"].unique())].copy()
    else:
        filtered_mut = mut_data.copy()
        
    if filtered_mut.empty:
        return _empty_plot("No protein position data available")
        
    col = "drug" if selected_virus == "HBV" and "drug" in filtered_mut.columns else ("gene" if "gene" in filtered_mut.columns else None)
    if protein_filter != "ALL" and col and col in filtered_mut.columns:
        filtered_mut = filtered_mut[filtered_mut[col] == protein_filter]
        
    if filtered_mut.empty:
        return _empty_plot(f"No mutations found for {protein_filter}")
        
    # Extract amino-acid position number
    def parse_pos(m):
        match = re.search(r'\d+', str(m))
        return int(match.group()) if match else None
        
    filtered_mut["pos"] = filtered_mut["mutation"].apply(parse_pos)
    filtered_mut = filtered_mut.dropna(subset=["pos"])
    
    if filtered_mut.empty:
        return _empty_plot("Could not parse amino-acid position numbers")
        
    total_seqs = max(df["ID"].nunique() if "ID" in df.columns else len(df), 1)
    pos_counts = filtered_mut.groupby(["pos", "mutation"]).agg(
        count=("ID", "nunique") if "ID" in filtered_mut.columns else ("mutation", "count")
    ).reset_index()
    pos_counts["freq_pct"] = (pos_counts["count"] / total_seqs) * 100
    # Filter out mutations with frequency less than 5%
    pos_counts = pos_counts[pos_counts["freq_pct"] >= 5.0]
    
    if pos_counts.empty:
        return _empty_plot("No mutations found with frequency ≥ 5%")
        
    pos_counts = pos_counts.sort_values("pos")
    
    fig = go.Figure()
    
    # Add vertical lollipop stems
    for _, row in pos_counts.iterrows():
        fig.add_shape(
            type="line",
            x0=row["pos"], y0=0,
            x1=row["pos"], y1=row["freq_pct"],
            line=dict(color="#64748B", width=1.5)
        )
        
    # Add lollipop head markers
    fig.add_trace(go.Scatter(
            x=pos_counts["pos"],
            y=pos_counts["freq_pct"],
            mode="markers",
            text=pos_counts["mutation"],
            marker=dict(
                size=np.clip(pos_counts["freq_pct"] * 2 + 8, 8, 24),
                color=pos_counts["freq_pct"],
                colorscale="Viridis",
                showscale=True,
                colorbar=dict(title=dict(text="Share %", side="top"))
            ),
            hovertemplate="<b>Mutation: %{text}</b><br>AA Position: <b>%{x}</b><br>Frequency Share: <b>%{y:.2f}%</b><extra></extra>"
        ))
        
    _label_declutter(fig, pos_counts["pos"].tolist(), pos_counts["freq_pct"].tolist(),
                      pos_counts["mutation"].tolist(), min_gap_frac=0.035)
    
    fig.update_layout(
        height=360, margin=dict(t=30, b=40, l=50, r=20),
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(title="Amino Acid Position", showgrid=True, gridcolor="#E2E8F0"),
        yaxis=dict(title="Frequency Share (%)", showgrid=True, gridcolor="#E2E8F0")
    )
    return fig


@callback(
    Output("mut-priority-table", "children"),
    Output("mut-genotype-heatmap-graph", "figure"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
)
def update_mut_priority_and_heatmap(filtered_json, virus):
    data = get_data_store()
    selected_virus = virus or "HBV"
    df = _df_from_json(filtered_json)
    mut_data = data.get("hbv_mut") if selected_virus == "HBV" else data.get("hcv_mut")
    
    if mut_data is None or mut_data.empty or df.empty or "mutation" not in mut_data.columns:
        return html.Div("No priority mutations available", className="text-muted p-3"), _empty_plot("No genotype association data")
        
    if "ID" in mut_data.columns and "ID" in df.columns:
        filtered_mut = mut_data[mut_data["ID"].isin(df["ID"].unique())].copy()
    else:
        filtered_mut = mut_data.copy()
        
    merged = _ensure_column(filtered_mut, df, "genotype")
        
    # 1. PRIORITY MUTATION TABLE
    total_seqs = max(df["ID"].nunique() if "ID" in df.columns else len(df), 1)
    if "type" in filtered_mut.columns:
        prio_df = filtered_mut[filtered_mut["type"].isin(["antiviral_resistance", "vaccine_escape"])].copy()
    else:
        prio_df = filtered_mut.copy()
        
    if not prio_df.empty:
        col_target = "drug" if selected_virus == "HBV" and "drug" in prio_df.columns else ("gene" if "gene" in prio_df.columns else "type")
        summary_prio = prio_df.groupby(["mutation", col_target]).agg(
            count=("ID", "nunique") if "ID" in prio_df.columns else ("mutation", "count")
        ).reset_index()
        summary_prio["Frequency %"] = (summary_prio["count"] / total_seqs * 100).round(1)
        summary_prio = summary_prio.sort_values("Frequency %", ascending=False).head(6)
        summary_prio["Significance"] = summary_prio[col_target].astype(str).str.title()
        
        table_prio = dash_table.DataTable(
            data=summary_prio[["mutation", "Significance", "Frequency %"]].to_dict("records"),
            columns=[{"name": col, "id": col} for col in ["mutation", "Significance", "Frequency %"]],
            style_table={"overflowX": "auto"},
            style_cell={"fontFamily": "IBM Plex Sans, sans-serif", "fontSize": "11px", "padding": "6px 8px", "textAlign": "left"},
            style_header={"backgroundColor": "#f8fafc", "fontWeight": "700", "fontSize": "11px"}
        )
    else:
        table_prio = html.Div("No priority resistance/escape variants detected", className="text-muted p-3")
        
    # 2. MUTATION X GENOTYPE HEATMAP
    if not merged.empty and "genotype" in merged.columns:
        merged_clean = merged[merged["genotype"].notna() & (merged["genotype"] != "Unknown")].copy()
        top_muts = merged_clean["mutation"].value_counts().head(10).index.tolist()
        merged_top = merged_clean[merged_clean["mutation"].isin(top_muts)]
        
        if not merged_top.empty:
            pivot = pd.crosstab(merged_top["mutation"], merged_top["genotype"])
            pivot_pct = pivot.div(pivot.sum(axis=0), axis=1) * 100
            
            fig_heat = go.Figure(go.Heatmap(
                z=pivot_pct.values,
                x=pivot_pct.columns.tolist(),
                y=pivot_pct.index.tolist(),
                colorscale="YlOrRd",
                colorbar=dict(title=dict(text="Share (%)", side="top")),
                hovertemplate="Mutation: <b>%{y}</b><br>Genotype: <b>%{x}</b><br>Prevalence Share: <b>%{z:.1f}%</b><extra></extra>"
            ))
            fig_heat.update_layout(
                height=350, margin=dict(t=20, b=40, l=60, r=20),
                paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                xaxis=dict(title="Viral Genotype"), yaxis=dict(title="Mutation")
            )
        else:
            fig_heat = _empty_plot("No genotype association data")
    else:
        fig_heat = _empty_plot("Genotype metadata not linked")
        
    return table_prio, fig_heat


@callback(
    Output("mut-cooccurrence-table", "children"),
    Output("mut-emerging-signals-graph", "figure"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
)
def update_mut_cooccurrence_and_signals(filtered_json, virus):
    data = get_data_store()
    selected_virus = virus or "HBV"
    df = _df_from_json(filtered_json)
    mut_data = data.get("hbv_mut") if selected_virus == "HBV" else data.get("hcv_mut")
    
    if mut_data is None or mut_data.empty or df.empty or "mutation" not in mut_data.columns:
        return html.Div("No co-occurrence data", className="text-muted p-3"), _empty_plot("No emerging signals data")
        
    if "ID" in mut_data.columns and "ID" in df.columns:
        filtered_mut = mut_data[mut_data["ID"].isin(df["ID"].unique())].copy()
    else:
        filtered_mut = mut_data.copy()
        
    # 1. CO-OCCURRING MUTATIONS TABLE
    if "ID" in filtered_mut.columns:
        grouped = filtered_mut.groupby("ID")["mutation"].apply(lambda s: sorted(list(set(s.dropna()))))
        pair_counts = {}
        for mut_list in grouped:
            if len(mut_list) >= 2:
                for pair in combinations(mut_list, 2):
                    pair_str = f"{pair[0]} + {pair[1]}"
                    pair_counts[pair_str] = pair_counts.get(pair_str, 0) + 1
                    
        if pair_counts:
            pair_df = pd.DataFrame(list(pair_counts.items()), columns=["Mutation Pair", "Genomes Co-Detected"]).sort_values("Genomes Co-Detected", ascending=False).head(6)
            table_co = dash_table.DataTable(
                data=pair_df.to_dict("records"),
                columns=[{"name": col, "id": col} for col in pair_df.columns],
                style_table={"overflowX": "auto"},
                style_cell={"fontFamily": "IBM Plex Sans, sans-serif", "fontSize": "11px", "padding": "6px 8px", "textAlign": "left"},
                style_header={"backgroundColor": "#f8fafc", "fontWeight": "700", "fontSize": "11px"}
            )
        else:
            table_co = html.Div("No co-occurring mutation pairs detected", className="text-muted p-3")
    else:
        table_co = html.Div("No sequence ID mapping available", className="text-muted p-3")
        
    # 2. EMERGING MUTATION SIGNALS SCATTER
    merged = _ensure_column(filtered_mut, df, "Year")
    if not merged.empty and "Year" in merged.columns and "ID" in merged.columns:
        merged["Year_num"] = pd.to_numeric(merged["Year"], errors="coerce")
        merged = merged.dropna(subset=["Year_num"])
        
        hist_df = merged[merged["Year_num"] < 2020]
        rec_df = merged[merged["Year_num"] >= 2020]
        
        hist_n = max(df[df["Year"] < 2020]["ID"].nunique(), 1) if ("Year" in df.columns and "ID" in df.columns) else 1
        rec_n = max(df[df["Year"] >= 2020]["ID"].nunique(), 1) if ("Year" in df.columns and "ID" in df.columns) else 1
        
        hist_props = (hist_df.groupby("mutation")["ID"].nunique() / hist_n * 100) if not hist_df.empty else pd.Series(dtype=float)
        rec_props = (rec_df.groupby("mutation")["ID"].nunique() / rec_n * 100) if not rec_df.empty else pd.Series(dtype=float)
        
        all_muts = sorted(list(set(hist_props.index).union(set(rec_props.index))))
        signals = []
        for m in all_muts:
            h_v = float(hist_props.get(m, 0.0))
            r_v = float(rec_props.get(m, 0.0))
            diff = r_v - h_v
            signals.append({"Mutation": m, "Historical %": round(h_v, 2), "Recent %": round(r_v, 2), "Shift %": round(diff, 2)})
            
        sig_df = pd.DataFrame(signals)
        fig_sig = go.Figure()
        max_val = max(sig_df["Historical %"].max(), sig_df["Recent %"].max(), 5.0) if not sig_df.empty else 10.0
        
        fig_sig.add_trace(go.Scatter(
            x=[0, max_val], y=[0, max_val],
            mode="lines", name="No Change (y=x)",
            line=dict(color="rgba(148,163,184,0.5)", dash="dash")
        ))
        
        if not sig_df.empty:
                    fig_sig.add_trace(go.Scatter(
                        x=sig_df["Historical %"],
                        y=sig_df["Recent %"],
                        mode="markers",
                        text=sig_df["Mutation"],
                        marker=dict(
                            size=12,
                            color=sig_df["Shift %"],
                            colorscale="Reds",
                            showscale=True,
                            colorbar=dict(title=dict(text="Shift %", side="top"))
                        )),
                        
        _label_declutter(
                        fig_sig,
                        sig_df["Historical %"].tolist(),
                        sig_df["Recent %"].tolist(),
                        sig_df["Mutation"].tolist(),
                        top_n=12,
                        min_gap_frac=0.05,
                        y_key=lambda i: abs(sig_df["Shift %"].iloc[i]),  # prioritize by biggest shift, not raw position
                    )
                )
            
        fig_sig.update_layout(
            height=360, margin=dict(t=30, b=40, l=50, r=20),
            paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
            xaxis=dict(title="Historical Frequency Share % (pre-2020)"),
            yaxis=dict(title="Recent Frequency Share % (>=2020)")
        )
    else:
        fig_sig = _empty_plot("No temporal sequence metadata available")
        
    return table_co, fig_sig