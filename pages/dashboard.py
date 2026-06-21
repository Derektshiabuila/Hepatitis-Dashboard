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
    register_theme, shade,
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
    if cause_filter not in ihme_df["cause"].values:
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

    # --- Defensive copy (CRITICAL to avoid UnboundLocalError) ---
    mut = mutation_df.copy()

    # --- Ensure we have an ID column ---
    if "ID" not in mut.columns:
        if "sample" in mut.columns:
            mut["ID"] = mut["sample"]
        elif "Taxa" in mut.columns:
            mut["ID"] = mut["Taxa"]
        else:
            raise KeyError(
                "Mutation DF has no ID/sample/Taxa to join on. Columns: "
                f"{mut.columns.tolist()}"
            )

    # --- Normalize ID just in case ---
    mut["ID"] = mut["ID"].astype(str).str.strip()

    # --- Ensure sequence_df has ID ---
    if "ID" not in sequence_df.columns:
        raise KeyError("Sequence DF has no ID column")

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

def create_world_map(
    country_data: pd.DataFrame, 
    country_genotype_counts: pd.DataFrame, 
    coord_lookup: dict[str, dict[str, float]], 
    virus_type: str = "HBV", 
    display_mode: str = "raw",  # "raw", "PerMillion", or "ihme"
    map_title: str = ""  # Add this to know what metric we're showing
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
            colorscale = coverage_colorscale
        elif virus_type == "HCV":
            colorscale = coverage_colorscale
        else:
            colorscale = coverage_colorscale
        
        # Determine colorbar title from map_title
        if "Prevalence" in map_title:
            colorbar_title = "Prevalence (Log10)"
        elif "Incidence" in map_title:
            colorbar_title = "Incidence (Log10)"
        elif "Deaths" in map_title:
            colorbar_title = "Deaths (Log10)"
        else:
            colorbar_title = "Value (Log10)"
        
        # Generate ticks for colorbar
        tick_min = int(np.floor(vmin))
        tick_max = int(np.ceil(vmax))
        tick_vals = list(range(tick_min, tick_max + 1))
        tick_text = [f"10^{x}" for x in tick_vals]
        
        fig.add_trace(
            go.Choropleth(
                locations=valid_nonzero["Country_standard"],
                locationmode="country names",
                z=valid_nonzero["log_value"],
                zmin=vmin,
                zmax=vmax,
                colorscale=colorscale,
                colorbar=dict(
                    title=colorbar_title,
                    len=0.6,
                    thickness=20,
                    tickvals=tick_vals,
                    ticktext=tick_text,
                ),
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
    
    # PER MILLION MODE
    elif display_mode == "PerMillion":
        # [Keep existing PerMillion code as is]
        z_vals = (
            valid["Metric_raw"].apply(lambda x: np.log10(x) if (pd.notna(x) and x > 0) else np.nan)
            .astype(float)
            .to_numpy()
        )
        if np.all(np.isnan(z_vals)):
            return _empty_world("No per‑million values available for current filters")
        vmin = float(np.nanmin(z_vals)) if np.isfinite(np.nanmin(z_vals)) else -5.0
        vmax = float(np.nanmax(z_vals)) if np.isfinite(np.nanmax(z_vals)) else 0.0


        if virus_type == "HBV":
            colorscale = coverage_colorscale
        elif virus_type == "HCV":
            colorscale = coverage_colorscale
        else:
            colorscale = coverage_colorscale

        fig.add_trace(
            go.Choropleth(
                locations=valid["Country_standard"],
                locationmode="country names",
                z=z_vals,
                zmin=vmin,
                zmax=vmax,
                colorscale=colorscale,
                colorbar_title="Log10 per million",
                marker_line_color="rgba(0,0,0,0.3)",
                marker_line_width=0.5,
                hovertemplate="<b>%{location}</b><br>Per million: %{customdata:.2f}<extra></extra>",
                customdata=valid["Metric_raw"].astype(float),
            )
        )
    
    # RAW COUNTS MODE (default)
    else:
        # [Keep existing raw counts code as is]
        bins = [0, 1, 5, 20, 100, 500, 2000, 4000, float("inf")]
        labels = ["0", "1–4", "5–19", "20–99", "100–499", "500–1,999", "2,000–3,999", "4,000+"]
        valid["bin"] = pd.cut(valid["Metric_raw"], bins=bins, labels=labels, include_lowest=True, right=False)
        bin_to_idx = {lab: i for i, lab in enumerate(labels)}
        valid["z_value"] = valid["bin"].map(bin_to_idx)

        nonzero = valid[valid["Metric_raw"] > 0].copy()
        driving = nonzero if not nonzero.empty else valid
        z_numeric = driving["z_value"].astype(float).to_numpy()
        locations = driving["Country_standard"]

        HCV_COLORS = ["#FFF7BC", "#FEE391", "#FEC44F", "#FE9929", "#EC7014", "#CC4C02", "#993404", "#662506"]
        HBV_COLORS = ["#F7FCF0", "#E0F3DB", "#A8DDB5", "#4EB3D3", "#2B8CBE", "#0868AC", "#084081", "#06214D"]
        HEV_COLORS = ["#F7FCF0", "#E0F3DB", "#A8DDB5", "#4EB3D3", "#2B8CBE", "#0868AC", "#084081", "#06214D"]
        
        colorscale = coverage_colorscale
        
        #colorscale = [[i / (len(labels) - 1), c] for i, c in enumerate(colors)]

        fig.add_trace(
            go.Choropleth(
                locations=locations,
                locationmode="country names",
                z=z_numeric,
                zmin=0,
                zmax=len(labels) - 1,
                colorscale=colorscale,
                showscale=True,
                colorbar=dict(
                    tickvals=list(range(len(labels))),
                    ticktext=labels,
                    title="Sequence count",
                    len=0.6,
                    thickness=20,
                ),
                marker_line_color="rgba(0,0,0,0.3)",
                marker_line_width=0.5,
                hovertext=driving.apply(
                    lambda r: f"<b>{r['Country_standard']}</b><br>Exact count: {float(r['Metric_raw']):.0f}<br>Range: {r['bin']}",
                    axis=1,
                ),
                hoverinfo="text",
            )
        )

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
        projection_type="natural earth",
        showcountries=True,
        countrycolor="rgba(0,0,0,0.2)",
        showsubunits=True,
        fitbounds="locations",
        domain=dict(x=[0, 1], y=[0.2, 1]),
    )
    fig.update_layout(height=700, margin=dict(t=30, b=30, l=10, r=10))
    return fig
    
def create_coverage_map(
    cov_df: pd.DataFrame,
    coord_lookup: dict[str, any],
    coords_df: pd.DataFrame | None = None,
    virus_type: str = "HBV",
    who_regions: list[str] | None = None,
    countries: list[str] | None = None,
) -> go.Figure:
    
    if cov_df is None or cov_df.empty:
        print("Warning: cov_df is empty")
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
    
    # Remove rows with zero or missing coverage
    plot_data = valid[valid["coverage_ratio"] > 0].copy()
    
    if plot_data.empty:
        return _empty_world("No coverage ratio data available for selected filters")
    
    # Calculate actual percentiles from the data
    coverage_values = plot_data["coverage_ratio"].values
    
    # Create meaningful bins based on actual distribution
    # Use unique quantiles that actually exist in the data
    unique_quantiles = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    bin_edges = np.quantile(coverage_values, unique_quantiles)
    
    # Ensure bin edges are unique (if data is highly clustered)
    bin_edges = np.unique(bin_edges)
    if len(bin_edges) < 3:
        # If data is extremely clustered, use linear spacing between min and max
        bin_edges = np.linspace(coverage_values.min(), coverage_values.max(), 5)
    
    # Create color mapping based on bins
    plot_data["coverage_bin"] = pd.cut(
        plot_data["coverage_ratio"], 
        bins=bin_edges, 
        include_lowest=True,
        labels=False
    )
    
    # Normalize bin indices to 0-1 for colorscale
    max_bin = plot_data["coverage_bin"].max()
    plot_data["coverage_display"] = plot_data["coverage_bin"] / max_bin if max_bin > 0 else 0
    
    # Create tick labels showing actual percentage ranges
    tick_positions = []
    tick_labels = []
    
    for i in range(len(bin_edges) - 1):
        low = bin_edges[i] * 100
        high = bin_edges[i + 1] * 100
        if i == len(bin_edges) - 2:
            tick_labels.append(f"{low:.2f}% – {high:.2f}%")
        else:
            tick_labels.append(f"{low:.2f}% – {high:.2f}%")
        tick_positions.append(i / max_bin if max_bin > 0 else 0)
    
    # Colorscale from red to green
    coverage_colorscale = [
        [0.0, "#1a9850"],   # Red (lowest coverage)
        [0.25, "#d9ef8b"],
        [0.5, "#fee08b"],
        [0.75, "#fc8d59"],
        [1.0, "#d73027"]    # Dark green (highest coverage)
    ]
    
    fig = go.Figure()
    
    fig.add_trace(go.Choropleth(
        locations=plot_data["Country_standard"],
        locationmode="country names",
        z=plot_data["coverage_display"],
        zmin=0,
        zmax=1,
        colorscale=coverage_colorscale,
        colorbar=dict(
            title="Coverage Ratio",
            len=0.6,
            thickness=20,
            tickvals=tick_positions,
            ticktext=tick_labels
        ),
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
            plot_data["coverage_bin"].apply(lambda x: f"Bin {int(x)+1}/{int(max_bin)+1}").values
        ], axis=1)
    ))
    
    fig.update_geos(
        projection_type="natural earth",
        showcountries=True,
        countrycolor="rgba(0,0,0,0.2)",
        showsubunits=True,
        fitbounds="locations",
        domain=dict(x=[0, 1], y=[0.2, 1]),
    )
    fig.update_layout(
        height=700, 
        margin=dict(t=30, b=30, l=10, r=10),
        title=dict(text=f"{virus_type} Sequencing Coverage Map (Relative Ranking)", font=dict(size=16))
    )
    
    return fig
    
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
        line_shape="spline",
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
    
    fig = px.scatter(timeline_data, 
                     x="Year", 
                     y="prevalence_pct",
                     color="mutation",
                     size="count",
                     hover_data={"count": True, "prevalence_pct": ":.2f"})
    
    # Add lines connecting the points
    for mutation in top_muts:
        mutation_data = timeline_data[timeline_data["mutation"] == mutation]
        fig.add_trace(go.Scatter(
            x=mutation_data["Year"],
            y=mutation_data["prevalence_pct"],
            mode='lines',
            line=dict(width=1, color='lightgray'),
            showlegend=False,
            hoverinfo='skip'
        ))
    
    fig.update_layout(
        yaxis_title="Prevalence (%)",
        yaxis=dict(range=[0, None], fixedrange=True, rangemode="nonnegative"),
        xaxis_title="Year",
        xaxis=dict(fixedrange=True),
        height=400,
        hovermode="closest"
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
            "font": {"size": 16}
        }],
        height=300
    )
    return fig

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

# Redundant sidebar helpers removed (now global in Full_Hepatitis_page.py)


def dashboard_page_header():
    """Top heading/actions area to the right of the fixed sidebar."""
    return html.Div(
        className="hep-page-header",
        children=[
            html.Div([
                html.H1(id="dashboard-page-title", className="hep-page-title"),
                html.Div(id="dashboard-page-subtitle", className="hep-page-subtitle"),
            ]),
            html.Div(
                className="hep-page-actions",
                children=[
                    dcc.Input(
                        id="global-search-input",
                        placeholder="Search by Sequence ID",
                        type="text",
                        className="hep-search-input",
                    ),
                    dbc.Button(
                        [html.I(className="fa fa-search me-2"), "Search"],
                        id="global-search-btn",
                        color="link",
                        className="hep-search-btn",
                    ),
                    dbc.DropdownMenu(
                        label="Download",
                        children=[
                            dbc.DropdownMenuItem("Data", id="btn-download-data"),
                            dbc.DropdownMenuItem("Reports", href="/about#report-section"),
                        ],
                        color="link",
                        className="hep-download-btn",
                    ),
                    dbc.Toast(
                        "Your download is starting…",
                        id="dl-toast",
                        header="Download",
                        is_open=False,
                        dismissable=True,
                        icon="success",
                        duration=3000,
                        className="position-fixed top-0 end-0 m-3",
                    ),
                ],
            ),
        ],
    )


def build_overview_summary_cards():
    """Horizontal summary cards shown above filters/plots."""
    def metric_card(icon_class, value_id, label):
        return dbc.Col(
            html.Div(
                className="hep-metric-card",
                children=[
                    html.Div(html.I(className=icon_class), className="hep-metric-icon"),
                    html.Div([
                        html.Div(id=value_id, className="hep-metric-value"),
                        html.Div(label, className="hep-metric-label"),
                    ], className="hep-metric-text"),
                ],
            ),
            xs=12,
            md=6,
            xl=3,
        )

    return dbc.Row(
        id="overview-summary-cards-row",
        children=[
            metric_card("fa-solid fa-dna", "indicator-total", "Sequences"),
            metric_card("fa-solid fa-globe", "indicator-countries", "Countries"),
            metric_card("fa-solid fa-code-branch", "indicator-genotypes", "Genotypes"),
            metric_card("fa-solid fa-triangle-exclamation", "indicator-mutations", "Mutations"),
        ],
        className="hep-summary-row g-4",
    )


def build_epi_summary_cards():
    """Horizontal summary cards shown above epidemiology filters/plots."""
    def metric_card(icon_class, value_id, label):
        return dbc.Col(
            html.Div(
                className="hep-metric-card",
                children=[
                    html.Div(html.I(className=icon_class), className="hep-metric-icon"),
                    html.Div([
                        html.Div(id=value_id, className="hep-metric-value"),
                        html.Div(label, className="hep-metric-label"),
                    ], className="hep-metric-text"),
                ],
            ),
            xs=12,
            md=6,
            xl=3,
        )

    return dbc.Row(
        id="epi-summary-cards-row",
        children=[
            metric_card("fa-solid fa-user-shield", "epi-card-livingwith", "Living with infection"),
            metric_card("fa-solid fa-virus-covid", "epi-card-newinfections", "New infections"),
            metric_card("fa-solid fa-skull-crossbones", "epi-card-deaths", "Deaths"),
            metric_card("fa-solid fa-chart-line", "epi-card-diag-treat", "Diagnosis / Treatment %"),
        ],
        className="hep-summary-row g-4",
        style={"display": "none"},
    )



def make_priority_table(priority_data):
    """Render priority ranking as a Dash DataTable instead of a Plotly go.Table."""
    if priority_data is None or priority_data.empty:
        return html.Div(
            "No priority data available for current filters",
            className="hep-empty-table",
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

#app = dash.Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)

# === APP LAYOUT ==============================================================
def create_dashboard_layout():
    stores = html.Div(
        [
            # selected-virus store is now global in Full_Hepatitis_page.py
            #dcc.Store(id="dashboard-active-page-store", data=True),
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

    return html.Div(
        id="hep-dashboard-shell-container",
        className="hep-dashboard-shell",
        children=[
            stores,
            html.Main(
                style={"width": "100%", "flex": "1"},
                children=[
                    html.Div(
                        className="hep-main-sticky-header",
                        children=[
                            dashboard_page_header(),
                        ]
                    ),
                    build_overview_summary_cards(),
                    build_epi_summary_cards(),

        # === COMMON FILTERS (ALWAYS VISIBLE) ===
        dbc.Row([
            dbc.Col([
                html.H6("Filters", className="hep-section-title"),
                dbc.Card([
                    dbc.CardBody([
                        dbc.Row([
                            # Year range slider (35-45% of filter row, i.e. width 5 of 12)
                            dbc.Col([
                                html.Div([
                                    html.Label("YEAR RANGE", className="fw-bold text-uppercase small mb-0", 
                                               style={"letterSpacing": "0.05em", "color": "rgba(255,255,255,0.6)"}),
                                    html.Span(id="year-range-label", className="fw-bold float-end text-primary", 
                                              style={"fontSize": "0.95rem"})
                                ], className="mb-2 clearfix"),
                                html.Div([
                                    dcc.RangeSlider(
                                        id="year-range-slider",
                                        min=1963,
                                        max=2024,
                                        value=[1963, 2024],
                                        step=1,
                                        marks={y: str(y) for y in [1963, 1970, 1980, 1990, 2000, 2010, 2020, 2024]},
                                        tooltip={"always_visible": True, "placement": "top"},
                                        className="hep-range-slider"
                                    )
                                ], style={"paddingLeft": "5px", "paddingRight": "5px"})
                            ], xs=12, lg=5, className="mb-3 mb-lg-0"),
                            
                            # Region dropdown
                            dbc.Col([
                                html.Label("REGION", className="fw-bold text-uppercase small mb-2", 
                                           style={"letterSpacing": "0.05em", "color": "rgba(255,255,255,0.6)"}),
                                dcc.Dropdown(
                                    id="continent-dropdown",
                                    multi=True,
                                    placeholder="All regions",
                                    className="hep-dropdown"
                                )
                            ], xs=12, lg=2, className="mb-3 mb-lg-0"),
                            
                            # Country dropdown
                            dbc.Col([
                                html.Label("COUNTRY", className="fw-bold text-uppercase small mb-2", 
                                           style={"letterSpacing": "0.05em", "color": "rgba(255,255,255,0.6)"}),
                                dcc.Dropdown(
                                    id="country-dropdown",
                                    multi=True,
                                    placeholder="All countries",
                                    className="hep-dropdown"
                                )
                            ], xs=12, lg=3, className="mb-3 mb-lg-0"),
                            
                            # Genotype dropdown
                            dbc.Col([
                                html.Label("GENOTYPE", className="fw-bold text-uppercase small mb-2", 
                                           style={"letterSpacing": "0.05em", "color": "rgba(255,255,255,0.6)"}),
                                dcc.Dropdown(
                                    id="genotype-dropdown",
                                    multi=True,
                                    placeholder="All genotypes",
                                    className="hep-dropdown"
                                )
                            ], xs=12, lg=2)
                        ], className="align-items-end")
                    ])
                ], className="mb-4 shadow-sm border-0")
            ], width=12)
        ], id="common-filters", style={"position": "relative", "zIndex": 100}),


        # === TAB 1: OVERVIEW CONTENT (DEFAULT) ===
        html.Div(id="overview-content", children=[
            # Display Toggles
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Display Mode:", className="fw-bold me-2"),
                                    dcc.RadioItems(
                                        id="display-mode",
                                        options=[
                                            {"label": "Raw Count", "value": "raw"},
                                            {"label": "Per Million", "value": "PerMillion"}
                                        ],
                                        value="raw",
                                        className="hep-radio-group"
                                    )
                                ], width=6),
                                
                                dbc.Col([
                                    html.Label("Map Mode:", className="fw-bold me-2"),
                                    dcc.RadioItems(
                                        id="map-mode",
                                        options=[
                                            {"label": "Sequences", "value": "sequences"},
                                            {"label": "Coverage", "value": "coverage"}, 
                                            {"label": "Epidemiology", "value": "epidemiology"},
                                        ],
                                        value="sequences",
                                        className="hep-radio-group"
                                    ),
                                ], width="auto"),
                            ], className="g-3 align-items-end mb-2"),
                        ])
                    ], className="mb-4 shadow-sm")
                ], width=12)
            ]),

            # ROW 1: Map Section
            dbc.Row([
                # Map Section
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H4(id="map-title-main", className="mb-1"),
                                html.Small(id="map-title-sub", className="text-muted")
                            ], className="mb-3"),
                            
                            dbc.Row([
                                dbc.Col([
                                    html.Label("IHME Metric:", className="fw-bold me-2"),
                                    dcc.Dropdown(
                                        id="ihme-metric-type",
                                        options=[
                                            {"label": "Deaths", "value": "Deaths|Number"},
                                            {"label": "Incidence", "value": "Incidence|Number"},
                                            {"label": "Prevalence", "value": "Prevalence|Number"},
                                        ],
                                        value="Prevalence|Number",
                                        clearable=False,
                                        style={"width": "200px"}
                                    )
                                ], width="auto", id="epidemiology-controls"),
                            ], className="mb-3 g-3"),
                            
                            dcc.Loading(
                                dcc.Graph(
                                    id="genotype-map", 
                                    config={'displayModeBar': True, 'displaylogo': False},
                                    style={"height": "100%"}
                                ),
                                type="circle"
                            )
                        ])
                    ], className="h-100")
                ], width=12)
            ], className="mb-4"),

            # ROW 1: Burden vs. Sequencing Correlation and Sequencing Priority
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Burden vs. Sequencing Correlation", className="mb-3"),
                            dcc.Loading(
                                dcc.Graph(
                                    id="burden-coverage-scatter",
                                    className="hep-graph",
                                    config={"displayModeBar": True, "displaylogo": False},
                                ),
                                type="circle"
                            )
                        ])
                    ], className="hep-card h-100 shadow-sm")
                ], width=7),

                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("Sequencing Priority Ranking", className="mb-0"),
                                dbc.Button(
                                    "Download CSV",
                                    id="priority-download-btn",
                                    color="link",
                                    className="hep-small-action-btn",
                                ),
                            ], className="hep-card-header-row"),
                            dcc.Loading(
                                html.Div(id="priority-ranking-table"),
                                type="circle"
                            )
                        ])
                    ], className="hep-card h-100 shadow-sm")
                ], width=5),
            ], className="g-4 mb-4"),
            
            # ROW 2: Time Series
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5(id="line-title-main", className="mb-3"),
                            dcc.Loading(
                                dcc.Graph(id="line-chart", style={"height": "100%"}),
                                type="circle"
                            )
                        ])
                    ], className="h-100 shadow-sm")
                ], width=12),           
            ], className="mb-4"),        
            
            # ROW 3: genotype and Country Distribution
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5(id="bar-title-main", className="mb-3"),
                            dcc.Loading(
                                dcc.Graph(id="genotype-bar-chart", style={"height": "100%"}),
                                type="circle"
                            )
                        ])
                    ], className="h-100 shadow-sm")
                ], width=6),
                
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Sequences by Country and genotype", className="mb-3"),
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Top N Countries:", className="fw-bold me-2"),
                                    dcc.Dropdown(
                                        id="top-countries-count",
                                        options=[{"label": str(i), "value": i} for i in [10, 15, 20, 25]],
                                        value=10,
                                        clearable=False,
                                        style={"width": "150px"}
                                    )
                                ], width="auto")
                            ], className="mb-3"),
                            dcc.Loading(
                                dcc.Graph(id="country-barchart", style={"height": "100%"}),
                                type="circle"
                            )
                        ])
                    ], className="h-100 shadow-sm")
                ], width=6)
            ], className="mb-4"),

            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Mutation Timeline", className="mb-3"),
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Top N Mutations:", className="fw-bold me-2"),
                                    dcc.Dropdown(
                                        id="top-mutations-count",
                                        options=[{"label": str(i), "value": i} for i in [5, 10, 15, 20]],
                                        value=10,
                                        clearable=False,
                                        style={"width": "150px"}
                                    )
                                ], width="auto")
                            ], className="mb-2"),
                            dcc.Loading(
                                dcc.Graph(id="mutation-timeline"),
                                type="circle"
                            )
                        ])
                    ], className="h-100 shadow-sm")
                ], width=12)
            ], className="mb-4"),

            
            # ROW 4: Epidemiology Summary
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Epidemiology Summary", className="mb-3"),
                            dbc.Row([
                                # Card 1: Total Prevalence
                                dbc.Col([
                                    dbc.Card([
                                        dbc.CardBody([
                                            html.H6("Total Prevalence", className="card-subtitle"),
                                            html.H4(id="epi-prevalence-total", className="card-title"),
                                            html.Small(id="epi-prevalence-trend", className="text-muted"),
                                            html.Div([
                                                html.I(className="bi bi-people-fill me-2"),
                                                html.Span("Global cases", className="small")
                                            ], className="mt-2")
                                        ])
                                    ], className="text-center h-100 border-start border-5 border-primary")
                                ], width=3),
                                
                                # Card 2: Incidence
                                dbc.Col([
                                    dbc.Card([
                                        dbc.CardBody([
                                            html.H6("Annual Incidence", className="card-subtitle"),
                                            html.H4(id="epi-incidence-total", className="card-title"),
                                            html.Small(id="epi-incidence-trend", className="text-muted"),
                                            html.Div([
                                                html.I(className="bi bi-graph-up-arrow me-2"),
                                                html.Span("New cases/year", className="small")
                                            ], className="mt-2")
                                        ])
                                    ], className="text-center h-100 border-start border-5 border-warning")
                                ], width=3),
                                
                                # Card 3: Deaths
                                dbc.Col([
                                    dbc.Card([
                                        dbc.CardBody([
                                            html.H6("Annual Deaths", className="card-subtitle"),
                                            html.H4(id="epi-deaths-total", className="card-title"),
                                            html.Small(id="epi-deaths-trend", className="text-muted"),
                                            html.Div([
                                                html.I(className="bi bi-heartbreak-fill me-2"),
                                                html.Span("Mortality", className="small")
                                            ], className="mt-2")
                                        ])
                                    ], className="text-center h-100 border-start border-5 border-danger")
                                ], width=3),
                                
                                # Card 4: Coverage Gap
                                dbc.Col([
                                    dbc.Card([
                                        dbc.CardBody([
                                            html.H6("Sequencing Coverage", className="card-subtitle"),
                                            html.H4(id="epi-coverage-percent", className="card-title"),
                                            html.Small(id="epi-coverage-status", className="text-muted"),
                                            html.Div([
                                                html.I(className="bi bi-clipboard-data-fill me-2"),
                                                html.Span("vs. estimated infections", className="small")
                                            ], className="mt-2")
                                        ])
                                    ], className="text-center h-100 border-start border-5 border-success")
                                ], width=3),
                            ], className="g-3 mb-3"),
                            
                            # Second row with more metrics
                            dbc.Row([
                                # Card 5: Sex Ratio
                                dbc.Col([
                                    dbc.Card([
                                        dbc.CardBody([
                                            html.H6("Male:Female Ratio", className="card-subtitle"),
                                            html.H4(id="epi-sex-ratio", className="card-title"),
                                            html.Small("Latest year", className="text-muted")
                                        ])
                                    ], className="text-center h-100")
                                ], width=2),
                                
                                # Card 6: Top Age Group
                                dbc.Col([
                                    dbc.Card([
                                        dbc.CardBody([
                                            html.H6("Most Affected Age", className="card-subtitle"),
                                            html.H4(id="epi-top-age-group", className="card-title"),
                                            html.Small(id="epi-age-percentage", className="text-muted")
                                        ])
                                    ], className="text-center h-100")
                                ], width=2),
                                
                                # Card 7: Top Region
                                dbc.Col([
                                    dbc.Card([
                                        dbc.CardBody([
                                            html.H6("Highest Burden Region", className="card-subtitle"),
                                            html.H4(id="epi-top-region", className="card-title"),
                                            html.Small(id="epi-region-percentage", className="text-muted")
                                        ])
                                    ], className="text-center h-100")
                                ], width=2),
                                
                                # Card 8: WHO 2030 Progress
                                dbc.Col([
                                    dbc.Card([
                                        dbc.CardBody([
                                            html.H6("WHO 2030 Target", className="card-subtitle"),
                                            html.H4(id="epi-2030-progress", className="card-title"),
                                            html.Small("Reduction needed", className="text-muted")
                                        ])
                                    ], className="text-center h-100")
                                ], width=3),
                                
                                # Card 9: Quick link to epidemiology tab
                                dbc.Col([
                                    dbc.Card([
                                        dbc.CardBody([
                                            html.Div([
                                                html.I(className="bi bi-clipboard2-pulse-fill", 
                                                       style={"fontSize": "2rem", "color": "#0d6efd"}),
                                            ], className="mb-2"),
                                            dbc.Button(
                                                "View Full Analysis",
                                                id="btn-go-to-epidemiology",
                                                color="primary",
                                                size="sm",
                                                className="w-100"
                                            )
                                        ], className="text-center")
                                    ], className="h-100")
                                ], width=3),
                            ], className="g-3"),
                        ])
                    ], className="shadow-sm mb-4")
                ], width=12)
            ], id="epidemiology-summary-row", className="mb-4"),
            
            # ROW 5: Mutation Summary
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5(id="mutation-section-title", className="mb-3"),  # Dynamic title
                            html.Div(id="mutation-section-content")  # Dynamic content
                        ])
                    ], className="shadow-sm")
                ], width=12)
            ], className="mb-4", id="mutation-summary-row"),
            
            # ROW 6: Quick Actions (New section)
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Quick Analysis Tools", className="mb-3"),
                            dbc.Row([
                                dbc.Col([
                                    dbc.Button(
                                        [html.I(className="bi bi-clipboard2-pulse me-2"), "View Burden Forecast"],
                                        id="btn-quick-forecast",
                                        color="outline-primary",
                                        className="w-100 mb-2 py-3"
                                    )
                                ], width=4),
                                dbc.Col([
                                    dbc.Button(
                                        [html.I(className="bi bi-sort-numeric-down me-2"), "Sequencing Priority"],
                                        id="btn-quick-priority",
                                        color="outline-success",
                                        className="w-100 mb-2 py-3"
                                    )
                                ], width=4),
                                dbc.Col([
                                    dbc.Button(
                                        [html.I(className="bi bi-clock-history me-2"), "Mutation Timeline"],
                                        id="btn-quick-timeline",
                                        color="outline-warning",
                                        className="w-100 mb-2 py-3"
                                    )
                                ], width=4),
                            ], className="g-3"),
                            html.Small("Click any button to jump to detailed analysis", className="text-muted mt-2 d-block text-center")
                        ])
                    ], className="shadow-sm")
                ], width=12)
            ], className="mb-4"),
        ]),

        # === TAB 2: MUTATIONS CONTENT (HIDDEN BY DEFAULT) ===
        html.Div(id="mutations-content", style={"display": "none"}, children=[
            # Mutation-specific filters
            dbc.Row(id="mutation-filters", children=[
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Mutation Type:", className="fw-bold me-2"),
                                    dcc.Dropdown(
                                        id="mutation-type-filter",
                                        options=[
                                            {"label": "All Types", "value": "all"},
                                            {"label": "Antiviral Resistance", "value": "antiviral_resistance"},
                                            {"label": "Vaccine Escape (HBV)", "value": "vaccine_escape"},
                                            {"label": "Substitutions of Interest", "value": "substitution_of_interest"},
                                            {"label": "No Resistance", "value": "no_resistance"}
                                        ],
                                        value="all",
                                        clearable=False,
                                        style={"width": "250px"}
                                    )
                                ], width=3),
                                
                                dbc.Col([
                                    html.Label("Gene/Drug:", className="fw-bold me-2"),
                                    dcc.Dropdown(
                                        id="mutation-category-filter",
                                        options=[],
                                        placeholder="Select category...",
                                        style={"width": "250px"}
                                    )
                                ], width=3),
                                
                                dbc.Col([
                                    html.Label("Show Top:", className="fw-bold me-2"),
                                    dcc.Dropdown(
                                        id="mutation-top-n",
                                        options=[{"label": str(i), "value": i} for i in [10, 20, 30, 50]],
                                        value=20,
                                        clearable=False,
                                        style={"width": "150px"}
                                    )
                                ], width=3),
                                
                                dbc.Col([
                                    html.Label("Group By:", className="fw-bold me-2"),
                                    dcc.RadioItems(
                                        id="mutation-group-by",
                                        options=[
                                            {"label": " Type", "value": "type"},
                                            {"label": " Drug", "value": "drug"},
                                            {"label": " Gene", "value": "gene"}
                                        ],
                                        value="type",
                                        inline=True,
                                        style={"marginTop": "8px"}
                                    )
                                ], width=3),
                            ], className="g-3 mb-3"),
                        ])
                    ], className="shadow-sm")
                ], width=12)
            ], className="mb-4"),
            
            # Mutation Visualizations
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Mutation Frequency", className="mb-3"),
                            dcc.Loading(
                                dcc.Graph(id="mutation-frequency-chart"),
                                type="circle"
                            )
                        ])
                    ], className="h-100 shadow-sm")
                ], width=8),
                
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Mutation Types", className="mb-3"),
                            dcc.Loading(
                                dcc.Graph(id="mutation-distribution-chart"),
                                type="circle"
                            )
                        ])
                    ], className="h-100 shadow-sm")
                ], width=4),
            ], className="mb-4"),
                        
            # Mutation Details Table
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("Mutation Details", className="d-inline mb-0"),
                                dbc.Button(
                                    "Back to Overview",
                                    id="btn-back-to-overview-from-mutations",
                                    color="primary",
                                    size="sm",
                                    className="float-end"
                                ),
                                dbc.Button(
                                    "Download CSV",
                                    id="btn-download-mutations",
                                    color="secondary",
                                    size="sm",
                                    className="float-end me-2"
                                )
                            ], className="mb-3"),
                            html.Div(id="mutation-details-table")
                        ])
                    ], className="shadow-sm")
                ], width=12)
            ], className="mb-4"),
        ]),
        # === TAB 2: MUTATIONS CONTENT (HIDDEN BY DEFAULT) ===
        html.Div(id="epidemiology-content", style={"display": "none"}, children=[
            # GHO Controls Row
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Select Metric:", className="fw-bold me-2"),
                                    dcc.Dropdown(
                                        id="epi-metric-dropdown",
                                        options=[
                                            {"label": "People living with infection", "value": "livingwith_num"},
                                            {"label": "New infections", "value": "new_infections_num"},
                                            {"label": "Deaths", "value": "deaths_num"},
                                            {"label": "Prevalence %", "value": "prevalence_pct"},
                                            {"label": "Diagnosis rate %", "value": "diagnosis_rate_pct"},
                                            {"label": "Treatment rate %", "value": "treatment_rate_diagnosed_pct"},
                                            {"label": "HBV vaccine coverage % (HepB3)", "value": "vaccine_hepb3_coverage_pct"}
                                        ],
                                        value="prevalence_pct",
                                        clearable=False,
                                        className="hep-dropdown"
                                    )
                                ], xs=12, md=6, lg=4),
                            ], className="g-3 align-items-center")
                        ])
                    ], className="mb-4 shadow-sm border-0")
                ], width=12)
            ], id="epi-controls-row"),

            # Row 1: Map
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Global Epidemiological Burden", className="mb-3"),
                            dcc.Loading(
                                dcc.Graph(id="gho-burden-map"),
                                type="circle"
                            )
                        ])
                    ], className="shadow-sm border-0 h-100")
                ], width=12),
            ], className="mb-4"),

            # Row 2: Cascade and Forecast Chart Side-by-Side
            dbc.Row([
                # Diagnosis & Treatment Cascade
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Diagnosis & Treatment Cascade", className="mb-3"),
                            dcc.Loading(
                                dcc.Graph(id="gho-cascade-chart"),
                                type="circle"
                            )
                        ])
                    ], className="shadow-sm border-0 h-100")
                ], xs=12, lg=6, className="mb-4 mb-lg-0"),
                
                # Forecast Chart
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Burden Forecast with Projections", className="mb-3"),
                            dcc.Loading(
                                dcc.Graph(
                                    id="forecast-chart",
                                    className="hep-graph",
                                    config={"displayModeBar": True, "displaylogo": False},
                                ),
                                type="circle"
                            )
                        ])
                    ], className="shadow-sm border-0 h-100")
                ], xs=12, lg=6)
            ], className="mb-4"),

            # Section: Program Response vs Burden
            dbc.Row([
                dbc.Col([
                    html.H4("Program Response vs Burden", className="text-white mt-4 mb-3 fw-bold")
                ], width=12)
            ]),

            # Row 3: HBV Indicators & Priority Table
            dbc.Row([
                # HBV Left Plot
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("HBV Vaccine Coverage vs Disease Burden", className="mb-3"),
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Y-Axis Burden Metric:", className="small fw-bold text-muted mb-1"),
                                    dcc.Dropdown(
                                        id="hbv-y-metric",
                                        options=[
                                            {"label": "Prevalence %", "value": "hbv_prevalence_pct"},
                                            {"label": "New Infections", "value": "hbv_new_infections_num"},
                                            {"label": "Living with HBV", "value": "hbv_livingwith_num"}
                                        ],
                                        value="hbv_prevalence_pct",
                                        clearable=False,
                                        className="hep-dropdown mb-3"
                                    )
                                ], width=12)
                            ]),
                            dcc.Loading(
                                dcc.Graph(id="hbv-epi-trend-plot"),
                                type="circle"
                            )
                        ])
                    ], className="shadow-sm border-0 h-100")
                ], xs=12, lg=6, className="mb-4 mb-lg-0"),
                
                # HBV Right Table
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("HBV Country Priority Ranking", className="mb-3"),
                            dcc.Loading(
                                html.Div(id="hbv-epi-priority-table"),
                                type="circle"
                            )
                        ])
                    ], className="shadow-sm border-0 h-100")
                ], xs=12, lg=6)
            ], id="hbv-epi-trends-row", className="mb-4", style={"display": "none"}),

            # Row 4: HCV Indicators & Priority Table
            dbc.Row([
                # HCV Left Plot
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("HCV Diagnosis/Treatment Access vs Disease Burden", className="mb-3"),
                            dbc.Row([
                                dbc.Col([
                                    html.Label("X-Axis Access Metric:", className="small fw-bold text-muted mb-1"),
                                    dcc.Dropdown(
                                        id="hcv-x-metric",
                                        options=[
                                            {"label": "Diagnosis Rate %", "value": "hcv_diagnosis_rate_pct"},
                                            {"label": "Treatment Rate %", "value": "hcv_treatment_rate_diagnosed_pct"}
                                        ],
                                        value="hcv_diagnosis_rate_pct",
                                        clearable=False,
                                        className="hep-dropdown mb-3"
                                    )
                                ], xs=12, md=6),
                                dbc.Col([
                                    html.Label("Y-Axis Burden Metric:", className="small fw-bold text-muted mb-1"),
                                    dcc.Dropdown(
                                        id="hcv-y-metric",
                                        options=[
                                            {"label": "Prevalence %", "value": "hcv_prevalence_pct"},
                                            {"label": "New Infections", "value": "hcv_new_infections_num"},
                                            {"label": "Living with HCV", "value": "hcv_livingwith_num"}
                                        ],
                                        value="hcv_prevalence_pct",
                                        clearable=False,
                                        className="hep-dropdown mb-3"
                                    )
                                ], xs=12, md=6)
                            ]),
                            dcc.Loading(
                                dcc.Graph(id="hcv-epi-trend-plot"),
                                type="circle"
                            )
                        ])
                    ], className="shadow-sm border-0 h-100")
                ], xs=12, lg=6, className="mb-4 mb-lg-0"),
                
                # HCV Right Table
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("HCV Country Priority Ranking", className="mb-3"),
                            dcc.Loading(
                                html.Div(id="hcv-epi-priority-table"),
                                type="circle"
                            )
                        ])
                    ], className="shadow-sm border-0 h-100")
                ], xs=12, lg=6)
            ], id="hcv-epi-trends-row", className="mb-4", style={"display": "none"}),

            # Row 5: Detailed DataTable
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("WHO GHO Country Profile Data", className="d-inline mb-0"),
                                html.Div([
                                    dbc.Button(
                                        "Back to Overview",
                                        id="btn-back-to-overview-from-epi",
                                        color="primary",
                                        size="sm",
                                        className="me-2"
                                    ),
                                    dbc.Button(
                                        "Download Data",
                                        id="btn-download-epi",
                                        color="secondary",
                                        size="sm",
                                    ),
                                ], className="d-flex align-items-center"),
                                dcc.Download(id="download-epi-data")
                            ], className="mb-3 d-flex justify-content-between align-items-center flex-wrap"),
                            dcc.Loading(
                                html.Div(id="gho-data-table"),
                                type="circle"
                            )
                        ])
                    ], className="shadow-sm border-0")
                ], width=12)
            ], className="mb-4")
        ]),

        # === TAB 4: USER SEQUENCE SUBMISSION ===
        user_seq_tab_content(),

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

                ],
            ),
        ],
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

    light = df[["Country_standard", "WHO_Regions", "Year", "genotype"]].copy()
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
    
# Overview Page Forecast Chart Callback
@callback(
    Output("forecast-chart", "figure"),
    Input("selected-virus", "data"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
)
def update_overview_forecast_chart(virus, regions, countries):
    data = get_data_store()
    ihme_df = data.get("ihme_df", pd.DataFrame())
    if ihme_df.empty:
        return _empty_plot("No historical GBD data available")
        
    fig = create_forecast_chart(
        ihme_df=ihme_df,
        selected_virus=virus,
        sex="Both",
        selected_regions=regions,
        selected_countries=countries
    )
    
    # Apply dark theme styling
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="rgba(255,255,255,0.7)"),
    )
    return fig

# — Indicator values —
@callback(
    Output("indicator-total", "children"),
    Output("indicator-countries", "children"),
    Output("indicator-genotypes", "children"),
    Output("indicator-mutations", "children"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
)
def update_indicators(filtered_json, selected_virus):
    df = _df_from_json(filtered_json)
    if df is None or df.empty:
        return "0", "0", "0", "0"

    total_genomes = len(df)
    unique_countries = df.get("Country_standard", pd.Series(dtype="object")).nunique()

    g = df.get("genotype", pd.Series(dtype="object")).astype("string").str.strip()
    recomb_mask = g.str.contains(r"recomb", case=False, na=False)
    base_genotype_count = g[~recomb_mask].dropna().nunique()

    data = get_data_store()
    selected_virus = selected_virus or "HBV"
    if selected_virus == "HBV":
        mut_df = data.get("hbv_mut", pd.DataFrame())
    elif selected_virus == "HCV":
        mut_df = data.get("hcv_mut", pd.DataFrame())
    elif selected_virus == "HEV":
        mut_df = data.get("hev_mut", pd.DataFrame())
    else:
        mut_df = pd.DataFrame()

    mutation_count = mut_df["mutation"].nunique() if not mut_df.empty and "mutation" in mut_df.columns else 0

    return (
        f"{total_genomes:,}",
        f"{unique_countries:,}",
        f"{base_genotype_count:,}",
        f"{mutation_count:,}",
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
        subtitle = f"{sequence_count:,} sequences · {mutation_count:,} mutations · Active surveillance"
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
    Output("tab-overview", "color"),
    Output("tab-mutations", "color"),
    Output("tab-epidemiology", "color"),
    Output("tab-user-seq", "color"),
    Output("tab-overview", "className"),
    Output("tab-mutations", "className"),
    Output("tab-epidemiology", "className"),
    Output("tab-user-seq", "className"),
    Output("tab-mutations", "disabled"),
    Output("active-tab-store", "data"),
    
    # Inputs
    Input("url", "pathname"),
    Input("tab-overview", "n_clicks"),
    Input("tab-mutations", "n_clicks"),
    Input("tab-epidemiology", "n_clicks"),
    Input("tab-user-seq", "n_clicks"),
    Input("selected-virus", "data"),
    
    # States to track initial load
    State("tab-overview", "n_clicks"),
    State("tab-mutations", "n_clicks"),
    State("tab-epidemiology", "n_clicks"),
)
def update_tab_highlights(pathname, overview_clicks, mutations_clicks, epidemiology_clicks,
                          user_seq_clicks, selected_virus,
                          overview_state, mutations_state, epidemiology_state):
    ctx = callback_context
    mutations_disabled = selected_virus == "HEV" if selected_virus else False

    if pathname not in ["/", "/dashboard"]:
        # Deactivate all dashboard tabs on non-dashboard pages
        print(f"DEBUG: update_tab_highlights called. pathname={pathname}. Return overview (inactive).")
        return ("link", "link", "link", "link", 
                "hep-nav-item", "hep-nav-item" + (" hep-nav-disabled" if mutations_disabled else ""), "hep-nav-item", "hep-nav-item", mutations_disabled, "overview")

    # Determine active tab using modern triggered_id
    trigger_id = ctx.triggered_id
    active_tab = "overview"
    if trigger_id == "tab-overview":
        active_tab = "overview"
    elif trigger_id == "tab-mutations" and not mutations_disabled:
        active_tab = "mutations"
    elif trigger_id == "tab-epidemiology":
        active_tab = "epidemiology"
    elif trigger_id == "tab-user-seq":
        active_tab = "user-seq"
    elif trigger_id == "selected-virus":
        active_tab = "overview"
    else:
        active_tab = "overview"

    colors = []
    classes = []
    for t in ["overview", "mutations", "epidemiology", "user-seq"]:
        colors.append("link")
        if t == "mutations" and mutations_disabled:
            classes.append("hep-nav-item hep-nav-disabled")
        elif t == active_tab:
            classes.append("hep-nav-item hep-nav-active")
        else:
            classes.append("hep-nav-item")

    print(f"DEBUG: update_tab_highlights finalized active_tab={active_tab}")
    return (*colors, *classes, mutations_disabled, active_tab)


# Client-side tab switching callback to avoid nonexistent DOM errors on page loading/switching
dash.clientside_callback(
    """
    function(activeTab, pathname, overviewId) {
        // If we are not on the dashboard page, do nothing to prevent nonexistent object errors
        if (pathname !== "/" && pathname !== "/dashboard") {
            return [
                window.dash_clientside.no_update, 
                window.dash_clientside.no_update, 
                window.dash_clientside.no_update, 
                window.dash_clientside.no_update, 
                window.dash_clientside.no_update,
                window.dash_clientside.no_update,
                window.dash_clientside.no_update
            ];
        }
        
        // Verify that the elements actually exist in the DOM before trying to update them
        const overview = document.getElementById("overview-content");
        if (!overview) {
            return [
                window.dash_clientside.no_update, 
                window.dash_clientside.no_update, 
                window.dash_clientside.no_update, 
                window.dash_clientside.no_update, 
                window.dash_clientside.no_update,
                window.dash_clientside.no_update,
                window.dash_clientside.no_update
            ];
        }
        
        const showBlock = {"display": "block"};
        const showFlex = {"display": "flex"};
        const hide = {"display": "none"};
        
        const active = activeTab || "overview";
        
        return [
            active === "overview" ? showBlock : hide,
            active === "mutations" ? showBlock : hide,
            active === "epidemiology" ? showBlock : hide,
            active === "user-seq" ? showBlock : hide,
            active === "user-seq" ? hide : showBlock,
            (active !== "epidemiology" && active !== "user-seq") ? showFlex : hide,
            active === "epidemiology" ? showFlex : hide
        ];
    }
    """,
    Output("overview-content", "style"),
    Output("mutations-content", "style"),
    Output("epidemiology-content", "style"),
    Output("user-seq-content", "style"),
    Output("common-filters", "style"),
    Output("overview-summary-cards-row", "style"),
    Output("epi-summary-cards-row", "style"),
    Input("active-tab-store", "data"),
    Input("url", "pathname"),
    Input("overview-content", "id"),
    prevent_initial_call=True
)


# Client-side callback to route quick/back buttons to global tab clicks
dash.clientside_callback(
    """
    function(back_mut, back_epi, go_epi, quick_forecast, quick_priority, quick_timeline) {
        const ctx = window.dash_clientside.callback_context;
        if (!ctx.triggered || ctx.triggered.length === 0) {
            return window.dash_clientside.no_update;
        }
        const trigger = ctx.triggered[0];
        if (!trigger || trigger.value === null || trigger.value === undefined || trigger.value === 0) {
            return window.dash_clientside.no_update;
        }
        const trigger_id = trigger.prop_id.split('.')[0];
        
        if (trigger_id === 'btn-back-to-overview-from-mutations' || trigger_id === 'btn-back-to-overview-from-epi') {
            const btn = document.getElementById('tab-overview');
            if (btn) btn.click();
        } else if (trigger_id === 'btn-go-to-epidemiology' || trigger_id === 'btn-quick-forecast' || trigger_id === 'btn-quick-priority') {
            const btn = document.getElementById('tab-epidemiology');
            if (btn) btn.click();
        } else if (trigger_id === 'btn-quick-timeline') {
            const btn = document.getElementById('tab-mutations');
            if (btn) btn.click();
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output("url", "id"),
    Input("btn-back-to-overview-from-mutations", "n_clicks"),
    Input("btn-back-to-overview-from-epi", "n_clicks"),
    Input("btn-go-to-epidemiology", "n_clicks"),
    Input("btn-quick-forecast", "n_clicks"),
    Input("btn-quick-priority", "n_clicks"),
    Input("btn-quick-timeline", "n_clicks"),
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
    Output("mutation-distribution-chart", "figure"),
    Input("filtered-store", "data"),
    Input("selected-virus", "data"),
    Input("mutation-type-filter", "value"),
)
def update_mutation_distribution_chart(filtered_json, virus, mutation_type):
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
    
    # Apply type filter
    if mutation_type != "all" and "type" in filtered_mutations.columns:
        filtered_mutations = filtered_mutations[filtered_mutations["type"] == mutation_type]
    
    if filtered_mutations.empty:
        return _empty_plot("No mutations in filtered data")
    
    # Group by type
    if "type" in filtered_mutations.columns:
        type_counts = filtered_mutations.groupby("type").size().reset_index(name="count")
        
        # Color mapping
        color_map = {
            "antiviral_resistance": "#e41a1c",
            "vaccine_escape": "#377eb8", 
            "substitution_of_interest": "#4daf4a",
            "no_resistance": "#999999"
        }
        
        # Create pie chart
        fig = px.pie(
            type_counts,
            names="type",
            values="count",
            hole=0.4,
            color="type",
            color_discrete_map=color_map,
            title=f"{selected_virus} Mutation Types Distribution"
        )
    else:
        # Fallback: group by mutation
        mutation_counts = filtered_mutations.groupby("mutation").size().reset_index(name="count")
        mutation_counts = mutation_counts.sort_values("count", ascending=False).head(10)
        
        fig = px.pie(
            mutation_counts,
            names="mutation",
            values="count",
            hole=0.4,
            title=f"{selected_virus} Top 10 Mutations"
        )
    
    fig.update_traces(
        textposition="inside",
        textinfo="percent+label",
        hovertemplate="<b>%{label}</b><br>Count: %{value}<br>Percentage: %{percent}<extra></extra>"
    )
    
    fig.update_layout(
        height=400,
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=-0.2,
            xanchor="center",
            x=0.5
        )
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
    
    # Create DataTable
    table = dash.dash_table.DataTable(
        data=filtered_mutations[columns_to_show].to_dict('records'),
        columns=[{"name": col.capitalize(), "id": col} for col in columns_to_show],
        page_size=10,
        style_table={'overflowX': 'auto'},
        style_cell={
            'textAlign': 'left',
            'padding': '10px',
            'overflow': 'hidden',
            'textOverflow': 'ellipsis',
            'maxWidth': 0,
        },
        style_header={
            'backgroundColor': 'rgb(230, 230, 230)',
            'fontWeight': 'bold'
        },
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
    Output("priority-ranking-table", "children"),
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
    
    return make_priority_table(priority_df), _df_to_json(priority_df)
    
# === EPIDEMIOLOGY CALLBACKS ===


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
    Input("selected-virus", "data"),
    State("epi-metric-dropdown", "value"),
)
def update_epi_metric_dropdown_options(virus, current_val):
    options = [
        {"label": "People living with infection", "value": "livingwith_num"},
        {"label": "New infections", "value": "new_infections_num"},
        {"label": "Deaths", "value": "deaths_num"},
        {"label": "Prevalence %", "value": "prevalence_pct"},
        {"label": "Diagnosis rate %", "value": "diagnosis_rate_pct"},
        {"label": "Treatment rate %", "value": "treatment_rate_diagnosed_pct"}
    ]
    
    if virus == "HBV":
        options.append({"label": "HBV vaccine coverage % (HepB3)", "value": "vaccine_hepb3_coverage_pct"})
        
    valid_vals = [opt["value"] for opt in options]
    val = current_val if current_val in valid_vals else "prevalence_pct"
    
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
    who_gho_df = data.get("who_gho_df", pd.DataFrame())
    
    if who_gho_df.empty:
        return _empty_plot("No WHO GHO data available")
        
    prefix = virus.lower() + "_"
    col_name = prefix + metric
    
    if col_name not in who_gho_df.columns:
        return _empty_plot(f"Metric '{metric}' not found for virus {virus}")
        
    if metric == "vaccine_hepb3_coverage_pct":
        non_null_years = who_gho_df[who_gho_df[col_name].notna()]["year"]
        map_year = int(non_null_years.max()) if not non_null_years.empty else 2022
    else:
        map_year = 2022
        
    df = who_gho_df[who_gho_df["year"] == map_year]
    
    if regions:
        df = df[df["WHO_Regions"].isin(regions)]
    if countries:
        df = df[df["Country_standard"].isin(countries)]
        
    if df.empty or df[col_name].isna().all():
        return _empty_plot(f"No data available for year {map_year} with current filters")
        
    metric_labels = {
        "livingwith_num": "People living with infection",
        "new_infections_num": "New infections",
        "deaths_num": "Deaths",
        "prevalence_pct": "Prevalence %",
        "diagnosis_rate_pct": "Diagnosis rate %",
        "treatment_rate_diagnosed_pct": "Treatment rate %",
        "vaccine_hepb3_coverage_pct": "HBV vaccine coverage % (HepB3)"
    }
    metric_label = metric_labels.get(metric, metric)
    
    virus_color = VIRUS_COLORS.get(virus, "#E84057")
    color_scale = sequential_scale(virus_color)
    
    fig = go.Figure(data=go.Choropleth(
        locations=df["country"],
        z=df[col_name],
        text=df["Country_standard"],
        locationmode="ISO-3",
        colorscale=color_scale,
        autocolorscale=False,
        reversescale=False,
        marker_line_color="rgba(255,255,255,0.15)",
        marker_line_width=0.5,
        colorbar_title=metric_label,
        colorbar_ticksuffix="%" if "pct" in metric or "rate" in metric else "",
        hovertemplate="<b>%{text}</b><br>" + metric_label + ": %{z:,.2f}<extra></extra>"
    ))
    
    fig.update_layout(
        title=dict(
            text=f"Global {metric_label} ({map_year})",
            font=dict(size=16, color="#ECEFF2"),
            x=0.05, y=0.95
        ),
        geo=dict(
            showframe=False,
            showcoastlines=True,
            projection_type='natural earth',
            landcolor='#101E2B',
            oceancolor='#07111A',
            showocean=True,
            showland=True,
            coastlinecolor="rgba(255,255,255,0.08)",
            countrycolor="rgba(255,255,255,0.06)",
        ),
        margin=dict(l=0, r=0, t=50, b=0),
        height=450,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
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
        "coverage_gap",
        "priority_score"
    ]
    
    available_cols = [c for c in columns_needed if c in display_df.columns]
    display_df = display_df[available_cols].copy()
    
    rename_map = {
        "rank": "Rank",
        "Country_standard": "Country",
        "burden": "Living with infection",
        diag_col: "Diagnosed %",
        treat_col: "Treatment %",
        "coverage_gap": "Estimated gap",
        "priority_score": "Priority score"
    }
    display_df = display_df.rename(columns=rename_map)
    
    # Format columns
    if "Priority score" in display_df.columns:
        display_df["Priority score"] = pd.to_numeric(display_df["Priority score"], errors="coerce").round(3)
        
    for col in ["Living with infection", "Estimated gap"]:
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
    
    return dash_table.DataTable(
        data=display_df.to_dict("records"),
        columns=[{"name": col, "id": col} for col in display_df.columns],
        page_size=6,
        sort_action="native",
        style_table={
            "overflowX": "auto",
            "backgroundColor": "transparent",
        },
        style_header=TABLE_HEADER_STYLE,
        style_cell=TABLE_CELL_STYLE,
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
        title=dict(
            text="HBV Vaccine Coverage vs Disease Burden (2022)",
            font=dict(color="white", size=14)
        ),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="rgba(255,255,255,0.7)"),
        xaxis=dict(
            gridcolor="rgba(255,255,255,0.08)",
            linecolor="rgba(255,255,255,0.25)",
            showgrid=True,
            zeroline=False,
        ),
        yaxis=dict(
            gridcolor="rgba(255,255,255,0.08)",
            linecolor="rgba(255,255,255,0.25)",
            showgrid=True,
            zeroline=False,
        ),
        height=380,
        margin=dict(l=40, r=40, t=40, b=40),
        hovermode="closest",
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
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
        title=dict(
            text="HCV Diagnosis/Treatment Access vs Disease Burden (2022)",
            font=dict(color="white", size=14)
        ),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="rgba(255,255,255,0.7)"),
        xaxis=dict(
            gridcolor="rgba(255,255,255,0.08)",
            linecolor="rgba(255,255,255,0.25)",
            showgrid=True,
            zeroline=False,
        ),
        yaxis=dict(
            gridcolor="rgba(255,255,255,0.08)",
            linecolor="rgba(255,255,255,0.25)",
            showgrid=True,
            zeroline=False,
        ),
        height=380,
        margin=dict(l=40, r=40, t=40, b=40),
        hovermode="closest",
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
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
    Output("hbv-epi-trends-row", "style"),
    Output("hcv-epi-trends-row", "style"),
    Input("selected-virus", "data"),
)
def toggle_epi_trends_rows(virus):
    if virus == "HBV":
        return {"display": "flex"}, {"display": "none"}
    elif virus == "HCV":
        return {"display": "none"}, {"display": "flex"}
    else:
        return {"display": "none"}, {"display": "none"}





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
    Output("year-range-slider", "value"),
    Output("year-range-slider", "marks"),
    Output('continent-dropdown', 'options'),
    Output('country-dropdown', 'options'),
    Output('genotype-dropdown', 'options'),
    Input("selected-virus", "data"),
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
    Output("year-range-label", "children"),
    Input("year-range-slider", "value")
)
def update_year_range_label(year_range):
    if not year_range or len(year_range) != 2:
        return ""
    return f"{year_range[0]} – {year_range[1]}"


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
    


# === SMALL UI TOGGLES ==============================================================
@callback(
    Output("selected-virus", "data"),
    Output("btn-hbv", "className"),
    Output("btn-hcv", "className"),
    Output("btn-hev", "className"),
    Input("btn-hbv", "n_clicks"),
    Input("btn-hcv", "n_clicks"),
    Input("btn-hev", "n_clicks"),
    prevent_initial_call=True
)
def update_virus(btn_hbv_clicks, btn_hcv_clicks, btn_hev_clicks):
    ctx = callback_context
    base = "hep-virus-item"
    active = "hep-virus-item hep-virus-active"

    if not ctx.triggered:
        return "HBV", active, base, base
    
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if triggered_id == "btn-hbv":
        return "HBV", active, base, base
    elif triggered_id == "btn-hcv":
        return "HCV", base, active, base
    elif triggered_id == "btn-hev":
        return "HEV", base, base, active
    
    return "HBV", active, base, base 


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

        