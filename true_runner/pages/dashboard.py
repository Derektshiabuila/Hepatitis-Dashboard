import dash
from dash import html, register_page, dcc, callback, Input, Output, State
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
pio.templates.default = "plotly_white"


register_page(__name__, path="/", name="Dashboard", order=0)
register_page(__name__, path="/dashboard", name="Dashboard")

# Import data loading functions
from data_loader import load_and_preprocess_data
try:
    from Full_Hepatitis_page import cache
except ImportError:
    class DummyCache:
        def memoize(self, *args, **kwargs):
            return lambda f: f
    cache = DummyCache()

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
    "HBV": ["hepadna", "hepatitis b", "hbsag", "hbcag", "hbv"],
    "HCV": ["flaviviri", "hepatitis c", "hcv", "ns5b", "ns3"],
    "HEV": ["hepeviri", "hepatitis e", "hev", "orf2"],
}


# Initialize data store (will be loaded on first access)
data_store = None

def get_data_store():
    """Helper function to access the global data store"""
    global data_store
    if data_store is not None:
        return data_store

    import flask
    try:
        if flask.has_app_context() and "DATA_STORE" in flask.current_app.config:
            data_store = flask.current_app.config["DATA_STORE"]
            if data_store is not None:
                print("Data store retrieved from flask.current_app.config")
                return data_store
    except Exception as e:
        print(f"Failed to access flask config: {e}")

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
HBV_GENOTYPE_COLORS = {
    "HBV-A": "#cdac02",
    "HBV-B": "#951de0",
    "HBV-C": "#016301",
    "HBV-D": "#1f23bb",
    "HBV-E": "#770104",
    "HBV-F": "#7f7340",
    "HBV-G": "#da760c",
    "HBV-H": "#dea8b1",
    "HBV-I": "#e9e905",
    "HBV-J": "#07e705",
    "Recombinant": "#e00603"
}

HCV_GENOTYPE_COLORS = {
    "HCV-1": "#cdac02",
    "HCV-2": "#951de0",
    "HCV-3": "#016301",
    "HCV-4": "#1f23bb",
    "HCV-5": "#770104",
    "HCV-6": "#7f7340",
    "HCV-7": "#da760c",
    "HCV-8": "#FF69B4",
    "Recombinant": "#e00603"
}

HEV_GENOTYPE_COLORS = {
    "HEV-1": "#cdac02",
    "HEV-2": "#951de0",
    "HEV-3": "#016301",
    "HEV-4": "#1f23bb",
    "HEV-5": "#770104",
    "HEV-6": "#7f7340",
    "HEV-7": "#da760c",
    "HEV-8": "#FF69B4",
    "Recombinant": "#e00603"
}

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

    burden = ihme_df[
        (ihme_df["cause"] == cause)
        & (ihme_df["measure"] == measure)
        & (ihme_df["metric"] == metric)
        & (ihme_df["sex"] == sex)
        & (ihme_df["age"] == "All ages")
    ].copy()

    # Filter by years
    if selected_years:
        y0, y1 = selected_years
        burden = burden[(burden["year"] >= y0) & (burden["year"] <= y1)]

    # Aggregate per country (latest year in range)
    if not burden.empty:
        burden = burden[burden["year"] == burden["year"].max()]
    burden = burden.groupby("Country_standard", as_index=False)["val"].sum().rename(columns={"val": "burden"})

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

    # Sum across age groups for each country
    if not latest.empty:
        # Group by country and sum values
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
        genotype_colors = HBV_GENOTYPE_COLORS
    elif virus_type == "HCV":
        genotype_colors = HCV_GENOTYPE_COLORS
    else:
        genotype_colors = HEV_GENOTYPE_COLORS
    
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
        HCV_COLOR_SCALE = [[0.0, "#FFF7BC"], [0.25, "#FEC44F"], [0.5, "#EC7014"], [0.75, "#993404"], [1.0, "#662506"]]
        HBV_COLOR_SCALE = [[0.0, "#F7FCF0"], [0.25, "#A8DDB5"], [0.5, "#2B8CBE"], [0.75, "#084081"], [1.0, "#06214D"]]
        HEV_COLOR_SCALE = [[0.0, "#F7FCF0"], [0.25, "#A8DDB5"], [0.5, "#2B8CBE"], [0.75, "#084081"], [1.0, "#06214D"]]
        
        if virus_type == "HBV":
           colorscale = HBV_COLOR_SCALE
        elif virus_type == "HCV":
            colorscale = HCV_COLOR_SCALE
        else:
            HEV_COLOR_SCALE
        
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

        HBV_COLOR_SCALE = [[0.0, "#deebf7"], [0.25, "#9ecae1"], [0.5, "#6baed6"], [0.75, "#3182bd"], [1.0, "#08519c"]]
        HCV_COLOR_SCALE = [[0.0, "#feedde"], [0.25, "#fdbe85"], [0.5, "#fd8d3c"], [0.75, "#e6550d"], [1.0, "#a63603"]]
        HEV_COLOR_SCALE = [[0.0, "#feedde"], [0.25, "#fdbe85"], [0.5, "#fd8d3c"], [0.75, "#e6550d"], [1.0, "#a63603"]]
        
        if virus_type == "HBV":
            colorscale = HBV_COLOR_SCALE
        elif virus_type == "HCV":
            colorscale = HCV_COLOR_SCALE
        else:
            HEV_COLOR_SCALE

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
        
        if virus_type == "HBV":
            colors = HBV_COLORS
        elif virus_type == "HCV":
            colors = HCV_COLORS
        else:
            colors = HEV_COLORS
        
        colorscale = [[i / (len(labels) - 1), c] for i, c in enumerate(colors)]

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
        [0.0, "#d73027"],   # Red (lowest coverage)
        [0.25, "#fc8d59"],
        [0.5, "#fee08b"],
        [0.75, "#d9ef8b"],
        [1.0, "#1a9850"]    # Dark green (highest coverage)
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
    
    # Group by year and genotype first
    line_data = (
        filtered_df.groupby(["Year", "genotype"])
        .size().reset_index(name="Genome Sequences")
    )
    
    # Apply rolling average for each genotype
    smoothed_data = []
    for genotype in line_data["genotype"].unique():
        genotype_df = line_data[line_data["genotype"] == genotype].copy()
        genotype_df = genotype_df.sort_values("Year")
        
        # Apply 3-year rolling average
        genotype_df["Smoothed_Sequences"] = (
            genotype_df["Genome Sequences"]
            .rolling(window=3, min_periods=1, center=True)
            .mean()
        )
        smoothed_data.append(genotype_df)
    
    smoothed_df = pd.concat(smoothed_data, ignore_index=True)
    
    # FIXED: Handle zeros for log scale
    smoothed_df["Smoothed_Sequences"] = smoothed_df["Smoothed_Sequences"].replace(0, 0.1)  # Avoid log(0)
    
    fig = px.line(
        smoothed_df,
        x="Year",
        y="Smoothed_Sequences",
        color="genotype",
        color_discrete_map=genotype_colors,
        markers=True,
        line_shape="spline"
    )
    
    fig.update_traces(mode='lines+markers', marker=dict(size=4))
    
    if not smoothed_df.empty:
        year_min = int(smoothed_df["Year"].min())
        year_max = int(smoothed_df["Year"].max())
        fig.update_layout(
            xaxis=dict(
                tickmode="array",
                tickvals=list(range(year_min, year_max+1, 4))
            )
        )
    
    # FIXED: Safe log range calculation
    y_pos = line_data["Genome Sequences"].replace(0, 0.1).dropna()  # Avoid zeros
    if len(y_pos):
        y_min = float(np.nanmax([0.1, np.nanmin(y_pos)]))  # Ensure positive
        y_max = float(np.nanmax(y_pos))
        lo = np.log10(max(0.1, y_min/1.5))
        hi = np.log10(y_max*1.5)
    else:
        lo, hi = 0, 2

    fig.update_layout(
        xaxis_title="Year",
        yaxis=dict(
            title="Number of Sequences (log scale)",
            type="log",
            autorange=False,
            range=[lo, hi],
            tickformat="~s",
            gridcolor="rgba(0,0,0,0.08)",
            linecolor="rgba(0,0,0,0.25)",
            zeroline=False
        ),
        height=400,
        margin=dict(t=100, b=0, l=0, r=0),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        plot_bgcolor="white",
        paper_bgcolor="white",
        meta={"y_full_log_range": [lo, hi]}
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
        ymin = float(vals[vals > 0].min()) if (vals > 0).any() else 0.1
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
        template="plotly_white",
    )

    # Hover text
    unit = " per million" if y_col == "PerMillion" else ""
    bar_fig.update_traces(
        hovertemplate="<b>%{x}</b><br>%{y:,}" + unit + " sequences<extra></extra>"
    )

    # Layout & legend
    # FIXED: Simplified y-axis configuration
    bar_fig.update_layout(
        title=f"Total {selected_virus.upper()} Sequences by genotype",
        xaxis_title="genotype",
        yaxis=dict(
            title=y_title + " (log scale)",
            type="log",
            # Let Plotly handle the autorange and ticks
            autorange=True,  # Change from False to True
            # Remove custom range to let Plotly calculate
            # Remove tickmode, tickvals, ticktext
            gridcolor="rgba(0,0,0,0.08)",
            zeroline=False,
            linecolor="rgba(0,0,0,0.25)",
            # Add this to ensure nice tick formatting
            tickformat=".0f",  # Format as integers
        ),
        bargap=0.25,
        plot_bgcolor="white",
        paper_bgcolor="white",
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
            height=450, plot_bgcolor="white", paper_bgcolor="white",
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
            height=450, plot_bgcolor="white", paper_bgcolor="white",
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
        template="plotly_white"
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
        plot_bgcolor="white", 
        paper_bgcolor="white",
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
            plot_bgcolor="white",
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
        template="plotly_white"
    )
    fig.update_layout(
        height=450,
        margin=dict(t=60, r=10, b=30, l=10),
        xaxis=dict(title=x_label, gridcolor="rgba(0,0,0,0.08)", linecolor="rgba(0,0,0,0.25)"),
        yaxis=dict(title=""),
        plot_bgcolor="white",
        paper_bgcolor="white",
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

@callback(
    Output("forecast-debug", "children"),
    Input("selected-virus", "data"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
    Input("correlation-sex", "value"),
)
def debug_forecast_data(virus, regions, countries, sex):
    data = get_data_store()
    ihme_df = data["ihme_df"]
    
    if ihme_df.empty:
        return html.Div("No IHME data loaded", style={"color": "red"})
    
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause_filter = cause_lookup.get((virus or "HBV").upper())
    
    # Check if cause exists
    if cause_filter not in ihme_df["cause"].values:
        return html.Div(f"Cause '{cause_filter}' not found in data", style={"color": "red"})
    
    # Filter data
    base = ihme_df[
        (ihme_df["sex"] == sex) &
        (ihme_df["cause"] == cause_filter) &
        (ihme_df["metric"] == "Number")
    ].copy()
    
    if regions:
        base = base[base["WHO_Regions"].isin(regions)]
    if countries:
        base = base[base["Country_standard"].isin(countries)]
    
    # Check what measures are available
    available_measures = base["measure"].unique()
    
    # Get sample data
    sample_data = base.head(3) if not base.empty else pd.DataFrame()
    
    return html.Div([
        html.P(f"Virus: {virus}, Cause: {cause_filter}"),
        html.P(f"Sex filter: {sex}"),
        html.P(f"Available measures: {list(available_measures)}"),
        html.P(f"Total rows: {len(base)}"),
        html.P("Sample data:"),
        html.Pre(sample_data[["year", "measure", "metric", "age", "sex", "val"]].to_string() if not sample_data.empty else "No data")
    ], style={"fontSize": "10px", "color": "#666", "padding": "10px", "backgroundColor": "#f0f0f0"})

#Time Series with Projections
def create_forecast_chart(ihme_df, selected_virus, sex, selected_regions=None, selected_countries=None):
    """Show historical trends with proper statistical forecasting"""
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause_filter = cause_lookup.get((selected_virus or "HBV").upper())
    
    if cause_filter not in ihme_df["cause"].values:
        return _empty_plot(f"No {selected_virus} burden data available for forecasting")
    
    # First, let's see what metrics are available for this cause
    cause_data = ihme_df[ihme_df["cause"] == cause_filter]
    available_metrics = cause_data["metric"].unique()
    
    # Try different metrics in order of preference
    metric_to_use = None
    for metric in ["Number", "Rate", "Percent"]:
        if metric in available_metrics:
            metric_to_use = metric
            break
    
    if not metric_to_use:
        return _empty_plot(f"No suitable metric found. Available: {list(available_metrics)}")
    
    # Handle "Both" sexes by summing Male and Female
    if sex == "Both":
        # Get Male data
        male_data = ihme_df[
            (ihme_df["cause"] == cause_filter) &
            (ihme_df["metric"] == metric_to_use) &
            (ihme_df["sex"] == "Male")
        ].copy()
        
        # Get Female data
        female_data = ihme_df[
            (ihme_df["cause"] == cause_filter) &
            (ihme_df["metric"] == metric_to_use) &
            (ihme_df["sex"] == "Female")
        ].copy()
        
        # Combine
        base = pd.concat([male_data, female_data])
    else:
        # Use the sex as-is
        base = ihme_df[
            (ihme_df["sex"] == sex) &
            (ihme_df["cause"] == cause_filter) &
            (ihme_df["metric"] == metric_to_use)
        ].copy()
    
    # Apply filters
    if selected_regions:
        base = base[base["WHO_Regions"].isin(selected_regions)]
    if selected_countries:
        base = base[base["Country_standard"].isin(selected_countries)]
    
    if base.empty:
        return _empty_plot(f"No data available for {selected_virus} with current filters")
    
    # Get the latest year
    if base["year"].notna().any():
        latest_data_year = int(base["year"].max())
    else:
        return _empty_plot("No valid year data")
    
    forecast_years = 8  # Forecast to 2030
    fig = go.Figure()
    
    # Colors for different measures
    measure_colors = {
        "Prevalence": "#1f77b4",
        "Incidence": "#ff7f0e", 
        "Deaths": "#d62728"
    }
    
    data_found = False
    
    for measure in ["Prevalence", "Incidence", "Deaths"]:
        measure_data = base[base["measure"] == measure].copy()
        
        if measure_data.empty:
            continue
            
        # Aggregate by year
        yearly_data = measure_data.groupby("year")["val"].sum().reset_index()
        yearly_data = yearly_data.sort_values("year")
        
        # Ensure we have positive values for log scale
        yearly_data = yearly_data[yearly_data["val"] > 0]
        
        if len(yearly_data) < 2:
            if not yearly_data.empty:
                # Show data even if insufficient for forecasting
                fig.add_trace(go.Scatter(
                    x=yearly_data["year"],
                    y=yearly_data["val"],
                    name=f"{measure} ({len(yearly_data)} points)",
                    line=dict(color=measure_colors[measure], width=2),
                    mode='markers',
                    marker=dict(size=8)
                ))
                data_found = True
            continue
        
        # Prepare data for forecasting
        years = yearly_data["year"].values.astype(float)
        values = yearly_data["val"].values.astype(float)
        
        # Try linear regression
        try:
            # Calculate linear regression
            n = len(years)
            sum_x = np.sum(years)
            sum_y = np.sum(values)
            sum_xy = np.sum(years * values)
            sum_x2 = np.sum(years * years)
            
            denominator = n * sum_x2 - sum_x * sum_x
            if denominator != 0:
                m = (n * sum_xy - sum_x * sum_y) / denominator
                b = (sum_y - m * sum_x) / n
                
                # Create forecast
                future_years = np.arange(latest_data_year + 1, latest_data_year + forecast_years + 1)
                forecast_pred = m * future_years + b
                forecast_pred = np.maximum(forecast_pred, 0)  # Ensure non-negative
                
                # Add historical data
                fig.add_trace(go.Scatter(
                    x=years,
                    y=values,
                    name=f"{measure} (Historical)",
                    line=dict(color=measure_colors[measure], width=3),
                    mode='lines+markers',
                    marker=dict(size=6),
                    hovertemplate=f"{measure}: %{{y:,.2f}}<extra></extra>"
                ))
                
                # Add forecast
                fig.add_trace(go.Scatter(
                    x=future_years,
                    y=forecast_pred,
                    name=f"{measure} (Forecast)",
                    line=dict(color=measure_colors[measure], width=2, dash='dash'),
                    mode='lines',
                    hovertemplate=f"{measure} Forecast: %{{y:,.2f}}<extra></extra>"
                ))
                
                data_found = True
                
            else:
                # If regression fails, just show historical
                fig.add_trace(go.Scatter(
                    x=years,
                    y=values,
                    name=f"{measure} (No Trend)",
                    line=dict(color=measure_colors[measure], width=2),
                    mode='lines+markers',
                    marker=dict(size=6)
                ))
                data_found = True
                
        except Exception as e:
            # Fallback to just showing data
            fig.add_trace(go.Scatter(
                x=years,
                y=values,
                name=f"{measure} (Error)",
                line=dict(color=measure_colors[measure], width=2),
                mode='markers',
                marker=dict(size=6)
            ))
            data_found = True
    
    if not data_found:
        # Check what measures ARE available
        available_measures = base["measure"].unique()
        return _empty_plot(
            f"No forecast data available for {selected_virus}. "
            f"Available measures: {list(available_measures)}. "
            f"Using metric: {metric_to_use}"
        )
    
    # Add reference lines
    if latest_data_year:
        fig.add_vline(
            x=latest_data_year, 
            line_dash="dot", 
            line_color="red",
            line_width=1.5,
            annotation_text="Data Limit"
        )
    
    # WHO 2030 target
    who_target_year = 2030
    fig.add_vline(
        x=who_target_year, 
        line_dash="dot", 
        line_color="green",
        line_width=1.5,
        annotation_text="WHO 2030"
    )
    
    # Set up layout
    fig.update_layout(
        xaxis_title="Year",
        yaxis_title=f"Value ({metric_to_use})",
        yaxis_type="log",
        height=400,
        hovermode="x unified",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="center",
            x=0.5
        ),
        margin=dict(t=40, b=40, l=60, r=20)
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
        xaxis_title="Year",
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
        yaxis_title="Number of Sequences",
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02,
            title="genotype"
        ),
        margin=dict(r=150),  # Add margin for legend
        xaxis={'categoryorder': 'total descending'}  # Sort by total sequences
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
    tickvals = MASTER_TICKS[mask].tolist()
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
            gridcolor="rgba(0,0,0,0.06)", showgrid=True
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
            zeroline=False
        ),
        plot_bgcolor="white",
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
        plot_bgcolor="white",
        paper_bgcolor="white"
    )
    
    return fig

# === HELPERS FOR LAYOUT =====================================================
def build_indicators(virus):
    color_class = "bg-primary" if virus == "HBV" else "bg-warning"

    return dbc.Col([
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.Row([
                        dbc.Col([
                            dbc.CardBody([
                                html.H6("Total Whole Genomes", className="card-subtitle"),
                                html.H4(id="indicator-total", className="card-title")
                            ])
                        ], width=10),
                        dbc.Col([
                            html.Div([
                                html.I(className="bi bi-bar-chart-fill", style={"fontSize": "2rem"})
                            ], className=f"d-flex align-items-center justify-content-center {color_class} text-white h-100")
                        ], width=2)
                    ], className="g-0")
                ], className="mb-4 shadow-sm")
            ], width=12),

            dbc.Col([
                dbc.Card([
                    dbc.Row([
                        dbc.Col([
                            dbc.CardBody([
                                html.H6("Countries"),
                                html.H4(id="indicator-countries")
                            ])
                        ], width=10),
                        dbc.Col([
                            html.Div([
                                html.I(className="bi bi-globe", style={"fontSize": "2rem"})
                            ], className=f"d-flex align-items-center justify-content-center {color_class} text-white h-100")
                        ], width=2)
                    ], className="g-0")
                ], className="mb-4 shadow-sm")
            ], width=12),

            dbc.Col([
                dbc.Card([
                    dbc.Row([
                        dbc.Col([
                            dbc.CardBody([
                                html.H6("Genotypes"),
                                html.H4(id="indicator-genotypes")
                            ])
                        ], width=10),
                        dbc.Col([
                            html.Div([
                                html.I(className="bi bi-diagram-3-fill", style={"fontSize": "2rem"})
                            ], className=f"d-flex align-items-center justify-content-center {color_class} text-white h-100")
                        ], width=2)
                    ], className="g-0")
                ], className="mb-4 shadow-sm")
            ], width=12),

            dbc.Col([
                dbc.Card([
                    dbc.Row([
                        dbc.Col([
                            dbc.CardBody([
                                html.H6("Years"),
                                html.H4(id="indicator-years")
                            ])
                        ], width=10),
                        dbc.Col([
                            html.Div([
                                html.I(className="bi bi-hourglass-split", style={"fontSize": "2rem"})
                            ], className=f"d-flex align-items-center justify-content-center {color_class} text-white h-100")
                        ], width=2)
                    ], className="g-0")
                ], className="mb-4 shadow-sm")
            ], width=12),

        ], className="mb-4 shadow-sm"),
    ], className="mb-4", width=3)

# === APP SETUP ==============================================================
external_stylesheets = [
    "https://cdn.jsdelivr.net/npm/bootstrap-icons@1.10.5/font/bootstrap-icons.css",
    dbc.themes.BOOTSTRAP
]

#app = dash.Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)

# === APP LAYOUT ==============================================================
def create_dashboard_layout():
    return dbc.Container([
        # --- page-level state stores ---
        dcc.Store(id="selected-virus", data="HBV"),
        dcc.Store(id="filtered-store"),
        dcc.Store(id="gap-store"),
        dcc.Store(id="computed-metrics-store"),
        dcc.Store(id="ihme-latest-store"),
        dcc.Store(id="priority-data-store"),
        dcc.Download(id="priority-download"),
        dcc.Download(id="download-mutations"),
        dcc.Download(id="download-mutation-report"),
        *USER_SEQ_STORES(),
        
        # Header Section
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.Img(src="/assets/logo.png", height="40px", className="me-2"),
                    html.H4("Hepatitis Virus Sequence Dashboard", className="mb-0")
                ], className="d-flex align-items-center")
            ], width=4),
            
            dbc.Col([
                html.Div([
                    dbc.Button("HBV", id="btn-hbv", color="primary", size="lg", 
                              className="me-2", outline=False, n_clicks=1),
                    dbc.Button("HCV", id="btn-hcv", color="warning", size="lg", 
                              className="me-2", outline=True, n_clicks=0),
                    dbc.Button("HEV", id="btn-hev", color="success", size="lg",
                              outline=True, n_clicks=0)
                ])
            ], width="auto"),
            
            dbc.Col([
                dbc.DropdownMenu(
                    label="Download",
                    children=[
                        dbc.DropdownMenuItem("Data", id="btn-download-data"),
                        dbc.DropdownMenuItem("Reports", href="/about#report-section"),
                    ],
                    color="success",
                    size="lg",
                    className="me-2"
                ),
                dbc.Toast("Your download is starting…", id="dl-toast", header="Download", is_open=False,
                  dismissable=True, icon="success", duration=3000, className="position-fixed top-0 end-0 m-3"),
                html.Div(id="download-trigger", style={"display": "none"}),
                dcc.Download(id="download-data"),
            ], width="auto"),
        ], className="mb-4 align-items-center", justify="between"),

        # === SIMPLE TABS ===
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H6("Navigation", className="text-muted mb-3"),
                        dbc.ButtonGroup([
                            dbc.Button("📊 Overview", id="tab-overview", color="primary", className="active", n_clicks=1),
                            dbc.Button("🧬 Mutations", id="tab-mutations", color="secondary", n_clicks=0),
                            dbc.Button("📈 Epidemiology", id="tab-epidemiology", color="secondary", n_clicks=0),
                            user_seq_tab_button(),
                        ], className="w-100")
                    ])
                ], className="mb-3 shadow-sm")
            ], width=12)
        ]),

        # === COMMON FILTERS (ALWAYS VISIBLE) ===
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H6("FILTERS", className="text-muted mb-3"),
                        
                        dbc.Row([
                            dbc.Col([
                                html.Label("Year Range", className="fw-bold mb-2"),
                                html.Div(id="selected-years-display", className="text-primary fw-bold mb-2"),
                                dcc.RangeSlider(
                                    id="year-slider",
                                    step=1,
                                    tooltip={"placement": "bottom", "always_visible": True},
                                    className="mb-3"
                                ),
                                dbc.Button("Reset Time Range", id="btn-reset-time", 
                                         color="secondary", outline=True, size="sm")
                            ], width=12)
                        ], className="mb-3"),
                        
                        dbc.Row([
                            dbc.Col([
                                html.Label("WHO Region(s)", className="fw-bold mb-2"),
                                dcc.Dropdown(id="continent-dropdown", multi=True, 
                                           placeholder="All Regions", className="mb-3")
                            ], width=4),
                            
                            dbc.Col([
                                html.Label("Country(s)", className="fw-bold mb-2"),
                                dcc.Dropdown(id="country-dropdown", multi=True, 
                                           placeholder="All Countries", className="mb-3")
                            ], width=4),
                            
                            dbc.Col([
                                html.Label("genotype(s)", className="fw-bold mb-2"),
                                dcc.Dropdown(id="genotype-dropdown", multi=True, 
                                           placeholder="Select virus first", className="mb-3")
                            ], width=4)
                        ])
                    ])
                ], className="mb-4 shadow-sm")
            ], width=12)
        ], id="common-filters"),

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
                                            {"label": " Raw Count", "value": "raw"},
                                            {"label": " Per Million", "value": "PerMillion"}
                                        ],
                                        value="raw",
                                        inline=True,
                                        className="g-3 align-items-end mb-2"
                                    )
                                ], width=6),
                                
                                dbc.Col([
                                    html.Label("Map Mode:", className="fw-bold me-2"),
                                    dcc.RadioItems(
                                        id="map-mode",
                                        options=[
                                            {"label": " Sequences", "value": "sequences"},
                                            {"label": " Coverage", "value": "coverage"}, 
                                            {"label": "Epidemiology", "value": "epidemiology"},
                                        ],
                                        value="sequences",
                                        inline=True,
                                        inputStyle={"marginRight": "6px", "marginLeft": "12px"}
                                    ),
                                ], width="auto"),
                            ], className="g-3 align-items-end mb-2"),
                        ])
                    ], className="mb-3 shadow-sm")
                ], width=12)
            ]),

            # ROW 1: Map + Indicators
            dbc.Row([
                # Indicators
                build_indicators("HBV"),
                
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
                ], width=9)
            ], className="mb-4"),
            
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
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Mutation Filters", className="mb-3"),
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

        # === TAB 3: EPIDEMIOLOGY CONTENT ===
        html.Div(id="epidemiology-content", style={"display": "none"}, children=[
            # Epidemiology Controls
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Epidemiology Analysis Settings", className="mb-3"),
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Burden Metric:", className="fw-bold me-2"),
                                    dcc.Dropdown(
                                        id="epi-burden-metric",
                                        options=[
                                            {"label": "Prevalence", "value": "Prevalence|Number"},
                                            {"label": "Incidence", "value": "Incidence|Number"},
                                            {"label": "Deaths", "value": "Deaths|Number"},
                                            {"label": "Prevalence Rate", "value": "Prevalence|Rate"},
                                            {"label": "Incidence Rate", "value": "Incidence|Rate"},
                                            {"label": "Death Rate", "value": "Deaths|Rate"},
                                        ],
                                        value="Prevalence|Number",
                                        clearable=False,
                                        style={"width": "200px"}
                                    )
                                ], width=3),
                                
                                dbc.Col([
                                    html.Label("Age Group:", className="fw-bold me-2"),
                                    dcc.Dropdown(
                                        id="epi-age-group",
                                        options=[],
                                        value="All ages",
                                        clearable=False,
                                        style={"width": "200px"}
                                    )
                                ], width=3),
                                
                                dbc.Col([
                                    html.Label("Sex:", className="fw-bold me-2"),
                                    dcc.RadioItems(
                                        id="epi-sex-filter",
                                        options=[
                                            {"label": " Both", "value": "Both"},
                                            {"label": " Male", "value": "Male"},
                                            {"label": " Female", "value": "Female"}
                                        ],
                                        value="Both",
                                        inline=True
                                    )
                                ], width=3),
                                
                                dbc.Col([
                                    html.Label("Display Mode:", className="fw-bold me-2"),
                                    dcc.RadioItems(
                                        id="epi-display-mode",
                                        options=[
                                            {"label": " Number", "value": "Number"},
                                            {"label": " Rate", "value": "Rate"},
                                            {"label": " Percent", "value": "Percent"}
                                        ],
                                        value="Number",
                                        inline=True
                                    )
                                ], width=3),
                            ], className="g-3")
                        ])
                    ], className="shadow-sm")
                ], width=12)
            ], className="mb-4"),
            
            # ROW 1: Global Burden Forecast and Top Countries
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Burden Forecast with Projections", className="mb-3"),
                            dcc.Loading(
                                dcc.Graph(id="forecast-chart"),
                                type="circle"
                            )
                        ])
                    ], className="h-100 shadow-sm")
                ], width=8),
                
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Top Countries by Burden", className="mb-3"),
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Show Top:", className="fw-bold me-2"),
                                    dcc.Dropdown(
                                        id="top-countries-n",
                                        options=[{"label": str(i), "value": i} for i in [5, 10, 15, 20]],
                                        value=10,
                                        clearable=False,
                                        style={"width": "150px"}
                                    )
                                ], width=12)
                            ], className="mb-2"),
                            dcc.Loading(
                                dcc.Graph(id="top-countries-burden"),
                                type="circle"
                            )
                        ])
                    ], className="h-100 shadow-sm")
                ], width=4),
            ], className="mb-4"),
            
            # ROW 2: Age and Sex Analysis
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Age Distribution", className="mb-3"),
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Year:", className="fw-bold me-2"),
                                    dcc.Dropdown(
                                        id="age-dist-year",
                                        options=[],
                                        placeholder="Select year...",
                                        style={"width": "150px"}
                                    )
                                ], width=6),
                                dbc.Col([
                                    html.Label("Sex:", className="fw-bold me-2"),
                                    dcc.RadioItems(
                                        id="age-dist-sex",
                                        options=[
                                            {"label": " Both", "value": "Both"},
                                            {"label": " Male", "value": "Male"},
                                            {"label": " Female", "value": "Female"}
                                        ],
                                        value="Both",
                                        inline=True
                                    )
                                ], width=6),
                            ], className="mb-2"),
                            dcc.Loading(
                                dcc.Graph(id="age-distribution-chart"),
                                type="circle"
                            )
                        ])
                    ], className="h-100 shadow-sm")
                ], width=6),
                
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Sex Ratio Over Time", className="mb-3"),
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Age Group:", className="fw-bold me-2"),
                                    dcc.Dropdown(
                                        id="sex-ratio-age",
                                        options=[],
                                        value="All ages",
                                        clearable=False,
                                        style={"width": "200px"}
                                    )
                                ], width=12)
                            ], className="mb-2"),
                            dcc.Loading(
                                dcc.Graph(id="sex-ratio-chart"),
                                type="circle"
                            )
                        ])
                    ], className="h-100 shadow-sm")
                ], width=6),
            ], className="mb-4"),

            # ROW 3: Regional Analysis
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Regional Burden Comparison", className="mb-3"),
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Compare:", className="fw-bold me-2"),
                                    dcc.RadioItems(
                                        id="region-compare-type",
                                        options=[
                                            {"label": " Total", "value": "total"},
                                            {"label": " Per Capita", "value": "per_capita"},
                                            {"label": " Age-Standardized", "value": "age_standardized"}
                                        ],
                                        value="total",
                                        inline=True
                                    )
                                ], width=12)
                            ], className="mb-2"),
                            dcc.Loading(
                                dcc.Graph(id="region-comparison-chart"),
                                type="circle"
                            )
                        ])
                    ], className="h-100 shadow-sm")
                ], width=6),
                
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Regional Age Patterns", className="mb-3"),
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Region:", className="fw-bold me-2"),
                                    dcc.Dropdown(
                                        id="region-age-pattern",
                                        options=[],
                                        placeholder="Select region...",
                                        style={"width": "200px"}
                                    )
                                ], width=12)
                            ], className="mb-2"),
                            dcc.Loading(
                                dcc.Graph(id="region-age-pattern-chart"),
                                type="circle"
                            )
                        ])
                    ], className="h-100 shadow-sm")
                ], width=6),
            ], className="mb-4"),
            
            # ROW 4: Sequencing Priority and Correlation
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Sequencing Priority Ranking", className="mb-3"),
                            dbc.Row([
                                dbc.Col([
                                    dbc.Button(
                                        "Download Priority Table (CSV)",
                                        id="priority-download-btn",
                                        color="secondary",
                                        className="ms-2"
                                    ),
                                    dcc.Download(id="priority-download")
                                ], width="auto"),
                            ], className="g-3 mb-2"),
                            dcc.Loading(
                                dcc.Graph(id="priority-ranking"),
                                type="circle"
                            )
                        ])
                    ], className="h-100 shadow-sm")
                ], width=6),
                
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Burden vs. Sequencing Correlation", className="mb-3"),
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Age Group:", className="fw-bold me-2"),
                                    dcc.Dropdown(
                                        id="correlation-age",
                                        options=[],
                                        value="All ages",
                                        clearable=False,
                                        style={"width": "200px"}
                                    )
                                ], width=4),
                                dbc.Col([
                                    html.Label("Sex:", className="fw-bold me-2"),
                                    dcc.RadioItems(
                                        id="correlation-sex",
                                        options=[
                                            {"label": " Male", "value": "Male"},
                                            {"label": " Female", "value": "Female"},
                                            {"label": " Both (summed)", "value": "Both"}
                                        ],
                                        value="Both",
                                        inline=True
                                    )
                                ], width=4),
                                dbc.Col([
                                    html.Label("Metric:", className="fw-bold me-2"),
                                    dcc.Dropdown(
                                        id="correlation-metric",
                                        options=[
                                            {"label": "Prevalence", "value": "Prevalence|Number"},
                                            {"label": "Incidence", "value": "Incidence|Number"},
                                            {"label": "Deaths", "value": "Deaths|Number"}
                                        ],
                                        value="Prevalence|Number",
                                        clearable=False,
                                        style={"width": "200px"}
                                    )
                                ], width=4),
                            ], className="mb-3"),
                            dcc.Loading(
                                dcc.Graph(id="burden-coverage-scatter"),
                                type="circle"
                            )
                        ])
                    ], className="h-100 shadow-sm")
                ], width=6),
            ], className="mb-4"),
            
            # Data Table
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.Div([
                                html.H5("Detailed Epidemiology Data", className="d-inline mb-0"),
                                dbc.Button(
                                    "Back to Overview",
                                    id="btn-back-to-overview-from-epi",
                                    color="primary",
                                    size="sm",
                                    className="float-end"
                                ),
                                dbc.Button(
                                    "Download Data",
                                    id="btn-download-epi",
                                    color="secondary",
                                    size="sm",
                                    className="float-end me-2"
                                ),
                                dcc.Download(id="download-epi-data")
                            ], className="mb-3"),
                            html.Div(id="epidemiology-data-table")
                        ])
                    ], className="shadow-sm")
                ], width=12)
            ], className="mb-4"),
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

    ], fluid=True, style={'backgroundColor': '#f8f9fa', 'minHeight': '100vh', 'padding': '20px'})
    
# === STATE STORES ==============================================================
# — Year bounds & dropdown option lists (fast) —
layout = create_dashboard_layout()  # contains the year-slider

# Callbacks from user_sequence_analysis register automatically at import time.

@callback(Output("dl-toast", "is_open"), Input("btn-download-data", "n_clicks"), prevent_initial_call=True)
def _show_toast(n): return True


@callback(
    Output("selected-years-display", "children"),
    Input("year-slider", "value"),
    prevent_initial_call=True,       # optional
)
def show_years(value):
    return f"{value[0]} – {value[1]}"

# — Filtered sequence dataframe store —
@callback(
    Output("filtered-store", "data"),
    Input("selected-virus", "data"),
    Input("year-slider", "value"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
    Input("genotype-dropdown", "value"),
)
def compute_filtered_store(virus, years, regions, countries, genotypes):
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

    # year range
    if years and len(years) == 2:
        y0, y1 = years
    else:
        y0, y1 = int(base["Year"].min()), int(base["Year"].max())

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
    Input("epi-burden-metric", "value"),
    Input("epi-sex-filter", "value"),   # ✅ Changed to epi-sex-filter
)
def update_gap_store(filtered_json, virus, ihme_metric_choice, sex):
    return compute_gap_from_filtered(
        filtered_json=filtered_json,
        virus=virus,
        ihme_metric_choice=ihme_metric_choice,
        sex=sex,
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
    Input("year-slider", "value"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
    Input("correlation-sex", "value"),  # Make sure this is the right sex selector
)
def compute_ihme_latest_store(virus, metric, years, regions, countries, sex):
    data = get_data_store()
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
    Input("year-slider", "value"),
)
def update_indicators(filtered_json, selected_years):
    df = _df_from_json(filtered_json)
    if df is None or df.empty:
        return "0", "0", "0", "N/A"

    # Totals
    total_genomes = len(df)
    unique_countries = df.get("Country_standard", pd.Series(dtype="object")).nunique()

    # Genotypes: exclude recombinants from the count, but flag presence
    g = df.get("genotype", pd.Series(dtype="object")).astype("string").str.strip()
    recomb_mask = g.str.contains(r"recomb", case=False, na=False)  # catches 'Recombinant', 'Recombinants', etc.
    base_genotype_count = g[~recomb_mask].dropna().nunique()
    has_recomb = bool(recomb_mask.any())

    if has_recomb:
        genotypes_text = f"{base_genotype_count} + Recombinants"
    else:
        genotypes_text = f"{base_genotype_count}"

    # Years label
    years_text = f"{selected_years[0]} - {selected_years[1]}" \
        if selected_years and len(selected_years) == 2 else "All years"

    return f"{total_genomes:,}", str(unique_countries), genotypes_text, years_text


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
    
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
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
    Input("year-slider", "value"),
)
def update_mutation_section(filtered_json, selected_virus, regions, countries, years):
    """Updates the mutation section based on selected virus and filters"""
    
    data = get_data_store()
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
        empty.update_layout(title="No data", xaxis={"visible": False}, yaxis={"visible": False})
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
    store = get_data_store()
    
    # FIXED: Properly handle all three viruses
    if virus == "HBV":
        df = store["hbv_grouped"]
    elif virus == "HCV":
        df = store["hcv_grouped"]
    elif virus == "HEV":
        df = store["hev_grouped"]  # ADDED: HEV data
    else:
        df = store["hbv_grouped"]
    
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
@cache.memoize(timeout=86400)
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
            title = f"{selected_virus} Whole Genome Sequence Map"
            subtitle = "Sequences per million population"
        else:
            title = f"{selected_virus} Whole Genome Sequence Map" 
            subtitle = "Sequences (count)"
        
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
@callback(
    Output("tab-overview", "color"),
    Output("tab-mutations", "color"),
    Output("tab-epidemiology", "color"),
    Output("tab-user-seq", "color"),
    Output("tab-overview", "className"),
    Output("tab-mutations", "className"),
    Output("tab-epidemiology", "className"),
    Output("tab-user-seq", "className"),
    Output("overview-content", "style"),
    Output("mutations-content", "style"),
    Output("epidemiology-content", "style"),
    Output("user-seq-content", "style"),
    Output("common-filters", "style"),
    Output("tab-mutations", "disabled"),
    
    # Inputs
    Input("tab-overview", "n_clicks"),
    Input("tab-mutations", "n_clicks"),
    Input("tab-epidemiology", "n_clicks"),
    Input("tab-user-seq", "n_clicks"),
    Input("btn-back-to-overview-from-mutations", "n_clicks"),
    Input("btn-back-to-overview-from-epi", "n_clicks"),
    Input("btn-go-to-epidemiology", "n_clicks"),
    Input("btn-quick-forecast", "n_clicks"),
    Input("btn-quick-priority", "n_clicks"),
    Input("btn-quick-timeline", "n_clicks"),
    Input("selected-virus", "data"),
    
    # States to track initial load
    State("tab-overview", "n_clicks"),
    State("tab-mutations", "n_clicks"),
    State("tab-epidemiology", "n_clicks"),
)
def switch_tabs(overview_clicks, mutations_clicks, epidemiology_clicks,
                user_seq_clicks,
                back_from_mutations_clicks, back_from_epi_clicks,
                go_epidemiology_clicks, quick_forecast_clicks, quick_priority_clicks,
                quick_timeline_clicks, selected_virus,
                overview_state, mutations_state, epidemiology_state):
    
    ctx = callback_context
    
    # Shared helpers
    show = {"display": "block"}
    hide = {"display": "none"}

    def make_return(active_tab, mutations_disabled):
        mutations_color = "light" if mutations_disabled else "secondary"
        tabs = ["overview", "mutations", "epidemiology", "user-seq"]
        colors    = []
        classes   = []
        styles    = []
        for t in tabs:
            is_active = (t == active_tab)
            if t == "mutations":
                colors.append("light" if mutations_disabled else ("primary" if is_active else "secondary"))
            else:
                colors.append("primary" if is_active else "secondary")
            classes.append("active" if is_active else "")
            styles.append(show if is_active else hide)
        
        filters_style = hide if active_tab == "user-seq" else show
        return (*colors, *classes, *styles, filters_style, mutations_disabled)

    mutations_disabled = selected_virus == "HEV" if selected_virus else False

    # Handle initial load
    if not ctx.triggered:
        return make_return("overview", mutations_disabled)
    
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    # Handle virus change — keep whichever tab is currently visible
    if button_id == "selected-virus":
        if mutations_clicks and mutations_clicks > 0 and not mutations_disabled:
            return make_return("mutations", mutations_disabled)
        elif epidemiology_clicks and epidemiology_clicks > 0:
            return make_return("epidemiology", mutations_disabled)
        elif user_seq_clicks and user_seq_clicks > 0:
            return make_return("user-seq", mutations_disabled)
        return make_return("overview", mutations_disabled)
    
    # Normal tab switching
    if button_id in ["tab-overview", "btn-back-to-overview-from-mutations", "btn-back-to-overview-from-epi"]:
        return make_return("overview", mutations_disabled)
    
    elif button_id in ["tab-mutations", "btn-quick-timeline"]:
        if mutations_disabled:
            return make_return("overview", mutations_disabled)
        return make_return("mutations", mutations_disabled)
    
    elif button_id in ["tab-epidemiology", "btn-go-to-epidemiology", "btn-quick-forecast", "btn-quick-priority"]:
        return make_return("epidemiology", mutations_disabled)

    elif button_id == "tab-user-seq":
        return make_return("user-seq", mutations_disabled)

    # Default
    return make_return("overview", mutations_disabled)

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
        plot_bgcolor="white",
        paper_bgcolor="white"
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
    State("year-slider", "value"),
    prevent_initial_call=True
)
def download_mutation_report(n_clicks, filtered_json, virus, regions, countries, years):
    if not n_clicks:
        return dash.no_update
    
    data = get_data_store()
    selected_virus = virus or "HBV"
    
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
    report_lines.append(f"Year Range: {years if years else 'All'}")
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

# Forecast Chart Callback
@callback(
    Output("forecast-chart", "figure"),
    Input("selected-virus", "data"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
    Input("correlation-sex", "value"),
)
def update_forecast_chart(virus, regions, countries, sex):
    data = get_data_store()
    return create_forecast_chart(
        data["ihme_df"],
        virus or "HBV",
        sex or "Both",
        regions, 
        countries
    )
    
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
    Output("priority-ranking", "figure"),
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
    
    return fig, _df_to_json(priority_df)
    
# === EPIDEMIOLOGY CALLBACKS ===
@callback(
    Output("correlation-age", "options"),
    Output("correlation-age", "value"),
    Input("selected-virus", "data"),
    Input("correlation-metric", "value"),
)
def update_correlation_age_options(virus, metric):
    data = get_data_store()
    ihme_df = data["ihme_df"]
    
    if ihme_df.empty:
        options = [{"label": "All ages", "value": "All ages"}]
        return options, "All ages"
    
    # Parse metric
    try:
        measure, metric_type = metric.split("|")
    except:
        measure, metric_type = "Prevalence", "Number"
    
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause = cause_lookup.get((virus or "HBV").upper())
    
    # Get available age groups for this cause and metric
    filtered = ihme_df[
        (ihme_df["cause"] == cause) &
        (ihme_df["measure"] == measure) &
        (ihme_df["metric"] == metric_type)
    ]
    
    if filtered.empty:
        # Fall back to all age groups for this cause
        filtered = ihme_df[ihme_df["cause"] == cause]
    
    # Get unique age groups and sort them logically
    age_groups = sorted(filtered["age"].dropna().unique())
    
    # Try to sort age groups in a logical order
    def sort_age_key(age):
        if not isinstance(age, str):
            return float('inf')
        
        age_lower = age.lower()
        
        # Handle months first
        if 'month' in age_lower:
            nums = re.findall(r'\d+', age)
            if nums:
                return int(nums[0]) - 1000
        
        # Handle years
        elif 'year' in age_lower:
            if '<' in age:
                nums = re.findall(r'\d+', age)
                if nums:
                    return int(nums[0]) - 500
            elif '-' in age:
                nums = re.findall(r'\d+', age)
                if nums:
                    return int(nums[0])
        
        return float('inf')
    
    try:
        age_groups.sort(key=sort_age_key)
    except:
        age_groups.sort()
    
    # Create options - include "All ages" first
    options = [{"label": "All ages (sum of all age groups)", "value": "All ages"}]
    
    # Add actual age groups
    for age in age_groups:
        if age != "All ages":
            options.append({"label": age, "value": age})
    
    # Default value
    default_value = "All ages"
    
    return options, default_value

@callback(
    Output("global-burden-timeline", "figure"),
    Input("selected-virus", "data"),
    Input("epi-burden-metric", "value"),
    Input("epi-age-group", "value"),
    Input("epi-sex-filter", "value"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
    Input("epi-display-mode", "value"),
)
def update_global_burden_timeline(virus, metric, age_group, sex, regions, countries, display_mode):
    data = get_data_store()
    ihme_df = data["ihme_df"]
    
    if ihme_df.empty:
        return _empty_plot("No IHME data available")
    
    # Parse metric
    try:
        measure, metric_type = metric.split("|")
        if display_mode in ["Rate", "Percent", "Number"]:
            metric_type = display_mode
    except:
        measure, metric_type = "Prevalence", display_mode or "Number"
    
    # Filter data
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause = cause_lookup.get((virus or "HBV").upper())
    
    # CORRECTED: Handle "Both" sexes by getting data for each sex separately
    if sex == "Both":
        # Get data for Male and Female separately
        if age_group == "All ages":
            male_data = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["sex"] == "Male")
            ].copy()
            female_data = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["sex"] == "Female")
            ].copy()
        else:
            male_data = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["age"] == age_group) &
                (ihme_df["sex"] == "Male")
            ].copy()
            female_data = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["age"] == age_group) &
                (ihme_df["sex"] == "Female")
            ].copy()
        
        # Combine the data (sum will happen in groupby)
        filtered = pd.concat([male_data, female_data])
    else:
        # Use the sex as-is
        if age_group == "All ages":
            filtered = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["sex"] == sex)
            ].copy()
        else:
            filtered = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["age"] == age_group) &
                (ihme_df["sex"] == sex)
            ].copy()
    
    # Apply region/country filters
    if regions:
        filtered = filtered[filtered["WHO_Regions"].isin(regions)]
    if countries:
        filtered = filtered[filtered["Country_standard"].isin(countries)]
    
    if filtered.empty:
        return _empty_plot(f"No {measure} data found for {virus} ({age_group}, {sex})")
    
    # Aggregate by year - this will sum Male and Female when sex="Both"
    yearly = filtered.groupby("year")["val"].sum().reset_index()
    yearly = yearly.sort_values("year")
    
    # Create figure
    fig = px.line(
        yearly,
        x="year",
        y="val",
        markers=True,
        title=f"{measure} Trend ({virus}, {age_group}, {sex})"
    )
    
    fig.update_traces(
        line=dict(width=3),
        hovertemplate="Year: %{x}<br>Value: %{y:,.2f}<extra></extra>"
    )
    
    fig.update_layout(
        height=400,
        xaxis_title="Year",
        yaxis_title=f"{measure} ({metric_type})",
        hovermode="x unified"
    )
    
    return fig


@callback(
    Output("top-countries-burden", "figure"),
    Input("selected-virus", "data"),
    Input("epi-burden-metric", "value"),
    Input("epi-age-group", "value"),
    Input("epi-sex-filter", "value"),
    Input("continent-dropdown", "value"),
    Input("top-countries-n", "value"),
)
def update_top_countries_burden(virus, metric, age_group, sex, regions, top_n):
    data = get_data_store()
    ihme_df = data["ihme_df"]
    
    if ihme_df.empty:
        return _empty_plot("No IHME data available")
    
    # Parse metric
    try:
        measure, metric_type = metric.split("|")
    except:
        measure, metric_type = "Prevalence", "Number"
    
    # Filter data
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause = cause_lookup.get((virus or "HBV").upper())
    
    # Handle "Both" sexes by summing Male and Female
    if sex == "Both":
        # Get Male data
        if age_group == "All ages":
            filtered_male = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["sex"] == "Male")
            ].copy()
        else:
            filtered_male = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["age"] == age_group) &
                (ihme_df["sex"] == "Male")
            ].copy()
        
        # Get Female data
        if age_group == "All ages":
            filtered_female = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["sex"] == "Female")
            ].copy()
        else:
            filtered_female = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["age"] == age_group) &
                (ihme_df["sex"] == "Female")
            ].copy()
        
        # Combine Male and Female data
        filtered = pd.concat([filtered_male, filtered_female])
    else:
        # Use the sex as-is
        if age_group == "All ages":
            filtered = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["sex"] == sex)
            ].copy()
        else:
            filtered = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["age"] == age_group) &
                (ihme_df["sex"] == sex)
            ].copy()
    
    # Apply region filter
    if regions:
        filtered = filtered[filtered["WHO_Regions"].isin(regions)]
    
    if filtered.empty:
        return _empty_plot(f"No {measure} data for {virus}")
    
    # Get latest year
    if filtered["year"].notna().any():
        latest_year = filtered["year"].max()
        latest_data = filtered[filtered["year"] == latest_year]
    else:
        return _empty_plot("No valid year data")
    
    if latest_data.empty:
        return _empty_plot(f"No data for year {latest_year}")
    
    # Aggregate by country
    country_data = latest_data.groupby("Country_standard")["val"].sum().reset_index()
    country_data = country_data.sort_values("val", ascending=False).head(top_n)
    
    # Create horizontal bar chart
    fig = px.bar(
        country_data,
        y="Country_standard",
        x="val",
        orientation="h",
        color="val",
        color_continuous_scale="Reds",
        title=f"Top {top_n} Countries by {measure} ({latest_year}, {age_group}, {sex})"
    )
    
    fig.update_traces(
        hovertemplate="<b>%{y}</b><br>%{x:,.2f}<extra></extra>"
    )
    
    fig.update_layout(
        height=400,
        xaxis_title=f"{measure} ({metric_type})",
        yaxis_title="Country",
        yaxis={'categoryorder': 'total ascending'},
        coloraxis_showscale=False
    )
    
    return fig


@callback(
    Output("age-dist-year", "options"),
    Output("age-dist-year", "value"),
    Input("selected-virus", "data"),
    Input("epi-burden-metric", "value"),
    Input("epi-sex-filter", "value"),
    Input("continent-dropdown", "value"),
)
def update_age_dist_year_options(virus, metric, sex, regions):
    data = get_data_store()
    ihme_df = data["ihme_df"]
    
    if ihme_df.empty:
        return [], None
    
    # Parse metric
    try:
        measure, metric_type = metric.split("|")
    except:
        measure, metric_type = "Prevalence", "Number"
    
    # Filter data
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause = cause_lookup.get((virus or "HBV").upper())
    
    # Get available years
    filtered = ihme_df[
        (ihme_df["cause"] == cause) &
        (ihme_df["measure"] == measure) &
        (ihme_df["metric"] == metric_type)
    ]
    
    # Apply region filter
    if regions:
        filtered = filtered[filtered["WHO_Regions"].isin(regions)]
    
    if filtered.empty:
        return [], None
    
    # Get available years
    years = sorted(filtered["year"].unique(), reverse=True)
    options = [{"label": str(year), "value": year} for year in years]
    
    return options, years[0] if years else None


@callback(
    Output("age-distribution-chart", "figure"),
    Input("selected-virus", "data"),
    Input("epi-burden-metric", "value"),
    Input("age-dist-year", "value"),
    Input("age-dist-sex", "value"),  # This should be "Both", "Male", or "Female"
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
)
def update_age_distribution_chart(virus, metric, year, sex, regions, countries):
    data = get_data_store()
    ihme_df = data["ihme_df"]
    
    if ihme_df.empty or not year:
        return _empty_plot("No data available")
    
    # Parse metric
    try:
        measure, metric_type = metric.split("|")
    except:
        measure, metric_type = "Prevalence", "Number"
    
    # Filter data
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause = cause_lookup.get((virus or "HBV").upper())
    
    # CORRECTED: Handle "Both" sexes
    if sex == "Both":
        # Get Male and Female data separately
        male_data = ihme_df[
            (ihme_df["cause"] == cause) &
            (ihme_df["measure"] == measure) &
            (ihme_df["metric"] == metric_type) &
            (ihme_df["year"] == year) &
            (ihme_df["sex"] == "Male")
        ].copy()
        
        female_data = ihme_df[
            (ihme_df["cause"] == cause) &
            (ihme_df["measure"] == measure) &
            (ihme_df["metric"] == metric_type) &
            (ihme_df["year"] == year) &
            (ihme_df["sex"] == "Female")
        ].copy()
        
        # Combine Male and Female data
        filtered = pd.concat([male_data, female_data])
    else:
        # Use the sex as-is
        filtered = ihme_df[
            (ihme_df["cause"] == cause) &
            (ihme_df["measure"] == measure) &
            (ihme_df["metric"] == metric_type) &
            (ihme_df["year"] == year) &
            (ihme_df["sex"] == sex)
        ].copy()
    
    # Exclude "All ages" for age distribution chart
    filtered = filtered[filtered["age"] != "All ages"]
    
    # Apply filters
    if regions:
        filtered = filtered[filtered["WHO_Regions"].isin(regions)]
    if countries:
        filtered = filtered[filtered["Country_standard"].isin(countries)]
    
    if filtered.empty:
        return _empty_plot(f"No age-specific data for {year} ({sex})")
    
    # Aggregate data by age - this will sum when sex="Both"
    age_data = filtered.groupby("age")["val"].sum().reset_index()
    
    # Sort age groups logically
    def sort_age_key(age):
        # Custom sorting logic for age groups
        if not isinstance(age, str):
            return float('inf')
        
        age_lower = age.lower()
        
        # Handle months
        if 'month' in age_lower:
            nums = re.findall(r'\d+', age)
            if nums:
                return int(nums[0]) - 1000
        
        # Handle years
        elif 'year' in age_lower:
            if '<' in age:
                nums = re.findall(r'\d+', age)
                if nums:
                    return int(nums[0]) - 500
            elif '-' in age:
                nums = re.findall(r'\d+', age)
                if nums:
                    return int(nums[0])
        
        return float('inf')
    
    try:
        age_data = age_data.sort_values("age", key=lambda x: x.map(sort_age_key))
    except:
        age_data = age_data.sort_values("age")
    
    # Create figure
    fig = px.bar(
        age_data,
        x="age",
        y="val",
        title=f"Age Distribution of {measure} ({virus}, {year}, {sex})"
    )
    
    fig.update_traces(
        hovertemplate="<b>%{x}</b><br>%{y:,.2f}<extra></extra>"
    )
    
    fig.update_layout(
        height=400,
        xaxis_title="Age Group",
        yaxis_title=f"{measure} ({metric_type})",
        xaxis_tickangle=-45
    )
    
    return fig


@callback(
    Output("sex-ratio-chart", "figure"),
    Input("selected-virus", "data"),
    Input("epi-burden-metric", "value"),
    Input("sex-ratio-age", "value"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
)
def update_sex_ratio_chart(virus, metric, age_group, regions, countries):
    data = get_data_store()
    ihme_df = data["ihme_df"]
    
    if ihme_df.empty:
        return _empty_plot("No IHME data available")
    
    # Parse metric
    try:
        measure, metric_type = metric.split("|")
    except:
        measure, metric_type = "Prevalence", "Number"
    
    # Filter data
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause = cause_lookup.get((virus or "HBV").upper())
    
    # Handle different age group selections
    if age_group == "All ages":
        # For "All ages", include all age groups
        all_data = ihme_df[
            (ihme_df["cause"] == cause) &
            (ihme_df["measure"] == measure) &
            (ihme_df["metric"] == metric_type) &
            (ihme_df["sex"].isin(["Male", "Female"]))
        ].copy()
    else:
        # For specific age groups, use exact match
        all_data = ihme_df[
            (ihme_df["cause"] == cause) &
            (ihme_df["measure"] == measure) &
            (ihme_df["metric"] == metric_type) &
            (ihme_df["age"] == age_group) &
            (ihme_df["sex"].isin(["Male", "Female"]))
        ].copy()
        
        # If no data found, try other metrics as fallback
        if all_data.empty:
            # Try with just Number metric (most common)
            all_data = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == "Number") &
                (ihme_df["age"] == age_group) &
                (ihme_df["sex"].isin(["Male", "Female"]))
            ].copy()
            
            # If still no data, try Prevalence if we weren't already using it
            if all_data.empty and measure != "Prevalence":
                all_data = ihme_df[
                    (ihme_df["cause"] == cause) &
                    (ihme_df["measure"] == "Prevalence") &
                    (ihme_df["metric"] == "Number") &
                    (ihme_df["age"] == age_group) &
                    (ihme_df["sex"].isin(["Male", "Female"]))
                ].copy()
    
    # Apply filters
    if regions:
        all_data = all_data[all_data["WHO_Regions"].isin(regions)]
    
    if countries:
        all_data = all_data[all_data["Country_standard"].isin(countries)]
    
    if all_data.empty:
        return _empty_plot(f"No data available for {age_group} age group")
    
    # Separate Male and Female data
    male_data = all_data[all_data["sex"] == "Male"].copy()
    female_data = all_data[all_data["sex"] == "Female"].copy()
    
    # Aggregate by year
    male_yearly = male_data.groupby("year")["val"].sum().reset_index()
    female_yearly = female_data.groupby("year")["val"].sum().reset_index()
    
    # Create complete year range
    all_years = sorted(set(male_yearly["year"]).union(set(female_yearly["year"])))
    
    # Create dataframes with all years
    male_all_years = pd.DataFrame({"year": all_years})
    female_all_years = pd.DataFrame({"year": all_years})
    
    # Merge with actual data
    male_all_years = male_all_years.merge(male_yearly, on="year", how="left")
    female_all_years = female_all_years.merge(female_yearly, on="year", how="left")
    
    # Fill missing values with interpolation
    male_all_years["val"] = male_all_years["val"].interpolate(method='linear')
    female_all_years["val"] = female_all_years["val"].interpolate(method='linear')
    
    # For edge years, use forward/backward fill
    male_all_years["val"] = male_all_years["val"].fillna(method='ffill').fillna(method='bfill')
    female_all_years["val"] = female_all_years["val"].fillna(method='ffill').fillna(method='bfill')
    
    # Merge and calculate ratio
    ratio_data = pd.merge(male_all_years, female_all_years, on="year", suffixes=("_male", "_female"))
    
    # Calculate ratio safely
    ratio_data["sex_ratio"] = np.where(
        ratio_data["val_female"] > 0,
        ratio_data["val_male"] / ratio_data["val_female"],
        np.nan
    )
    
    ratio_data = ratio_data.dropna(subset=["sex_ratio"])
    
    if ratio_data.empty:
        return _empty_plot("No valid ratio data")
    
    # Create figure
    fig = px.line(
        ratio_data,
        x="year",
        y="sex_ratio",
        markers=True,
        title=f"Male-to-Female Ratio of {measure} ({virus}, {age_group})"
    )
    
    fig.add_hline(
        y=1.0,
        line_dash="dash",
        line_color="gray",
        annotation_text="Equal (1:1)",
        annotation_position="bottom right"
    )
    
    fig.update_traces(
        line=dict(width=3),
        hovertemplate="Year: %{x}<br>Male:Female Ratio: %{y:.2f}<extra></extra>"
    )
    
    # Add note about interpolation if needed
    missing_male = len(male_yearly) < len(all_years)
    missing_female = len(female_yearly) < len(all_years)
    
    if missing_male or missing_female:
        fig.add_annotation(
            text="Note: Some years interpolated (alternating Male/Female data)",
            xref="paper", yref="paper",
            x=0.02, y=0.02,
            showarrow=False,
            font=dict(size=10, color="gray"),
            bgcolor="rgba(255,255,255,0.8)"
        )
    
    fig.update_layout(
        height=400,
        xaxis_title="Year",
        yaxis_title="Male-to-Female Ratio",
        hovermode="x unified"
    )
    
    return fig

@callback(
    Output("region-age-pattern", "options"),
    Input("selected-virus", "data"),
)
def update_region_age_pattern_options(virus):
    data = get_data_store()
    ihme_df = data["ihme_df"]
    
    if ihme_df.empty:
        return []
    
    # Get unique regions
    regions = sorted(ihme_df["WHO_Regions"].dropna().unique())
    options = [{"label": region, "value": region} for region in regions]
    
    return options


@callback(
    Output("region-age-pattern-chart", "figure"),
    Input("selected-virus", "data"),
    Input("epi-burden-metric", "value"),
    Input("region-age-pattern", "value"),
    Input("epi-sex-filter", "value"),
)
def update_region_age_pattern_chart(virus, metric, region, sex):
    data = get_data_store()
    ihme_df = data["ihme_df"]
    
    if ihme_df.empty or not region:
        return _empty_plot("Select a region")
    
    # Parse metric
    try:
        measure, metric_type = metric.split("|")
    except:
        measure, metric_type = "Prevalence", "Number"
    
    # Filter data
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause = cause_lookup.get((virus or "HBV").upper())
    
    # Handle "Both" sexes by summing Male and Female
    if sex == "Both":
        # Get Male data
        filtered_male = ihme_df[
            (ihme_df["cause"] == cause) &
            (ihme_df["measure"] == measure) &
            (ihme_df["metric"] == metric_type) &
            (ihme_df["WHO_Regions"] == region) &
            (ihme_df["sex"] == "Male")
        ].copy()
        
        # Get Female data
        filtered_female = ihme_df[
            (ihme_df["cause"] == cause) &
            (ihme_df["measure"] == measure) &
            (ihme_df["metric"] == metric_type) &
            (ihme_df["WHO_Regions"] == region) &
            (ihme_df["sex"] == "Female")
        ].copy()
        
        # Combine Male and Female data
        filtered = pd.concat([filtered_male, filtered_female])
    else:
        filtered = ihme_df[
            (ihme_df["cause"] == cause) &
            (ihme_df["measure"] == measure) &
            (ihme_df["metric"] == metric_type) &
            (ihme_df["WHO_Regions"] == region) &
            (ihme_df["sex"] == sex)
        ].copy()
    
    # Exclude "All ages"
    filtered = filtered[filtered["age"] != "All ages"]
    
    if filtered.empty:
        return _empty_plot(f"No age-specific data for {region}")
    
    # Get latest year
    latest_year = filtered["year"].max()
    latest_data = filtered[filtered["year"] == latest_year]
    
    if latest_data.empty:
        return _empty_plot(f"No data for {latest_year}")
    
    # Aggregate by age
    age_data = latest_data.groupby("age")["val"].mean().reset_index()
    age_data = age_data.sort_values("age")
    
    # Create figure
    fig = px.bar(
        age_data,
        x="age",
        y="val",
        title=f"Age Pattern in {region} ({virus}, {latest_year}, {sex})"
    )
    
    fig.update_traces(
        hovertemplate="<b>%{x}</b><br>%{y:,.2f}<extra></extra>"
    )
    
    fig.update_layout(
        height=400,
        xaxis_title="Age Group",
        yaxis_title=f"{measure} ({metric_type})",
        xaxis_tickangle=-45
    )
    
    return fig


@callback(
    Output("region-comparison-chart", "figure"),
    Input("selected-virus", "data"),
    Input("epi-burden-metric", "value"),
    Input("epi-age-group", "value"),
    Input("epi-sex-filter", "value"),
    Input("region-compare-type", "value"),
)
def update_region_comparison(virus, metric, age_group, sex, compare_by):
    data = get_data_store()
    ihme_df = data["ihme_df"]
    
    if ihme_df.empty:
        return _empty_plot("No IHME data available")
    
    # Parse metric
    try:
        measure, metric_type = metric.split("|")
    except:
        measure, metric_type = "Prevalence", "Number"
    
    # Filter data
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause = cause_lookup.get((virus or "HBV").upper())
    
    # Handle "Both" sexes by summing Male and Female
    if sex == "Both":
        # Get Male data
        if age_group == "All ages":
            filtered_male = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["sex"] == "Male")
            ].copy()
        else:
            filtered_male = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["age"] == age_group) &
                (ihme_df["sex"] == "Male")
            ].copy()
        
        # Get Female data
        if age_group == "All ages":
            filtered_female = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["sex"] == "Female")
            ].copy()
        else:
            filtered_female = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["age"] == age_group) &
                (ihme_df["sex"] == "Female")
            ].copy()
        
        # Combine Male and Female data
        filtered = pd.concat([filtered_male, filtered_female])
    else:
        # Use the sex as-is
        if age_group == "All ages":
            filtered = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["sex"] == sex)
            ].copy()
        else:
            filtered = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["age"] == age_group) &
                (ihme_df["sex"] == sex)
            ].copy()
    
    if filtered.empty:
        return _empty_plot("No data for selected filters")
    
    # Get latest year
    latest_year = filtered["year"].max()
    latest_data = filtered[filtered["year"] == latest_year]
    
    # Simple aggregation by region
    region_data = latest_data.groupby("WHO_Regions")["val"].sum().reset_index()
    region_data = region_data.sort_values("val", ascending=False)
    
    # Create bar chart
    fig = px.bar(
        region_data,
        x="WHO_Regions",
        y="val",
        color="WHO_Regions",
        title=f"{measure} by WHO Region ({latest_year})"
    )
    
    fig.update_traces(
        hovertemplate="<b>%{x}</b><br>Value: %{y:,.0f}<extra></extra>"
    )
    
    fig.update_layout(
        height=400,
        xaxis_title="WHO Region",
        yaxis_title=f"{measure} ({metric_type})",
        showlegend=False,
        xaxis_tickangle=-45
    )
    
    return fig


@callback(
    Output("burden-coverage-scatter", "figure"),
    Input("selected-virus", "data"),
    Input("correlation-metric", "value"),
    Input("correlation-age", "value"),
    Input("correlation-sex", "value"),
    Input("continent-dropdown", "value"),
)
def update_burden_coverage_scatter(virus, metric, age_group, sex, regions):
    data = get_data_store()
    ihme_df = data["ihme_df"]
    seq_data = data["hbv_data"] if (virus or "HBV") == "HBV" else data["hcv_data"]
    population_df = data["population_df"]
    
    if ihme_df.empty or seq_data.empty:
        return _empty_plot("Insufficient data for correlation analysis")
    
    # Parse x metric
    try:
        measure, metric_type = metric.split("|")
    except:
        measure, metric_type = "Prevalence", "Number"
    
    # Get burden data
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause = cause_lookup.get((virus or "HBV").upper())
    
    # Handle "Both" sexes by summing Male and Female
    if sex == "Both":
        # Get Male data
        if age_group == "All ages":
            male_data = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["sex"] == "Male")
            ].copy()
        else:
            male_data = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["age"] == age_group) &
                (ihme_df["sex"] == "Male")
            ].copy()
        
        # Get Female data
        if age_group == "All ages":
            female_data = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["sex"] == "Female")
            ].copy()
        else:
            female_data = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["age"] == age_group) &
                (ihme_df["sex"] == "Female")
            ].copy()
        
        # Combine Male and Female data
        burden_data = pd.concat([male_data, female_data])
    else:
        # Use the sex as-is
        if age_group == "All ages":
            burden_data = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["sex"] == sex)
            ].copy()
        else:
            burden_data = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["age"] == age_group) &
                (ihme_df["sex"] == sex)
            ].copy()
    
    if burden_data.empty:
        return _empty_plot(f"No burden data for {virus} ({age_group}, {sex})")
    
    # Get latest year
    latest_year = burden_data["year"].max()
    burden_latest = burden_data[burden_data["year"] == latest_year]
    
    if burden_latest.empty:
        return _empty_plot(f"No burden data for year {latest_year}")
    
    # Aggregate burden by country
    country_burden = burden_latest.groupby("Country_standard")["val"].sum().reset_index()
    
    # Get sequence counts by country
    country_sequences = seq_data.groupby("Country_standard").size().reset_index(name="sequence_count")
    
    # Merge data
    merged = pd.merge(country_burden, country_sequences, on="Country_standard", how="inner")
    
    # Apply region filter
    if regions:
        merged = merged[merged["Country_standard"].isin(
            ihme_df[ihme_df["WHO_Regions"].isin(regions)]["Country_standard"].unique()
        )]
    
    if merged.empty:
        return _empty_plot("No overlapping data between burden and sequences")
    
    # Add population data for per-capita calculations (optional)
    if not population_df.empty:
        pop_latest = population_df[population_df["Year"] == latest_year]
        if not pop_latest.empty:
            merged = merged.merge(pop_latest[["Country_standard", "Population"]], 
                                on="Country_standard", how="left")
            merged["sequences_per_million"] = (merged["sequence_count"] / merged["Population"]) * 1_000_000
            # You could use sequences_per_million instead of sequence_count for y-axis
    
    # Create scatter plot
    fig = px.scatter(
        merged,
        x="val",
        y="sequence_count",
        size="sequence_count",
        hover_name="Country_standard",
        log_x=True,
        log_y=True,
        title=f"Burden vs. Sequences Correlation ({virus}, {latest_year}, {age_group}, {sex})",
        hover_data={"Country_standard": True, "val": ":,.0f", "sequence_count": True}
    )
    
    # Add trend line
    if len(merged) > 1:
        # Calculate linear regression for trend line
        x = np.log(merged["val"] + 1)  # Add 1 to avoid log(0)
        y = np.log(merged["sequence_count"] + 1)
        
        # Fit line
        coefficients = np.polyfit(x, y, 1)
        polynomial = np.poly1d(coefficients)
        
        # Generate points for trend line
        x_line = np.linspace(x.min(), x.max(), 100)
        y_line = polynomial(x_line)
        
        # Add trend line to plot
        fig.add_trace(go.Scatter(
            x=np.exp(x_line),
            y=np.exp(y_line),
            mode='lines',
            name='Trend Line',
            line=dict(color='red', dash='dash'),
            hovertemplate='Trend Line<extra></extra>'
        ))
    
    fig.update_layout(
        height=400,
        xaxis_title=f"{measure} ({metric_type})",
        yaxis_title="Number of Sequences",
        hovermode="closest",
        showlegend=True
    )
    
    return fig


@callback(
    Output("epidemiology-data-table", "children"),
    Input("selected-virus", "data"),
    Input("epi-burden-metric", "value"),
    Input("epi-age-group", "value"),
    Input("epi-sex-filter", "value"),
    Input("continent-dropdown", "value"),
    Input("country-dropdown", "value"),
)
def update_epidemiology_table(virus, metric, age_group, sex, regions, countries):
    data = get_data_store()
    ihme_df = data["ihme_df"]
    
    if ihme_df.empty:
        return html.P("No epidemiology data available", className="text-muted")
    
    # Parse metric
    try:
        measure, metric_type = metric.split("|")
    except:
        measure, metric_type = "Prevalence", "Number"
    
    # Filter data
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause = cause_lookup.get((virus or "HBV").upper())
    
    # Handle "Both" sexes by summing Male and Female
    if sex == "Both":
        # Get Male data
        if age_group == "All ages":
            filtered_male = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["sex"] == "Male")
            ].copy()
        else:
            filtered_male = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["age"] == age_group) &
                (ihme_df["sex"] == "Male")
            ].copy()
        
        # Get Female data
        if age_group == "All ages":
            filtered_female = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["sex"] == "Female")
            ].copy()
        else:
            filtered_female = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["age"] == age_group) &
                (ihme_df["sex"] == "Female")
            ].copy()
        
        # Combine Male and Female data
        filtered = pd.concat([filtered_male, filtered_female])
    else:
        # Use the sex as-is
        if age_group == "All ages":
            filtered = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["sex"] == sex)
            ].copy()
        else:
            filtered = ihme_df[
                (ihme_df["cause"] == cause) &
                (ihme_df["measure"] == measure) &
                (ihme_df["metric"] == metric_type) &
                (ihme_df["age"] == age_group) &
                (ihme_df["sex"] == sex)
            ].copy()
    
    # Apply filters
    if regions:
        filtered = filtered[filtered["WHO_Regions"].isin(regions)]
    if countries:
        filtered = filtered[filtered["Country_standard"].isin(countries)]
    
    if filtered.empty:
        return html.P("No data for selected filters", className="text-muted")
    
    # Select and rename columns
    display_cols = ["Country_standard", "WHO_Regions", "year", "sex", "age", "val"]
    display_cols = [col for col in display_cols if col in filtered.columns]
    
    display_df = filtered[display_cols].copy()
    display_df = display_df.sort_values(["year", "Country_standard"])
    
    # Create DataTable
    table = dash.dash_table.DataTable(
        data=display_df.to_dict('records'),
        columns=[{"name": col.replace("_", " ").title(), "id": col} for col in display_cols],
        page_size=10,
        style_table={'overflowX': 'auto'},
        style_cell={
            'textAlign': 'left',
            'padding': '8px',
            'overflow': 'hidden',
            'textOverflow': 'ellipsis',
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
    Output("epi-age-group", "options"),
    Output("epi-age-group", "value"),
    Output("sex-ratio-age", "options"),
    Output("sex-ratio-age", "value"),
    Input("selected-virus", "data"),
    Input("epi-burden-metric", "value"),
)
def update_age_group_options(virus, metric):
    data = get_data_store()
    ihme_df = data["ihme_df"]
    
    if ihme_df.empty:
        options = [{"label": "All ages", "value": "All ages"}]
        return options, "All ages", options, "All ages"
    
    # Parse metric
    try:
        measure, metric_type = metric.split("|")
    except:
        measure, metric_type = "Prevalence", "Number"
    
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E",
    }
    cause = cause_lookup.get((virus or "HBV").upper())
    
    # First try with the selected metric
    filtered = ihme_df[
        (ihme_df["cause"] == cause) &
        (ihme_df["measure"] == measure) &
        (ihme_df["metric"] == metric_type)
    ]
    
    # If no data with that metric, try with just the measure and any metric
    if filtered.empty:
        filtered = ihme_df[
            (ihme_df["cause"] == cause) &
            (ihme_df["measure"] == measure)
        ]
    
    # If still no data, try with any measure for this cause
    if filtered.empty:
        filtered = ihme_df[ihme_df["cause"] == cause]
    
    # Get unique age groups
    age_groups = sorted(filtered["age"].dropna().unique())
    
    # Try to sort age groups logically
    def sort_age_key(age):
        if not isinstance(age, str):
            return float('inf')
        
        age_lower = age.lower()
        
        # Handle months first
        if 'month' in age_lower:
            nums = re.findall(r'\d+', age)
            if nums:
                return int(nums[0]) - 1000
        
        # Handle years
        elif 'year' in age_lower:
            if '<' in age:
                nums = re.findall(r'\d+', age)
                if nums:
                    return int(nums[0]) - 500
            elif '-' in age:
                nums = re.findall(r'\d+', age)
                if nums:
                    return int(nums[0])
        
        return float('inf')
    
    try:
        age_groups.sort(key=sort_age_key)
    except:
        age_groups.sort()
    
    # Create options - include "All ages" first
    options = [{"label": "All ages (sum of all age groups)", "value": "All ages"}]
    
    # Add actual age groups
    for age in age_groups:
        if age != "All ages":  # Don't add if it's already there
            options.append({"label": age, "value": age})
    
    # Default value
    default_value = "All ages"
    
    return options, default_value, options, default_value

# Priority Download Callback
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
    Output("year-slider", "min"),
    Output("year-slider", "max"),
    Output("year-slider", "value"),
    Output("year-slider", "marks"),
    Output('continent-dropdown', 'options'),
    Output('country-dropdown', 'options'),
    Output('genotype-dropdown', 'options'),
    Input("selected-virus", "data"),
)
def init_controls(virus):
    data = get_data_store()
    
    # FIXED: Properly handle all three viruses
    if virus == "HBV":
        base = data['hbv_data']
    elif virus == "HCV":
        base = data['hcv_data']
    elif virus == "HEV":
        base = data['hev_data']
    else:
        base = data['hbv_data']
    
    if base.empty:
        return 2000, 2023, [2000, 2023], {}, [], [], []
    
    y0, y1 = int(base["Year"].min()), int(base["Year"].max())
    marks = {y: (str(y) if y % 5 == 0 else "") for y in range(y0, y1 + 1)}
    cont_opts = [{"label": r, "value": r} for r in sorted(base["WHO_Regions"].dropna().unique()) if r!="Unknown"]
    country_opts = [{"label": c, "value": c} for c in sorted(base["Country_standard"].dropna().unique()) if c!="Unknown"]
    geno_opts = [{"label": g, "value": g} for g in sorted(base["genotype"].dropna().unique())]
    return y0, y1, [y0, y1], marks, cont_opts, country_opts, geno_opts

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

@callback(
    Output("btn-go-to-epidemiology", "n_clicks"),  # Add this output
    Input("btn-go-to-epidemiology", "n_clicks"),
    prevent_initial_call=True
)
def go_to_epidemiology_tab(n_clicks):
    # This triggers the tab switch via the existing tab navigation callback
    return n_clicks

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
    State("year-slider", "value"),
    State("continent-dropdown", "value"),
    State("country-dropdown", "value"),
    State("genotype-dropdown", "value"),
    prevent_initial_call=True,
)
def download_main_data_with_taxa(n_clicks, filtered_json, virus, years, regions, countries, genotypes):
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
    Output("btn-hbv", "outline"),
    Output("btn-hcv", "outline"),
    Output("btn-hev", "outline"),
    Input("btn-hbv", "n_clicks"),
    Input("btn-hcv", "n_clicks"),
    Input("btn-hev", "n_clicks"),
    prevent_initial_call=True
)

def update_virus(btn_hbv_clicks, btn_hcv_clicks, btn_hev_clicks):
    ctx = callback_context
    if not ctx.triggered:
        return dash.no_update, dash.no_update, dash.no_update
    
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if triggered_id == "btn-hbv":
        return "HBV", False, True, True
    elif triggered_id == "btn-hcv":
        return "HCV", True, False, True
    elif triggered_id == "btn-hev":
        return "HEV", True, True, False
    
    return "HBV", False, True, True 

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

        
# - Year reset -button
@callback(
    Output("year-slider", "value", allow_duplicate=True),
    Input("btn-reset-time", "n_clicks"),
    State("selected-virus", "data"),
    prevent_initial_call=True
)
def reset_time_range(n_clicks, selected_virus):
    if n_clicks is None or n_clicks == 0:
        raise PreventUpdate
        
    data = get_data_store()
    df = data['hbv_data'] if selected_virus == "HBV" else data['hcv_data']
    min_year = int(df["Year"].min())
    max_year = int(df["Year"].max())
    
    return [min_year, max_year]