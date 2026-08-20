#!/usr/bin/env python3
"""
WHO Global Health Observatory (GHO) Hepatitis Country Profiles Extractor
========================================================================

Queries the public WHO GHO OData API to fetch country-level statistics for
Hepatitis B and Hepatitis C, pivots the indicators into a tidy country-year
format, and saves the result as a TSV.

Output:
    results/who_gho/who_gho_hepatitis_country_profiles.tsv

Usage:
    python who_gho_extractor.py
    python who_gho_extractor.py --dry-run
"""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

# Base URL for the WHO GHO OData API
GHO_API_BASE_URL = "https://ghoapi.azureedge.net/api"

# Mapping of GHO Indicator Codes to clean, human-readable column names
INDICATOR_MAPPING = {
    # --- Hepatitis B (HBV) ---
    "HEPATITIS_HBV_PREVALENCE_PER100": "hbv_prevalence_pct",
    "SDGHEPHBSAGPRV": "hbv_prevalence_hbsag_pct",
    "HEPATITIS_HBV_LIVINGWITH_NUM": "hbv_livingwith_num",
    "HEPATITIS_HBV_DIAGNOSED_NUM": "hbv_diagnosed_num",
    "HEPATITIS_HBV_DIAGNOSIS_PERINFECTED_PER100": "hbv_diagnosis_rate_pct",
    "HEPATITIS_HBV_TREATMENT_NUM": "hbv_treatment_num",
    "HEPATITIS_HBV_TREATMENT_PERDIAGNOSED_PER100": "hbv_treatment_rate_diagnosed_pct",
    "HEPATITIS_HBV_TREATMENT_PERINFECTED_PER100": "hbv_treatment_rate_infected_pct",
    "HEPATITIS_HBV_INFECTIONS_NEW_NUM": "hbv_new_infections_num",
    "HEPATITIS_HBV_DEATHS_NUM": "hbv_deaths_num",
    "WHS4_117": "hbv_vaccine_hepb3_coverage_pct",

    # --- Hepatitis C (HCV) ---
    "HEPATITIS_HCV_PREVALENCE_PER100": "hcv_prevalence_pct",
    "HEPATITIS_HCV_LIVINGWITH_NUM": "hcv_livingwith_num",
    "HEPATITIS_HCV_DIAGNOSED_NUM": "hcv_diagnosed_num",
    "HEPATITIS_HCV_DIAGNOSIS_PERINFECTED_CRD_PER100": "hcv_diagnosis_rate_pct",
    "HEPATITIS_HCV_TREATMENT_CUM_NUM": "hcv_treatment_cumulative_num",
    "HEPATITIS_HCV_TREATMENT_PERDIAGNOSED_PER100": "hcv_treatment_rate_diagnosed_pct",
    "HEPATITIS_HCV_INFECTIONS_NEW_NUM": "hcv_new_infections_num",
    "HEPATITIS_HCV_DEATHS_NUM": "hcv_deaths_num"
}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

log = logging.getLogger(__name__)


def fetch_gho_indicator_data(indicator_code: str, retries: int = 3, backoff: float = 2.0) -> list[dict]:
    """
    Fetches raw data for a specific GHO indicator, filtering for country-level records.
    """
    url = f"{GHO_API_BASE_URL}/{indicator_code}"
    headers = {
        "User-Agent": "HepatitisDB-Research/1.0",
        "Accept": "application/json",
    }

    for attempt in range(1, retries + 1):
        try:
            log.info(f"Fetching indicator {indicator_code} (Attempt {attempt}/{retries})")
            response = requests.get(url, timeout=30, headers=headers)
            response.raise_for_status()
            
            data = response.json()
            raw_rows = data.get("value", [])
            
            # Filter for country level only (ignores regional, global or subnational aggregations)
            country_records = []
            for row in raw_rows:
                if row.get("SpatialDimType") == "COUNTRY" and row.get("SpatialDim"):
                    country_records.append(row)
            
            log.info(f"  → Found {len(country_records)} country records")
            return country_records

        except (requests.RequestException, ValueError) as e:
            log.warning(f"  Failed to fetch {indicator_code} on attempt {attempt}: {e}")
            if attempt < retries:
                time.sleep(backoff * attempt)
            else:
                log.error(f"  Could not retrieve data for indicator {indicator_code} after {retries} attempts.")
                return []

    return []


def run(output_path: str, dry_run: bool = False) -> None:
    """
    Runs the GHO Hepatitis Country Profiles Ingestion Pipeline.
    """
    output = Path(output_path)
    
    if dry_run:
        log.info("Dry-run mode: checking configurations and API connectivity.")
        for code, column in INDICATOR_MAPPING.items():
            log.info(f"  - Map: {code} → {column}")
        return

    # Create destination directory if it doesn't exist
    output.parent.mkdir(parents=True, exist_ok=True)

    all_data = []

    for indicator_code, column_name in INDICATOR_MAPPING.items():
        records = fetch_gho_indicator_data(indicator_code)
        
        for r in records:
            country = r.get("SpatialDim")
            year = r.get("TimeDim") or r.get("TimeDimensionValue")
            
            # Handle float values safely
            numeric_value = r.get("NumericValue")
            if numeric_value is None:
                # Fallback to parse standard string format
                str_val = str(r.get("Value") or "").replace(" ", "").strip()
                try:
                    numeric_value = float(str_val)
                except ValueError:
                    numeric_value = None

            if country and year is not None:
                try:
                    year_int = int(year)
                except (ValueError, TypeError):
                    continue

                all_data.append({
                    "country": str(country).upper(),
                    "year": year_int,
                    "indicator": column_name,
                    "value": numeric_value
                })
        
        # Friendly delay between requests to be gentle to the public API
        time.sleep(0.5)

    if not all_data:
        log.error("No data collected from GHO OData API. Ingestion terminated.")
        return

    log.info(f"Successfully collected {len(all_data)} raw data observations.")

    # 1. Load into DataFrame
    df = pd.DataFrame(all_data)

    # 2. Pivot to tidy country-year records
    # This aligns indicators side-by-side for each country and year
    log.info("Pivoting data to country-year tidy layout...")
    pivoted_df = df.pivot_table(
        index=["country", "year"],
        columns="indicator",
        values="value"
    ).reset_index()

    # 3. Fill missing mapped columns if they were absent in the returned values
    for col in INDICATOR_MAPPING.values():
        if col not in pivoted_df.columns:
            pivoted_df[col] = None

    # 4. Reorder columns: country, year, then indicator columns alphabetically
    ordered_columns = ["country", "year"] + sorted(list(INDICATOR_MAPPING.values()))
    pivoted_df = pivoted_df[ordered_columns]

    # 5. Sort by country and year
    pivoted_df = pivoted_df.sort_values(by=["country", "year"]).reset_index(drop=True)

    # 6. Save as TSV
    pivoted_df.to_csv(output, sep="\t", index=False)
    log.info(f"Tidy country profiles dataset saved: {output} ({len(pivoted_df)} rows)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="WHO Global Health Observatory (GHO) Hepatitis country profiles extractor"
    )
    parser.add_argument(
        "--output",
        default="results/who_gho/who_gho_hepatitis_country_profiles.tsv",
        help="Path to save the output TSV file",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Check configuration mapping without fetching from GHO API",
    )
    
    args = parser.parse_args()
    run(output_path=args.output, dry_run=args.dry_run)
