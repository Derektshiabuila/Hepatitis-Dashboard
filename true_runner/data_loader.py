import pandas as pd
import numpy as np
import re
from functools import lru_cache
import pycountry
import os

# Configuration
ENCODINGS = ['utf-8-sig', 'utf-8', 'latin1', 'iso-8859-1', 'cp1252', 'windows-1252', 'mac_roman']

SPECIAL_CASES = {
    'USA': 'United States',
    'US': 'United States',
    'CapeVerde': 'Cape Verde',
    'Cabo Verde': 'Cape Verde',
    'Democratic Republic of the Congo': 'Congo [DRC]',
    'Congo, Dem. Rep.': 'Congo [DRC]',
    'Congo, Rep.': 'Congo [Republic]',
    'Congo': 'Congo [Republic]',
    'Federal Republic of Nigeria': 'Nigeria',
    'Republic of Niger': 'Niger',
    'Republic of Serbia': 'Serbia',
    'Myanmar': 'Myanmar [Burma]',
    'VietNam': 'Viet Nam',
    'Vietnam': 'Viet Nam',
    'Iran, Islamic Republic of': 'Iran',
    'Iran (Islamic Republic of)': 'Iran',
    'Iran, Islamic Rep.': 'Iran',
    'Iran': 'Iran',
    'Korea, Republic of': 'South Korea',
    'South Korea': 'South Korea',
    'North Korea': 'North Korea',
    "Lao People's Democratic Republic": 'Laos',
    'Lao PDR': 'Laos',
    'Laos': 'Laos',
    'Côte d\'Ivoire': 'Ivory Coast',
    'Cote d\'Ivoire': 'Ivory Coast',
    'Bolivia (Plurinational State of)': 'Bolivia, Plurinational State of',
    'Bolivia': 'Bolivia, Plurinational State of',
    'Venezuela (Bolivarian Republic of)': 'Venezuela, Bolivarian Republic of',
    'Venezuela': 'Venezuela, Bolivarian Republic of',
    'Tanzania, United Republic of': 'Tanzania',
    'United Republic of Tanzania': 'Tanzania',
    'Tanzania': 'Tanzania',
    'Swaziland': 'Eswatini',
    'Turkey': 'Turkey',
    'St. Lucia': 'Saint Lucia',
    'St. Kitts and Nevis': 'Saint Kitts and Nevis',
    'St. Vincent and the Grenadines': 'Saint Vincent and the Grenadines',
    'Hong Kong': 'China',
    'Hong-Kong': 'China',
    'HongKong': 'China',
    'Hong Kong SAR': 'China',
    'Hong Kong SAR, China': 'China',
    'Hong Kong, SAR China': 'China',
    'Hong Kong Special Administrative Region': 'China',
    'Hong Kong S.A.R.': 'China',
    'HKSAR': 'China',
    'HK': 'China',     # sometimes appears as a short code in metadata
    'HKG': 'China',    # ISO3 code occasionally found in CSVs
    'Venezuela (Bolivarian Republic of)': 'Venezuela',
    'Democratic People\'s Republic of Korea': 'North Korea',
    'Republic of Korea': 'South Korea',
    'Taiwain, Province of China': 'Taiwan, Province of China',
    'Taiwan (Province of China)': 'Taiwan, Province of China',
    'Taiwan': 'Taiwan, Province of China',
    'Micronesia (Federated States of)': 'Micronesia, Federated States of',
    'Micronesia, Federated States of': 'Micronesia, Federated States of',
    'Micronesia': 'Micronesia, Federated States of',
    'Namibia': 'Namibia',  # This should already map to itself
    'Venezuela, Bolivarian Republic of': 'Venezuela',
    'Venezuela (Bolivarian Republic of)': 'Venezuela',
    'United States Virgin Islands': 'U.S. Virgin Islands',
    'Sao Tome and Principe': 'Sao Tome and Principe',
    'United States Virgin Islands': 'U.S. Virgin Islands',
    'Venezuela': 'Venezuela',
    'North Macedonia': 'Macedonia [FYROM]',
    'Macedonia': 'Macedonia [FYROM]',
    'Republic of North Macedonia': 'Macedonia [FYROM]',
    'Uknown': None,
    'unknown': None,
    'Not provided': None,
    'not provided': None,
}

# === DATA LOADING ============================================================
# Pre-compile regex patterns
ACCESSION_PATTERN = re.compile(r"([A-Z0-9]+\.\d+)")
WHITESPACE_PATTERN = re.compile(r"\s+")

# Get the data directory path
def get_data_path(file_name):
    """Get absolute path to data file"""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(current_dir, '.')
    return os.path.join(data_dir, file_name)

def load_with_encoding(file_name, sep='\t'):
    """Load file with encoding fallback and better error handling"""
    file_path = get_data_path(file_name)
    
    # Check if file exists first
    if not os.path.exists(file_path):
        print(f"⚠️  File not found: {file_path}")
        return pd.DataFrame()  # Return empty DataFrame instead of crashing
    
    for encoding in ENCODINGS:
        try:
            print(f"Trying encoding {encoding} for {file_name}...")
            df = pd.read_csv(file_path, sep=sep, encoding=encoding, 
                           on_bad_lines='skip', low_memory=False)
            print(f"✅ Successfully loaded {file_name} with {encoding} encoding")
            return df
        except Exception as e:
            print(f"Failed with encoding {encoding} for {file_name}: {e}")
            continue
    
    print(f"❌ Could not read {file_path} with any encoding")
    return pd.DataFrame()  # Return empty DataFrame instead of crashing
    
def compute_burden_adjusted_coverage_from_sequences(
    seq_df: pd.DataFrame,
    ihme_df: pd.DataFrame,
    virus: str
) -> pd.DataFrame:
    cause_lookup = {
        "HBV": "Total burden related to hepatitis B",
        "HCV": "Total burden related to hepatitis C",
        "HEV": "Total burden related to hepatitis E"
    }

    cause = cause_lookup[virus]
    
    # Check if we have data
    if seq_df.empty or ihme_df.empty:
        print(f"⚠️  Empty data for {virus} coverage calculation")
        return pd.DataFrame()

    # --- Filter IHME to prevalence numbers - FIXED: Handle Male/Female separately ---
    # Get Male data
    ihme_male = ihme_df[
        (ihme_df["cause"] == cause) &
        (ihme_df["measure"] == "Prevalence") &
        (ihme_df["metric"] == "Number") &
        (ihme_df["sex"] == "Male") &
        (ihme_df["age"] == "All ages")
    ].copy()
    
    # Get Female data
    ihme_female = ihme_df[
        (ihme_df["cause"] == cause) &
        (ihme_df["measure"] == "Prevalence") &
        (ihme_df["metric"] == "Number") &
        (ihme_df["sex"] == "Female") &
        (ihme_df["age"] == "All ages")
    ].copy()
    
    # Combine Male and Female (sum them for total burden)
    if not ihme_male.empty and not ihme_female.empty:
        # Sum Male and Female values per country/year
        ihme_combined = pd.concat([ihme_male, ihme_female])
        ihme_prev = ihme_combined.groupby(["Country_standard", "year"], as_index=False)["val"].sum()
    elif not ihme_male.empty:
        ihme_prev = ihme_male[["Country_standard", "year", "val"]].copy()
    elif not ihme_female.empty:
        ihme_prev = ihme_female[["Country_standard", "year", "val"]].copy()
    else:
        # Try without age filter as fallback
        ihme_male = ihme_df[
            (ihme_df["cause"] == cause) &
            (ihme_df["measure"] == "Prevalence") &
            (ihme_df["metric"] == "Number") &
            (ihme_df["sex"] == "Male")
        ].copy()
        
        ihme_female = ihme_df[
            (ihme_df["cause"] == cause) &
            (ihme_df["measure"] == "Prevalence") &
            (ihme_df["metric"] == "Number") &
            (ihme_df["sex"] == "Female")
        ].copy()
        
        if not ihme_male.empty or not ihme_female.empty:
            ihme_combined = pd.concat([ihme_male, ihme_female])
            if not ihme_combined.empty:
                ihme_prev = ihme_combined.groupby(["Country_standard", "year"], as_index=False)["val"].sum()
            else:
                print(f"⚠️  No prevalence data found for {virus}")
                return pd.DataFrame()
        else:
            print(f"⚠️  No prevalence data found for {virus}")
            return pd.DataFrame()

    if ihme_prev.empty:
        print(f"⚠️  IHME prevalence data empty for {virus}")
        return pd.DataFrame()

    # Use latest IHME year per country
    ihme_latest = (
        ihme_prev.sort_values("year")
        .groupby("Country_standard", as_index=False)
        .tail(1)
        .rename(columns={"val": "Prevalence_abs"})
        [["Country_standard", "Prevalence_abs"]]
    )
    
    print(f"📊 {virus} coverage: {len(ihme_latest)} countries with prevalence data")

    # --- Sequence counts by country & genotype ---
    seq_counts = (
        seq_df
        .dropna(subset=["Country_standard", "genotype"])
        .groupby(["Country_standard", "genotype"])
        .size()
        .reset_index(name="Seq_count")
    )
    
    if seq_counts.empty:
        print(f"⚠️  No sequence counts for {virus}")
        return pd.DataFrame()

    totals = (
        seq_counts
        .groupby("Country_standard")["Seq_count"]
        .sum()
        .reset_index(name="Total_seq")
    )

    seq_counts = seq_counts.merge(totals, on="Country_standard", how="left")
    seq_counts["Genotype_prop"] = seq_counts["Seq_count"] / seq_counts["Total_seq"]

    # --- Merge prevalence ---
    out = seq_counts.merge(ihme_latest, on="Country_standard", how="left")
    
    out = out.dropna(subset=["Prevalence_abs"])
    
    if out.empty:
        print(f"⚠️  No matches between sequence countries and IHME countries for {virus}")
        return pd.DataFrame()

    # --- Estimate infections & coverage ---
    out["Est_infections_genotype"] = out["Prevalence_abs"] * out["Genotype_prop"]
    out["Coverage_ratio"] = np.where(
        out["Est_infections_genotype"] > 0,
        out["Seq_count"] / out["Est_infections_genotype"],
        0
    )

    out["Virus"] = virus

    print(f"✅ {virus} coverage computed: {len(out)} genotype-country combinations")
    
    return out[
        [
            "Country_standard",
            "genotype",
            "Seq_count",
            "Total_seq",
            "Genotype_prop",
            "Prevalence_abs",
            "Est_infections_genotype",
            "Coverage_ratio",
            "Virus",
        ]
    ]

def load_recombinant_ids(file_name):
    file_path = get_data_path(file_name)
    
    # Check if file exists
    if not os.path.exists(file_path):
        print(f"⚠️  Recombinant IDs file not found: {file_path}")
        return set()  # Return empty set
    
    try:
        with open(file_path) as f:
            return set(line.strip() for line in f)
    except Exception as e:
        print(f"❌ Error loading recombinant IDs: {e}")
        return set()

def load_csv_file(file_name, **kwargs):
    """Load CSV file with error handling"""
    file_path = get_data_path(file_name)
    
    if not os.path.exists(file_path):
        print(f"⚠️  File not found: {file_path}")
        return pd.DataFrame()
    
    try:
        return pd.read_csv(file_path, **kwargs)
    except Exception as e:
        print(f"❌ Error loading {file_name}: {e}")
        return pd.DataFrame()

def standardize_country_column(df, column_name):
    """Standardize country column in a dataframe"""
    df['Country_standard'] = df[column_name].apply(standardize_country_cached)
    return df
    
@lru_cache(maxsize=5000)

def standardize_country_cached(name):
    if not isinstance(name, str) or not name.strip():
        return None
    
    name = name.strip()
    
    # First use special case mapping
    if name in SPECIAL_CASES:
        return SPECIAL_CASES[name]
    
    # Then try fuzzy match
    try:
        match = pycountry.countries.search_fuzzy(name)
        if match:
            return match[0].name
    except:
        pass
    
    # If still no match, return original name
    return name

def merge_population_nearest(counts_df: pd.DataFrame, pop_df: pd.DataFrame, tol_years: int = 3) -> pd.DataFrame:
    counts = counts_df.copy()
    pop = pop_df.copy()

    # Precompute standardized country names
    for df in (counts, pop):
        df["Country_standard"] = df["Country_standard"].astype(str)
        df["Country_standard"] = df["Country_standard"].str.strip()
        df["Country_standard"] = df["Country_standard"].str.replace(WHITESPACE_PATTERN, " ", regex=True)
        df["Year"] = pd.to_numeric(df["Year"], errors="coerce")

    pop["Population"] = pd.to_numeric(pop["Population"], errors="coerce")

    counts = counts.dropna(subset=["Country_standard", "Year"])
    pop = pop.dropna(subset=["Country_standard", "Year", "Population"])

    counts["Year"] = counts["Year"].astype("int32")
    pop["Year"] = pop["Year"].astype("int32")

    out = []
    # Use vectorized operations where possible
    for country, g in counts.groupby("Country_standard", sort=False, observed=True):
        gp = pop[pop["Country_standard"] == country]
        if gp.empty:
            g["Population"] = np.nan
            out.append(g)
            continue
        
        # Sort once and use merge_asof
        g_sorted = g.sort_values("Year", kind='stable')
        gp_sorted = gp[["Year", "Population"]].sort_values("Year", kind='stable')
        merged = pd.merge_asof(g_sorted, gp_sorted, on="Year", direction="nearest", tolerance=tol_years)
        out.append(merged)

    return pd.concat(out, ignore_index=True)

def merge_population_nearest_two_pass(counts_df, pop_df, tol_years_first=3, tol_years_wide=50):
    first = merge_population_nearest(counts_df, pop_df, tol_years=tol_years_first)
    if first["Population"].isna().any():
        widened = merge_population_nearest(counts_df, pop_df, tol_years=tol_years_wide)
        first["Population"] = first["Population"].fillna(widened["Population"])
    return first

def merge_who_region(df, coords_df):
    """Optimized WHO region merging"""
    # Create lookup dictionary for faster merging
    region_lookup = coords_df.set_index('Country_standard')['WHO'].rename('WHO_Regions')
    df = df.join(region_lookup, on='Country_standard', how='left')
    df["WHO_Regions"] = df["WHO_Regions"].fillna("Unknown").astype('category')
    return df

def standardize_country_column(df, column_name):
    """Standardize country column in a dataframe"""
    df['Country_standard'] = df[column_name].apply(standardize_country_cached)
    return df

def get_countries_missing_who_region(data_store):
    hbv_data = data_store['hbv_data']
    hcv_data = data_store['hcv_data']
    hev_data = data_store['hev_data']
    ihme_df = data_store['ihme_df']
    
    # Get unique unknown countries from each dataset
    hbv_unknown = set(hbv_data[hbv_data["WHO_Regions"] == "Unknown"]["Country_standard"].unique())
    hcv_unknown = set(hcv_data[hcv_data["WHO_Regions"] == "Unknown"]["Country_standard"].unique())
    hev_unknown = set(hev_data[hev_data["WHO_Regions"] == "Unknown"]["Country_standard"].unique())
    ihme_unknown = set(ihme_df[ihme_df["WHO_Regions"] == "Unknown"]["Country_standard"].unique())
    
    # Combine all
    all_unknown = hbv_unknown | hcv_unknown | hev_unknown | ihme_unknown
    
    return sorted(all_unknown)

def optimize_dataframe_dtypes(df):
    for col in df.columns:
        if df[col].dtype == 'object':
            # Convert to category if low cardinality
            if df[col].nunique() / len(df[col]) < 0.5:
                df[col] = df[col].astype('category')
        elif df[col].dtype in ['int64', 'float64']:
            # Downcast numeric types
            if pd.api.types.is_integer_dtype(df[col]):
                df[col] = pd.to_numeric(df[col], downcast='integer')
            else:
                df[col] = pd.to_numeric(df[col], downcast='float')
    return df
    
def get_mutation_summary_stats(mutation_df, sequence_df, virus_type):
    if mutation_df.empty or sequence_df.empty:
        return {
            'total_samples': 0,
            'total_mutations': 0,
            'antiviral_resistance': 0,
            'vaccine_escape': 0,
            'substitutions_of_interest': 0,
            'top_mutations': pd.DataFrame(),
            'mutation_types': pd.DataFrame()
        }
    
    # Calculate basic stats
    total_samples = sequence_df['ID'].nunique() if 'ID' in sequence_df.columns else 0
    total_mutations = mutation_df['ID'].nunique() if 'ID' in mutation_df.columns else 0
    
    # Count by mutation type
    if 'type' in mutation_df.columns:
        mutation_types = mutation_df.groupby('type').size().reset_index(name='count')
        
        # Get counts for specific types
        antiviral_resistance = mutation_types[mutation_types['type'] == 'antiviral_resistance']['count'].sum()
        vaccine_escape = mutation_types[mutation_types['type'] == 'vaccine_escape']['count'].sum()
        substitutions_of_interest = mutation_types[mutation_types['type'] == 'substitution_of_interest']['count'].sum()
        
        # Get top mutations
        top_mutations = (mutation_df.groupby(['mutation', 'type'])
                        .size()
                        .reset_index(name='count')
                        .sort_values('count', ascending=False)
                        .head(10))
    else:
        mutation_types = pd.DataFrame()
        antiviral_resistance = vaccine_escape = substitutions_of_interest = 0
        top_mutations = pd.DataFrame()
    
    return {
        'total_samples': total_samples,
        'total_mutations': total_mutations,
        'antiviral_resistance': antiviral_resistance,
        'vaccine_escape': vaccine_escape,
        'substitutions_of_interest': substitutions_of_interest,
        'top_mutations': top_mutations,
        'mutation_types': mutation_types
    }
    
def normalize_accession_id(s):
    if not isinstance(s, str):
        return s
    return s.strip().upper()

def normalize_genotype_label(genotype: str, virus: str) -> str:
    if not isinstance(genotype, str):
        return genotype

    g = genotype.strip()

    # Handle recombinants explicitly
    if g.lower().startswith("recombinant"):
        return "Recombinant"

    if virus.upper() == "HBV":
        # "Genotype A" -> "HBV-A"
        if g.lower().startswith("genotype"):
            return "HBV-" + g.split()[-1].upper()

    if virus.upper() == "HCV":
        # "Genotype 1" -> "HCV-1"
        if g.lower().startswith("genotype"):
            return "HCV-" + g.split()[-1]
    
    if virus.upper() == "HEV":
        if g.lower().startswith("al_"):
            genotype = g.split()[-1]
            return "HEV-" + genotype.split('_')[-1]

    return g

def load_and_preprocess_data():
    cache_file = get_data_path("results/preprocessed_data_store.pkl")
    
    # List of source files to check for changes
    source_files = [
        "results/hbv/hbv_metadata.tsv",
        "results/hcv/hcv_metadata.tsv",
        "results/hev/hev_metadata.tsv",
        "results/hbv/final_resistance.tsv",
        "results/hcv/final_resistance.tsv",
        "results/hev/final_resistance.tsv",
        "data/population_by_country_year.csv",
        "data/IHME-GBD_2021_DATA-9e7ec2c0-1.csv",
        "data/WHO_regions_countries_coordinates.txt"
    ]
    
    # Check if cache is valid (exists and newer than all existing source files)
    cache_valid = False
    if os.path.exists(cache_file):
        try:
            cache_mtime = os.path.getmtime(cache_file)
            cache_valid = True
            for sf in source_files:
                sf_path = get_data_path(sf)
                if os.path.exists(sf_path) and os.path.getmtime(sf_path) > cache_mtime:
                    cache_valid = False
                    print(f"🔄 Cache invalidated: {sf} has been updated.")
                    break
        except Exception as e:
            print(f"⚠️ Error checking cache mtimes: {e}")
            cache_valid = False
                
    if cache_valid:
        print("🚀 Loading preprocessed data from cache...")
        try:
            import pickle
            with open(cache_file, 'rb') as f:
                data_store = pickle.load(f)
                print("✅ Cache loaded successfully.")
                return data_store
        except Exception as e:
            print(f"⚠️ Failed to load cache: {e}, falling back to full preprocessing...")

    print("📥 Loading data...")
    
    try:
        # Load main data files
        hbv_data = load_with_encoding("results/hbv/hbv_metadata.tsv")
        hcv_data = load_with_encoding("results/hcv/hcv_metadata.tsv")
        hev_data = load_with_encoding("results/hev/hev_metadata.tsv")
        coords = load_with_encoding("data/WHO_regions_countries_coordinates.txt").dropna()

        # --- Normalize metadata columns ---
        def normalize_metadata(df: pd.DataFrame) -> pd.DataFrame:
            if df.empty:
                return df
        
            # accession_id -> ID
            if "accession_id" in df.columns:
                df = df.rename(columns={"accession_id": "ID"})
        
            # submitter_country -> Country (what your standardize function expects)
            if "submitter_country" in df.columns:
                df = df.rename(columns={"submitter_country": "Country"})
        
            # Extract year from isolate_collection_date (handles "2019", "2024/04/16", etc.)
            if "isolate_collection_date" in df.columns:
                df["Year"] = (
                    df["isolate_collection_date"]
                    .astype(str)
                    .str.extract(r"(\d{4})")
                )
                df["Year"] = pd.to_numeric(df["Year"], errors="coerce").astype("Int32")
        
            return df

        hbv_data = normalize_metadata(hbv_data)
        hbv_data["ID"] = hbv_data["ID"].apply(normalize_accession_id)
        hcv_data = normalize_metadata(hcv_data)
        hcv_data["ID"] = hcv_data["ID"].apply(normalize_accession_id)
        hev_data = normalize_metadata(hev_data)
        hev_data["ID"] = hev_data["ID"].apply(normalize_accession_id)

        # Load additional data with better error handling
        population_df = load_csv_file("data/population_by_country_year.csv")
        ihme_df = load_csv_file("data/IHME-GBD_2021_DATA-9e7ec2c0-1.csv")

        # Load mutation data (handle missing files gracefully)
        hbv_mut = load_with_encoding("results/hbv/final_resistance.tsv")
        hcv_mut = load_with_encoding("results/hcv/final_resistance.tsv")
        hev_mut = load_with_encoding("results/hev/final_resistance.tsv")
        
        # --- Normalize GLUE mutations ---
        for df in (hbv_mut, hcv_mut, hev_mut):
            if not df.empty and "sample" in df.columns:
                df.rename(columns={"sample": "ID"}, inplace=True)
        
        if not hbv_mut.empty:
            hbv_mut["ID"] = hbv_mut["ID"].apply(normalize_accession_id)
            
        if not hcv_mut.empty:
            hcv_mut["ID"] = hcv_mut["ID"].apply(normalize_accession_id)
        
        if not hev_mut.empty:
            hev_mut["ID"] = hev_mut["ID"].apply(normalize_accession_id)
        
        # --- Extract genotype/subgenotype from GLUE ---
        hbv_geno = (
            hbv_mut[["ID", "genotype", "subgenotype"]]
            .dropna(subset=["ID"])
            .drop_duplicates()
        ) if not hbv_mut.empty else pd.DataFrame(columns=["ID", "genotype", "subgenotype"])
        
        hcv_geno = (
            hcv_mut[["ID", "genotype", "subgenotype"]]
            .dropna(subset=["ID"])
            .drop_duplicates()
        ) if not hcv_mut.empty else pd.DataFrame(columns=["ID", "genotype", "subgenotype"])
        
        hev_geno = (
            hev_mut[["ID", "genotype", "subgenotype"]]
            .dropna(subset=["ID"])
            .drop_duplicates()
        ) if not hev_mut.empty else pd.DataFrame(columns=["ID", "genotype", "subgenotype"])

        
        hbv_data = hbv_data.merge(hbv_geno, on="ID", how="left")
        hcv_data = hcv_data.merge(hcv_geno, on="ID", how="left")
        hev_data = hev_data.merge(hev_geno, on="ID", how="left")
        
        # --- Normalize genotype labels to dashboard canonical form ---
        if "genotype" in hbv_data.columns:
            hbv_data["genotype"] = hbv_data["genotype"].apply(
                lambda g: normalize_genotype_label(g, virus="HBV")
            )
        
        if "genotype" in hcv_data.columns:
            hcv_data["genotype"] = hcv_data["genotype"].apply(
                lambda g: normalize_genotype_label(g, virus="HCV")
            )

        if "genotype" in hev_data.columns:
            hev_data["genotype"] = hev_data["genotype"].apply(
                lambda g: normalize_genotype_label(g, virus="HEV")
            )

        
        # Only process non-empty dataframes
        hbv_cov = pd.DataFrame()
        hcv_cov = pd.DataFrame()
        hev_cov = pd.DataFrame()
        
        if not hbv_data.empty:
            hbv_data = standardize_country_column(hbv_data, 'Country')
        
        if not hcv_data.empty:
            hcv_data = standardize_country_column(hcv_data, 'Country')

        if not hev_data.empty:
            hev_data = standardize_country_column(hev_data, 'Country')

        
        if not coords.empty:
            coords = standardize_country_column(coords, 'name')
        
        if not population_df.empty:
            population_df = standardize_country_column(population_df, 'Country_standard')
        
        if not ihme_df.empty:
            ihme_df = standardize_country_column(ihme_df, 'location')
            
        if not ihme_df.empty and "Country_standard" in ihme_df.columns:
            ihme_df["Country_standard"] = (
                ihme_df["Country_standard"]
                .astype(str)
                .str.strip()
                .apply(standardize_country_cached)
            )
        
        if not ihme_df.empty:
            for c in ["measure", "sex", "age", "cause", "metric"]:
                if c in ihme_df.columns:
                    ihme_df[c] = ihme_df[c].astype(str).str.strip()
            if "year" in ihme_df.columns:
                ihme_df["year"] = pd.to_numeric(ihme_df["year"], errors="coerce")
            for c in ["val", "upper", "lower"]:
                if c in ihme_df.columns:
                    ihme_df[c] = pd.to_numeric(ihme_df[c], errors="coerce")

        # Sanitize strings efficiently for non-empty dataframes
        for df in [hbv_data, hcv_data, hev_data, coords]:
            if not df.empty:
                df["Country_standard"] = df["Country_standard"].astype(str).str.strip()
    
        # Drop any coords with missing standardization
        if not coords.empty:
            coords = coords.dropna(subset=["Country_standard"])
    
        print("📌 Merging WHO Regions...")
        
        # Only merge if we have data
        if not hbv_data.empty and not coords.empty:
            hbv_data = merge_who_region(hbv_data, coords)
            
        if not hcv_data.empty and not coords.empty:
            hcv_data = merge_who_region(hcv_data, coords)

        if not hev_data.empty and not coords.empty:
            hev_data = merge_who_region(hev_data, coords)

        if not ihme_df.empty and not coords.empty:
            ihme_df = merge_who_region(ihme_df, coords)
    
        # Fix population data if available
        if not population_df.empty:
            population_df["Year"] = pd.to_numeric(population_df["Year"], errors="coerce").astype("Int32")
            population_df["Population"] = pd.to_numeric(population_df["Population"], errors="coerce")
            population_df = (
                population_df
                .dropna(subset=["Year", "Population"])
                .sort_values(["Country_standard", "Year"])
                .drop_duplicates(["Country_standard", "Year"], keep="last")
            )
            population_df.loc[population_df["Population"] <= 0, "Population"] = np.nan
    
        print("📊 Precomputing groupings...")
        # Only group if we have data
        hbv_grouped = pd.DataFrame()
        hcv_grouped = pd.DataFrame()
        hev_grouped = pd.DataFrame()
        
        if not hbv_data.empty:
            grouping_cols = ["Year", "genotype", "WHO_Regions", "Country_standard"]
            hbv_grouped = (hbv_data.groupby(grouping_cols, observed=True)
                           .size()
                           .reset_index(name="Count"))
            
            # Merge with population if available
            if not population_df.empty:
                hbv_grouped = pd.merge(hbv_grouped, population_df, on=["Country_standard", "Year"], how="left")
                hbv_grouped["PerMillion"] = (hbv_grouped["Count"] / hbv_grouped["Population"]) * 1_000_000
        
        if not hcv_data.empty:
            grouping_cols = ["Year", "genotype", "WHO_Regions", "Country_standard"]
            hcv_grouped = (hcv_data.groupby(grouping_cols, observed=True)
                           .size()
                           .reset_index(name="Count"))
            
            if not population_df.empty:
                hcv_grouped = pd.merge(hcv_grouped, population_df, on=["Country_standard", "Year"], how="left")
                hcv_grouped["PerMillion"] = (hcv_grouped["Count"] / hcv_grouped["Population"]) * 1_000_000
                
        if not hev_data.empty:
            grouping_cols = ["Year", "genotype", "WHO_Regions", "Country_standard"]
            hev_grouped = (hev_data.groupby(grouping_cols, observed=True)
                           .size()
                           .reset_index(name="Count"))
            
            if not population_df.empty:
                hev_grouped = pd.merge(hev_grouped, population_df, on=["Country_standard", "Year"], how="left")
                hev_grouped["PerMillion"] = (hev_grouped["Count"] / hev_grouped["Population"]) * 1_000_000

        
        hbv_cov = compute_burden_adjusted_coverage_from_sequences(
            hbv_data, ihme_df, virus="HBV"
        )
        
        hcv_cov = compute_burden_adjusted_coverage_from_sequences(
            hcv_data, ihme_df, virus="HCV"
        )

        hev_cov = compute_burden_adjusted_coverage_from_sequences(
            hev_data, ihme_df, virus="HEV"
        )

        print("📍 Creating coordinate lookup...")
        coord_lookup = {}
        if not coords.empty:
            coord_lookup = (
                coords.drop_duplicates(subset='Country_standard')
                      .set_index('Country_standard')[['latitude', 'longitude']]
                      .to_dict('index')
            )
            
            # Add manual overrides
            coord_lookup.update({
                "Nigeria": {"latitude": 9.0820, "longitude": 7.49508},
                "New Caledonia": {"latitude": -21.450553, "longitude": 165.505710}
            })

        print("🧬 Processing mutation data...")
        hbv_mut_stats = pd.DataFrame()
        hcv_mut_stats = pd.DataFrame()
        hev_mut_stats = pd.DataFrame()

        # Merge mutations onto metadata+genotype table (hbv_data/hcv_data)
        if not hbv_mut.empty and not hbv_data.empty:
            hbv_mut = hbv_mut.merge(
                hbv_data[["ID", "Country_standard", "Year", "genotype", "subgenotype", "WHO_Regions"]].drop_duplicates("ID"),
                on="ID",
                how="left"
            )
            hbv_mut = hbv_mut.dropna(subset=["Country_standard", "Year"])

            hbv_mut_stats = (
                hbv_mut.groupby(["Country_standard", "type", "drug", "mutation", "gene"], observed=True)
                .size()
                .reset_index(name="Count")
            )

        if not hcv_mut.empty and not hcv_data.empty:
            hcv_mut = hcv_mut.merge(
                hcv_data[["ID", "Country_standard", "Year", "genotype", "subgenotype", "WHO_Regions"]].drop_duplicates("ID"),
                on="ID",
                how="left"
            )
            hcv_mut = hcv_mut.dropna(subset=["Country_standard", "Year"])

            hcv_mut_stats = (
                hcv_mut.groupby(["Country_standard", "type", "drug", "mutation", "gene"], observed=True)
                .size()
                .reset_index(name="Count")
            )

        if not hev_mut.empty and not hev_data.empty:
            hev_mut = hev_mut.merge(
                hev_data[["ID", "Country_standard", "Year", "genotype", "subgenotype", "WHO_Regions"]].drop_duplicates("ID"),
                on="ID",
                how="left"
            )
            hev_mut = hev_mut.dropna(subset=["Country_standard", "Year"])

            hev_mut_stats = (
                hev_mut.groupby(["Country_standard", "type", "drug", "mutation", "gene"], observed=True)
                .size()
                .reset_index(name="Count")
            )

        print("✅ Done loading data.")
        
        data_store = {
            'hbv_data': hbv_data,
            'hcv_data': hcv_data,
            'hev_data': hev_data,
            'coords': coords,
            'hbv_grouped': hbv_grouped,
            'hcv_grouped': hcv_grouped,
            'hev_grouped': hev_grouped,
            'coord_lookup': coord_lookup,
            'population_df': population_df,
            'ihme_df': ihme_df,
            'hbv_mut': hbv_mut,
            'hcv_mut': hcv_mut,
            'hev_mut': hev_mut,
            'hbv_mut_stats': hbv_mut_stats,
            'hcv_mut_stats': hcv_mut_stats,
            'hev_mut_stats': hev_mut_stats,
            'cov_hbv': hbv_cov,
            'cov_hcv': hcv_cov,
            'cov_hev': hev_cov,
            'hbv_summary_raw': hbv_data.copy(),
            'hcv_summary_raw': hcv_data.copy(),
            'hev_summary_raw': hev_data.copy()
        }
        
        # Save to cache
        print("💾 Saving preprocessed data to cache...")
        try:
            import pickle
            os.makedirs(os.path.dirname(cache_file), exist_ok=True)
            with open(cache_file, 'wb') as f:
                pickle.dump(data_store, f, protocol=pickle.HIGHEST_PROTOCOL)
            print("✅ Cache saved successfully.")
        except Exception as e:
            print(f"⚠️ Failed to save cache: {e}")
            
        return data_store
        
    except Exception as e:
        print(f"❌ Critical error loading data: {e}")
        # Return empty data structure instead of crashing
        return {
            'hbv_data': pd.DataFrame(),
            'hcv_data': pd.DataFrame(),
            'hev_data': pd.DataFrame(),
            'coords': pd.DataFrame(),
            'hbv_grouped': pd.DataFrame(),
            'hcv_grouped': pd.DataFrame(),
            'hev_grouped': pd.DataFrame(),
            'coord_lookup': {},
            'population_df': pd.DataFrame(),
            'ihme_df': pd.DataFrame(),
            'hbv_mut': pd.DataFrame(),
            'hcv_mut': pd.DataFrame(),
            'hev_mut': pd.DataFrame(),
            'hbv_mut_stats': pd.DataFrame(),
            'hcv_mut_stats': pd.DataFrame(),
            'hev_mut_stats': pd.DataFrame(),
            'cov_hbv': pd.DataFrame(),
            'cov_hcv': pd.DataFrame(),
            'cov_hev': pd.DataFrame(),
            'hbv_summary_raw': pd.DataFrame(),
            'hcv_summary_raw': pd.DataFrame(),
            'hcv_summary_raw': pd.DataFrame()
        }