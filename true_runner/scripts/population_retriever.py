import requests
import pandas as pd
import tempfile
import zipfile
import glob
import os

WORLD_BANK_URL = "http://api.worldbank.org/v2/en/indicator/SP.POP.TOTL?downloadformat=csv"
OUTPUT_PATH = "data/population_by_country_year.csv"


def fetch_population_data_all_years() -> pd.DataFrame:
    print("🌍 Downloading population data from World Bank...")

    with tempfile.TemporaryDirectory() as tmpdir:
        zip_path = os.path.join(tmpdir, "population_data.zip")

        response = requests.get(WORLD_BANK_URL)
        response.raise_for_status()

        with open(zip_path, "wb") as f:
            f.write(response.content)

        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            zip_ref.extractall(tmpdir)

        matches = glob.glob(os.path.join(tmpdir, "API_SP.POP.TOTL*.csv"))
        if not matches:
            raise FileNotFoundError("World Bank population CSV not found in ZIP")

        data_file = matches[0]

        df = pd.read_csv(data_file, skiprows=4)

        df_long = df.melt(
            id_vars=["Country Name", "Country Code"],
            var_name="Year",
            value_name="Population"
        )

        df_long["Year"] = pd.to_numeric(df_long["Year"], errors="coerce")
        df_long["Population"] = pd.to_numeric(df_long["Population"], errors="coerce")

        df_long = df_long.dropna(subset=["Year", "Population"])

        df_long = df_long.rename(columns={
            "Country Name": "Country_standard"
        })

        return df_long[["Country_standard", "Year", "Population"]]


if __name__ == "__main__":
    df = fetch_population_data_all_years()

    os.makedirs("data", exist_ok=True)
    df.to_csv(OUTPUT_PATH, index=False)

    print(f"✅ Population data written to {OUTPUT_PATH}")
