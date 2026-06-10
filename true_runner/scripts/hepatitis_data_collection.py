import time
import pandas as pd
import requests as api
from Bio import Entrez, SeqIO
import ssl

from typing import List, Optional
from pydantic import BaseModel

ncbi_api_url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"

# ensure SSL compatibility
ssl._create_default_https_context = ssl._create_unverified_context

###########################################################
# Pydantic Models (metadata structures)
###########################################################

class VirusIsolate(BaseModel):
    name: Optional[str] = None
    source: Optional[str] = None
    collection_date: Optional[str] = None

class Lineage(BaseModel):
    tax_id: Optional[int] = None
    name: Optional[str] = None

class VirusHost(BaseModel):
    tax_id: Optional[int] = None
    sci_name: Optional[str] = None
    organism_name: Optional[str] = None
    common_name: Optional[str] = None
    lineage: Optional[List[Lineage]] = None
    strain: Optional[str] = None

class Virus(BaseModel):
    tax_id: Optional[int] = None
    sci_name: Optional[str] = None
    organism_name: Optional[str] = None
    common_name: Optional[str] = None
    lineage: Optional[List[Lineage]] = None
    strain: Optional[str] = None

class ReportLocation(BaseModel):
    geographic_location: Optional[str] = None
    geographic_region: Optional[str] = None

class ReportNucleotide(BaseModel):
    seq_id: Optional[str] = None
    accession_version: Optional[str] = None
    title: Optional[str] = None

class ReportSubmitter(BaseModel):
    names: Optional[List[str]] = None
    affiliation: Optional[str] = None
    country: Optional[str] = None

###########################################################
# Response Structures
###########################################################

class VirusReport(BaseModel):
    accession: Optional[str] = None
    is_complete: Optional[bool] = None
    is_annotated: Optional[bool] = None
    isolate: Optional[VirusIsolate] = None
    host: Optional[VirusHost] = None
    virus: Optional[Virus] = None
    location: Optional[ReportLocation] = None
    nucleotide: Optional[ReportNucleotide] = None
    length: Optional[int] = None
    gene_count: Optional[int] = None
    submitter: Optional[ReportSubmitter] = None
    update_date: Optional[str] = None
    release_date: Optional[str] = None
    source_database: Optional[str] = None

class VirusMetadataResponse(BaseModel):
    reports: Optional[List[VirusReport]] = None
    total_count: Optional[int] = None
    next_page_token: Optional[str] = None

class VirusFilter(BaseModel):
    taxon: Optional[str] = None
    accessions: Optional[List[str]] = None
    complete_only: Optional[bool] = None
    annotated_only: Optional[bool] = None
    host: Optional[str] = None
    geo_location: Optional[str] = None

class VirusRequestPayload(BaseModel):
    filter: Optional[VirusFilter] = None
    page_size: Optional[int] = None
    page_token: Optional[str] = None

###########################################################
# Download FASTA sequences
###########################################################

def fetch_sequences_batch(accession_list: List[str], batch_size=200, output_file=None, sleep_time=0.1):
    """
    Saves FASTA sequences to output_file.
    Snakemake supplies output_file, so no hard-coded paths are allowed.
    """
    if output_file is None:
        raise ValueError("output_file must be supplied")

    with open(output_file, "w") as out_handle:
        for i in range(0, len(accession_list), batch_size):
            batch = accession_list[i:i + batch_size]

            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=",".join(batch),
                    rettype="fasta",
                    retmode="text"
                )
                SeqIO.write(SeqIO.parse(handle, "fasta"), out_handle, "fasta")
                handle.close()
                time.sleep(sleep_time)

            except Exception as e:
                print(f"[ERROR] batch {batch}: {e}")

    return output_file

###########################################################
# Download metadata from NCBI Datasets API
###########################################################

def fetch_virus_metadata(payload: VirusRequestPayload, max_retries=5, backoff_factor=2) -> VirusMetadataResponse:
    url = f"{ncbi_api_url}/virus"
    data = payload.model_dump(exclude_none=True)
    
    for attempt in range(max_retries):
        try:
            response = api.post(url, json=data)
            
            # Handle rate limiting (429) and temporary server issues (5xx)
            if response.status_code == 429 or response.status_code >= 500:
                sleep_time = (backoff_factor ** attempt) + 1
                status_desc = "Rate limited (429)" if response.status_code == 429 else f"Server error ({response.status_code})"
                print(f"[WARNING] {status_desc} from NCBI. Retrying in {sleep_time}s (attempt {attempt + 1}/{max_retries})...")
                time.sleep(sleep_time)
                continue
                
            response.raise_for_status()
            
            try:
                response_json = response.json()
                return VirusMetadataResponse(**response_json)
            except Exception as e:
                # If decoding failed, it might have been an HTML error response under status code 200
                print(f"[ERROR] Failed to parse JSON from NCBI response. Status: {response.status_code}. Text: {response.text[:500]}")
                if attempt < max_retries - 1:
                    sleep_time = (backoff_factor ** attempt) + 1
                    print(f"Retrying in {sleep_time}s...")
                    time.sleep(sleep_time)
                    continue
                raise e
                
        except Exception as e:
            if attempt == max_retries - 1:
                print(f"[FATAL] NCBI datasets request failed after {max_retries} attempts.")
                raise e
            sleep_time = (backoff_factor ** attempt) + 1
            print(f"[WARNING] Request failed: {e}. Retrying in {sleep_time}s (attempt {attempt + 1}/{max_retries})...")
            time.sleep(sleep_time)

def fetch_all_virus_metadata(filters: VirusFilter, page_size=200):
    reports = []
    token = None

    while True:
        payload = VirusRequestPayload(
            filter=filters,
            page_size=page_size,
            page_token=token
        )
        response = fetch_virus_metadata(payload)

        if not response.reports:
            break

        reports.extend(response.reports)

        if not response.next_page_token:
            break

        token = response.next_page_token
        # Sleep for 0.6s to respect NCBI API rate limit of 5 requests per second across parallel jobs
        time.sleep(0.6)

    return reports

###########################################################
# Convert nested metadata → flat dictionaries
###########################################################

def parse_response_to_dictionaries(records: List[VirusReport]) -> List[dict]:
    out = []
    for r in records:
        out.append({
            "accession_id": r.accession,
            "length": r.length,
            "is_complete": r.is_complete,
            "is_annotated": r.is_annotated,
            "isolate_collection_date": r.isolate.collection_date if r.isolate else None,
            "location": r.location.geographic_location if r.location else None,
            "region": r.location.geographic_region if r.location else None,
            "submitter_country": r.submitter.country if r.submitter else None,
            "release_date": r.release_date,
            "update_date": r.update_date,
        })
    return out
