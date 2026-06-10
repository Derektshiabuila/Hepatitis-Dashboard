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

def fetch_virus_metadata(payload: VirusRequestPayload) -> VirusMetadataResponse:
    response = api.post(
        f"{ncbi_api_url}/virus",
        json=payload.model_dump(exclude_none=True)
    ).json()

    return VirusMetadataResponse(**response)

def fetch_all_virus_metadata(filters: VirusFilter, page_size=1000):
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
        time.sleep(0.1)

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
