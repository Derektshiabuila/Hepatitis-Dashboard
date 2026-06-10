import sys
import xml.etree.ElementTree as ET
import pandas as pd

def extract_hev_genotype(root):
    """
    Extract genotype and subgenotype from HEV genotypeCommandResult XML.
    Returns (genotype, subgenotype)
    """
    genotype = "NA"
    subgenotype = "NA"

    row = root.find(".//row")
    if row is not None:
        values = row.findall("value")
        if len(values) >= 2:
            genotype = values[1].text or "NA"
        if len(values) >= 3:
            subgenotype = values[2].text or "NA"

    return genotype, subgenotype

def write_empty(sample, virus, out_tsv, reason):
    df = pd.DataFrame([{
        "sample": sample,
        "virus": virus.upper(),
        "genotype": "NA",
        "subgenotype": "NA",
        "gene": "NA",
        "drug": "NA",
        "mutation": "NA",
        "detected": "NA",
        "type": reason
    }])
    df.to_csv(out_tsv, sep="\t", index=False)

def main():
    xml_file = sys.argv[1]
    out_tsv = sys.argv[2]

    sample = xml_file.split("/")[-1].replace(".xml", "")

    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
    except ET.ParseError as e:
        write_empty(sample, "hcv", out_tsv, "xml_parse_error")
        sys.exit(0)

    genotype, subgenotype = extract_hev_genotype(root)

    # Single row output to match HBV schema
    rows = [{
        "sample": sample,
        "virus": "HEV",
        "genotype": genotype,
        "subgenotype": subgenotype,
        "gene": "NA",
        "drug": "NA",
        "mutation": "NA",
        "detected": "NA",
        "type": "genotype"
    }]

    df = pd.DataFrame(rows)
    df.to_csv(out_tsv, sep="\t", index=False)

    print(f"Extracted HEV genotype info")
    print(f"Genotype: {genotype}, Subgenotype: {subgenotype}")


if __name__ == "__main__":
    main()

