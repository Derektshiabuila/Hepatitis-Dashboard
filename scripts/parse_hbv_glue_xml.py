import sys
import xml.etree.ElementTree as ET
import pandas as pd

def extract_genotype_info(root):
    genotype = "NA"
    subgenotype = "NA"

    genotype_elem = root.find(".//genotypeCategoryResult")
    if genotype_elem is not None:
        final_clade = genotype_elem.find("finalCladeRenderedName")
        if final_clade is not None and final_clade.text:
            genotype = final_clade.text

    subgenotype_elem = root.find(".//subgenotypeCategoryResult")
    if subgenotype_elem is not None:
        final_clade = subgenotype_elem.find("finalCladeRenderedName")
        if final_clade is not None and final_clade.text:
            subgenotype = final_clade.text
        else:
            short_name = subgenotype_elem.find("shortRenderedName")
            if short_name is not None and short_name.text:
                subgenotype = short_name.text

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
        write_empty(sample, "hbv", out_tsv, "xml_parse_error")
        sys.exit(0)

    # Handle GLUE failure XML
    if root.tag == "glueFailure":
        df = pd.DataFrame([{
            "sample": sample,
            "virus": "HBV",
            "genotype": "NA",
            "subgenotype": "NA",
            "gene": "NA",
            "drug": "NA",
            "mutation": "NA",
            "detected": "NA",
            "type": "glue_failure"
        }])
        df.to_csv(out_tsv, sep="\t", index=False)
        return

    genotype, subgenotype = extract_genotype_info(root)
    rows = []

    for res in root.findall(".//antiviralResistance"):
        rows.append({
            "sample": sample,
            "virus": "HBV",
            "genotype": genotype,
            "subgenotype": subgenotype,
            "gene": res.findtext("virusDomain", default="RT"),
            "drug": res.findtext("drug/name", default="NA"),
            "mutation": res.findtext("description", default="NA"),
            "detected": res.findtext("detectedPattern", default="NA"),
            "type": "antiviral_resistance"
        })

    for vac in root.findall(".//vaccineEscape"):
        rows.append({
            "sample": sample,
            "virus": "HBV",
            "genotype": genotype,
            "subgenotype": subgenotype,
            "gene": vac.findtext("virusDomain", default="S"),
            "drug": "Vaccine",
            "mutation": vac.findtext("description", default="NA"),
            "detected": vac.findtext("detectedPattern", default="NA"),
            "type": "vaccine_escape"
        })

    if not rows:
        rows.append({
            "sample": sample,
            "virus": "HBV",
            "genotype": genotype,
            "subgenotype": subgenotype,
            "gene": "POL",
            "drug": "NA",
            "mutation": "NA",
            "detected": "NA",
            "type": "none"
        })

    pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
