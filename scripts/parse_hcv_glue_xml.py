#!/usr/bin/env python3
"""
Parse HCV GLUE XML output and extract resistance information in HBV-like format.
"""
import sys
import xml.etree.ElementTree as ET
import pandas as pd

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


def parse_hcv_glue_xml(xml_file, out_tsv):
    """
    Parse HCV GLUE XML file and extract resistance mutations and drug assessments.
    """
    # Extract sample name from file path
    sample = xml_file.split("/")[-1].replace(".xml", "")
    
    # Parse XML
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
    except ET.ParseError as e:
        write_empty(sample, "hcv", out_tsv, "xml_parse_error")
        sys.exit(0)
    
    rows = []
    
    # Get genotype and subtype information
    genotype = "Unknown"
    subtype = "Unknown"
    
    # Extract genotype and subtype from genotypingResult
    genotyping_result = root.find(".//genotypingResult")
    if genotyping_result is not None:
        # Get genotype - try multiple approaches
        genotype_elem = genotyping_result.find(".//categoryName[.='genotype']/..")
        if genotype_elem is not None:
            # Try finalCladeRenderedName first
            final_clade_elem = genotype_elem.find("finalCladeRenderedName")
            if final_clade_elem is not None and final_clade_elem.text:
                genotype = final_clade_elem.text.replace("HCV ", "")
            else:
                # Try shortRenderedName as fallback
                short_name_elem = genotype_elem.find("shortRenderedName")
                if short_name_elem is not None and short_name_elem.text:
                    genotype = f"Genotype {short_name_elem.text}"
        
        # Get subtype - try multiple approaches
        subtype_elem = genotyping_result.find(".//categoryName[.='subtype']/..")
        if subtype_elem is not None:
            # Try finalCladeRenderedName first
            final_clade_elem = subtype_elem.find("finalCladeRenderedName")
            if final_clade_elem is not None and final_clade_elem.text:
                subtype = final_clade_elem.text.replace("HCV ", "")
            else:
                # Try shortRenderedName as fallback
                short_name_elem = subtype_elem.find("shortRenderedName")
                if short_name_elem is not None and short_name_elem.text:
                    subtype = f"Subtype {short_name_elem.text}"
                else:
                    # Try finalClade as last resort
                    final_clade_elem = subtype_elem.find("finalClade")
                    if final_clade_elem is not None and final_clade_elem.text:
                        # Try to parse clade name like "AL_1b" to "1b"
                        clade = final_clade_elem.text
                        if clade.startswith("AL_"):
                            subtype = f"Subtype {clade[3:]}"
    
    # ====================================================================
    # 1. Look for detected RAS mutations in rasScanResults
    # ====================================================================
    ras_scan_sections = root.findall(".//rasScanResults")
    for ras_scan in ras_scan_sections:
        # Check if this RAS is present in the sample
        present_elem = ras_scan.find("present")
        if present_elem is None or present_elem.text != "true":
            continue
            
        # Get basic mutation info
        feature_name_elem = ras_scan.find("featureName")  # e.g., "NS3", "NS5B"
        variation_name_elem = ras_scan.find("variationName")  # e.g., "phdr_ras:NS3:54S"
        
        if feature_name_elem is None or variation_name_elem is None:
            continue
            
        gene = feature_name_elem.text
        
        # Extract mutation from variation name
        # Format: "phdr_ras:NS3:54S" -> mutation is "54S"
        variation_name = variation_name_elem.text
        if ":" in variation_name:
            mutation = variation_name.split(":")[-1]
        else:
            mutation = variation_name
        
        # Get detailed RAS information
        ras_details = ras_scan.find("rasDetails")
        if ras_details is not None:
            # Look for drug resistance info in alignmentRasDrug elements
            alignment_ras_drugs = ras_details.findall(".//alignmentRasDrug")
            
            if alignment_ras_drugs:
                for drug_elem in alignment_ras_drugs:
                    drug_name_elem = drug_elem.find("drug")
                    resistance_cat_elem = drug_elem.find("resistanceCategory")
                    
                    if drug_name_elem is not None:
                        drug = drug_name_elem.text
                        
                        # Determine resistance type
                        resistance_category = "Unknown"
                        if resistance_cat_elem is not None:
                            resistance_category = resistance_cat_elem.text
                        
                        # Check if this is significant resistance
                        if resistance_category in ["category_I", "category_II", "category_III", "I", "II", "III"]:
                            detected = "Detected"
                            resistance_type = "antiviral_resistance"
                        elif resistance_category in ["insignificant", "-"]:
                            detected = "Not significant"
                            resistance_type = "no_resistance"
                        else:
                            detected = "Detected" if resistance_category != "insignificant" else "Not significant"
                            resistance_type = "antiviral_resistance" if resistance_category not in ["insignificant", "-"] else "no_resistance"
                        
                        rows.append({
                            "sample": sample,
                            "virus": "HCV",
                            "genotype": genotype,
                            "subgenotype": subtype,
                            "gene": gene,
                            "drug": drug,
                            "mutation": mutation,
                            "detected": detected,
                            "type": resistance_type
                        })
            else:
                # Just record the mutation without specific drug info
                rows.append({
                    "sample": sample,
                    "virus": "HCV",
                    "genotype": genotype,
                    "subgenotype": subtype,
                    "gene": gene,
                    "drug": "NA",
                    "mutation": mutation,
                    "detected": "Detected",
                    "type": "mutation_detected"
                })
    
    # ====================================================================
    # 2. Look at drugScores for comprehensive resistance summary
    # ====================================================================
    drug_score_categories = root.findall(".//drugScores")
    for category_elem in drug_score_categories:
        category_name_elem = category_elem.find("category")
        category_name = category_name_elem.text if category_name_elem is not None else "Unknown"
        
        drug_assessments = category_elem.findall(".//drugAssessments")
        for assessment in drug_assessments:
            # Get drug name
            drug_elem = assessment.find(".//drug/id")
            if drug_elem is None:
                continue
            drug = drug_elem.text
            
            # Get resistance assessment
            drug_score_elem = assessment.find("drugScoreDisplay")
            if drug_score_elem is None:
                continue
            resistance_text = drug_score_elem.text
            
            # Extract specific RAS mutations causing resistance
            mutations_found = []
            
            # Look for RAS scores in different categories
            for ras_prefix in ["rasScores_category_I", "rasScores_category_II", "rasScores_category_III"]:
                ras_scores = assessment.findall(f".//{ras_prefix}")
                for ras in ras_scores:
                    gene_elem = ras.find("gene")
                    structure_elem = ras.find("structure")
                    display_elem = ras.find("displayStructure")
                    
                    if gene_elem is not None and structure_elem is not None:
                        mutation = structure_elem.text
                        display_mutation = display_elem.text if display_elem is not None else mutation
                        
                        # Determine resistance type
                        category = ras_prefix.replace("rasScores_category_", "")
                        if category in ["I", "II", "III"]:
                            detected = "Detected"
                            resistance_type = "antiviral_resistance"
                        else:
                            detected = "Not significant" if category == "insignificant" else "Detected"
                            resistance_type = "no_resistance" if category == "insignificant" else "antiviral_resistance"
                        
                        rows.append({
                            "sample": sample,
                            "virus": "HCV",
                            "genotype": genotype,
                            "subgenotype": subtype,
                            "gene": gene_elem.text,
                            "drug": drug,
                            "mutation": mutation,
                            "detected": detected,
                            "type": resistance_type
                        })
            
            # If no specific mutations but drug resistance is indicated
            if not mutations_found and any(word in resistance_text.lower() for word in ["resistance", "probable"]):
                rows.append({
                    "sample": sample,
                    "virus": "HCV",
                    "genotype": genotype,
                    "subgenotype": subtype,
                    "gene": "NA",
                    "drug": drug,
                    "mutation": "NA",
                    "detected": "Detected" if "resistance" in resistance_text.lower() else "Potential",
                    "type": "antiviral_resistance" if "resistance" in resistance_text.lower() else "potential_resistance"
                })
    
    # ====================================================================
    # 3. Check for substitutions of interest
    # ====================================================================
    substitutions = root.findall(".//substitutionsOfInterest")
    for sub in substitutions:
        protein_elem = sub.find("virusProtein")
        display_elem = sub.find("displayStructure")
        reason_elem = sub.find("displayReasonForInterest")
        
        if protein_elem is not None and display_elem is not None:
            mutation = display_elem.text
            reason = reason_elem.text if reason_elem is not None else "Substitution of interest"
            
            rows.append({
                "sample": sample,
                "virus": "HCV",
                "genotype": genotype,
                "subgenotype": subtype,
                "gene": protein_elem.text,
                "drug": "NA",
                "mutation": mutation,
                "detected": "Of interest",
                "type": "substitution_of_interest"
            })
    
    # ====================================================================
    # 4. Create DataFrame and handle duplicates
    # ====================================================================
    if not rows:
        # No resistance detected
        rows.append({
            "sample": sample,
            "virus": "HCV",
            "genotype": genotype,
            "subgenotype": subtype,
            "gene": "NA",
            "drug": "NA",
            "mutation": "NA",
            "detected": "None",
            "type": "no_resistance"
        })
    
    df = pd.DataFrame(rows)
    
    # Remove exact duplicates
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.drop_duplicates()
    
    # Sort by gene, then drug
    if not df.empty:
        df = df.sort_values(["gene", "drug", "detected"])
    
    # Save to TSV
    df.to_csv(out_tsv, sep="\t", index=False)
    
    # Print summary
    num_resistant = len(df[df["type"].str.contains("antiviral_resistance|potential", na=False)])
    print(f"[INFO] Processed {sample}")
    print(f"[INFO] Genotype: {genotype}, Subtype: {subtype}")
    print(f"[INFO] Found {len(df)} total records, {num_resistant} with resistance")
    print(f"[INFO] Saved to {out_tsv}")
    
    return df

def main():
    if len(sys.argv) != 3:
        print("Usage: python parse_hcv_glue_xml.py <input.xml> <output.tsv>")
        sys.exit(1)
    
    xml_file = sys.argv[1]
    out_tsv = sys.argv[2]
    
    try:
        df = parse_hcv_glue_xml(xml_file, out_tsv)
        
        # Display a preview
        if not df.empty:
            print("\nPreview of extracted data:")
            print(df.head().to_string())
    except Exception as e:
        print(f"[ERROR] Failed to parse {xml_file}: {e}")
        # Create an empty output file with just headers to avoid breaking the pipeline
        with open(out_tsv, 'w') as f:
            f.write("sample\tvirus\tgenotype\tsubgenotype\tgene\tdrug\tmutation\tdetected\ttype\n")
        print(f"[INFO] Created empty output file {out_tsv}")
        sys.exit(0)  # Exit with 0 to not break the snakemake pipeline

if __name__ == "__main__":
    main()