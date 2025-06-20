import argparse
import re
import os
import csv
import pandas as pd
from Bio import SeqIO

def read_fasta(file: str, part_index: int, separator: str = " ") -> dict:
    """Reads a FASTA file and returns a dictionary of SeqRecord objects, keyed by the specified part of the header."""
    seq_records = {}
    print(f"Reading FASTA file: {file}")
    with open(file, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            header_parts = record.description.split(separator)
            if 1 <= part_index <= len(header_parts):
                part = header_parts[part_index - 1].strip()
                seq_records[part] = record
                # print("Header:", record.description)  # Print the header
    print(f"Processed {len(seq_records)} sequences from FASTA file")
    return seq_records

def read_gtf(file: str, column_type: str, pattern: str, postfix: str = "") -> dict:
    """Reads a GTF file and returns a dictionary of gene information, keyed by specified ID type."""
    gtf_data = {}
    regex = re.compile(pattern)
    print(f"Reading GTF file: {file}")
    with open(file, "r") as gtf_file:
        for line in gtf_file:
            if not line.startswith("#"):
                columns = line.strip().split("\t")
                if columns[2] == column_type:
                    match = regex.search(columns[8])
                    if match is not None:
                        identifier = match.group(1) + postfix
                        gtf_data[identifier] = {"chr": columns[0], "start": int(columns[3]), "end": int(columns[4])}
    print(f"Processed {len(gtf_data)} entries from GTF file")
    return gtf_data

def read_csv(file: str) -> dict:
    """Reads a CSV file and returns a dictionary of gene descriptions, keyed by protein ID."""
    descriptions = {}
    print(f"Reading descriptions CSV file: {file}")
    with open(file, "r") as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            target_organisms = row["Target Organism"].strip('"').split(", ")
            for organism in target_organisms:
                protein_id = organism
                gene_info = {"Locus": row["Locus"], "Description": row["Description"]}
                descriptions[protein_id] = gene_info
    print(f"Processed {len(descriptions)} entries from descriptions file")
    return descriptions

def read_model_descriptions(file: str) -> dict:
    """Reads model gene descriptions from CSV file."""
    model_descriptions = {}
    with open(file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            locus = row["locus"].strip()
            model_descriptions[locus] = {
                "Description": row["description"].split("Description: ")[-1].split(", Other Name:")[0],
                "Other Name": row["description"].split("Other Name: ")[-1].split(", Keywords:")[0],
                "Keywords": row["description"].split("Keywords: ")[-1].strip()
            }
    return model_descriptions

def read_target_mappings(file: str) -> dict:
    """Reads target to orthogroup/model gene mappings."""
    target_mappings = {}
    with open(file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            target_genes = row["Target Organism"].strip('"').split(", ")
            for target in target_genes:
                if target not in target_mappings:
                    target_mappings[target] = []
                target_mappings[target].append({
                    "Orthogroup": row["Orthogroup"],
                    "Model Gene": row["Locus"]
                })
    return target_mappings

def genes_in_range(fasta_file: str, gtf_file: str, chromosome: str, start=None, end=None, descriptions_directory=None,
                   fasta_part_index=1, fasta_separator=" ", gtf_column_type="CDS", gtf_pattern="", postfix=""):
    """Finds all proteins located within a specified range on a specified chromosome."""
    print("\n--- STEP 1: Loading data files ---")
    seq_records = read_fasta(fasta_file, fasta_part_index, fasta_separator)
    gtf_data = read_gtf(gtf_file, gtf_column_type, gtf_pattern, postfix)
    
    # Path to the descriptions file in the specified directory
    descriptions_file = None
    if descriptions_directory:
        descriptions_file = os.path.join(descriptions_directory, "target_matched_to_model_orthologs.csv")
    
    descriptions = read_csv(descriptions_file) if descriptions_file and os.path.exists(descriptions_file) else {}

    # Print contents of seq_records
    # print("Contents of seq_records:")
    # for key, value in seq_records.items():
    #     print(key, value)

    # Print contents of gtf_data
    # print("\nContents of gtf_data:")
    # for key, value in gtf_data.items():
    #     print(key, value)

    # Print contents of descriptions
    # print("\nContents of descriptions:")
    # for key, value in descriptions.items():
    #     print(key, value)

    print(f"\n--- STEP 2: Finding orthologs in chromosome {chromosome} {f'from position {start} to {end}' if start and end else ''} ---")
    num_orthologs_found = 0
    results = []

    for protein_id, record in seq_records.items():
        gene_info = gtf_data.get(protein_id)
        if gene_info:
            if gene_info["chr"] == chromosome:
                if start is None or gene_info["start"] >= start:
                    if end is None or gene_info["end"] <= end:
                        description = descriptions.get(protein_id, 'N/A') if descriptions else "N/A"
                        # message = f"Protein {protein_id} ({record.description}) is located within the specified range at position {gene_info['start']} to {gene_info['end']}."
                        # print(message)
                        # Store result for later parsing
                        results.append({
                            'protein_id': protein_id,
                            'description': record.description,
                            'start': gene_info['start'],
                            'end': gene_info['end'],
                            'model_description': description
                        })
                        num_orthologs_found += 1
        # else:
        #     message = f"Protein {protein_id} not found in the gtf file."
        #     print(message)

    if num_orthologs_found == 0:
        print("No orthologs were found in the specified region of interest.")
        return None
    
    print(f"Found {num_orthologs_found} orthologs in the specified region")
    return results

def find_multi_matches(results: list) -> dict:
    """Find genes with multiple Arabidopsis matches."""
    gene_matches = {}
    
    # Group results by target gene
    for result in results:
        protein_id = result['Protein']
        model_gene = result['Arabidopsis Gene Model']
        orthogroup = result['Orthogroup']
        
        if protein_id not in gene_matches:
            gene_matches[protein_id] = {
                'models': set(),
                'orthogroup': orthogroup
            }
        gene_matches[protein_id]['models'].add(model_gene)
    
    # Filter for only multi-matches
    multi_matches = {
        gene: info for gene, info in gene_matches.items() 
        if len(info['models']) > 1
    }
    
    return multi_matches

def write_warning_file(multi_matches: dict, output_dir: str):
    """Write warning file for genes with multiple matches."""
    warning_file = os.path.join(output_dir, "genes_with_warnings.csv")
    
    with open(warning_file, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['Target_Gene', 'Orthogroup', 'Arabidopsis_Models'])
        
        for gene, info in multi_matches.items():
            writer.writerow([
                gene,
                info['orthogroup'],
                ';'.join(sorted(info['models']))
            ])

def parse_protein_data(results, descriptions_directory, output_csv):
    """Parse protein data directly from the results of genes_in_range and merge with locus data."""
    if not results:
        print("No results to parse.")
        return
        
    print("\n--- STEP 3: Parsing ortholog data ---")
    
    # Load model descriptions and target mappings
    model_desc_file = os.path.join(descriptions_directory, "model_locus_matched_to_orthogroup_output.csv")
    target_map_file = os.path.join(descriptions_directory, "target_matched_to_model_orthologs.csv")
    
    model_descriptions = read_model_descriptions(model_desc_file)
    target_mappings = read_target_mappings(target_map_file)
    
    # Initialize output data
    output_rows = []
    
    # Process each result
    for result in results:
        protein_id = result['protein_id']
        
        # Get all model genes mapped to this target
        if protein_id in target_mappings:
            for mapping in target_mappings[protein_id]:
                model_gene = mapping["Model Gene"]
                orthogroup = mapping["Orthogroup"]
                
                # Get description for this specific model gene
                desc_data = model_descriptions.get(model_gene, {
                    "Description": "N/A",
                    "Other Name": "N/A",
                    "Keywords": "N/A"
                })
                
                output_rows.append({
                    "Protein": protein_id,
                    "Orthogroup": orthogroup,
                    "Arabidopsis Gene Model": model_gene,
                    "Start": result['start'],
                    "Stop": result['end'],
                    "Description": desc_data["Description"],
                    "Other Name": desc_data["Other Name"],
                    "Keywords": desc_data["Keywords"]
                })
    
    # Create and save DataFrame
    df = pd.DataFrame(output_rows)
    df_sorted = df.sort_values(by=["Start", "Stop"])
    df_sorted.to_csv(output_csv, index=False)
    
    # Find and report multiple matches
    multi_matches = find_multi_matches(output_rows)
    
    if multi_matches:
        print("\nWARNING! These genes have multiple possible Arabidopsis thaliana orthologs.")
        print("View the gene trees for these orthogroups to assess these genes further.\n")
        
        print("Multiple match summary:")
        for gene, info in multi_matches.items():
            print(f"Target gene: {gene}")
            print(f"Orthogroup: {info['orthogroup']}")
            print(f"Arabidopsis models: {', '.join(sorted(info['models']))}")
            print()
            
        # Write warning file
        output_dir = os.path.dirname(output_csv)
        write_warning_file(multi_matches, output_dir)
        print(f"\nWarning details written to: {os.path.join(output_dir, 'genes_with_warnings.csv')}")
    
    print(f"\nSuccessfully saved {len(df_sorted)} entries to {output_csv}")

def main():
    parser = argparse.ArgumentParser(description="Find orthologs in a specific region and parse the results.")
    
    # Combined arguments
    parser.add_argument("-f", "--fasta", required=True, help="Peptide fasta file.")
    parser.add_argument("-g", "--gtf", required=True, help="GTF file.")
    parser.add_argument("-c", "--chromosome", required=True, help="Chromosome.")
    parser.add_argument("-s", "--start", type=int, default=None, help="Start position.")
    parser.add_argument("-e", "--end", type=int, default=None, help="End position.")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file for the final parsed results.")
    parser.add_argument("-v", "--verbose_output", required=True, 
                        help="Directory containing 'target_matched_to_model_orthologs.csv' and 'model_locus_matched_to_orthogroup_output.csv'.")
    parser.add_argument("--fasta_part_index", type=int, default=2, help="Index of the part in the FASTA header (1-based index). Default is 2.")
    parser.add_argument("--fasta_separator", type=str, default=" ", help="Separator used in the FASTA header. Default is space (' ').")
    parser.add_argument("--gtf_column_type", type=str, default="transcript", help="Type of column to extract information from in the GTF file. Default is 'transcript'.")
    parser.add_argument("--gtf_pattern", type=str, default="(.*)", help="Regular expression pattern to extract the ID from the GTF file. Default searches anything in Column 8 (0-based index) of the gtf.")
    parser.add_argument("--postfix", type=str, default="", help="Custom post-fix to add to identifiers extracted from GTF file.")
    
    args = parser.parse_args()
    
    if args.start is None and args.end is not None:
        raise ValueError("Start position must be specified if end position is specified.")
    elif args.end is None and args.start is not None:
        raise ValueError("End position must be specified if start position is specified.")
    
    print("\n=== ORTHOLOG FINDER AND PARSER ===")
    print(f"Input FASTA file: {args.fasta}")
    print(f"Input GTF file: {args.gtf}")
    print(f"Chromosome: {args.chromosome}")
    print(f"Region: {f'From {args.start} to {args.end}' if args.start and args.end else 'All positions'}")
    print(f"Output file: {args.output}")
    print(f"Descriptions directory: {args.verbose_output}")
    
    # Check if verbose_output directory exists
    if not os.path.isdir(args.verbose_output):
        print(f"Warning: Directory {args.verbose_output} does not exist. Creating it...")
        os.makedirs(args.verbose_output, exist_ok=True)
    
    # Run find_orthologs_in_region function without writing to a file
    results = genes_in_range(
        args.fasta, args.gtf, args.chromosome, args.start, args.end, 
        args.verbose_output, args.fasta_part_index, args.fasta_separator, 
        args.gtf_column_type, args.gtf_pattern, args.postfix
    )
    
    # If results found, directly parse them and write to the output file
    if results:
        parse_protein_data(results, args.verbose_output, args.output)
        print("\n=== PROCESS COMPLETED SUCCESSFULLY ===")
    else:
        print("\n=== PROCESS COMPLETED - NO ORTHOLOGS FOUND ===")

if __name__ == "__main__":
    main()