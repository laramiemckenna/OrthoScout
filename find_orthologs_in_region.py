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

def parse_protein_data(results, descriptions_directory, output_csv):
    """Parse protein data directly from the results of genes_in_range and merge with locus data."""
    if not results:
        print("No results to parse.")
        return
        
    print("\n--- STEP 3: Parsing ortholog data ---")
    
    # Initialize lists to store parsed data
    proteins = []
    orthogroups = []
    starts = []
    stops = []
    descriptions = []
    other_names = []
    keywords = []
    
    # Regular expression patterns for extracting data
    orthogroup_pattern = re.compile(r"(OG\d+)")
    description_pattern = re.compile(r"Description: (.*?)(?:, Other Name|, Keywords)")
    other_name_pattern = re.compile(r"Other Name: (.*?)(?:, Keywords|})")
    keywords_pattern = re.compile(r"Keywords: (.*?)}")
    
    print("Extracting ortholog information from results...")
    # Process each result
    for result in results:
        protein_id = result['protein_id']
        full_desc = result['description']
        
        # Extract orthogroup
        orthogroup_match = orthogroup_pattern.search(full_desc)
        orthogroup = orthogroup_match.group(1) if orthogroup_match else "N/A"
        
        # Extract other metadata from the description
        description_match = description_pattern.search(str(result['model_description']))
        other_name_match = other_name_pattern.search(str(result['model_description']))
        keywords_match = keywords_pattern.search(str(result['model_description']))
        
        # Append data
        proteins.append(protein_id)
        orthogroups.append(orthogroup)
        starts.append(result['start'])
        stops.append(result['end'])
        descriptions.append(description_match.group(1) if description_match else "N/A")
        other_names.append(other_name_match.group(1) if other_name_match else "N/A")
        keywords.append(keywords_match.group(1) if keywords_match else "N/A")
    
    # Create a DataFrame
    print("Creating dataframe with ortholog information...")
    df = pd.DataFrame({
        "Protein": proteins,
        "Orthogroup": orthogroups,
        "Start": starts,
        "Stop": stops,
        "Description": descriptions,
        "Other Name": other_names,
        "Keywords": keywords
    })
    
    # Sort by Start position
    df_sorted = df.sort_values(by="Start")
    
    # Path to the locus data file in the specified directory
    locus_file = os.path.join(descriptions_directory, "model_locus_matched_to_orthogroup_output.csv")
    
    print("\n--- STEP 4: Merging with locus information ---")
    # Load the locus data if file exists
    if os.path.exists(locus_file):
        print(f"Reading locus file: {locus_file}")
        locus_df = pd.read_csv(locus_file)
        print(f"Found {len(locus_df)} locus entries")
        
        # Merge with locus data
        print("Merging ortholog data with locus information...")
        merged_df = df_sorted.merge(locus_df, on="Orthogroup", how="left")
        
        # Rename columns
        merged_df.rename(columns={"locus": "Arabidopsis Gene Model"}, inplace=True)
        
        # Reorder columns
        merged_df = merged_df[["Protein", "Orthogroup", "Arabidopsis Gene Model", "Start", "Stop", "Description", "Other Name", "Keywords"]]
    else:
        print(f"Warning: Locus file not found at {locus_file}. Proceeding without locus information.")
        merged_df = df_sorted
        # Add empty column for Arabidopsis Gene Model
        merged_df["Arabidopsis Gene Model"] = "N/A"
        # Reorder columns
        merged_df = merged_df[["Protein", "Orthogroup", "Arabidopsis Gene Model", "Start", "Stop", "Description", "Other Name", "Keywords"]]
    
    # Save to CSV
    print(f"\n--- STEP 5: Saving results to {output_csv} ---")
    merged_df.to_csv(output_csv, index=False)
    print(f"Successfully saved {len(merged_df)} orthologs to {output_csv}")

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