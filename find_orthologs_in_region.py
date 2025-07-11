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
            if separator == " ":
                header_parts = record.description.split()
            else:
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

def genes_in_range(fasta_file: str, gtf_file: str, chromosome: str, start_coords=None, end_coords=None, complement_on=False, full_chr_on=False, descriptions_directory=None,
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

    # Check if regions were provided
    sorted_coords_list = []
    if start_coords is not None:       # Main function checks that starting and ending coordinates are provided in reasonable pairs
        # Collect non-overlapping region coordinates into a list, sorted from low to high coordinates
        coords_list = []
        n = 0
        for start in start_coords:
            coords = [start, end_coords[n]]
            coords_list.append(coords)
            n += 1
        sorted_coords_list = sorted(coords_list, key=lambda x: x[0])    # List of non-overlapping regions sorted by position on chromosome from lowest to highest coordinates
    
    print(f"\n--- STEP 2: Finding orthologs in chromosome {chromosome} {f'in the following region(s): {sorted_coords_list}' if start_coords and end_coords else ''} ---")
    results_dict = {}

    if start_coords is not None:
        results_region = []
        num_orthologs_found_in_region = 0
    if complement_on == True:
        num_orthologs_found_in_complement = 0
        results_complement = []
    if full_chr_on == True or start_coords is None:
        num_orthologs_found_on_chr = 0
        results_chr = []

    for protein_id, record in seq_records.items():
        gene_info = gtf_data.get(protein_id)
        if gene_info:
            if gene_info["chr"] == chromosome:
                found = False
                if start_coords is not None:
                    for [start, end] in sorted_coords_list:
                        if gene_info["start"] >= start and gene_info["end"] <= end:
                            found = True
                            description = descriptions.get(protein_id, 'N/A') if descriptions else "N/A"
                            # message = f"Protein {protein_id} ({record.description}) is located within the specified range at position {gene_info['start']} to {gene_info['end']}."
                            # print(message)
                            # Store result for later parsing
                            results_region.append({
                                'protein_id': protein_id,
                                'description': record.description,
                                'start': gene_info['start'],
                                'end': gene_info['end'],
                                'model_description': description
                            })
                            num_orthologs_found_in_region += 1
                            break
                if complement_on == True and found == False:
                    description = descriptions.get(protein_id, 'N/A') if descriptions else "N/A"
                    # message = f"Protein {protein_id} ({record.description}) is located within the specified range at position {gene_info['start']} to {gene_info['end']}."
                    # print(message)
                    # Store result for later parsing
                    results_complement.append({
                        'protein_id': protein_id,
                        'description': record.description,
                        'start': gene_info['start'],
                        'end': gene_info['end'],
                        'model_description': description
                    })
                    num_orthologs_found_in_complement += 1
                if start_coords is None or full_chr_on == True:
                    description = descriptions.get(protein_id, 'N/A') if descriptions else "N/A"
                    # message = f"Protein {protein_id} ({record.description}) is located within the specified range at position {gene_info['start']} to {gene_info['end']}."
                    # print(message)
                    # Store result for later parsing
                    results_chr.append({
                        'protein_id': protein_id,
                        'description': record.description,
                        'start': gene_info['start'],
                        'end': gene_info['end'],
                        'model_description': description
                    })
                    num_orthologs_found_on_chr += 1
        # else:
        #     message = f"Protein {protein_id} not found in the gtf file."
        #     print(message)

    if start_coords is None:
        if num_orthologs_found_on_chr == 0:
            print("No orthologs were found on the chromosome of interest.")
            return None
    else:
        if num_orthologs_found_in_region == 0:
            print("No orthologs were found in the specified region of interest.")

    # Add results for each type of region searched to results dictionary
    if start_coords is not None:
        results_dict['region'] = results_region
        print(f"Found {num_orthologs_found_in_region} orthologs in the specified region.")
    if complement_on == True:
        results_dict['complement'] = results_complement
        print(f"Found {num_orthologs_found_in_complement} orthologs outside specified region.")
    if start_coords is None or full_chr_on == True:
        results_dict['full_chr'] = results_chr
        print(f"Found {num_orthologs_found_on_chr} orthologs on the chromosome")

    return results_dict

def find_multi_matches(results_dict: dict) -> dict:
    """Find genes with multiple Arabidopsis matches."""
    gene_matches = {}
    
    # Group results by target gene
    if 'full_chr' in results_dict:
        results_all = results_dict['full_chr']
    elif 'complement' in results_dict:
        results_all = results_dict['region'] + results_dict['complement']
    else:
        results_all = results_dict['region']
    
    for _, result in results_all.iterrows():
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

def write_warning_file(multi_matches: dict, warning_file_name: str):
    """Write warning file for genes with multiple matches."""
    
    with open(warning_file_name, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['Target_Gene', 'Orthogroup', 'Arabidopsis_Models'])
        
        for gene, info in multi_matches.items():
            writer.writerow([
                gene,
                info['orthogroup'],
                ';'.join(sorted(info['models']))
            ])

def parse_protein_data(results_dict, descriptions_directory, input_file: str, chromosome: str, start_coords=None, label="BEG", region_name="target_region", complement_name="complement", chr_name=""):
    """Parse protein data directly from the results of genes_in_range and merge with locus data."""
    if start_coords is not None:
        if not results_dict['region']:
            print("No results to parse.")
            return
    else:
        if not results_dict['full_chr']:
            print("No results to parse.")
            return
        
    print("\n--- STEP 3: Parsing ortholog data ---")
    
    # Load model descriptions and target mappings
    model_desc_file = os.path.join(descriptions_directory, "model_locus_matched_to_orthogroup_output.csv")
    target_map_file = os.path.join(descriptions_directory, "target_matched_to_model_orthologs.csv")
    
    model_descriptions = read_model_descriptions(model_desc_file)
    target_mappings = read_target_mappings(target_map_file)

    # Create output names list for naming files after extracting orthologs
    output_names = {}
    if 'region' in results_dict.keys():
        output_names['region'] = region_name
    if 'complement' in results_dict.keys():
        output_names['complement'] = complement_name
    
    # Extracting ortholog information from results
    results_df_dict = {}
    multi_matched_dict = {}
    for type in results_dict.keys():
        results = results_dict[type]
        output_rows = []    # Initialize output data
        
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
        results_df_dict[type] = df_sorted

        # Save to CSV
        # Collect directory for output from input file
        output_csv_file_path = f"{os.path.dirname(input_file)}"

        # Create basename for region file
        #if start_coords != None and (type == 'region' or type == 'complement'):
        if label == "BEG":
            if type != "full_chr":
                output = f"{output_names[type]}-"
            else:
                output = ""
            if chr_name != "":
                name = f"{output}{chr_name}-{chromosome}-{os.path.splitext(os.path.basename(input_file))[0]}"
            else:
                name = f"{output}{chromosome}-{os.path.splitext(os.path.basename(input_file))[0]}"
        else:
            if output != "":
                output = f"-{output_names[type]}"
            else:
                output = ""
            if chr_name != "":
                name = f"{os.path.splitext(os.path.basename(input_file))[0].replace('-orthologs', '')}-{chromosome}-{chr_name}{output}-orthologs"
            else:
                name = f"{os.path.splitext(os.path.basename(input_file))[0].replace('-orthologs', '')}-{chromosome}{output}-orthologs"
            
        file = os.path.join(output_csv_file_path, f"{name}.csv")
        df_sorted.to_csv(os.path.join(output_csv_file_path, f"{name}.csv"), index=False)
        print(f"\nSuccessfully saved {len(df_sorted)} entries to {file}")
    
    # Find and report multiple matches
    multi_matches = find_multi_matches(results_df_dict)

    if multi_matches:
        print("\nWARNING! These genes have multiple possible Arabidopsis thaliana orthologs.")
        print("View the gene trees for these orthogroups to assess these genes further.\n")
        
        print("Multiple match summary:")
        for gene, info in multi_matches.items():
            print(f"Target gene: {gene}")
            print(f"Orthogroup: {info['orthogroup']}")
            print(f"Arabidopsis models: {', '.join(sorted(info['models']))}")
            print()

        multi_matched_dict[type] = multi_matches

def main():
    parser = argparse.ArgumentParser(description="Find orthologs in a specific region and parse the results.")
    
    # Combined arguments
    parser.add_argument("-f", "--fasta", required=True, help="Peptide fasta file.")
    parser.add_argument("-g", "--gtf", required=True, help="GTF file.")
    parser.add_argument("-c", "--chromosome", required=True, help="Chromosome.")
    parser.add_argument("-s", "--start", type=int, nargs="+", default=None, help="Start position(s).")
    parser.add_argument("-e", "--end", type=int, nargs="+", default=None, help="End position(s).")
    parser.add_argument("-v", "--verbose_output", required=True, 
                        help="Directory containing 'target_matched_to_model_orthologs.csv' and 'model_locus_matched_to_orthogroup_output.csv'.")
    parser.add_argument("--region_name", type=str, default="target_region", help="Replaces 'target_region' with custom label in region output file names. For example: SDR (for sex determining region)")
    parser.add_argument("--complement_on", action='store_true', help="Generates an additional file of all orthologs located outside of region, when specified.")
    parser.add_argument("--complement_name", type=str, default="complement", help="Replaces 'complement' with custom label in complement output file names. For example: PAR (for pseudoautosomal region)")
    parser.add_argument("--full_chr_on", action='store_true', help="Generates an additional file of all orthologs located on the chromosome *when region is specified.*")
    parser.add_argument("--chr_name", type=str, default="", help="Adds custom chromosome name in all output file names. For example: chr12")
    parser.add_argument("-l", "--label", default="BEG", help="Assigns whether labels go at beginning or end of output file names using two options: BEG or END")
    parser.add_argument("--fasta_part_index", type=int, default=2, help="Index of the part in the FASTA header (1-based index). Default is 2.")
    parser.add_argument("--fasta_separator", type=str, default=" ", help="Separator used in the FASTA header. Default is space (' ').")
    parser.add_argument("--gtf_column_type", type=str, default="transcript", help="Type of column to extract information from in the GTF file. Default is 'transcript'.")
    parser.add_argument("--gtf_pattern", type=str, default="(.*)", help="Regular expression pattern to extract the ID from the GTF file. Default searches anything in Column 8 (0-based index) of the gtf.")
    parser.add_argument("--postfix", type=str, default="", help="Custom post-fix to add to identifiers extracted from GTF file.")
    
    args = parser.parse_args()
    
    # Check whether start and end coordinates are provided
    if args.start == None or args.end == None:
        if not args.start == None and args.end == None:
            parser.error("Error: Starting and ending coordinates must be provided in pairs if they are used.\
                                 However, coordinates are not required. Do not use the optional --start START or --end END flags if you wish to check across an entire chromosome.")
        else:
            # Check that flags that are supposed to be used only when regions are provided are not on
            # Turn flags off and print warnings if flags are turned on
            if args.region_name != "target_region":
                print("WARNING: The --region_name flag should only be used when a region is provided. No region was provided, so the entire chromosome was searched.")
            if args.complement_on == True:
                args.complement_on = False
                print("Warning: The --complement_on flag should only be used when a region is provided. No region was provided.")
            if args.full_chr_on == True:
                args.full_chr_on = False
                print("Warning: The --full_chr_on flag should only be used when a region is provided. No region was provided.")
    # If start and end coordinates are provided, check that they are reasonable.
    else:
        # Check that there are an equal number of start and end coordinates provided
        if not (len(args.start) == len(args.end)):
            parser.error("The same number of start and end positions must be provided.")
            print(len(args.start) + " start positions and " + len(args.end) + " end positions were indicated.")
        # Check that the bounds are appropriate
        else:
            # Check that starting coordinate is less than ending coordinate
            wrong_coord = False     # Used to raise an error for the coordinates and stop the program
            coords_list = []
            n = 0
            for start in args.start:
                if start < args.end[n]:
                    coords = [start, args.end[n]]
                    coords_list.append(coords)
                else:
                    wrong_coord = True
                    print(f"Error: Starting coordinate {start} is not less than the ending coordinate {args.end[n]}.")
                    print(f"       Please adjust the bounds so that the starting coordinate is less than the ending coordinate.")
                n += 1
            sorted_coords_list = sorted(coords_list, key=lambda x: x[0])
            
            if wrong_coord == True:
                parser.error("Some of the regions have matching start and end bounds. See the above warnings.\nPlease adjust your starting and ending coordinates so that they do not match.")

            # Check that the bounds are not overlapping
            overlaps = False     # Used to raise an error if overlapping regions are detected and stop the program
            n = 0
            for n in range(len(sorted_coords_list)-1):
                start_1 = sorted_coords_list[n][0]
                end_1 = sorted_coords_list[n][1]
                start_2 = sorted_coords_list[n+1][0]
                end_2 = sorted_coords_list[n+1][1]

                wrong_start = False

                # Check starting position
                if start_2 < end_1 or start_2 == start_1:
                    wrong_start = True
                    print(f"Error: Regions {start_1} to {end_1} and {start_2} to {end_2} are overlapping.")
                    print(f"       The start coordinate {start_2} is within the range of {start_1} to {end_1}.")
                    overlaps = True
            if overlaps == True:
                parser.error("Some of the regions provided are overlapping. See the above warnings.\nPlease adjust your starting and ending coordinates so that they do not overlap.")
            # Check whether additional region flags (--complement_on and --full_chr_on) are on if naming flags are used
            # If so, print warning message and turn on associated region flag(s).
            else:
                if args.complement_name != "complement" and args.complement_on == False:
                    print("Warning: The complementary region was named, but the complement flag was not turned on by user.\n\
                                    The --complement_on flag was automatically turned on.")
                    args.complement_on = True

    if args.label != "BEG" and args.label != "END":
        raise ValueError("Option --label can only have one of two values: BEG or END (all caps).")
    
    print("\n=== ORTHOLOG FINDER AND PARSER ===")
    print(f"Input FASTA file: {args.fasta}")
    print(f"Input GTF file: {args.gtf}")
    print(f"Chromosome: {args.chromosome}")
    print(f"Region(s): {f'{sorted_coords_list}' if args.start and args.end else 'All positions'}")
    print(f"Descriptions directory: {args.verbose_output}")
    
    # Check if verbose_output directory exists
    if not os.path.isdir(args.verbose_output):
        print(f"Warning: Directory {args.verbose_output} does not exist. Creating it...")
        os.makedirs(args.verbose_output, exist_ok=True)
    
    # Run find_orthologs_in_region function without writing to a file
    results = genes_in_range(
        args.fasta, args.gtf, args.chromosome, args.start, args.end, 
        args.complement_on, args.full_chr_on,
        args.verbose_output, args.fasta_part_index, args.fasta_separator, 
        args.gtf_column_type, args.gtf_pattern, args.postfix
    )
    
    # If results found, directly parse them and write to the output file
    if results:
        parse_protein_data(results, args.verbose_output, args.fasta, args.chromosome,
                           args.start, args.label, args.region_name, args.complement_name, args.chr_name)
        print("\n=== PROCESS COMPLETED SUCCESSFULLY ===")
    else:
        print("\n=== PROCESS COMPLETED - NO ORTHOLOGS FOUND ===")

if __name__ == "__main__":
    main()