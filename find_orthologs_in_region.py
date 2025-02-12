import argparse
import re
import os
import csv
from Bio import SeqIO

def read_fasta(file: str, part_index: int, separator: str = " ") -> dict:
    """Reads a FASTA file and returns a dictionary of SeqRecord objects, keyed by the specified part of the header."""
    seq_records = {}
    with open(file, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            header_parts = record.description.split(separator)
            if 1 <= part_index <= len(header_parts):
                part = header_parts[part_index - 1].strip()
                seq_records[part] = record
                print("Header:", record.description)  # Print the header
    return seq_records

def read_gtf(file: str, column_type: str, pattern: str, postfix: str = "") -> dict:
    """Reads a GTF file and returns a dictionary of gene information, keyed by specified ID type."""
    gtf_data = {}
    regex = re.compile(pattern)
    with open(file, "r") as gtf_file:
        for line in gtf_file:
            if not line.startswith("#"):
                columns = line.strip().split("\t")
                if columns[2] == column_type:
                    match = regex.search(columns[8])
                    if match is not None:
                        identifier = match.group(1) + postfix
                        gtf_data[identifier] = {"chr": columns[0], "start": int(columns[3]), "end": int(columns[4])}
    return gtf_data

def read_csv(file: str) -> dict:
    """Reads a CSV file and returns a dictionary of gene descriptions, keyed by protein ID."""
    descriptions = {}
    with open(file, "r") as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            target_organisms = row["Target Organism"].strip('"').split(", ")
            for organism in target_organisms:
                protein_id = organism
                gene_info = {"Locus": row["Locus"], "Description": row["Description"]}
                descriptions[protein_id] = gene_info
    return descriptions

def genes_in_range(fasta_file: str, gtf_file: str, chromosome: str, start=None, end=None, output_file=None, csv_file=None,
                   fasta_part_index=1, fasta_separator=" ", gtf_column_type="CDS", gtf_pattern="", postfix=""):
    """Finds all proteins located within a specified range on a specified chromosome."""
    seq_records = read_fasta(fasta_file, fasta_part_index, fasta_separator)
    gtf_data = read_gtf(gtf_file, gtf_column_type, gtf_pattern, postfix)
    descriptions = read_csv(csv_file) if csv_file else {}

    # Print contents of seq_records
    print("Contents of seq_records:")
    for key, value in seq_records.items():
        print(key, value)

    # Print contents of gtf_data
    print("\nContents of gtf_data:")
    for key, value in gtf_data.items():
        print(key, value)

    # Print contents of descriptions
    print("\nContents of descriptions:")
    for key, value in descriptions.items():
        print(key, value)

    num_orthologs_found = 0

    with open(output_file, "w") as output:
        for protein_id, record in seq_records.items():
            gene_info = gtf_data.get(protein_id)
            if gene_info:
                if gene_info["chr"] == chromosome:
                    if start is None or gene_info["start"] >= start:
                        if end is None or gene_info["end"] <= end:
                            description = f" (Description of Model Ortholog: {descriptions.get(protein_id, 'N/A')})" if csv_file else ""
                            message = f"Protein {protein_id} ({record.description}) is located within the specified range at position {gene_info['start']} to {gene_info['end']}{description}.\n"
                            print(message)
                            output.write(message)
                            num_orthologs_found += 1
            else:
                message = f"Protein {protein_id} not found in the gtf file.\n"
                print(message)
                output.write(message)

    if num_orthologs_found == 0:
        print("No orthologs were found in the specified region of interest.")
        if output_file:
            os.remove(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find proteins in a specific range on a specific chromosome.")
    parser.add_argument("-f", "--fasta", required=True, help="Peptide fasta file.")
    parser.add_argument("-g", "--gtf", required=True, help="GTF file.")
    parser.add_argument("-c", "--chromosome", required=True, help="Chromosome.")
    parser.add_argument("-s", "--start", type=int, default=None, help="Start position.")
    parser.add_argument("-e", "--end", type=int, default=None, help="End position.")
    parser.add_argument("-o", "--output", required=True, help="Output text file.")
    parser.add_argument("--find_descriptions", type=str, help="CSV file with gene descriptions.", default=None)
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

    genes_in_range(args.fasta, args.gtf, args.chromosome, args.start, args.end, args.output, args.find_descriptions,
                   args.fasta_part_index, args.fasta_separator, args.gtf_column_type, args.gtf_pattern, args.postfix)
