import argparse
import subprocess

def run_blast(genomic_fasta: str, query_peptide_fasta: str, output_file: str, subject_num: str, seq_range: str):
    """Runs BLAST on the specified input files and saves relevant sequences to an output file.
    
    Args:
        genomic_fasta (str): The path to the input genomic fasta file.
        query_peptide_fasta (str): The path to the input query peptide fasta file.
        output_file (str): The path to the output file.
        subject_num (str): The subject/chromosome identifier.
        seq_range (str): The sequence range in the format "start-end", where "start" and "end" are integers.
    """
    
    # Create a directory for the database
    db_dir = "db_files"
    subprocess.run(["mkdir", db_dir])
    
    # Run makeblastdb to create a database
    db_name = output_file + "_db"
    db_path = db_dir + "/" + db_name
    makeblastdb_command = ["makeblastdb", "-in", genomic_fasta, "-dbtype", "nucl", "-out", db_path, "-parse_seqids"]
    subprocess.run(makeblastdb_command)
    
    # Run tblastn to search the database
    output_file = output_file + ".txt"
    tblastn_command = ["tblastn", "-query", query_peptide_fasta, "-db", db_path, "-out", output_file, "-outfmt", "7"]
    subprocess.run(tblastn_command)
    
    # Parse the BLAST output file and save relevant sequences
    if seq_range:
        start_pos, end_pos = map(int, seq_range.split('-'))
    else:
        start_pos, end_pos = 0, float('inf')
    output_name = subject_num + "_" + seq_range + ".txt"
    with open(output_name, 'w') as out_file, open(output_file) as blast_output:
        for line in blast_output:
            if line.startswith('#'):
                continue
            fields = line.strip().split("\t")
            if subject_num == fields[1] and int(fields[8]) >= start_pos and int(fields[9]) <= end_pos:
                out_file.write(line)
    
    print("Done! Result saved in " + output_name)


if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Run BLAST and parse output file.')
    parser.add_argument('-g', '--genomic_fasta', help='Path to the input genomic fasta file', required=True)
    parser.add_argument('-q', '--query_peptide_fasta', help='Path to the input query peptide fasta file', required=True)
    parser.add_argument('-o', '--output_file', help='Path to the output file', required=True)
    parser.add_argument('-s', '--subject_num', help='Subject/chromosome identifier', required=True)
    parser.add_argument('-r', '--range', help='Sequence range in the format "start-end"', required=True)
    args = parser.parse_args()

    # Run BLAST
    run_blast(args.genomic_fasta, args.query_peptide_fasta, args.output_file, args.subject_num, args.range)
