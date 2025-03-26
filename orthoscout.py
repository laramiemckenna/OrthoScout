import time
import argparse
import os
from typing import List, Dict

from utils import (
    read_input_files,
    match_locus,
    read_tsv,
    match_orthogroups,
    write_to_csv,
    extract_gene_ids,
    match_sequences,
    write_fasta,
)

"""
Description of orthoscount.py:

This script was built to match the orthologs of the model species genes with those of the target species --
it is the "funnel" of the pipeline that utilizes orthology inference via OrthoFinder2.
It utilizes the OrthoFinder TSV file, the subsetted model species pattern CSV file,
and information about your target species to give you a final target species peptide file for further analysis.

"""

def parse_arguments() -> argparse.Namespace:
    """
    This function allows for command-line inputs by utilizing argparse to parse command line arguments 
    and returning them as a Namespace object.

    Returns:
        Namespace: An object containing the parsed command line arguments.
    """
    # Create an ArgumentParser object to parse command line arguments
    parser = argparse.ArgumentParser(
        description="Process OrthoFinder TSV, model species gene information CSV, and target protein FASTA files."
    )
    # Add a required argument to specify the input folder path
    parser.add_argument(
    "-i",
    "--input-folder",
    required=True,
    help="Path to the input folder containing the OrthoFinder TSV, model species gene information CSV, and target protein FASTA files",
    dest="input_folder"
    )   
    # Add optional argument to specify the output file path
    parser.add_argument(
        "-o",
        "--output",
        default="target_species_genes.pep.fa",
        help="Output FASTA file path",
    )
    # Add required argument to specify the target column name in the TSV file
    parser.add_argument(
        "-t", "--target-column", required=True, help="Target column name in the TSV file"
    )
    # Add optional argument to specify whether to save intermediate output files in a verbose_output folder
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Save intermediate output files in the 'verbose_output' folder"
    )

    # Parse the command line arguments and return them as a Namespace object
    return parser.parse_args()


def find_files(input_folder: str, extensions: List[List[str]]) -> Dict[str, str]:
    """
    This functions searches for the specified input folder for files with the specified extensions and 
    returns a dictionary containing the file paths for each extension group.

    Args:
        input_folder (str): The path to the input folder to search.
        extensions (list): A list of lists, with each inner list containing the file extensions to search for.

    Returns:
        dict: A dictionary containing the file paths for each extension group.
    """
    # Initialize an empty dictionary to store found files
    found_files = {}

    # Iterate over each group of file extensions
    for ext_group in extensions:

        # Iterate over each file extension in the group until a file is found
        for ext in ext_group:
            for file in os.listdir(input_folder):

                # If a file with the extension is found, add its path to the dictionary and break out of the loop
                if file.endswith(ext):
                    found_files[ext] = os.path.join(input_folder, file)
                    break

            # If a file with the extension was found, break out of the loop for the current group of extensions
            if ext in found_files:
                break

        # If no file was found for any of the extensions in the current group, raise a FileNotFoundError
        else:
            exts = ', '.join(ext_group)
            raise FileNotFoundError(
                f"No file with any of the following extensions was found in the specified folder: {exts}"
            )

    # Return the dictionary containing the file paths for each extension group
    return found_files


def main() -> None:
    """
    This function is the central functionality of the pipeline, 
    which helps users match the orthogroups of their target species genes 
    with those of the model species genes from the mutant phenotype database.
    It parses the command-line arguments, matches the locus tags of the target species 
    with their corresponding Orthogroup, matches your target species gene IDs with those Orthogroups,
    then, matches those gene IDs to a peptide fasta file and generates a sub-setted output for further analysis.

    Args
        None.

    Returns:
        None.

    If the --verbose flag is specified in the command line, 
    then this the function will save intermediate outputs as CSV files in the 'verbose_output' directory 
    that can later be used for debugging or double-checking.
    """
    # parse command line arguments
    args = parse_arguments()

    # get input folder path from arguments
    input_folder = args.input_folder

    # check if the input folder exists and is a directory
    if not os.path.isdir(input_folder):
        raise NotADirectoryError(
            f"The specified input folder does not exist or is not a directory: {input_folder}"
        )

    # get output file path, target species column name, and verbosity setting from arguments
    output_file = args.output
    target_species_column_name = args.target_column
    verbose = args.verbose

    # create the 'verbose_output' directory if --verbose is specified
    if verbose:
        os.makedirs("verbose_output", exist_ok=True)

    # find input files with the specified extensions in the input folder
    input_files = find_files(input_folder, [['.tsv'], ['.csv'], ['pep.fa', 'pep.faa', 'pep.fasta']])
    orthofinder_tsv = input_files['.tsv']
    model_gene_info_csv = input_files['.csv']
    input_fasta_ext = next(ext for ext in ['pep.fa', 'pep.faa', 'pep.fasta'] if ext in input_files)
    input_fasta_file = input_files[input_fasta_ext]

    # read input files as data frames
    start_time = time.perf_counter()
    orthofinder_df, model_gene_info_df = read_input_files(orthofinder_tsv, model_gene_info_csv)
    end_time = time.perf_counter()
    tot_time = end_time - start_time
    print(f"Time spent reading input files as data frames: {tot_time} seconds.")

    # log a message to indicate successful file reading
    print("Files read successfully.")

    # match the locus tags of the model species genes to those in the OrthoFinder TSV file
    start_time = time.perf_counter()
    matched_df = match_locus(orthofinder_df, model_gene_info_df)
    end_time = time.perf_counter()
    tot_time = end_time - start_time
    print(f"Time spent matching locus tags of the model species to those in the OrthoFinder TSV file: {tot_time} seconds.")
    print("Locus matching completed.")

    # determine the path for the matched output file and save the matched data frame as CSV
    if verbose:
        # matched_df_path = "verbose_output/model_locus_matched_to_orthogroup_output.csv"
        path = "verbose_output/model_locus_matched_to_orthogroup_output.csv"
        matched_df.to_csv(path, index=False, mode='w')
        print(f"Matched output saved to {path}.")


    # read the OrthoFinder TSV file as a dictionary
    start_time = time.perf_counter()
    orthofinder_data = read_tsv(orthofinder_tsv)
    end_time = time.perf_counter()
    tot_time = end_time - start_time
    print(f"Time spent reading OrthoFinder TSV as a dictionary: {tot_time} seconds.")

    # Convert the matched DataFrame to a dictionary
    start_time = time.perf_counter()
    matched_model_locus_data = matched_df.to_dict()
    end_time = time.perf_counter()
    tot_time = end_time - start_time
    print(f"Time spent converting the matched data frame to a dictionary: {tot_time} seconds.")

    # match the orthogroups of the target species genes with those of the model species genes
    print("target_species_column_name:", target_species_column_name)
    start_time = time.perf_counter()
    target_matched_results = match_orthogroups(orthofinder_data, matched_model_locus_data, target_species_column_name)
    end_time = time.perf_counter()
    tot_time = end_time - start_time
    print(f"Time spent matching orthogroups: {tot_time} seconds.")
    print("Orthogroups matched.")

    # save the target match data frame as a CSV file if verbosity is enabled
    if verbose:
        path ="verbose_output/target_matched_to_model_orthologs.csv"
        write_to_csv(path, target_matched_results)
        print(f"Target match saved to {path}.")

    # extract the gene IDs from the target match data frame
    start_time = time.perf_counter()
    gene_ids = extract_gene_ids(target_matched_results)
    end_time = time.perf_counter()
    tot_time = end_time - start_time
    print(f"Time spent extracting gene IDs from the target match data frame: {tot_time} seconds.")

    # match the gene IDs with the sequences in the input FASTA file
    start_time = time.perf_counter()
    matched_target_sequences = match_sequences(input_fasta_file, gene_ids)
    end_time = time.perf_counter()
    tot_time = end_time - start_time
    print(f"Time spent matching the gene IDs with the FASTA sequences: {tot_time} seconds.")

    # check if any sequences were matched from the input FASTA file
    if not matched_target_sequences:
        raise ValueError("No sequences were matched from the input FASTA file.")

    # write the final output to a FASTA file
    write_fasta(output_file, matched_target_sequences)
    print(f"Final output saved to {output_file}.")

if __name__ == "__main__":
    main()
