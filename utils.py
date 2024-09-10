import pandas as pd
import csv
from typing import List, Dict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

"""
Description of utils.py:

This file contains functions that perform various operations on data relevant to the orthoscout.py script.
It includes functions to read and match data from OrthoFinder TSV file to the model gene information CSV file,
read and write the input TSV and CSV files, match specific target gene IDs in a CSV file to sequences
in a target peptide FASTA file, and write these matched sequences to a new FASTA file for downstream analysis.
It uses the Pandas and BioPython libraries to handle data and perform sequence analysis.
"""

def read_input_files(orthofinder_tsv: str, model_gene_info_csv: str) -> tuple:
    """
    This function reads the OrthoFinder TSV file and model species gene information CSV file generated bu search_and_extract.py into DataFrames.

    Args:
        orthofinder_tsv (str): The path to the OrthoFinder TSV file.
        model_gene_info_csv (str): The path to the model species gene information CSV file.

    Returns:
        tuple: A tuple containing the OrthoFinder pandas DataFrame and model species gene information pandas DataFrame.
    """
    # Load the OrthoFinder TSV file into a pandas DataFrame
    orthofinder_df = pd.read_csv(orthofinder_tsv, sep="\t", dtype={i: str for i in range(1, 7)})

    # Load the model species gene information CSV file into a DataFrame
    model_gene_info_df = pd.read_csv(model_gene_info_csv)

    # Return a tuple containing the OrthoFinder DataFrame and model species gene information DataFrame
    return orthofinder_df, model_gene_info_df


def match_locus(orthofinder_df: pd.DataFrame, model_gene_info_df: pd.DataFrame) -> pd.DataFrame:
    """
    This function searches for matching "locus" values from the model species gene information pandas DataFrame in in the OrthoFinder pandas DataFrame.

    Args:
        orthofinder_df (pd.DataFrame): The OrthoFinder DataFrame.
        model_gene_info_df (pd.DataFrame): The model species gene information DataFrame.

    Returns:
        pd.DataFrame: A new DataFrame containing the Orthogroup, original model locus ID, and the description of the phenotype associated with that locus..
    """
    # Create an empty list to store the matched rows
    matched_rows = []

    # Iterate over each row in the model species gene information DataFrame
    for index, row in model_gene_info_df.iterrows():
        locus = row['Gene_model']
        description = row['Description']

        # Iterate over each row in the OrthoFinder DataFrame
        for _, ortho_row in orthofinder_df.iterrows():
            orthogroup = ortho_row['Orthogroup']
            proteins = ', '.join([str(x) for x in ortho_row[1:] if not pd.isnull(x)])

            # Check if the locus value is found anywhere in the protein list
            # If a match is found, add a new row to the matched_rows list
            if locus in proteins:
                matched_rows.append({'Orthogroup': orthogroup, 'locus': locus, 'description': description})
                break

    # Create a new DataFrame from the matched_rows list
    matched_df = pd.DataFrame(matched_rows)

    # Return the new DataFrame
    return matched_df


def read_tsv(file_path: str) -> List:
    """
    This function reads the OrthoFinder TSV file and returns its contents (specifically the rows) as a list of lists.

    Args:
        file_path (str): The path to the TSV file.

    Returns:
        list: A list of lists, where each inner list is a row from the TSV file.
    """
    # Open the TSV file and read the contents
    with open(file_path, "r") as f:
        # Add the list of values to the rows list
        rows = [line.strip().split("\t") for line in f.readlines()]

    # Return the list of rows
    return rows


def match_orthogroups(tsv_data: List[List[str]], model_species_locus_tag_data: Dict, target_column: str) -> List[List[str]]:
    """
    This function matches IDs in a TSV file with corresponding information in a dictionary.

    Args:
        tsv_data (list): A list of lists representing a TSV file.
        model_species_locus_tag_data (dict): A dictionary representing model species locus tag data with 'Orthogroup', 'locus', and 'description' keys.
        target_column (str): The name of the column in the TSV file to match with the dictionary.

    Returns:
        list: A list of lists, where each inner list contains the matched information for each matched key in the form of (orthogroup_id, locus_id, description, target_column_value).
    """
    # Create a dictionary to store the information for each orthogroup in the CSV file
    orthogroups = {}
    for data_idx in model_species_locus_tag_data['Orthogroup']:
        group = model_species_locus_tag_data['Orthogroup'][data_idx]

        # create locus, description tuple
        locus = model_species_locus_tag_data['locus'][data_idx]
        description = model_species_locus_tag_data['description'][data_idx]
        locus_description = (locus, description)

        if group not in orthogroups:
            # If the orthogroup ID doesn't exist in the dictionary, add it and initialize it with a list containing the locus information
            orthogroups[group] = [locus_description]
        else:
            orthogroups[group].append(locus_description)

    # Get the index of the target column in the TSV data
    if target_column not in tsv_data[0]:
        print(f"Error: '{target_column}' not found in tsv_data headers.")
        return []

    target_index = tsv_data[0].index(target_column)

    # Create a dictionary to store the matched results
    matched_results = {}

    # Iterate over each row in the TSV data
    for row_idx, row in enumerate(tsv_data[1:], start=1):  # Start from index 1 to skip the header row
        orthogroup_id = row[0]
        if orthogroup_id in orthogroups:
            # Check if target_index is valid for this row
            if target_index < len(row):
                # For each matching orthogroup, iterate over its locus information in the CSV data
                for locus_description in orthogroups[orthogroup_id]:
                    # Create a key for the matched result using the orthogroup ID and the target column value
                    key = (orthogroup_id, row[target_index])
                    value = locus_description[0]
                    # If the key doesn't exist in the matched_results dictionary, add a new entry
                    if key not in matched_results:
                        matched_results[key] = [orthogroup_id] + list(locus_description) + [row[target_index]]
                    # If the key already exists, append the new locus ID to the existing entry
                    else:
                        matched_results[key][1] += f"; {value}"
            else:
                print(f"Warning: target_index ({target_index}) is out of range for row {row_idx}. Skipping this row.")
                continue  # Skip the current row and proceed to the next row

    # Convert the matched_results dictionary to a list of lists and return it
    return list(matched_results.values())



def write_to_csv(file_name: str, data: List[List[str]]) -> None:
    """
    This function writes data to a CSV file.

    Args:
        file_name (str): The name of the output CSV file.
        data (list): A list of lists containing the data to write to the CSV file.

    Returns:
        None
    """
    # Open the output CSV file in write mode
    # Set newline='' to prevent blank lines from being inserted between rows
    with open(file_name, 'w', newline='') as csv_file:
        # Create a CSV writer object
        csv_writer = csv.writer(csv_file)

        # Write the header row to the CSV file
        csv_writer.writerow(['Orthogroup', 'Locus', 'Description', 'Target Organism'])

        # Iterate over each row in the data list and write it to the CSV file
        for row in data:
            csv_writer.writerow(row)


def extract_gene_ids(csv_data: List[List[str]]) -> Dict[str, str]:
    """
    This function extracts gene IDs from the CSV data and returns a dictionary with gene IDs as keys
    and corresponding Orthogroups as values.

    Args:
        csv_data (list): A list of lists representing a CSV file.

    Returns:
        dict: A dictionary with gene IDs as keys and corresponding Orthogroups as values.
    """
    # Create an empty dictionary to store the gene IDs and corresponding Orthogroups
    gene_ids = {}

    # Iterate over each row in the CSV data, starting from the second row
    for row in csv_data[1:]:
        orthogroup = row[0]
        target_organisms = row[3].split(", ")
        # For each target organism in the row, extract the gene ID and add it to the gene_ids dictionary
        for organism in target_organisms:
            gene_id = organism.split("_")[-1]
            gene_ids[gene_id] = orthogroup

    # Return the gene_ids dictionary
    return gene_ids


def match_sequences(fasta_file: str, gene_ids: Dict[str, str]) -> List[SeqRecord]:
    """
    This function searches the input FASTA for sequences with matching gene IDs and returns a list
    of matched sequences with modified headers containing the corresponding Orthogroup information.

    Args:
        fasta_file (str): The path to the input FASTA file.
        gene_ids (dict): A dictionary with gene IDs as keys and corresponding Orthogroups as values.

    Returns:
        list: A list of matched sequences with modified headers.
    """
    # Create an empty list to store the matched sequences
    matched_sequences = []

    # Iterate over each sequence record in the input FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_id = record.id
        # If the gene ID matches one of the gene IDs in the gene_ids dictionary,
        # modify the record description to include the corresponding Orthogroup information
        if gene_id in gene_ids:
            record.description = f"{gene_ids[gene_id]} {record.description}"
            matched_sequences.append(record)

    # Return the list of matched sequences
    return matched_sequences


def write_fasta(file_name: str, sequences: List[SeqRecord]) -> None:
    """
    This function writes the matched sequences to a FASTA file in the standard format with an empty line after each sequence.

    Args:
        file_name (str): The name of the output FASTA file.
        sequences (list): A list of SeqRecord objects containing the matched sequences.

    Returns:
        None
    """
    # Open the output FASTA file in write mode
    with open(file_name, "w") as output_file:
        # Iterate over each sequence in the sequences list
        for seq in sequences:
            # Format the FASTA record header
            fasta_record = f">{seq.description}\n"
            # Split the sequence into lines of 80 characters or less
            sequence_parts = [str(seq.seq)[i:i+80] for i in range(0, len(seq.seq), 80)]
            # Join the sequence lines together with newline characters
            sequence_formatted = "\n".join(sequence_parts) + "\n"
            # Write the FASTA record to the output file, with an empty line after each sequence
            output_file.write(fasta_record + sequence_formatted + "\n")

