import os
import sqlite3
import csv
import re
from typing import List, Tuple, Optional
import argparse

"""
Description of search_and_extract.py

This script is built to search an SQLite3 database 
and provide a subset of that database by matching a pattern
and column specified in the commandline.
"""

def connect_to_database(database_path: str) -> Optional[sqlite3.Connection]:
    """
    This function establishes a connection to the SQLite database and returns the connection object.

    Args:
        database_path (str): The path to the SQLite database file.

    Returns:
        A connection to the SQLite database.
    """
    try:
        # Connect to the SQLite database
        conn = sqlite3.connect(database_path)

        # Return the connection object
        return conn
    except sqlite3.Error as e:
        print(f"An error occurred while connecting to the database: {e}")
        return None

def find_pattern_rows(conn: sqlite3.Connection, pattern: str, column_name: str) -> List[Tuple]:
    """
    This function will search for and find rows in an SQLite database that match a pattern in a specified column.
    It is called later in the find_and_write_rows function.

    Args:
        conn (sqlite3.Connection): A connection to the SQLite database.
        pattern (str): The pattern to search for, enclosed in quotes - can be a word or exact phrase.
        column_name (str): The name of the column to search in.

    Returns:
        A list of tuples representing the matching rows.
    """
    try:
        # Create a cursor object from the connection
        cur = conn.cursor()

        # Check if the pattern contains special characters
        if re.search(r'[^\w\s]', pattern):
            print("Error: Pattern should only contain letters, numbers, and spaces.")
            return []

        # Remove the quotes from the pattern
        pattern = pattern.replace('"', '')

        # Construct the query string to search for rows that match the given pattern in the specified column
        query = 'SELECT * FROM genes WHERE {} LIKE ?'.format(column_name)
        query += ' AND {} LIKE "%{}%"'.format(column_name, pattern.replace('"', '""'))


        # Execute the query with the pattern and retrieve the matching rows
        cur.execute(query, ('%',))
        matching_rows = cur.fetchall()

        # Return the list of tuples representing the matching rows
        return matching_rows
    except sqlite3.Error as e:
        print(f"An error occurred while searching for rows: {e}")
        return []
    except Exception as e:
        print(f"An unexpected error occurred while searching for rows: {e}")
        return []
    finally:
        # Close the cursor
        cur.close()


def write_pattern_rows_to_csv(filename: str, rows: List[Tuple], headers: List[str]) -> None:
    """
    This function writes a list of rows to a CSV file with the given filename.
    It is called later in the find_and_write_rows function.

    Args:
        filename (str): The name of the output CSV file.
        rows (list): A list of tuples representing the rows to write.
        headers (list): A list of column names to use as headers.

    Returns:
        None
    """
    # Open the output file in write mode and create a CSV writer object
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)

        # Write the header row to the CSV file
        writer.writerow(headers)

        # Write the rows to the CSV file
        writer.writerows(rows)

def find_and_write_rows(args: argparse.Namespace) -> None:
    """
    This function finds rows in the SQLite database specified by the command line arguments that match a pattern in the
    specified column, and writes the matching rows to a new CSV file.

    Args:
        args (argparse.Namespace): The parsed command-line arguments.

    Returns:
        None

    The output CSV file will have the following format:
    - The first row will contain the column headers.
    - Each subsequent row will represent a matching row in the database table, with each column separated by a comma.
    - The filename of the output CSV file will be <pattern>_rows.csv, and it will be written to the same directory as the
    input SQLite database file.


    If no rows in the database match the specified pattern, an empty CSV file will be created with the same format as
    described above.
    """
    # Connect to the SQLite database
    conn = connect_to_database(args.database_path)
    if conn is None:
        return

    # Validate column_name argument
    cur = conn.cursor()
    cur.execute("PRAGMA table_info(genes)")
    columns = [column[1] for column in cur.fetchall()]
    if args.column_name not in columns:
        print(f"Error: {args.column_name} is not a valid column in the database table.")
        conn.close()
        return

    # Find matching rows
    rows = find_pattern_rows(conn, args.pattern, args.column_name)

    # Get column headers
    headers = [header.capitalize() for header in columns]

    # Write matching rows to CSV file
    filename = os.path.join(os.path.dirname(args.database_path), f"{args.pattern}_rows.csv")
    write_pattern_rows_to_csv(filename, rows, headers)

    # Print number of matching rows
    print(f"{len(rows)} rows matched the pattern '{args.pattern}'")

    # Close the connection
    conn.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Search an SQLite database for rows that match a pattern in a specified column, and write the matching rows to a CSV file.')
    parser.add_argument('--pattern', '-p', type=str, required=True, help='The pattern to search for, enclosed in quotes - can be a word or exact phrase.')
    parser.add_argument('--column-name', '-c', type=str, required=True, help='The name of the column to search in.')
    parser.add_argument('--database-path', '--db', type=str, required=True, help='The path to the SQLite database file.')
    args = parser.parse_args()

    find_and_write_rows(args)
