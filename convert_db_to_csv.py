import sqlite3
import csv
import argparse
import os

def convert_db_to_csv(db_file: str) -> None:
    """Converts an SQLite database to a single CSV file.
    
    Args:
        db_file (str): The path to the SQLite database file.
    """
    
    # Connect to the SQLite database
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()

    # Get a list of all tables in the database
    cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = [table[0] for table in cur.fetchall()]

    # Get the output CSV file name and path
    csv_file = os.path.splitext(os.path.basename(db_file))[0] + '.csv'
    output_dir = os.path.dirname(db_file)
    csv_file_path = os.path.join(output_dir, csv_file)

    # Create and open the output CSV file
    with open(csv_file_path, 'w', newline='', encoding='utf-8') as csvfile:
        csv_writer = csv.writer(csvfile)
        
        # Loop through all the tables and write their records to the CSV file
        for table_name in tables:
            # Query to select all records from the table
            query = f'SELECT * FROM {table_name}'
            
            # Execute the query
            cur.execute(query)
            
            # Get the column names from the table
            column_names = [desc[0].capitalize() for desc in cur.description]  # Capitalize column names
            
            # Write the header if it's the first table
            if table_name == tables[0]:
                csv_writer.writerow(column_names)
            
            # Write the records to the CSV file
            csv_writer.writerows(cur.fetchall())

    # Close the database connection
    conn.close()

if __name__ == '__main__':
    # Set up the command line argument parser
    parser = argparse.ArgumentParser(description='Convert an SQLite database to a single CSV file.')
    parser.add_argument('--db', '-d', required=True, help='Path to the SQLite database file')

    # Parse the command line arguments
    args = parser.parse_args()

    # Get the database file from the parsed arguments
    db_file = args.db
    
    # Convert the database to CSV
    convert_db_to_csv(db_file)
