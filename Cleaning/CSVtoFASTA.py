## Author: Nolan Top
## Date: February 8, 2024
## Description: This script will convert CSV files to FASTA,
##          as long as they follow the format described below


import csv
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

# Get file paths from user
csv_file = input("Enter the path to the CSV file: ")
fasta_file = input("Enter the path for your output FASTA file: ")

# Read the CSV and write the FASTA
with open(csv_file, 'r') as csv_handle, open(fasta_file, 'w') as fasta_handle:
    reader = csv.reader(csv_handle)
    for row in reader:
        # In my csv the first column is the scientific name, and the second is the sequence
        scientific_name, sequence = row[0], row[1]
        # Replace spaces with underscores in the scientific name for ID
        sanitized_name = scientific_name.replace(" ", "_")
        # Create a SeqRecord object
        record = SeqRecord(Seq(sequence), id=sanitized_name, description="")
        # Write the SeqRecord to the FASTA file
        SeqIO.write(record, fasta_handle, "fasta")
