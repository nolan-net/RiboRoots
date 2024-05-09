## Author: Nolan Top
## Date: February 8, 2024
## Description: This script will convert FASTA files to CSV,
##          as long as they follow the format described below


from Bio import SeqIO
import csv

# Get file paths from user
fasta_file = input("Enter the path to your FASTA file: ")
csv_file = input("Enter the path for the output CSV file: ")

# Open the FASTA file and the CSV file
with open(fasta_file, mode='r') as fasta, open(csv_file, mode='w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    # Write the header
    csv_writer.writerow(['Sequence ID', 'Description', 'Sequence'])

    # Iterate over the FASTA file and write to the CSV file
    for record in SeqIO.parse(fasta, 'fasta'):
        csv_writer.writerow([record.id, record.description, str(record.seq)])

print(f"Conversion complete. CSV file saved as '{csv_file}'")

