# Author: Nolan Top

# Description: This is our cleaning file. It will ensure completely unique dna strands,
# Then it will assign unique id's to the remaining (assuming we have multiple bacteria of the same species)
# Then it will output it to a new csv

# Assumption: The FASTA is in csv form and has been cleaned to have a 'Description' column containing the name
#              followed by a 'Sequence' column with the sequence.

import pandas as pd
from Bio.Seq import Seq
import re

# Function to clean special characters from descriptions
def clean_description(desc):
    # Replace problematic characters with underscores or remove them
    desc = re.sub(r'[,\(\)\[\]:;]', '_', desc)  # Characters to be replaced
    desc = re.sub(r'\s+', '_', desc)  # Replace whitespace with '_'
    desc = re.sub(r'_{2,}', '_', desc)  # Don't use multiple underscores
    return desc.strip('_')  # Remove leading or trailing underscores

# Input name of the csv file
input_file = input("Enter the path to your CSV file: ")

# Input name of the output csv file
output_file = input("Enter the path for the output CSV file: ")

# Load the data from the CSV file
data = pd.read_csv(input_file)

# This dictionary will store the unique forward sequences for each Bacteria with updated descriptions
unique_sequences = {}

# To track IDs and ensure uniqueness
id_counter = {}

# Let's go through the dataframe
for index, row in data.iterrows():

    # Extract the sequence as a biopython sequence and check if it needs to be reversed
    sequence = Seq(row['Sequence'])
    if row['Description'].endswith('-'):
        sequence = sequence.reverse_complement()

    # Convert the sequence back to a string
    sequence_str = str(sequence)

    # Remove trailing '+' and '-' and the whitespace before
    key = clean_description(row['Description'][:-2])

    # Update the ID counter for the key and create a new unique description
    id_count = id_counter.get(key, 0) + 1
    id_counter[key] = id_count
    new_description = f"{key}_{id_count}"

    # Check if the sequence is added
    unique_sequences[new_description] = sequence_str

# Convert the back into dataframe
unique_data = [[desc, seq] for desc, seq in unique_sequences.items()]
unique_df = pd.DataFrame(unique_data, columns=['Description', 'Sequence'])

# Remove duplicate sequences based on the 'Sequence' column
unique_df = unique_df.drop_duplicates(subset='Sequence', keep='first')

# Save to new csv
unique_df.to_csv(output_file, index=False)

print("Done!")
