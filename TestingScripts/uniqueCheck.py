## Author: Nolan Top
## Description: Quick script to make sure only unique entries are found

import pandas as pd

# Path to the processed data
csv_file_path = input("Enter the path to the CSV file: ")

# Load the processed data
data = pd.read_csv(csv_file_path)

# This dictionary will store sequences for each unique description (minus the last character)
unique_sequences_check = {}

# Variables to track potential duplicates
duplicates_found = False
duplicate_details = []

# Iterate over the dataframe
for index, row in data.iterrows():
    # Create a key from the 'Description' column excluding the last character
    key = row['Description'][:-2]  # Exclude last character which is '+' or '-'

    # Sequence as string for comparison
    sequence_str = row['Sequence']

    # Check if the sequence is already in the list for this key
    if key in unique_sequences_check:
        if sequence_str in unique_sequences_check[key]:
            # Duplicate found
            duplicates_found = True
            duplicate_details.append((key, sequence_str))
        else:
            # Add the sequence to the list for this key
            unique_sequences_check[key].append(sequence_str)
    else:
        # If the key is not in the dictionary, add it with this sequence in a new list
        unique_sequences_check[key] = [sequence_str]

# Check if any duplicates were found
if duplicates_found:
    print("Duplicates found:")
    for dup in duplicate_details:
        print(f"Duplicate group: {dup[0]}, Sequence: {dup[1]}")
else:
    print("No duplicates found. Each sequence is unique for its group.")
