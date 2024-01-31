# This script will clean our Michigan data by making sure all dna sequences are unique
# While the initial data having several entries for the same sequence is not a problem
# for our analysis we don't need to have the same sequence multiple times in our data

import csv

input_csv_file_path = input("Input the path to the CSV file to be cleaned: ")
output_csv_file_path = input("Output filename: ")


def main():
    # Open the input and output files   
    with open(input_csv_file_path, 'r') as input_csv_file:
        with open(output_csv_file_path, 'w') as output_csv_file:
            # Create our reader and writer
            csv_reader = csv.reader(input_csv_file)
            csv_writer = csv.writer(output_csv_file)
            # Create a set to store the unique sequences
            unique_sequences = set()
            # Iterate through the rows of the CSV file
            for row in csv_reader:
                # Check if the sequence (3rd column) is not in the set
                if row[2] not in unique_sequences:
                    # Add the sequence to the set
                    unique_sequences.add(row[2])
                    
                    # Extract the scientific name and the forward/reverse indicator
                    parts = row[1].split('|')
                    scientific_name = parts[0].strip()  # Get the first part and strip whitespace
                    forward_reverse_indicator = row[1][-1]  # Get the last character of the second column
                    
                    # Prepare the row for writing, replacing the second column with scientific name and forward/reverse indicator
                    row_to_write = row[:1] + [scientific_name + ' ' + forward_reverse_indicator] + row[2:]
                    
                    # Write the modified row to the output file
                    csv_writer.writerow(row_to_write)


if __name__ == '__main__':
    main()