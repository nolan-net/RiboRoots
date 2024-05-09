## Author: Nolan Top
## Description: Find all genera and order them most entries to least
## Use Case: This helped me figure out what genus had a small enough number to
##           of entries to test quickly, but not too little to show off the programs uses.

from Bio import SeqIO
from collections import Counter

# File path for the original FASTA
input_fasta = 'Final.fasta'  # Update this as necessary

# Load the original sequences
original_sequences = list(SeqIO.parse(input_fasta, "fasta"))

# Extract the genus from each sequence description
genera = [record.description.split('_')[0] for record in original_sequences if '_' in record.description]

# Count the occurrences of each genus
genus_counts = Counter(genera)

# Sort the genera by their frequency
sorted_genera = sorted(genus_counts.items(), key=lambda x: x[1], reverse=True)

# Print out the sorted genera and their counts
print("Most popular genera sorted by frequency:")
for genus, count in sorted_genera:
    print(f"{genus}: {count}")
