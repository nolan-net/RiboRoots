## Author: Nolan Top
## Date: March 12, 2024
## Description: This program takes a dataset of sequences, optionally filters and aligns them,
##              and constructs a phylogenetic tree using the selected method (UPGMA or NJ).
##              The tree is outputted as a .nwk file for phylogenetic analysis.
##              The user can also choose to display the tree graphically using ete3.
##              Please ensure your data is in FASTA format, with each entry formatted as a
##              sequence name followed by the DNA sequence.

import subprocess
import os
import sys
from importlib.metadata import distribution, PackageNotFoundError
from ete3 import Tree, TreeStyle
from Bio import AlignIO, SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# Automatic installation of packages if necessary
def install(package):
    confirm = input(f"Package {package} not found. Would you like to install it? (yes/no): ")
    if confirm.lower() == 'yes':
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
            print(f"{package} installed successfully.")
        except subprocess.CalledProcessError:
            print(f"Failed to install {package}. Please try installing it manually.")
    else:
        print("Installation aborted. The package is required for the program to run.")

required_packages = ['biopython', 'ete3', 'numpy']

for package in required_packages:
    try:
        distribution(package)
    except PackageNotFoundError:
        install(package)

def cleanup_temp_files():
    # List of temporary files that might not get deleted if quit/crash unexpectedly
    print("Cleaning up before starting...")
    temp_files = ['temp_sequences.fasta', 'Aligned.fasta']
    for file in temp_files:
        if os.path.exists(file):
            os.remove(file)

def display_tree(newick_file):
    # Load the tree using format=1 to handle internal node names
    tree = Tree(newick_file, format=1)  

    # Set up the tree style
    ts = TreeStyle()
    ts.show_branch_length = True
    tree.show(tree_style=ts)

# Input wrapper to allow safe exit()
def safe_input(prompt):
    user_input = input(prompt)
    if user_input.lower() == 'quit':
        print("Exiting program.")
        exit()  
    return user_input

def show_help():
    help_text = """
    Help - Phylogenetic Tree Constructor
    Usage:
        - Enter the full path to your FASTA file: Provide a path to a .fasta format file.
        - Do the sequences need to be aligned? (yes/no): Indicate whether the tool should align the sequences.
        - Choose the tree construction method (UPGMA/NJ): Select either 'UPGMA' or 'NJ' for constructing the tree.
        - Enter comma-separated terms to filter sequences: Filter input sequences by these terms.
        - Progress will be shown during data processing.
        - You can exit or restart at any decision point.
    Necessities:
        - A FASTA file that contains two columns: [sequence name/id], [gene sequence]
        - Clean your dataset if necessary, i.e: remove duplicate sequences, make sure all sequences are forward or backward
        - Example scripts to help with this are in the github repo: https://github.com/nolan-net/RiboRoots
    """
    print(help_text)

    while True:
        decision = safe_input("Do you want to continue with the execution or quit? (continue/quit): \n")
        if decision.lower() == 'continue':
            print("Continuing with execution...")
            break  
        else:
            print("Invalid input. Please enter 'continue' or 'quit'.")
    


def main():
    print("\nWelcome to the Phylogenetic Tree Constructor!")
    print("You may type 'quit' during any requested input to exit.")
    cleanup_temp_files()  # Clean up if necessary

    
    while True:
        user_choice = safe_input("Would you like to see the help menu or start creating your tree? (help/start):\n").lower()
        if user_choice == 'help':
            show_help()
            break  
        elif user_choice == 'start':
            break 
        else:
            print("Invalid input. Please type 'help' to view the help menu, 'start' to begin tree creation, or 'quit to exit.")

    print("Ensure your data is in FASTA format, with each entry formatted as a sequence name followed by the DNA sequence.")
    print("The data should be cleaned and ready to process.")

    # User input for dataset path + check file exists
    while True:
        input_file = safe_input("Enter the path to your FASTA file: \n")
        if os.path.exists(input_file):
            break 
        else:
            print("File not found. Please enter a valid path to a FASTA file.")

    # Sequence alignment necessity check
    while True:
        align_needed = safe_input("Do the sequences need to be aligned? (yes/no):\n").lower()
        if align_needed in ['yes', 'no']:
            break
        print("Invalid input. Please enter 'yes' or 'no'.")

    input_terms = safe_input("Enter comma-separated terms to filter sequences by (e.g., 'Salmonella, Leclercia, Enterobacter'):\n")
    terms = [term.strip().lower() for term in input_terms.split(',')]

    temp_fasta = 'temp_sequences.fasta'
    output_aligned = 'Aligned.fasta'
    output_tree = 'phylogenetic_tree.nwk'

    original_sequences = list(SeqIO.parse(input_file, "fasta"))
    filtered_sequences = [record for record in original_sequences if any(term in record.description.lower() for term in terms)]

    if align_needed == 'yes':
        print("\nAligning the data with Clustal...")
        SeqIO.write(filtered_sequences, temp_fasta, "fasta")
        subprocess.run(['clustalo', '-i', temp_fasta, '-o', output_aligned, '--auto', '-v'], check=True)
        alignment = AlignIO.read(output_aligned, "fasta")
        os.remove(temp_fasta)
    else:
        alignment = filtered_sequences

    print("\nCalculating distance matrix...")
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    # Tree construction method selection
    while True:
        tree_method = safe_input("Choose the tree construction method (UPGMA/NJ):\n").lower()
        constructor = DistanceTreeConstructor()
        if tree_method == 'nj':
            tree = constructor.nj(distance_matrix)
            break
        elif tree_method == 'upgma':
            tree = constructor.upgma(distance_matrix)
            break
        print("Invalid input. Please choose either 'UPGMA' or 'NJ'.")

    print("\nConstructing the tree...")
    
    Phylo.write(tree, output_tree, "newick")
    if align_needed == 'yes':
        os.remove(output_aligned)

    while True:
        view_graphically = safe_input("Would you like to visually display the phylogenetic tree using ete3? The tree file will be saved either way. (yes/no):\n").lower()
        if view_graphically == 'yes':
            display_tree(output_tree)
            print("Tree displayed (If you're having trouble viewing, change device display settings to 'light mode').\nThe Newick file is also saved at:", output_tree)
            break
        elif view_graphically == 'no':
            print("Tree visualization skipped. The Newick file is available at:", output_tree)
            break
        else:
            print("Invalid input. Please enter 'yes' or 'no'.")


    print("\nPhylogenetic tree construction complete. Output saved to:", output_tree)

if __name__ == "__main__":
    main()