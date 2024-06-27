import os

def read_fasta(file_path):
    """ Read a FASTA file and return a dictionary of unique sequences with their headers. """
    sequences = {}
    with open(file_path, 'r') as file:
        while True:
            header = file.readline().strip()
            if not header:
                break
            sequence = file.readline().strip()
            if sequence not in sequences:  # Check for uniqueness
                sequences[sequence] = header
    return sequences

def combine_fastas(directory, output_file):
    all_sequences = {}
    for filename in os.listdir(directory):
        if filename.endswith(".fa"):  # Adjust if your files have a different extension
            file_path = os.path.join(directory, filename)
            sequences = read_fasta(file_path)
            all_sequences.update(sequences)  # Update with the new sequences from each file
    
    # Write the unique sequences to the output file
    with open(output_file, 'w') as f:
        for seq, header in all_sequences.items():
            f.write(f"{header}\n{seq}\n")

# Directory containing FASTA files and the name for the output file
fasta_directory = '/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/reference_genome/24nt_unzip/'
output_file = '/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/reference_genome/combined_unique_pl24nt_siRNA.fasta'
combine_fastas(fasta_directory, output_file)
