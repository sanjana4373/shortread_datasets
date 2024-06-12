from Bio import SeqIO
import pandas as pd

# Load sequence names from the CSV file
csv_file = '/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/up_regulated_sequence_names.csv'  # Change to your CSV file path
df = pd.read_csv(csv_file)
sequence_names = set(df['sequence_name'].values)  # Assuming 'SequenceName' is the column header

# Path to your large FASTA file
fasta_file = '/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/filtered_siRNAs_unique_sequences.fasta'  # Change to your FASTA file path
output_file = '/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/up_sequences.fasta'  # Output file path

# Read FASTA and write matching sequences to a new file
with open(output_file, 'w') as out_fasta:
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id in sequence_names:
            SeqIO.write(record, out_fasta, 'fasta')

print("Extraction complete.")
