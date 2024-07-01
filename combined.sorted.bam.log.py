import pandas as pd
import os

# Define the directory where your files are stored
directory = "/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_loci/idxstats/"

# Initialize an empty list to hold data from each file
dataframes = []

# Define the filenames explicitly
filenames = [
    "Col0SR_1.sorted.bam.log", "Col0SR_2.sorted.bam.log", "Col0SR_3.sorted.bam.log",
    "u1c_1.sorted.bam.log", "u1c_2.sorted.bam.log", "u1c_3.sorted.bam.log",
    "u170k_1.sorted.bam.log", "u170k_2.sorted.bam.log", "u170k_3.sorted.bam.log"
]

# Loop through each filename
for filename in filenames:
    file_path = os.path.join(directory, filename)
    # Read only the first and third columns (0-indexed)
    df = pd.read_csv(file_path, sep="\t", usecols=[0, 2], header=None, names=['reference sequence name', 'mapped read-segments'])
    
    # Rename the 'mapped read-segments' column to be unique based on the filename
    df.rename(columns={'mapped read-segments': filename.replace('.sorted.bam.log', '')}, inplace=True)
    
    # Append to the list of dataframes
    dataframes.append(df)

# Merge all dataframes on the 'reference sequence name' column
combined_df = dataframes[0]
for df in dataframes[1:]:
    combined_df = pd.merge(combined_df, df, on="reference sequence name", how='outer')

# Replace NaN with 0 where there are no data for some files
combined_df.fillna(0, inplace=True)

# Save the combined data to a new CSV file
combined_df.to_csv("/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_loci/idxstats/combined_siRNAs_reads.tsv", sep='\t', index=False)

