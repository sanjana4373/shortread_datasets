import os
import gzip
import shutil

def decompress_and_move(source_folder, target_folder):
    # Ensure the target folder exists
    os.makedirs(target_folder, exist_ok=True)

    # Process each file in the source folder
    for filename in os.listdir(source_folder):
        if filename.endswith('.fa.gz'):
            # Full path to the compressed file
            compressed_file = os.path.join(source_folder, filename)
            # Define the path for the decompressed file
            decompressed_file = os.path.join(target_folder, filename[:-3])  # Remove '.gz' from the filename

            # Decompress the file
            with gzip.open(compressed_file, 'rb') as f_in:
                with open(decompressed_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

# Paths to the source and target folders
source_folder = '/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/reference_genome/filter_data'
target_folder = '/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/reference_genome/24nt_unzip'

# Function call to start the process
decompress_and_move(source_folder, target_folder)

