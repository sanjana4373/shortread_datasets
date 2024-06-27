def split_fasta(/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/reference_genome/pln24_siRNA_sorted_copy.fasta, 25000):
    """
    Splits a FASTA file into multiple smaller files, each containing a specified number of genes.
    """
    try:
        with open(filename, 'r') as fasta:
            current_file = None
            gene_count = 0
            file_index = 1

            for line in fasta:
                if line.startswith('>'):  # This is the start of a new gene
                    if gene_count % genes_per_file == 0:
                        if current_file:
                            current_file.close()  # Close previous file if open
                        # Open a new file for the next batch of genes
                        current_file = open(f"{filename.split('.')[0]}_part_{file_index}.fasta", 'w')
                        file_index += 1
                    gene_count += 1

                if current_file:
                    current_file.write(line)  # Write the line to the active file

            if current_file:
                current_file.close()  # Make sure the last file is closed

    except IOError as e:
        print(f"Error reading file: {e}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python split_fasta.py <filename.fasta> <genes_per_file>")
    else:
        filename = sys.argv[1]
        genes_per_file = int(sys.argv[2])
        split_fasta(filename, genes_per_file)
