from scipy.stats import hypergeom
import pandas as pd

def hypergeometric_test(file1, file2, regulation_type='UP'):
    """
    Perform a hypergeometric test for overlap of differentially expressed genes.
    
    Parameters:
    - file1: str, path to the first CSV/TSV file containing gene expression data.
    - file2: str, path to the second CSV/TSV file containing gene expression data.
    - regulation_type: str, 'UP' or 'DOWN', type of regulation to test for overlap.
    
    Returns:
    - p_value: float, p-value for the overlap of specified regulation type.
    - overlap_genes_details: DataFrame, details of overlapping genes.
    """
    
    # Load the datasets
    if file1.endswith('.csv'):
        df1 = pd.read_csv(file1)
    elif file1.endswith('.tsv'):
        df1 = pd.read_csv(file1, sep='\t')
    else:
        raise ValueError("Unsupported file format for file1. Please provide a CSV or TSV file.")
    
    if file2.endswith('.csv'):
        df2 = pd.read_csv(file2)
    elif file2.endswith('.tsv'):
        df2 = pd.read_csv(file2, sep='\t')
    else:
        raise ValueError("Unsupported file format for file2. Please provide a CSV or TSV file.")
    
    # Total number of genes in each dataset
    total_genes = len(df1)
    
    # Number of specified regulation type genes in each dataset
    reg_genes1 = df1[df1['diffexpressed'] == regulation_type].shape[0]
    reg_genes2 = df2[df2['diffexpressed'] == regulation_type].shape[0]
    
    # Identifying overlapping regulation type genes
    genes1 = set(df1[df1['diffexpressed'] == regulation_type]['sequence_name'])
    genes2 = set(df2[df2['diffexpressed'] == regulation_type]['sequence_name'])
    
    overlap_genes = genes1.intersection(genes2)
    overlap_count = len(overlap_genes)
    
    # Calculate the p-value using hypergeometric test
    M = total_genes  # Population size
    n = reg_genes1  # Number of successes in the population
    N = reg_genes2  # Number of draws
    x = overlap_count  # Number of observed successes
    
    # Hypergeometric test
    p_value = hypergeom.sf(x-1, M, n, N)
    
    # Extract details of overlapping genes
    overlap_genes_details = df1[df1['sequence_name'].isin(overlap_genes)]
    
    return p_value, overlap_genes_details

# Example usage:
file1 = '/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_loci/idxstats/volcano_siRNAs_u170k.csv'
file2 = '/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_loci/idxstats/de_u1c_siRNAs.tsv'
p_value, overlap_genes_details = hypergeometric_test(file1, file2, regulation_type='UP')

print(f"P-value for UP-regulated gene overlap: {p_value}")
print(overlap_genes_details)
