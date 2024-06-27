import pandas as pd

# Load the files
de_u1c_siRNAs = pd.read_csv('/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_loci/idxstats/de_u1c_siRNAs.tsv', sep='\t')
volcano_siRNAs_u170k = pd.read_csv('/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_loci/idxstats/volcano_siRNAs_u170k.csv')

# Display the first few rows of each file to understand their structure
print(de_u1c_siRNAs.head())
print(volcano_siRNAs_u170k.head())

# Extract UP-regulated clusters from both datasets
up_de_u1c_siRNAs = de_u1c_siRNAs[de_u1c_siRNAs['diffexpressed'] == 'UP']
up_volcano_siRNAs_u170k = volcano_siRNAs_u170k[volcano_siRNAs_u170k['diffexpressed'] == 'UP']

# Extract DOWN-regulated clusters from both datasets
down_de_u1c_siRNAs = de_u1c_siRNAs[de_u1c_siRNAs['diffexpressed'] == 'DOWN']
down_volcano_siRNAs_u170k = volcano_siRNAs_u170k[volcano_siRNAs_u170k['diffexpressed'] == 'DOWN']

# Find the overlapping UP-regulated clusters
overlap_up = pd.merge(up_de_u1c_siRNAs, up_volcano_siRNAs_u170k, on='sequence_name')

# Find the overlapping DOWN-regulated clusters
overlap_down = pd.merge(down_de_u1c_siRNAs, down_volcano_siRNAs_u170k, on='sequence_name')

print("Overlapping UP-regulated Clusters")
print(overlap_up)
print("Overlapping DOWN-regulated Clusters")
print(overlap_down)

# Save the overlapping dataframes to CSV files for record
overlap_up.to_csv('/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_loci/idxstats/overlapping_up_clusters.csv', index=False)
overlap_down.to_csv('/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_loci/idxstats/overlapping_down_clusters.csv', index=False)
