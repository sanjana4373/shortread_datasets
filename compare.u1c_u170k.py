import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the files with proper delimiters
U170K = pd.read_csv('/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_unique/idxstats/volcano_siRNAs_u170k.csv')
U1C = pd.read_csv('/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_unique/idxstats/de_u1c_siRNAs.tsv', delimiter='\t')

# Merge the DataFrames on sequence_name
merged_df = pd.merge(U170K[['sequence_name', 'diffexpressed']], U1C[['sequence_name', 'diffexpressed']], on='sequence_name', suffixes=('_file1', '_file2'))

# Determine if the expression status matches
merged_df['comparison'] = merged_df.apply(lambda x: 'Match' if x['diffexpressed_file1'] == x['diffexpressed_file2'] else 'Mismatch', axis=1)

# Summary of matches and mismatches
summary_counts = merged_df['comparison'].value_counts()
print(summary_counts)

# Detailed breakdown of mismatches
mismatch_detail = merged_df[merged_df['comparison'] == 'Mismatch'].pivot_table(index='diffexpressed_file1', columns='diffexpressed_file2', aggfunc='size', fill_value=0)
print(mismatch_detail)
mismatch_detail.to_csv('/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_unique/idxstats/mismatches_detail.csv')
# Save the detailed mismatches to a CSV file
mismatches = merged_df[merged_df['comparison'] == 'Mismatch']
mismatches.to_csv('/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_unique/idxstats/mismatches.csv', index=False)

print("Detailed Mismatches:")
print(mismatches)

# Heatmap of matches and mismatches
plt.figure(figsize=(10, 8))
sns.heatmap(mismatch_detail, annot=True, cmap="Blues", fmt="d")
plt.title('Mismatched Expression Comparisons')
plt.xlabel('Expressions in U1C')
plt.ylabel('Expressions in U170K')
plt.savefig('mismatch_heatmap.png')
plt.show()

