import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Load the files
U170K = pd.read_csv('/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_loci/idxstats/volcano_siRNAs_u170k.csv')
U1C = pd.read_csv('/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_loci/idxstats/de_u1c_siRNAs.tsv', delimiter='\t')

# Assuming 'diffexpressed' contains the categories 'UP', 'DOWN', 'NO'
up_set1 = set(U170K[U170K['diffexpressed'] == 'UP']['sequence_name'])
up_set2 = set(U1C[U1C['diffexpressed'] == 'UP']['sequence_name'])

down_set1 = set(U170K[U170K['diffexpressed'] == 'DOWN']['sequence_name'])
down_set2 = set(U1C[U1C['diffexpressed'] == 'DOWN']['sequence_name'])

no_set1 = set(U170K[U170K['diffexpressed'] == 'NO']['sequence_name'])
no_set2 = set(U1C[U1C['diffexpressed'] == 'NO']['sequence_name'])

# Create a Venn diagram for 'UP'
plt.figure(figsize=(8, 8))
venn2([up_set1, up_set2], ('U170K UP', 'U1C UP'))
plt.title("Venn Diagram of UP Expression")
plt.show()
plt.annotate('Clusters : 9026\nP-value: 6.007e-55', xy=(0.5, -0.1), xycoords='axes fraction', ha='center')
plt.savefig('/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_loci/idxstats/UP.png')

# Similarly for 'DOWN'
plt.figure(figsize=(8, 8))
venn2([down_set1, down_set2], ('U170K DOWN', 'U1C DOWN'))
plt.title("Venn Diagram of DOWN Expression")
plt.show()
plt.annotate('Clusters : 9026\nP-value: 2.3e-32', xy=(0.5, -0.1), xycoords='axes fraction', ha='center')
plt.savefig('/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_loci/idxstats/DOWN.png')

# And for 'NO'
plt.figure(figsize=(8, 8))
venn2([no_set1, no_set2], ('U170K NO', 'U1C NO'))
plt.title("Venn Diagram of NO Expression")
plt.show()
plt.savefig('/vol/data/smallRNAseq_12102023/X204SC23083269-Z01-F001/MasterThesis/results_siRNAs_loci/idxstats/No.png')
