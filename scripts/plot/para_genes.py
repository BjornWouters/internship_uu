import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

absent_file = '../../results/plot_data/absent_genes.csv'
absent_df = pd.read_csv(absent_file, header=0, index_col='mrna_id')

sb.clustermap(absent_df)
plt.show()
