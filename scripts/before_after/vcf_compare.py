from Bio import SeqIO
import vcf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import re

vcf_file = '../data/extra/combined_vcf/Pfs1_clean.g.vcf'

annotation_information_file = '../../data/annotation_information.csv'


reader = vcf.Reader(filename=vcf_file)
annotation_information_df = pd.read_csv(annotation_information_file, header=0,
                                        index_col='TranscriptId')

# effector = annotation_information_df.Effector
# annotation_information_df = annotation_information_df.loc[list(effector[effector == 'Yes'].index)]
#
# signalp = annotation_information_df.Signalp
# annotation_information_df = annotation_information_df.loc[list(signalp[signalp.notna()].index)]
#
# repeat = annotation_information_df.repeatregion
# annotation_information_df = annotation_information_df.loc[list(repeat[repeat.isna()].index)]

for vcf_call in vcf:
    print(vcf_call)




#
# lut = dict(zip(effector.unique(), 'br'))
# row_colors = effector.map(lut)
# effector_list = list(effector[effector == 'Yes'].index)
# cl_map = sns.clustermap(count_df.astype(float), figsize=(15, 15), row_colors=row_colors)
# cl_map.fig.set_tight_layout(False)
# for tick_label in cl_map.ax_heatmap.axes.get_yticklabels():
#     tick_text = tick_label.get_text()
#     if tick_text in effector_list:
#         tick_label.set_color('red')
#
# plt.show()
