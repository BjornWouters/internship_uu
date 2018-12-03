import os

import pandas as pd

coverage_folder = '../data/gene_cov/'
output_df = pd.DataFrame()

for filename in os.listdir(coverage_folder):
    file_path = os.path.join(coverage_folder, filename)
    coverage_file = pd.read_csv(file_path, header=0, index_col='mrna_id')
    coverage = coverage_file['coverage']

    isolate_name = filename.split('_')[1].split('.')[0].title()
    output_df[isolate_name] = coverage

# Normalize around the mean of the isolates
output_df = output_df / output_df.mean()

para_offset = output_df.mean() + output_df.std()*2

output_df[output_df > para_offset] = 10
output_df[output_df <= para_offset] = 0
output_df[output_df == 10] = 1

annotation_information_file = '../data/annotation_information.csv'
annotation_information_df = pd.read_csv(annotation_information_file, header=0,
                                        index_col='TranscriptId')

# FILTERS #

# annotation_information_df['paralogous'] = False
# annotation_information_df.loc[list(coverage_df.index), 'paralogous'] = True
# annotation_information_df.to_csv('../data/anno_info.csv', index=True)

# effector = annotation_information_df.Effector
# annotation_information_df = annotation_information_df.loc[list(effector[effector == 'Yes'].index)].dropna()
# coverage_df = output_df.loc[list(effector.dropna().index)]

# repeat = annotation_information_df.repeatregion
# coverage_df = coverage_df.loc[list(repeat[repeat.isna()].index)].dropna()

# signalp = annotation_information_df.Signalp
# coverage_df = coverage_df.loc[list(signalp[signalp.notna()].index)].dropna()


sum_df = output_df.transpose().sum()
filter_df1 = sum_df[(sum_df > 0) & (sum_df < 24)]
output_df = output_df.loc[filter_df1.index]

output_df.to_csv('../results/plot_data/para_genes.csv', index=True)
