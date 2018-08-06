import os

import pandas as pd
import matplotlib.pyplot as plt

high_impact_offset = 5
moderate_impact_offset = 80
prepared_cov_file = '../data/prepared_dataset.csv'
snp_folder = '../data/snps/'
output_df = pd.DataFrame()

for filename in os.listdir(snp_folder):
    file_path = os.path.join(snp_folder, filename)
    snp_file = pd.read_csv(file_path, header=0, index_col='TranscriptId', sep='\t')

    high_impact_df = snp_file['variants_impact_HIGH']
    high_impact_df[high_impact_df < high_impact_offset] = 0
    high_impact_df[high_impact_df >= high_impact_offset] = 1

    moderate_impact_df = snp_file['variants_impact_MODERATE']
    moderate_impact_df[moderate_impact_df < high_impact_offset] = 0
    moderate_impact_df[moderate_impact_df >= high_impact_offset] = 1

    isolate_name = filename.split('_')[0]

    if isolate_name.startswith('Pfs'):
        output_df[isolate_name.lower()] = high_impact_df + moderate_impact_df

    # Plot high impact mutations`
    # plt.plot(range(len(snp_file['variants_impact_HIGH'].sort_values())), snp_file['variants_impact_HIGH'].sort_values())
    # plt.show()

prepared_cov_file = pd.read_csv(prepared_cov_file, header=0, index_col='mrna_id')
valid_mrna_ids = prepared_cov_file.index
output_df = output_df.loc[list(valid_mrna_ids)]
output_df = output_df.fillna(0)

output_df.index.name = 'mrna_id'

pres_abs_file = '../data/prepared_dataset.csv'
pres_abs_df = pd.read_csv(pres_abs_file, header=0, index_col='mrna_id')
output_df = output_df.join(pres_abs_df, lsuffix='_snp', rsuffix='_cov')

annotation_information_file = '../data/annotation_information.tsv'
annotation_information_df = pd.read_csv(annotation_information_file, header=0,
                                        index_col='TranscriptId', sep='\t')
signalp = annotation_information_df.Signalp
output_df = output_df.loc[list(signalp[signalp.notna()].index)].dropna()

output_df.to_csv('../data/parsed_snps.csv', index=True)
