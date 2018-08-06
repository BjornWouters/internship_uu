import os
import pandas as pd
import numpy as np

min_contig_length = 1000
identifier_file = '../data/gene_identifiers.csv'
coverage_folder = '../data/all_coverage'

id_file = pd.read_csv(identifier_file, index_col=False)
id_file = id_file.set_index('mrna_id')
output_df = pd.DataFrame()
test_df = pd.DataFrame()

for filename in os.listdir(coverage_folder):
    file_path = os.path.join(coverage_folder, filename)
    coverage_file = pd.read_csv(file_path, sep='\t', index_col=False, header=None, names=[
        'name', 'start', 'end', 'mrna_id', 'dot', 'strand', 'dot2', 'type', 'dot3', 'description',
        'coverage'
    ])
    coverage_file = coverage_file.set_index('mrna_id')

    temp_df = pd.DataFrame(id_file.loc[list(coverage_file.index)])

    temp_df = temp_df[temp_df.contig_length > min_contig_length]

    temp_df['coverage'] = coverage_file.coverage
    # temp_df = temp_df.set_index('gene_name')

    name = filename.strip('_coverage.txt')
    if filename.startswith('pfs'):
        output_df[name] = temp_df['coverage']
    else:
        test_df[name] = temp_df['coverage']

output_df = output_df.loc[output_df.index.notna()]
test_df = test_df.loc[test_df.index.notna()]
# output_df.index = output_df.index.str.replace('pfs1\|', '')
# output_df.index = pd.to_numeric(output_df.index)
output_df = output_df.sort_index()
test_df = test_df.sort_index()
# output_df.index = 'pfs1|' + output_df.index.astype(str)
output_df.to_csv('../data/gene_cov/read_cov.csv')
test_df.to_csv('../data/gene_cov/read_test_cov.csv')
