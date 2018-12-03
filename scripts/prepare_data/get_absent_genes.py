import pandas as pd


absent_df = pd.DataFrame()

coverage_file = '../../data/all_coverage/combined.csv'
coverage_df = pd.read_csv(coverage_file, index_col='mrna_id', header=0)

isolates = coverage_df.columns

for isolate in isolates:
    if isolate == 'Pfs1':
        continue

    print(isolate)

    coverage_df = pd.read_csv(coverage_file, index_col='mrna_id', header=0)

    # Exclude genes with Pfs1 coverage of 0
    coverage_df = coverage_df[coverage_df.Pfs1 > 0]

    # Normalize around the mean of the isolates
    coverage_df = coverage_df / coverage_df.mean()

    coverage_df = pd.DataFrame(coverage_df[isolate])
    coverage_df = coverage_df.reset_index()
    coverage_df.columns = ['mrna_id', 'coverage']

    # Add pres/abs gene column
    coverage_df.index = coverage_df.mrna_id
    coverage_df['pres_abs'] = 1

    annotation_information_file = '../../data/annotation_information.csv'
    annotation_information_df = pd.read_csv(annotation_information_file, header=0)

    annotation_information_df.index = annotation_information_df.TranscriptId

    # Correct for badly mapped reads Pfs1
    read_cov_file = '../../data/all_coverage/read_cov.csv'
    read_cov_df = pd.read_csv(read_cov_file, index_col='mrna_id')

    absent = read_cov_df[read_cov_df[isolate] == 0]
    absent = absent[absent.max(axis=1) > 0]
    absent_index = coverage_df.loc[absent.index].dropna().index
    coverage_df.loc[absent_index, 'pres_abs'] = 0

    exclude_gene_df = read_cov_df[read_cov_df.Pfs1 < 1][read_cov_df.max(axis=1) < 1]
    annotation_information_df = annotation_information_df.loc[
        ~annotation_information_df.TranscriptId.isin(exclude_gene_df.index)]

    # Filter out repeats
    # repeat = annotation_information_df.repeatregion
    # annotation_information_df = annotation_information_df.loc[list(repeat[repeat.isna()].index)]

    # Only take effectors into account
    # effector = annotation_information_df.Effector
    # annotation_information_df = annotation_information_df.loc[list(effector[effector == 'Yes'].index)]

    chromosome_file = '../../data/extra/pfs1_chromosomes.tsv'
    chromosome_df = pd.read_csv(chromosome_file, header=0,
                                sep='\t')

    merged_anno_df = pd.merge(annotation_information_df, chromosome_df,
                              left_on='Contigname', right_on='contig')
    merged_cov_file = pd.merge(merged_anno_df, coverage_df,
                               left_on='TranscriptId', right_on='mrna_id')

    length_cov_df = merged_cov_file[['length', 'coverage', 'pres_abs', 'mrna_id']]
    length_cov_df = length_cov_df[length_cov_df.length > 1000]

    absent = length_cov_df[length_cov_df.pres_abs == 0]
    absent.index = absent.mrna_id
    absent.coverage = 1
    temp_absent = pd.DataFrame(absent.coverage)
    temp_absent.columns = [isolate]
    absent_df = pd.concat([absent_df, temp_absent], axis=1)

absent_df.index.name = 'mrna_id'
absent_df.to_csv('../../results/plot_data/absent_genes.csv', na_rep=0)
