import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

absent_df = pd.DataFrame()


coverage_file = '../data/all_coverage/combined.csv'
# coverage_file = '../data/gene_cov/coverage_pfs12.csv'
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

    annotation_information_file = '../data/annotation_information.csv'
    annotation_information_df = pd.read_csv(annotation_information_file, header=0)

    annotation_information_df.index = annotation_information_df.TranscriptId

    # Correct for badly mapped reads Pfs1
    read_cov_file = '../data/all_coverage/read_cov.csv'
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

    chromosome_file = '../data/extra/pfs1_chromosomes.tsv'
    chromosome_df = pd.read_csv(chromosome_file, header=0,
                                sep='\t')

    merged_anno_df = pd.merge(annotation_information_df, chromosome_df, left_on='Contigname', right_on='contig')
    merged_cov_file = pd.merge(merged_anno_df, coverage_df, left_on='TranscriptId', right_on='mrna_id')

    length_cov_df = merged_cov_file[['length', 'coverage', 'repeatregion', 'pres_abs', 'mrna_id']]
    length_cov_df.loc[length_cov_df.repeatregion.notna(), 'repeatregion'] = 'yes'
    length_cov_df.loc[length_cov_df.repeatregion.isna(), 'repeatregion'] = 'no'
    length_cov_df = length_cov_df[length_cov_df.length > 1000]

    bins = np.linspace(length_cov_df.length.min()-1, length_cov_df.length.max(), 300)
    groups = length_cov_df.groupby(pd.cut(length_cov_df.length, bins))
    bin_std = groups.std().coverage
    bin_mean = groups.mean().coverage
    boundary = length_cov_df.coverage.mean() - groups.std().coverage
    length_cov_df['std'] = bin_std.loc[list(length_cov_df.length)].reset_index().coverage
    length_cov_df['mean'] = bin_mean.loc[list(length_cov_df.length)].reset_index().coverage

    absent = length_cov_df[length_cov_df.pres_abs == 0]
    absent.index = absent.mrna_id
    absent.coverage = 1
    temp_absent = pd.DataFrame(absent.coverage)
    temp_absent.columns = [isolate]
    # absent_df[isolate] = absent.coverage
    absent_df = pd.concat([absent_df, temp_absent], axis=1)

    # length_cov_df.sort_values(by='coverage')
    # sns.pairplot(x_vars=["length"], y_vars=["coverage"], data=length_cov_df, hue="repeatregion", size=5, plot_kws=dict(alpha=.3), markers="+")
    # # length_cov_df.plot.scatter(x='length', y='coverage', c='repeatregion')
    # plt.axhline(y=1.8)
    # plt.ylim(0, 3.5)
    # plt.title(isolate)
    # plt.savefig('../results/coverage_effectors/coverage_effector_{}.png'.format(isolate))
    # plt.show()

# absent_df.to_csv('../results/absent_genes.csv', na_rep=0)