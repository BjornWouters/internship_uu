import pandas as pd
from sklearn import preprocessing
import matplotlib.pyplot as plt


def import_data():
    # columns = [
    #     'genename', 'Pfs1RZ', 'Pfs2RZ', 'Pfs3RZ_C', 'Pfs4PV_C', 'Pfs5RZ_C', 'Pfs6RZ', 'Pfs7PV',
    #     'Pfs8RZ_C', 'Pfs9RZ', 'Pfs10RZ', 'Pfs11PV', 'Pfs12RZ', 'Pfs13PV', 'Pfs14RZ',
    #     'Pfs15RZ', 'Pfs16PV_C'
    # ]
    train_data = '../data/gene_cov/read_cov.csv'
    gene_data = '../data/gene_cov/read_test_cov.csv'
    coverage_df = pd.read_csv(gene_data, index_col='mrna_id')
    train_df = pd.read_csv(train_data, index_col='mrna_id')
    return coverage_df, train_df


def prepare_data(df):
    offset = 0.1
    offset_para = 3

    # Lose all the genes that have a mean coverage lower than threshold
    mean_thr = 1
    df = df[df.mean(axis=1) > mean_thr]

    # Normalize the mean and standard deviation between the isolates
    # normalized_df = preprocessing.scale(df, axis=1)

    # Divide all the genes by the mean over all isolates
    corrected_df = df.div(df.mean(axis=1), axis=0)
    # corrected_df.plot.density()
    # plt.savefig('../results/density_plot.png', dpi=900)

    corrected_df[corrected_df <= offset] = int(0)
    corrected_df[(corrected_df >= offset) & (corrected_df <= offset_para)] = 1
    corrected_df[corrected_df > offset_para] = int(2)
    return corrected_df


def main():
    coverage_df, train_df = import_data()
    prepared_df = prepare_data(coverage_df)
    train_df = prepare_data(train_df)

    # Correction for equal features only use with test set
    # valid_rows = prepared_df.loc[prepared_df['ES1314'].notnull()].index.values
    # prepared_df = prepared_df.loc[valid_rows]
    # train_df = train_df.loc[valid_rows]

    prepared_df.to_csv('../data/prepared_test_dataset.csv')
    train_df.to_csv('../data/prepared_dataset.csv')


if __name__ == '__main__':
    main()
