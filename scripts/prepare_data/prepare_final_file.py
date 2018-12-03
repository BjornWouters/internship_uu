import pandas as pd
import numpy as np
from sklearn import preprocessing
import matplotlib.pyplot as plt


def import_data():
    # columns = [
    #     'genename', 'Pfs1RZ', 'Pfs2RZ', 'Pfs3RZ_C', 'Pfs4PV_C', 'Pfs5RZ_C', 'Pfs6RZ', 'Pfs7PV',
    #     'Pfs8RZ_C', 'Pfs9RZ', 'Pfs10RZ', 'Pfs11PV', 'Pfs12RZ', 'Pfs13PV', 'Pfs14RZ',
    #     'Pfs15RZ', 'Pfs16PV_C'
    # ]
    train_data = '../data/all_coverage/read_cov.csv'
    # gene_data = '../data/gene_cov/read_test_cov.csv'
    # coverage_df = pd.read_csv(gene_data, index_col='mrna_id')
    train_df = pd.read_csv(train_data, index_col='mrna_id')
    return train_df


def prepare_data(df):
    offset = 1
    offset_para = 0.7

    # Lose all the genes that have a mean coverage lower than threshold
    mean_thr = 1
    # df = df[df.mean(axis=1) > mean_thr]

    plt.plot(range(0, len(df.Pfs1)), df.Pfs1.sort_values())
    plt.show()
    import sys; sys.exit()
    # Normalize the mean and standard deviation between the isolates
    # normalized_df = preprocessing.scale(df, axis=1)



    # Divide all the genes by the mean over all isolates
    # corrected_df = df.div(df.mean(axis=1), axis=0)
    # corrected_df.plot.density()
    # plt.savefig('../results/density_plot.png', dpi=900)

    df[df < offset_para] = int(0)
    # df[df >= offset_para] = int(1)
    # df[df < offset] = np.NaN
    df[(df < offset) & (df >= offset_para)] = np.NaN
    df = df.dropna()
    # corrected_df[corrected_df > offset_para] = int(2)
    return df


def main():
    train_df = import_data()
    #prepared_df = prepare_data(coverage_df)
    train_df = prepare_data(train_df)

    # Correction for equal features only use with test set
    # valid_rows = prepared_df.loc[prepared_df['ES1314'].notnull()].index.values
    # prepared_df = prepared_df.loc[valid_rows]
    # train_df = train_df.loc[valid_rows]

    # prepared_df.to_csv('../data/prepared_test_dataset.csv')
    train_df.to_csv('../data/prepared_dataset.csv')


if __name__ == '__main__':
    main()
