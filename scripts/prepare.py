import pandas as pd


def import_data():
    columns = [
        'genename', 'Pfs1RZ', 'Pfs2RZ', 'Pfs3RZ_C', 'Pfs4PV_C', 'Pfs5RZ_C', 'Pfs6RZ', 'Pfs7PV',
        'Pfs8RZ_C', 'Pfs9RZ', 'Pfs10RZ', 'Pfs11PV', 'Pfs12RZ', 'Pfs13PV', 'Pfs14RZ',
        'Pfs15RZ', 'Pfs16PV_C'
    ]
    gene_data = '../data/gene_cov/Pfs_genecov_all.tsv'
    coverage_df = pd.read_csv(gene_data, sep='\t', usecols=columns)
    return coverage_df


def prepare_data(df):
    df = df.set_index('genename')
    df[df > 0.02] = int(1)
    df[df <= 0.02] = int(0)
    return df


def main():
    coverage_df = import_data()
    prepared_df = prepare_data(coverage_df)
    prepared_df.to_csv('../data/prepared_dataset.csv')


if __name__ == '__main__':
    main()
