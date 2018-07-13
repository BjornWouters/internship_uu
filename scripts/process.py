import pandas as pd
import numpy as np
from sklearn.model_selection import LeaveOneOut
from sklearn.utils import resample
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn import linear_model


def import_data():
    gene_data = '../data/prepared_dataset.csv'
    resistance_data = '../data/resistance_matrix.csv'
    coverage_df = pd.read_csv(gene_data, index_col=0)
    resistance_df = pd.read_csv(resistance_data, index_col=0)
    return resistance_df, coverage_df


def process_rows(differential, coverage_df):
    gene_df = pd.DataFrame()
    for i, name in enumerate(differential.index):
        gene_series = coverage_df[name]
        resistance = pd.Series(differential[i])
        resistance.index = ['resistance']
        gene_series = gene_series.append(resistance)
        gene_series.name = name
        gene_df = gene_df.append(gene_series)

    return gene_df


def resample_dataset(x_train, y_train, n_samples, balance_value):
    zero_results = list(y_train[y_train == balance_value].index)
    zero_df_train = x_train.loc[zero_results]
    zero_df_test = y_train.loc[zero_results]

    df_train_minority_upsampled = resample(zero_df_train,
                                           replace=True,
                                           n_samples=n_samples,
                                           random_state=123)
    df_test_minority_upsampled = resample(zero_df_test,
                                          replace=True,
                                          n_samples=n_samples,
                                          random_state=123)

    x_train = pd.concat([x_train, df_train_minority_upsampled])
    y_train = pd.concat([y_train, df_test_minority_upsampled])
    return x_train, y_train


def train_differntial_model(df, name, randomized):
    loo = LeaveOneOut()
    y = df['resistance']
    x = df.drop('resistance', 1)
    # clf = RandomForestClassifier(
    #     max_depth=4, random_state=0, criterion='gini',
    #     n_estimators=10)
    clf = GradientBoostingClassifier(
        random_state=0,
        max_depth=3,
        n_estimators=100,
        min_samples_leaf=6,
        learning_rate=0.3,
    )
    # clf = linear_model.SGDClassifier(loss='log', penalty='elasticnet', random_state=0)
    count = 0
    correct = 0

    for train, test in loo.split(df):
        y_train = y[train]
        y_test = y[test]
        x_train = x.iloc[train]
        x_test = x.iloc[test]

        try:
            one_count = y_train.value_counts()[1]
            zero_count = y_train.value_counts()[0]
        except KeyError:
            continue

        # Resample if unequal class distribution
        if one_count > zero_count:
            x_train, y_train = resample_dataset(x_train, y_train, one_count, 0)
        elif one_count < zero_count:
            x_train, y_train = resample_dataset(x_train, y_train, zero_count, 1)

        clf.fit(x_train, y_train)
        if clf.predict(x_test)[0] == y_test[0]:
            correct += 1
        count += 1

    score = correct/count*100

    if not randomized:
        top_5_features_index = np.argsort(clf.feature_importances_)[-5:]
        print('Partial score of differential ' + name + ': ' + str(score) + '%')
        print(top_5_features_index)
        print('Probabilities: ' + ' '.join([str(clf.feature_importances_[i])
                                            for i in top_5_features_index]))

    return score


def main():
    resistance_df, coverage_df = import_data()
    total_score = 0
    iterations = 5
    randomize = False
    for i in range(iterations):
        if i >= 1:
            randomize = True

        for differential in resistance_df.iterrows():
            try:
                row = differential[1]

                if randomize:
                    row = differential[1].sample(frac=1)

                # Check for valid rows
                count_one = row.value_counts()[1]
                count_zero = row.value_counts()[0]

                row_df = process_rows(row, coverage_df)
                score = train_differntial_model(row_df, differential[0], randomize)
                total_score += score
            except KeyError:
                continue

        if not randomize:
            print('\nTotal accuracy: ' + str(total_score/len(resistance_df[1:])) + '%')
        else:
            print('\nTotal accuracy (random): ' + str(total_score/len(resistance_df[1:])) + '%')

        total_score = 0


if __name__ == '__main__':
    main()
