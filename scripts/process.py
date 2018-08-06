# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Version: 1.0.0
# Author: Björn Wouters
# Date: 16-7-2018
# Function:
# Learning and testing susceptibility of differentials based on the genomic data of the different
# races of Pfs through a machine learning algorithm.
# Known bugs: None found yet
# # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import auc
from sklearn.model_selection import LeaveOneOut
from sklearn.utils import resample


class Roc:
    """
    Roc class
    Attributes and function for calculating Roc plot features
    """
    def __init__(self):
        self.tp = 0
        self.fp = 0
        self.tpr_list = list()
        self.fpr_list = list()
        self.probabilities = list()

    def calc_roc(self):
        """
        Calculates tpr and fpr based on the sorted probabilities.
        """
        counted_values = [value[1] for value in self.probabilities]
        scores = [value[0] for value in self.probabilities]
        min_score = min(scores)
        max_score = max(scores)
        # Divides the probabilities into 500 equal distributed values
        thr = np.linspace(min_score, max_score, 500)
        total_positives = sum(counted_values)
        total_negatives = len(counted_values) - total_positives

        for T in thr:
            for i in range(0, len(scores)):
                if scores[i] > T:
                    if counted_values[i] == 1:
                        self.tp += 1
                    else:
                        self.fp += 1

            self.tpr_list.append(self.tp/float(total_positives))
            self.fpr_list.append(self.fp/float(total_negatives))
            self.tp = 0
            self.fp = 0


def import_data():
    """
    Data import function, converts csv to pandas dataframes
    """
    gene_data = '../data/parsed_snps.csv'
    resistance_data = '../data/resistance_matrix.csv'
    test_data = '../data/prepared_test_dataset.csv'
    coverage_df = pd.read_csv(gene_data, index_col=0)
    resistance_df = pd.read_csv(resistance_data, index_col=0)
    test_df = pd.read_csv(test_data, index_col=0).transpose()
    # Commend this line if the test set has to be predicted
    test_df = pd.DataFrame()
    return resistance_df, coverage_df, test_df


def process_rows(differential, coverage_df):
    """
    Calculates the resistance outcomes for testing purposes
    :param differential: Differential series of the resistance matrix
    :param coverage_df: abs/pres dataset pre-calculated
    :return: Matrix with classes for ML training purpose
    """
    gene_df = pd.DataFrame()
    variants = ['_snp']

    for variant in variants:
        for i, name in enumerate(differential.index):
            gene_series = coverage_df[name + variant]
            resistance = pd.Series(differential[i])
            resistance.index = ['resistance']
            gene_series = gene_series.append(resistance)
            gene_series.name = name + variant
            gene_df[name + variant] = gene_series

    return gene_df.transpose()


def resample_dataset(x_train, y_train, n_samples, balance_value):
    """
    :param x_train: All the training values of the x trainingset
    :param y_train:  All the classes (outcomes) of the trainingset
    :param n_samples: Number of samples to be resampled
    :param balance_value: Either the 0 or 1 to be resampled
    :return: Resampled x and y dataset.
    """
    zero_results = list(y_train[y_train == balance_value].index)
    zero_df_train = x_train.loc[zero_results]
    zero_df_test = y_train.loc[zero_results]

    # Both the x and y dataset has to be resampled
    df_train_minority_upsampled = resample(zero_df_train,
                                           replace=True,
                                           n_samples=n_samples,
                                           random_state=123)
    df_test_minority_upsampled = resample(zero_df_test,
                                          replace=True,
                                          n_samples=n_samples,
                                          random_state=123)

    # Adds the resampled instances to the actual dataframe
    x_train = pd.concat([x_train, df_train_minority_upsampled])
    y_train = pd.concat([y_train, df_test_minority_upsampled])
    return x_train, y_train


def fit_classifier(clf, x_train, y_train):
    """
    Train the given classifier
    :param clf: Classifier to be fitted
    :param x_train: Values training set
    :param y_train: Prediction classes training set
    """
    one_count = y_train.value_counts()[1]
    zero_count = y_train.value_counts()[0]

    # Resample if unequal class distribution
    if one_count > zero_count:
        x_train, y_train = resample_dataset(x_train, y_train, one_count, 0)
    elif one_count < zero_count:
        x_train, y_train = resample_dataset(x_train, y_train, zero_count, 1)

    clf.fit(x_train, y_train)


def train_differntial_model(df, predict_df, name, roc, randomized, test_data=pd.DataFrame()):
    """
    :param df: Input dataframe; the genes i.c.w. their outcome (1/0)
    :param predict_df: Dataframe to save the outcomes of the predictions
    :param name: Name of the current differential
    :param roc: Roc class declared in the main function
    :param randomized: Boolean of there are multiple iterations (i > 1 = True)
    :param test_data: Dataframe of the non-pfs isolates (optional)
    :return:
    """
    loo = LeaveOneOut()
    y = df['resistance']
    x = df.drop('resistance', 1)
    clf = GradientBoostingClassifier(
        random_state=0,
        max_depth=2,
        n_estimators=200,
        min_samples_leaf=6,
        learning_rate=0.2,
        subsample=0.8
    )
    # clf = RandomForestClassifier(
    #     random_state=0,
    #     n_estimators=100,
    # )
    count = 0
    correct = 0

    # Only run this a test set of non-pfs isolates is given
    if not test_data.empty:
        fit_classifier(clf, x, y)
        prediction = pd.Series(clf.predict(test_data))
        prediction.index = test_data.index
        prediction.name = name
        predict_df[name] = prediction
        return 0, roc

    # Iterates while leaving out one class at the time
    for train, test in loo.split(df):
        y_train = y[train]
        y_test = y[test]
        x_train = x.iloc[train]
        x_test = x.iloc[test]

        try:
            fit_classifier(clf, x_train, y_train)
        except KeyError:
            # Only possible if there is only one class that can be predicted (skip)
            continue

        prediction = clf.predict(x_test)[0]
        actual = y_test[0]
        probability = clf.predict_proba(x_test).max()

        # True positive
        if prediction == actual:
            correct += 1
            roc.probabilities.append([probability, 1])
        # False positive
        else:
            roc.probabilities.append([probability, 0])

        # Set prediction dataframe
        gene_name = y_test.index[0]

        # Only first non-random iteration
        if not randomized:
            # True positive
            if actual == 1 and prediction == 1:
                predict_df.loc[name, gene_name] = 2
            # False positive
            elif actual == 0 and prediction == 1:
                predict_df.loc[name, gene_name] = 1
            # False negative
            elif actual == 1 and prediction == 0:
                predict_df.loc[name, gene_name] = -1
            # True positive
            elif actual == 0 and prediction == 0:
                predict_df.loc[name, gene_name] = -2

            # predict_df.loc[name, gene_name] = probability

        count += 1

    score = correct/count*100

    if not randomized:
        top_5_features_index = np.argsort(clf.feature_importances_)[-5:]
        print('Partial score of differential ' + name + ': ' + str(score) + '%')
        print(x_train.iloc[:, top_5_features_index].columns.values)
        print('Probabilities: ' + ' '.join([str(clf.feature_importances_[i])
                                            for i in top_5_features_index]))

    return score, roc


def plot_auc(roc):
    roc.calc_roc()
    tpr = roc.tpr_list
    fpr = roc.fpr_list
    roc_auc = auc(fpr, tpr)
    plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
             lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic')
    plt.legend(loc="lower right")
    plt.show()


def main():
    """
    Main function
    """
    # Number of iterations. > 1 means you compare it against permuted resistance matrices
    iterations = 1

    resistance_df, coverage_df, test_df = import_data()
    resistance_predicted_df = pd.DataFrame()
    total_score = 0
    roc = Roc()
    randomize = False

    for i in range(iterations):
        # Don't permute the first iteration
        if i >= 1:
            randomize = True

        # Start from the second row, viroflay isn't interesting
        for differential in resistance_df[1:].iterrows():
            row = differential[1]
            if randomize:
                row = differential[1].sample(frac=1)

            processed_coverage = process_rows(row, coverage_df)

            score, roc = train_differntial_model(processed_coverage, resistance_predicted_df,
                                                 differential[0], roc, randomize, test_df)

            total_score += score

        if not randomize:
            print('\nTotal accuracy: ' + str(total_score/len(resistance_df[1:])) + '%')
            # Show the predicted differentials
            if not test_df.empty:
                sns.heatmap(resistance_predicted_df.transpose())
                plt.show()
                sys.exit(0)
            sns.heatmap(resistance_predicted_df)
            plt.title('Predicted resistance matrix')
            # plt.show()
            plot_auc(roc)
        else:
            print('\nTotal accuracy (random): ' + str(total_score/len(resistance_df[1:])) + '%')

        total_score = 0


if __name__ == '__main__':
    main()
