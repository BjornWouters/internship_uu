# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Version: 1.0.0
# Author: Björn Wouters
# Date: 16-7-2018
# Function:
# Learning and testing susceptibility of differentials based on the genomic data of the different
# races of Pfs through a machine learning algorithm.
# Known bugs: None found yet
# # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import csv
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pprint
import seaborn as sns
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier, \
    RandomForestRegressor
from sklearn.metrics import auc
from sklearn.model_selection import LeaveOneOut
from sklearn.neural_network import MLPClassifier
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
        # Divides the probabilities into n(scores) equal distributed values
        thr = np.linspace(min_score, max_score, len(scores))
        total_positives = sum(counted_values)
        total_negatives = len(counted_values) - total_positives

        for T in thr:
            for i in range(0, len(scores)):
                if scores[i] > T:
                    if counted_values[i] == 1:
                        self.tp += 1
                    else:
                        self.fp += 1

            self.tpr_list.append(self.tp / float(total_positives))
            self.fpr_list.append(self.fp / float(total_negatives))
            self.tp = 0
            self.fp = 0


class Results:
    def __init__(self):
        self.final_results = dict()

    def update_results(self, results, differential):
        if differential not in self.final_results:
            self.final_results.update({differential: dict()})

        for gene, prob in results:
            if gene in self.final_results[differential]:
                self.final_results[differential][gene] += prob
            else:
                self.final_results[differential].update({gene: prob})

    def write_results(self, file):
        processed_genes = set()
        with open(file, 'w') as csvfile:
            differentials = self.final_results.keys()
            header = ['gene'] + list(differentials)
            writer = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_NONE, quotechar='',
                                escapechar='\\')
            writer.writerow(header)
            for differential in differentials:
                for gene in self.final_results[differential]:
                    if gene not in processed_genes:
                        output = [0] * len(differentials)
                        for i, differential_check in enumerate(differentials):
                            if gene in self.final_results[differential_check]:
                                output[i] = self.final_results[differential_check][gene]

                        processed_genes.add(gene)
                        writer.writerow([gene] + output)


def import_data():
    """
    Data import function, converts csv to pandas dataframes
    """
    gene_data = '../data/all_features.csv'
    resistance_data = '../data/resistance_matrix.csv'
    weight_data = '../data/reliability_scores.csv'
    coverage_df = pd.read_csv(gene_data, index_col=0)
    resistance_df = pd.read_csv(resistance_data, index_col=0)
    test_df = coverage_df[[
        'A-03', 'Es-13', 'F-05', 'Nl-05', 'Us-11', 'Us-13a', 'Us-13b', 'Us-15'
    ]]
    return resistance_df, coverage_df, test_df, weight_data


def process_rows(differential, coverage_df):
    """
    Calculates the resistance outcomes for testing purposes
    :param differential: Differential series of the resistance matrix
    :param coverage_df: abs/pres dataset pre-calculated
    :return: Matrix with classes for ML training purpose
    """
    gene_df = pd.DataFrame()
    variants = ['']

    for variant in variants:
        for i, name in enumerate(differential.index):
            gene_series = coverage_df[name.title() + variant]
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


def train_differntial_model(df, predict_df, name, roc, randomized,
                            test_data=pd.DataFrame(), weight=None):
    """
    :param df: Input Dataframe; the genes i.c.w. their outcome (1/0)
    :param predict_df: Dataframe to save the outcomes of the predictions
    :param name: Name of the current differential
    :param roc: Roc class declared in the main function
    :param randomized: Boolean of there are multiple iterations (i > 1 = True)
    :param test_data: Dataframe of the non-pfs isolates (optional)
    :return:
    """
    feature_dict = dict()
    loo = LeaveOneOut()
    y = df['resistance']
    x = df.drop('resistance', 1)
    # clf = GradientBoostingClassifier(
    #     random_state=0,
    #     n_estimators=100,
    # )
    clf = RandomForestClassifier(
        random_state=0,
        n_estimators=100,
        max_leaf_nodes=4,
        max_features=10,
    )
    # clf = MLPClassifier(hidden_layer_sizes=(15, 4), random_state=0)
    count = 0
    correct = 0

    # Only run this a test set of non-pfs isolates is given
    if not test_data.empty:
        fit_classifier(clf, x, y)
        prediction = pd.Series(clf.predict(test_data))
        prediction.index = test_data.index
        prediction.name = name
        if weight == 0.95:
            predict_df[name] = prediction
        else:
            predict_df[name] += prediction
        return 0, roc, None

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

        # Set prediction dataframe                row = differential[1].sample(frac=1)

        gene_name = y_test.index[0]

        # Only first non-random iteration
        if not randomized:

            if weight == 0.95:
                # True positive
                if prediction == 1:
                    predict_df.loc[name, gene_name] = 1
                else:
                    predict_df.loc[name, gene_name] = 0
            else:
                # True positive
                if prediction == 1:
                    predict_df.loc[name, gene_name] += 1
                else:
                    predict_df.loc[name, gene_name] += 0
            # # True positive
            # if actual == 1 and prediction == 1:
            #     predict_df.loc[name, gene_name] = 2
            # # False positive
            # elif actual == 0 and prediction == 1:
            #     predict_df.loc[name, gene_name] = 1
            # # False negative
            # elif actual == 1 and prediction == 0:
            #     predict_df.loc[name, gene_name] = -1
            # # True positive
            # elif actual == 0 and prediction == 0:
            #     predict_df.loc[name, gene_name] = -2

            # predict_df.loc[name, gene_name] = probability

        count += 1

        top_10_features_index = np.argsort(clf.feature_importances_)[-10:]
        top_10_genes = x_train.iloc[:, top_10_features_index].columns.values
        probabilities = [clf.feature_importances_[i] for i in top_10_features_index]
        for i, gene in enumerate(top_10_genes):
            if gene not in feature_dict:
                feature_dict.update({gene: probabilities[i]})
            else:
                feature_dict[gene] += probabilities[i]

    score = correct / count * 100

    if not randomized:
        # top_5_features_index = np.argsort(clf.feature_importances_)[-5:]
        results = sorted(feature_dict.items(), key=lambda kv: kv[1])
        # print('Partial score of differential ' + name + ': ' + str(score) + '%')
        # pprint.pprint(sorted(feature_dict.items(), key=lambda kv: kv[1]))
        #
        # print(x_train.iloc[:, top_5_features_index].columns.values)
        # print('Probabilities: ' + ' '.join([str(clf.feature_importances_[i])
        #                                     for i in top_5_features_index]))

    return score, roc, results


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
    results_output = '../results/gene_predictions_test.csv'
    resistance_df, coverage_df, test_df, weight_data = import_data()
    weight_df = pd.read_csv(weight_data, index_col='mrna_id', header=0)
    resistance_predicted_df = pd.DataFrame()
    total_score = 0
    roc = Roc()
    results = Results()
    randomize = False

    for weight in np.linspace(0.95, 0.50, 10):
        temp_weight_df = weight_df[weight_df.mean(axis=1) > weight]
        temp_coverage_df = coverage_df.loc[temp_weight_df.index].dropna()
        temp_test_df = test_df.loc[temp_weight_df.index].dropna().transpose()
        temp_test_df = pd.DataFrame()

        for i in range(iterations):
            # Don't permute the first iteration
            if i >= 1:
                randomize = True

            # Start from the second row, viroflay isn't interesting
            for differential in resistance_df[1:].iterrows():
                row = differential[1]
                if randomize:
                    row = differential[1].sample(frac=1)

                processed_coverage = process_rows(row, temp_coverage_df)

                score, roc, clf_results = train_differntial_model(
                    processed_coverage,
                    resistance_predicted_df,
                    differential[0],
                    roc,
                    randomize,
                    temp_test_df,
                    weight
                )

                if temp_test_df.empty:
                    results.update_results(clf_results, row.name)

                total_score += score

    if not randomize:
        print('\nTotal accuracy: ' + str(total_score / len(resistance_df[1:])) + '%')
        # Show the predicted differentials
        if not temp_test_df.empty:
            sns.heatmap(resistance_predicted_df.transpose())
            plt.show()
            sys.exit(0)

        results.write_results(results_output)
        sns.heatmap(resistance_predicted_df)
        plt.title('Predicted resistance matrix')
        plot_auc(roc)
    else:
        print('\nTotal accuracy (random): ' + str(total_score / len(resistance_df[1:])) + '%')

    total_score = 0


if __name__ == '__main__':
    main()
