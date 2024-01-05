#reference: https://towardsdatascience.com/stop-using-0-5-as-the-threshold-for-your-binary-classifier-8d7168290d44
#/home1/jialh/anaconda3/envs/pyg/bin/python
#https://github.com/ploomber/posts/blob/master/threshold/fit.ipynb
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn import datasets
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
from functools import partial

def metric_at_threshold(metric, y_test, y_scores, threshold):
    return metric(y_test, y_scores >= threshold)

def flagged(y_scores, threshold):
    return (y_scores >= threshold).sum()

def get_the_best_threshold(train_study, test_study, stage, model, feat_type, fs_method, y_test, y_prob, thresholds, outprefix):
    #accuracy_score_ = partial(accuracy_score, zero_division=1)
    precision_score_ = partial(precision_score, zero_division=1)
    recall_score_ = partial(recall_score, zero_division=1)
    f1_score_ = partial(f1_score, zero_division=1)

    accuracy = [metric_at_threshold(accuracy_score, y_test, y_prob, threshold)
                 for threshold in thresholds]
    precision = [metric_at_threshold(precision_score_, y_test, y_prob, threshold)
                 for threshold in thresholds]

    recall = [metric_at_threshold(recall_score_, y_test, y_prob, threshold)
              for threshold in thresholds]

    f1 = [metric_at_threshold(f1_score_, y_test, y_prob, threshold)
          for threshold in thresholds]

    n_flagged = [flagged(y_prob, threshold) for threshold in thresholds]

    df = pd.DataFrame({'threshold': thresholds,
                       'accuracy': accuracy,
                       'precision': precision,
                       'recall': recall,
                       'f1': f1,
                       'n_flagged': n_flagged})
    df.to_csv('{}_predict_metrics.csv'.format(outprefix))
    plt.plot(thresholds, accuracy, label='accuracy', color="red")
    plt.plot(thresholds, precision, label='precision', color="orange")
    plt.plot(thresholds, recall, label='recall', color="green")
    plt.plot(thresholds, f1, label='f1', color="blue")

    plt.xlabel('Threshold')
    plt.legend()
    plt.grid()
    plt.ylabel('Metric value')

    ax = plt.twinx()
    plt.ylabel('Flagged')
    ax.plot(thresholds, n_flagged, label='flagged', color='grey')
    ax.legend(loc=5)
    plt.savefig('{}_predict_metrics.png'.format(outprefix), dpi=300)
    plt.savefig('{}_predict_metrics.pdf'.format(outprefix), dpi=300)
    plt.clf()
    idx = np.argmax(f1)
    best_threshold = thresholds[idx]
    best_accuracy = accuracy[idx]
    best_f1 = f1[idx]
    best_precision = precision[idx]
    best_recall = recall[idx]
    result = [train_study, test_study, stage, model, feat_type, fs_method, best_threshold, best_accuracy, best_f1, best_precision, best_recall]
    return result

def main():
    workdir = "/home1/jialh/brain/01meta/multikingdom/06classification_20231115"
    train_studies = ['CHN+CHN2']
    #taxons = ['species', 'KOs', 'pathways', 'species+ConQuR', 'KOs+ConQuR', 'pathways+ConQuR', 'species+KOs+pathways', 'species+KOs+pathways+ConQuR']
    taxons = ['species+ConQuR', 'KOs+ConQuR', 'pathways+ConQuR', 'GMMs+ConQuR', 'GBMs+ConQuR']
    fs_methods = ['RandomForest+RFECV', 'RandomForest+RFECV_X2', 'RandomForest+RFECV_X4', 'all+RFECV',
                  'all+Boruta', 'MMUPHin+Boruta', 'MMUPHin+RFECV', 'MMUPHin_stage+RFECV']
    models = ['RandomForest']
    stages = ['AD', 'MCI', 'SCD', 'SCS']
    #thresholds = np.arange(0, 1, step=0.02)
    result_list = []
    performance_merge_df = pd.DataFrame()
    for train_study in train_studies:
        for taxon in taxons:
            if taxon == "species":
                feat_types = ['ABFV']
            elif taxon == 'species+ConQuR':
                #feat_types = ['ABFV+ConQuR']
                feat_types = ['A', 'B', 'F', 'V', 'AB', 'AF', 'AV', 'BF', 'BV', 'FV', 'ABF', 'ABV', 'AFV', 'BFV','ABFV']
            elif taxon in ['pathways', 'KOs', 'pathways+ConQuR', 'KOs+ConQuR', 'GMMs+ConQuR', 'GBMs+ConQuR']:
                feat_types = [taxon.strip().split('+')[0]]
            elif taxon == 'species+KOs+pathways':
                feat_types = ['ABFV+KOs+pathways']
            elif taxon == 'species+KOs+pathways+ConQuR':
                feat_types = ['ABFV+KOs+pathways+ConQuR']
            else:
                feat_types = [taxon]
            for feat_type in feat_types:
                for model in models:
                    inputdir = os.path.join(workdir, "{}_{}".format(train_study, taxon), model)
                    # thresholds_stat = pd.DataFrame(columns=['best_threshold', 'accuracy', 'f1', 'precision', 'recall'])
                    # result_list = []
                    for stage in stages:
                        for fs_method in fs_methods:
                            ####
                            result_prefix = "01results_{}_{}_{}_{}_{}_None".format(train_study, stage, model, feat_type, fs_method)
                            performance_file =  os.path.join(inputdir, "{}.csv".format(result_prefix))
                            if os.path.exists(performance_file):
                                performance_df = pd.read_csv(performance_file, header=0)
                                performance_df['taxon'] = taxon
                                if performance_merge_df.shape[0] == 0:
                                    performance_merge_df = performance_df
                                else:
                                    performance_merge_df = pd.concat([performance_merge_df, performance_df], axis=0)
                            else:
                                print("performance_file: {} does not exist.".format(performance_file))

                            # ###<------------------select the best threshold for precision and recall------------------------->###
                            # for test_study in ['GER', 'JPN']:
                            #     transfer_prefix = "transfer_from{}to{}_{}_{}_{}_{}_None".format(train_study, test_study, stage, model, feat_type, fs_method)
                            #     infile = os.path.join(inputdir, "{}_predict_metrix.csv".format(transfer_prefix))
                            #     if os.path.exists(infile):
                            #         print("deal with {}".format(infile))
                            #         data = pd.read_csv(infile, header=0, index_col=0)
                            #         y_test = data['label'].to_numpy()
                            #         y_prob = data['mean_proba'].to_numpy()
                            #         labels = np.unique(y_test)
                            #         outprefix = os.path.join(inputdir, transfer_prefix)
                            #         result = get_the_best_threshold(
                            #             train_study, test_study, stage, model, feat_type, fs_method, y_test, y_prob, thresholds, outprefix)
                            #         result_list.append(result)
    # result_df = pd.DataFrame(result_list, columns=[
    #     'train_study', 'test_study', 'stage', 'model', 'feat_type', 'fs_method',
    #     'best_threshold', 'accuracy', 'f1', 'precision', 'recall'])
    # result_df.to_csv(os.path.join(workdir, "thresholds_stat.csv"), header=True, index=False)
    performance_merge_df.to_csv(os.path.join(workdir, "performance_stat.csv"), header=True, index=False)

if __name__ == '__main__':
    main()