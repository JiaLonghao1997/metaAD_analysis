import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score,precision_score,recall_score,roc_auc_score,auc
from sklearn.metrics import RocCurveDisplay
import matplotlib.pyplot as plt

import os
import sys
sys.path.append("/scripts")
from prepare_data import data_preprocess


def study_to_study_transfer(train_data, train_label, train_study, test_data, test_label, test_study, stage, feat_type, model, fs, outdir, repeats=1, k_folds=5):
    print('Study to Study Transfer: training study({}), testing study:{}'.format(train_study, test_study))
    accuracy = []
    precision = []
    recall = []
    auroc = []
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    fig, ax = plt.subplots(figsize=(6, 6))
    # predict_matrix = pd.DataFrame(np.zeros(shape=(data.shape[0], 1)), index=sample_list,
    #                               columns=[stage])
    k = 0

    train_sample_id = train_data.index.values.tolist()
    x = train_data.values
    y = train_label.values
    # train_pre_matrix = pd.DataFrame(np.zeros(shape=(train_data.shape[0], 10)),
    #                                 index=train_sample_id, columns=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'])
    # right_num_study = 0
    # right_num_positive = 0
    # parameters = {'C':[0.1,1,5,10,100,1000]}
    for i in range(repeats):
        skf = StratifiedKFold(n_splits=k_folds)
        predict_result = np.zeros((len(test_label), k_folds))
        predict_prob = np.zeros((len(test_label), k_folds))
        for train_index, test_index in skf.split(x, y):
            train_x = x[train_index]
            train_y = y[train_index]
            test_x = x[test_index]
            test_y = y[test_index]

            """
                  remove features with std=0
                  reference: https://www.cnblogs.com/pinard/p/9032759.html
            """
            log_n0 = 1e-5
            sd_min_q = 10
            std_preprocess = np.std(train_x, axis=0)

            train_x = train_x[:, std_preprocess != 0]
            train_x = np.log10(train_x + log_n0)
            mean = np.mean(train_x, axis=0)
            std = np.std(train_x, axis=0)
            q = np.percentile(std, sd_min_q)
            train_x = (train_x - mean) / (std + q)

            test_x = test_x[:, std_preprocess != 0]
            test_x = np.log10(test_x + log_n0)
            test_x = (test_x - mean) / (std + q)

            clf = LogisticRegression(penalty='l1',solver='liblinear')
            param_grid = {'C': [0.01, 0.1, 1, 10, 100, 1000]}
            grid = GridSearchCV(clf, param_grid, cv=5, n_jobs=-1, pre_dispatch=5)
            grid.fit(train_x,train_y)
            param = grid.best_params_
            print("best param: ", param)

            clf = LogisticRegression(penalty='l1',solver='liblinear',C=param['C'])
            clf.fit(train_x,train_y)

            predict_proba = clf.predict_proba(test_x)[:, 1]
            predict_y = clf.predict(test_x)

            ###<--------------------------外部测试集---------------------------->###
            test_stst_index = test_data._stat_axis.values.tolist()
            test_stst_y = test_label.values

            test_stst_x = test_data.values #test_stst_x = test_stst_x.values
            test_stst_x = test_stst_x[:, std_preprocess != 0]
            test_stst_x = np.log10(test_stst_x + log_n0)
            test_stst_x = (test_stst_x - mean) / (std + q)

            viz = RocCurveDisplay.from_estimator(
                clf,
                test_stst_x,
                test_stst_y,
                name="ROC fold {}".format(k),
                alpha=0.3,
                lw=1,
                ax=ax,
            )
            interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)
            aucs.append(viz.roc_auc)

            test_stst_predict_y = clf.predict(test_stst_x)
            test_stst_y_proba = clf.predict_proba(test_stst_x)[:, 1]  ##(nsample, nfeatures)
            #testdata = data_df.iloc[[19, 20], :]
            #test_samples = data_df.iloc[test_index,].index.values.tolist()
            #predict_matrix.loc[data_df.iloc[test_index,].index.values.tolist(), stage] += y_proba / (cv * repeat)
            #predict_result[:, k] = np.array(predict_y).reshape(len(predict_y), 1)
            #predict_prob[:, k] = np.array(y_proba).reshape(len(predict_y), 1)
            k = k + 1

            #test_proba_y = clf.predict_proba(test_stst_x)[:, 1]
            accuracy.append(accuracy_score(test_stst_y, test_stst_predict_y))
            precision.append(precision_score(test_stst_y, test_stst_predict_y))
            recall.append(recall_score(test_stst_y, test_stst_predict_y))
            auroc.append(roc_auc_score(test_stst_y, test_stst_y_proba))

    accuracy = np.array(accuracy).mean()
    precision = np.array(precision).mean()
    recall = np.array(recall).mean()
    auroc = np.array(auroc).mean()
    ###<-------------------------------------------->###'
    ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(
        mean_fpr,
        mean_tpr,
        color="b",
        label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
        lw=2,
        alpha=0.8,
    )

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(
        mean_fpr,
        tprs_lower,
        tprs_upper,
        color="grey",
        alpha=0.2,
        label=r"$\pm$ 1 std. dev.",
    )

    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        title="Receiver operating characteristic example",
    )
    ax.legend(loc="lower right")

    plt.savefig(os.path.join(outdir, "auroc_{}_{}_{}.png".format(stage, model, fs)), dpi=300)
    plt.show()
    result = [stage, train_study, test_study, feat_type, model, fs, accuracy, precision, recall, auroc]
    return result


def train_transfer(train_data,train_label,test_data,test_label,repeat=10,type='function',model='lasso',is_TCA=True):
    print('Transfer training using {}'.format(model))
    if model == 'svm':
        svm = SVC(class_weight='balanced')
        param_grid = {'C': [0.0001, 0.001, 0.1, 1, 10, 100, 1000], 'gamma': [0.001, 0.0001],
                      'kernel': ['rbf', 'linear']}
        grid = GridSearchCV(svm, param_grid, cv=10, scoring='accuracy')
        grid.fit(train_data, train_label)
        param = grid.best_params_
        print(param)
        accuracy = []
        precision = []
        recall = []
        for i in range(repeat):
            svm = SVC(kernel=param['kernel'], C=param['C'], gamma=param['gamma'],class_weight='balanced')
            svm.fit(train_data, train_label)
            predict_y = svm.predict(test_data)
            # accuracy_result.append(accuracy_score(test_y,predict_y))
            # print(accuracy_result)
            accuracy.append(accuracy_score(test_label, predict_y))
            precision.append(precision_score(test_label, predict_y, average=None))
            recall.append(recall_score(test_label, predict_y, average=None))
        print("accuracy:", np.array(accuracy).mean())
        print("precision:", np.array(precision).mean(axis=0))
        print("recall:", np.array(recall).mean(axis=0))

    if model == 'randomforest':
        rf = RandomForestClassifier(n_jobs=-1,class_weight='balanced')
        param_grid = {'n_estimators': [10, 50, 100,500,1000]}
        grid = GridSearchCV(rf, param_grid, cv=10, scoring='accuracy')
        grid.fit(train_data,train_label)
        param = grid.best_params_
        print(param)
        accuracy = []
        precision=[]
        recall=[]
        for i in range(repeat):
            rf = RandomForestClassifier(n_estimators =param['n_estimators'],class_weight='balanced')
            rf.fit(train_data,train_label)
            predict_y = rf.predict(test_data)
            accuracy.append(accuracy_score(test_label,predict_y))
            precision.append(precision_score(test_label,predict_y,average=None))
            recall.append(recall_score(test_label,predict_y,average=None))
        print("accuracy:",np.array(accuracy).mean())
        print("precision:" ,np.array(precision).mean(axis=0))
        print("recall:" , np.array(recall).mean(axis=0))

    if model == 'lasso':
        lr = LogisticRegression(penalty='l1', solver='liblinear', n_jobs=-1,class_weight='balanced')
        param_grid = {'C':[ 0.0001,0.001, 0.1, 1, 10, 100, 1000]}
        grid = GridSearchCV(lr, param_grid, cv=10, scoring='accuracy')
        grid.fit(train_data,train_label)
        param = grid.best_params_
        print(param)
        accuracy = []
        precision=[]
        recall=[]
        for i in range(repeat):
            lr = LogisticRegression(penalty='l1', solver='liblinear', n_jobs=-1,C=param['C'],class_weight='balanced')
            lr.fit(train_data,train_label)
            predict_y = lr.predict(test_data)
            accuracy.append(accuracy_score(test_label,predict_y))
            precision.append(precision_score(test_label,predict_y,average=None))
            recall.append(recall_score(test_label,predict_y,average=None))
        print("accuracy:",np.array(accuracy).mean())
        print("precision:" ,np.array(precision).mean(axis=0))
        print("recall:" , np.array(recall).mean(axis=0))