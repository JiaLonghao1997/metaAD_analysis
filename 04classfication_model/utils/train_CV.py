import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score,precision_score,recall_score, roc_auc_score
#from script.utils.feature_selection import feature_selection
from sklearn.metrics import RocCurveDisplay
from sklearn.metrics import auc, roc_curve
import matplotlib.pyplot as plt
import os
from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import train_test_split
import sys

## plot auroc: https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html
## plot auroc: https://machinelearningmastery.com/roc-curves-and-precision-recall-curves-for-classification-in-python
## keep probabilities for the positive outcome only

def plot_auroc_curve(study, stage, model, feat_type, fs_method, auroc_mean, auroc_std, mean_tpr, mean_fpr, outdir):
    #mean_fpr = np.linspace(0, 1, 100)
    fig, ax = plt.subplots(figsize=(6, 6))

    ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)
    ##<---------------------如何保存绘制AUROC曲线所需数据。----------------------------->##
    ax.plot(
        mean_fpr,
        mean_tpr,
        color="b",
        label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (auroc_mean, auroc_std),
        lw=2,
        alpha=0.8,
    )

    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        title="AUROC ({} {} {} {} {})".format(study, stage, model, feat_type, fs_method),
    )
    ax.legend(loc="lower right")

    plt.savefig(os.path.join(outdir, "auroc_{}_{}_{}_{}_{}.png".format(study, stage, model, feat_type, fs_method)),
                dpi=300)
    plt.show()


def train_cv(data, label, model, stage, study, sample_list, feature_list, feat_type, fs_method, n_features_in, n_features, selected_features, outdir, repeats, kfold):
    """
    model:'svm','randomforest','lasso'
    :return:
    """
    print('CrossValidation using {}'.format(model))
    accuracy_list = []
    precision_list = []
    recall_list = []
    auroc_list = []
    tpr_list = []
    #aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    predict_matrix = pd.DataFrame(data = np.zeros(shape=(data.shape[0], repeats)),
                                  index = sample_list,
                                  columns= ['repeat'+str(i) for i in range(0, repeats)])
    #print("predict_matrix columns: {}".format(predict_matrix.columns.values.tolist()))
    predict_matrix['label'] = label.tolist()
    #print("predict_matrix columns: {}".format(predict_matrix.columns.values.tolist()))

    ###<-----------记录特征的重要性--------->###
    feature_importance_df = pd.DataFrame(data=np.zeros(shape=(len(feature_list), repeats * kfold)),
                                         index=feature_list,
                                         columns=[str(i) for i in range(0, repeats * kfold)])

    #print("feature_importance columns: {}".format(feature_importance_df.columns.values.tolist()))


    ###参考：https://github.com/SegataLab/metaml/blob/master/classification.py
    for repeat in range(repeats):
        skf = StratifiedKFold(n_splits=kfold, shuffle=True, random_state=2022)
        for fold, (train_index, test_index) in enumerate(skf.split(data, label)):
            #print("{} repeat, {} fold, {} train".format(repeat + 1, fold + 1, repeat * kfold + fold + 1))
            train_X = data[train_index]
            train_y = label[train_index]
            test_X = data[test_index]
            test_y = label[test_index]
            clf = ''
            ##参数调优参考: https://www.kaggle.com/code/ldfreeman3/a-data-science-framework-to-achieve-99-accuracy/notebook#5.12-Tune-Model-with-Hyper-Parameters
            if model == 'SVM':
                svc = SVC(probability=True)
                param_grid = {'C': [0.0001, 0.001, 0.1, 1, 10, 100, 1000], 'gamma': [0.001, 0.0001],
                              'kernel': ['rbf', 'linear']}
                grid = GridSearchCV(estimator=svc, param_grid=param_grid, cv=kfold, scoring='roc_auc', n_jobs=-1)
                grid.fit(train_X, train_y)
                param = grid.best_params_
                clf = SVC(probability=True, C=param['C'], gamma=param['gamma'], kernel=param['kernel'])
                clf.fit(train_X, train_y)

            if model == 'RandomForest':
                rf = RandomForestClassifier(n_jobs=-1)
                param_grid = {'n_estimators': [100, 300, 500]}
                grid = GridSearchCV(estimator=rf, param_grid=param_grid, cv=kfold, scoring='roc_auc', n_jobs=-1)
                grid.fit(train_X, train_y)
                param = grid.best_params_
                clf = RandomForestClassifier(n_estimators=param['n_estimators'], n_jobs=-1)
                clf.fit(train_X, train_y)
                feature_importance_df[str(repeat * kfold + fold)] = clf.feature_importances_

            if model == 'Lasso':
                lr = LogisticRegression(penalty='l1', solver='liblinear', max_iter=1000)
                param_grid = {'C': [0.01, 0.1, 1, 10, 100, 1000]}
                grid = GridSearchCV(estimator=lr, param_grid=param_grid, cv=kfold, scoring='roc_auc', n_jobs=-1)
                grid.fit(train_X, train_y)
                param = grid.best_params_
                clf = LogisticRegression(penalty='l1', max_iter=1000, solver='liblinear', C=param['C'])
                clf.fit(train_X, train_y)
                ##Feature importance: https://machinelearningmastery.com/calculate-feature-importance-with-python/
                feature_importance_df[str(repeat * kfold + fold)] = clf.coef_.ravel().tolist()
            ###<------------------------交叉验证中，每一次训练时都画一条曲线---------------------------->###
            # feature_importance_df[str(repeat*kfold + fold)] = clf.feature_importances_

            predict_y = clf.predict(test_X)
            y_proba = clf.predict_proba(test_X)[:, 1] ##(nsample, nfeatures)
            predict_matrix.loc[np.array(sample_list)[test_index], 'repeat'+str(repeat)] = y_proba
            fpr, tpr, thresholds = roc_curve(test_y, y_proba)
            interp_tpr = np.interp(mean_fpr, fpr, tpr)  ##np.interp()用于线性插值。

            accuracy_list.append(accuracy_score(test_y, predict_y))
            precision_list.append(precision_score(test_y, predict_y))
            recall_list.append(recall_score(test_y, predict_y))
            auroc_list.append(roc_auc_score(test_y, y_proba))
            tpr_list.append(interp_tpr)

    predict_matrix['mean_proba'] = predict_matrix.iloc[:,0:repeats].mean(axis=1)
    #print(predict_matrix.columns.values.tolist())
    #print("predict_matrix columns: {}".format(predict_matrix.columns.values.tolist()))
    outprefix="CV_{}_{}_{}_{}_{}".format(study, stage, model, feat_type, fs_method)
    predict_matrix.to_csv(os.path.join(outdir, '{}_predict_metrix.csv'.format(outprefix)), index=False)
    ##reference: https://cmdlinetips.com/2018/01/how-to-create-pandas-dataframe-from-multiple-lists/
    classification_matrix = pd.DataFrame({'accuracy': accuracy_list,
                                          'precision': precision_list,
                                          'recall': recall_list,
                                          'auroc': auroc_list})   ##True Positive Rate
    classification_matrix.to_csv(os.path.join(outdir, '{}_classification_metrix.csv'.format(outprefix)), index=False)
    np.save(os.path.join(outdir, '{}_tpr.npy'.format(outprefix)), np.array(tpr_list))
    ##
    feature_importance_df.to_csv(os.path.join(outdir, '{}_feature_importance.csv'.format(outprefix)),index=True, index_label="features")

    ###<=======================================================>###
    accuracy = np.array(accuracy_list).mean()
    precision = np.array(precision_list).mean()
    recall = np.array(recall_list).mean()
    auroc_mean = np.array(auroc_list).mean()
    auroc_std = np.array(auroc_list).std()
    mean_tpr = np.mean(tpr_list, axis=0)   ##np.mean(axis=0),压缩的是行。
    mean_tpr[0] = 0
    mean_tpr[-1] = 1.0

    plot_auroc_curve(study, stage, model, feat_type, fs_method, auroc_mean, auroc_std, mean_tpr, mean_fpr, outdir)
    result = [stage, study, model, feat_type, fs_method, accuracy, precision, recall, auroc_mean, n_features_in, n_features, selected_features]
    result_df = pd.DataFrame([result],
                             columns=['stage', 'study', 'model', 'feat_type', 'fs_method',
                                      'accuracy', 'precision', 'recall', 'auroc',
                                      'n_features_in', 'n_features', 'selected_features'])
    outfile = os.path.join(outdir, '01CV_{}_{}_{}_{}_{}.csv'.format(study, stage, model, feat_type, fs_method))
    result_df.to_csv(outfile, index=False)
    return result_df

def main():
    workdir = "D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom"
    print("workdir: {}".format(workdir))
    X, y = load_breast_cancer(return_X_y=True, as_frame=True)
    feature_list = X.columns.values.tolist()
    sample_list = X.index.values.tolist()
    X = X.values
    y = y.to_numpy()
    #X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
    outdir = os.path.join(workdir, "testwork")
    stage = "BC"
    study = "sklearn"
    feat_type = "30features"
    model = "Lasso"
    fs_method = "all"
    repeats = 3
    kfold = 5
    n_features_in = len(feature_list)
    n_features = len(feature_list)
    selected_features = ",".join(feature_list)
    ##train_cv(data, label, model, stage, study, sample_list, feature_list, feat_type, fs_method, n_features_in, n_features, selected_features, outdir, repeats, kfold)
    result = train_cv(X, y, model, stage, study, sample_list, feature_list, feat_type, fs_method,
                      n_features_in, n_features, selected_features, outdir, repeats, kfold)
    print(result)


if __name__ == '__main__':
    main()