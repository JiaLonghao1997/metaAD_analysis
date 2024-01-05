from sklearn.feature_selection import RFECV
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import ElasticNet
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import GradientBoostingClassifier
from imblearn.ensemble import BalancedRandomForestClassifier
from sklearn.model_selection import GridSearchCV
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score,precision_score,recall_score, roc_auc_score, f1_score
#from script.utils.feature_selection import feature_selection
from sklearn.metrics import RocCurveDisplay
from sklearn.metrics import auc, roc_curve

from sklearn.feature_selection import RFE
from sklearn.metrics import roc_auc_score
import os
import numpy as np
import pandas as pd
import warnings
import matplotlib.pyplot as plt
from collections import Counter
warnings.filterwarnings('ignore')

def plot_auroc_curve(train_study, stage, model, feat_type, fs_method, sampler, auroc_curve_dic, outdir):
    #mean_fpr = np.linspace(0, 1, 100)
    fig, ax = plt.subplots(figsize=(6, 6))

    ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)
    AD_color_dic = {'CHN(CV)': '#000000', 'CHN2(CV)': '#000000', 'CHN+CHN2(CV)': '#000000',
                    'CHN+CHN2+GER(CV)': '#000000', 'CHN+CHN2+JPN(CV)': '#000000',
                 'CHN':'#E03837', 'CHN2':'#E03837',
                 'JPN':'#6AB42D', 'GER':'#03ACDA'}
    CRC_color_dic = {'AT-CRC(CV)': '#000000',
                     'CN-CRC(CV)': '#000000',
                     'DE-CRC(CV)': '#000000',
                     'FR-CRC(CV)': '#000000',
                     'US-CRC(CV)': '#000000',
                     'AT-CRC':'#ea716d',
                     'CN-CRC':'#9f9c23',
                     'DE-CRC':'#37b178',
                     'FR-CRC':'#49a1d0',
                     'US-CRC':'#b376b1'}
    ##<---------------------如何保存绘制AUROC曲线所需数据。----------------------------->##
    ##auroc_mean, auroc_std, mean_tpr, mean_fpr
    for key, value in auroc_curve_dic.items():
        study = key
        color_dic = {}
        if study in AD_color_dic.keys():
            color_dic = AD_color_dic
        elif study in CRC_color_dic.keys():
            color_dic = CRC_color_dic
        color = color_dic[study]
        auroc_mean = value[0]
        auroc_std = value[1]
        mean_tpr = value[2]
        mean_fpr = value[3]
        ax.plot(
            mean_fpr,
            mean_tpr,
            color=color,
            label=r"%s (AUC = %0.2f $\pm$ %0.2f)" % (study, auroc_mean, auroc_std),
            lw=2,
            alpha=0.8,
        )

    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        title="AUROC ({} {} {} {} {})".format(train_study, stage, model, feat_type, fs_method),
    )
    ax.legend(loc="lower right")

    plt.savefig(os.path.join(outdir, "auroc_{}_{}_{}_{}_{}_{}.png".format(
        train_study, stage, model, feat_type, fs_method, sampler)), dpi=300)
    plt.savefig(os.path.join(outdir,"auroc_{}_{}_{}_{}_{}_{}.pdf".format(
        train_study, stage, model, feat_type, fs_method, sampler)), dpi=300)
    #plt.show()

##pool.apply_async(func=RFECV_clf, args=(X, y, feature_list, outdir, stage, study, feat_type, model, repeats=1, fold=5))
##train_X, train_y, test_study_dic, feature_list, fs_method, outdir, stage, train_study, feat_type, model, repeats, fold
def CV_transfer(data, label, test_study_dic, sample_list, feature_list, fs_method, sampler,
                outdir, stage, train_study, test_studies, feat_type,
                n_features_in, n_features, selected_features,
                model, repeats, kfold):
    """
        model:'svm','randomforest','lasso'
        :return:
        """
    print('CrossValidation using {}'.format(model))
    ##https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html
    accuracy_list = []
    precision_list = []
    recall_list = []
    f1_list = []
    auroc_list = []
    tpr_list = []

    transfer_accuracy_list = {}
    transfer_precision_list = {}
    transfer_recall_list = {}
    transfer_f1_list = {}
    transfer_auroc_list = {}
    transfer_tpr_list = {}
    transfer_predict_matrix_dic = {}

    for test_study, value in test_study_dic.items():
        test_X = value['test_X']
        test_y = value['test_y']
        test_sample_list = value['test_sample_list']
        test_feature_list = value['test_feature_list']
        transfer_accuracy_list[test_study] = []
        transfer_precision_list[test_study] = []
        transfer_recall_list[test_study] = []
        transfer_f1_list[test_study] = []
        transfer_auroc_list[test_study] = []
        transfer_tpr_list[test_study] = []
        ###<-------------------预测矩阵--------------------->###
        predict_matrix = pd.DataFrame(data=np.zeros(shape=(test_X.shape[0], repeats)),
                                      index=test_sample_list,
                                      columns=['repeat' + str(i) for i in range(0, repeats)])

        # print("predict_matrix columns: {}".format(predict_matrix.columns.values.tolist()))
        predict_matrix['label'] = test_y.tolist()
        transfer_predict_matrix_dic[test_study] = predict_matrix

    # aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    predict_matrix = pd.DataFrame(data=np.zeros(shape=(data.shape[0], repeats)),
                                  index=['sample' + str(i) for i in range(0, data.shape[0])],
                                  columns=['repeat' + str(i) for i in range(0, repeats)])
    sample_list = predict_matrix.index.values.tolist()
    # print("predict_matrix columns: {}".format(predict_matrix.columns.values.tolist()))
    predict_matrix['label'] = label.tolist()
    # print("predict_matrix columns: {}".format(predict_matrix.columns.values.tolist()))

    ###<-----------记录特征的重要性--------->###
    feature_importance_df = pd.DataFrame(data=np.zeros(shape=(len(feature_list), repeats * kfold)),
                                         index=feature_list,
                                         columns=[str(i) for i in range(0, repeats * kfold)])

    # print("feature_importance columns: {}".format(feature_importance_df.columns.values.tolist()))

    ###参考：https://github.com/SegataLab/metaml/blob/master/classification.py
    for repeat in range(repeats):
        skf = StratifiedKFold(n_splits=kfold, shuffle=True, random_state=2022)
        ###<------------------------------------------------------------------>###
        for fold, (train_index, test_index) in enumerate(skf.split(data, label)):
            # print("{} repeat, {} fold, {} train".format(repeat + 1, fold + 1, repeat * kfold + fold + 1))
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
                rf = RandomForestClassifier(n_jobs=-1, class_weight="balanced")
                param_grid = {'n_estimators': [500]}
                grid = GridSearchCV(estimator=rf, param_grid=param_grid, cv=kfold, scoring='roc_auc', n_jobs=-1)
                grid.fit(train_X, train_y)
                param = grid.best_params_
                clf = RandomForestClassifier(n_estimators=param['n_estimators'], class_weight="balanced", n_jobs=-1)
                clf.fit(train_X, train_y)
                feature_importance_df[str(repeat * kfold + fold)] = clf.feature_importances_
            if model == "BalancedRandomForest":
                rf = BalancedRandomForestClassifier(n_jobs=-1, class_weight="balanced")
                param_grid = {'n_estimators': [500]}
                grid = GridSearchCV(estimator=rf, param_grid=param_grid, cv=kfold, scoring='roc_auc', n_jobs=-1)
                grid.fit(train_X, train_y)
                param = grid.best_params_
                clf = BalancedRandomForestClassifier(n_estimators=param['n_estimators'], class_weight="balanced", n_jobs=-1)
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
            y_proba = clf.predict_proba(test_X)[:, 1]  ##(nsample, nfeatures)
            predict_matrix.loc[np.array(sample_list)[test_index], 'repeat' + str(repeat)] = y_proba
            fpr, tpr, thresholds = roc_curve(test_y, y_proba)
            interp_tpr = np.interp(mean_fpr, fpr, tpr)  ##np.interp()用于线性插值。
            interp_tpr[0] = 0.0
            accuracy_list.append(accuracy_score(test_y, predict_y))
            precision_list.append(precision_score(test_y, predict_y))
            recall_list.append(recall_score(test_y, predict_y))
            f1_list.append(f1_score(test_y, predict_y))
            auroc_list.append(roc_auc_score(test_y, y_proba))
            tpr_list.append(interp_tpr)

        ##<--------------------------------数据集迁移-------------------------------------------->####
        transfer_clf = RandomForestClassifier(n_jobs=-1, class_weight="balanced")
        if model == 'SVM':
            svc = SVC(probability=True)
            param_grid = {'C': [0.0001, 0.001, 0.1, 1, 10, 100, 1000], 'gamma': [0.001, 0.0001],
                          'kernel': ['rbf', 'linear']}
            grid = GridSearchCV(estimator=svc, param_grid=param_grid, cv=kfold, scoring='roc_auc', n_jobs=-1)
            grid.fit(data, label)
            param = grid.best_params_
            transfer_clf = SVC(probability=True, C=param['C'], gamma=param['gamma'], kernel=param['kernel'])
            transfer_clf.fit(data, label)

        if model == 'RandomForest':
            rf = RandomForestClassifier(n_jobs=-1, class_weight="balanced")
            param_grid = {'n_estimators': [500]}
            grid = GridSearchCV(estimator=rf, param_grid=param_grid, cv=kfold, scoring='roc_auc', n_jobs=-1)
            grid.fit(data, label)
            param = grid.best_params_
            transfer_clf = RandomForestClassifier(n_estimators=param['n_estimators'], class_weight="balanced", n_jobs=-1)
            transfer_clf.fit(data, label)
            #feature_importance_df[str(repeat * kfold + fold)] = transfer_clf.feature_importances_
        if model == "BalancedRandomForest":
            rf = BalancedRandomForestClassifier(n_jobs=-1, class_weight="balanced")
            param_grid = {'n_estimators': [500]}
            grid = GridSearchCV(estimator=rf, param_grid=param_grid, cv=kfold, scoring='roc_auc', n_jobs=-1)
            grid.fit(data, label)
            param = grid.best_params_
            transfer_clf = BalancedRandomForestClassifier(n_estimators=param['n_estimators'], class_weight="balanced", n_jobs=-1)
            transfer_clf.fit(data, label)
            #feature_importance_df[str(repeat * kfold + fold)] = clf.feature_importances_

        if model == 'Lasso':
            lr = LogisticRegression(penalty='l1', solver='liblinear', max_iter=1000)
            param_grid = {'C': [0.01, 0.1, 1, 10, 100, 1000]}
            grid = GridSearchCV(estimator=lr, param_grid=param_grid, cv=kfold, scoring='roc_auc', n_jobs=-1)
            grid.fit(data, label)
            param = grid.best_params_
            transfer_clf = LogisticRegression(penalty='l1', max_iter=1000, solver='liblinear', C=param['C'])
            transfer_clf.fit(data, label)
            ##Feature importance: https://machinelearningmastery.com/calculate-feature-importance-with-python/
            #feature_importance_df[str(repeat * kfold + fold)] = transfer_clf.coef_.ravel().tolist()


        ####<---------------------迁移到其他数据集的预测结果-------------------->####
        for test_study, value in test_study_dic.items():
            test_X = value['test_X']
            test_y = value['test_y']
            test_sample_list = value['test_sample_list']
            test_feature_list = value['test_feature_list']
            if len(Counter(test_y)) > 1:
                proba_threshold = 0.5
                if test_study == "JPN":
                    proba_threshold = 0.7
                predict_y = transfer_clf.predict(test_X)
                test_y_proba = transfer_clf.predict_proba(test_X)[:, 1]
                transfer_predict_matrix_dic[test_study].loc[
                    np.array(test_sample_list), 'repeat' + str(repeat)] = test_y_proba
                ####
                fpr, tpr, thresholds = roc_curve(test_y, test_y_proba)
                interp_tpr = np.interp(mean_fpr, fpr, tpr)  ##np.interp()用于线性插值。
                interp_tpr[0] = 0.0
                transfer_accuracy_list[test_study].append(accuracy_score(test_y, test_y_proba>proba_threshold))
                transfer_precision_list[test_study].append(precision_score(test_y, test_y_proba>proba_threshold))
                transfer_recall_list[test_study].append(recall_score(test_y, test_y_proba>proba_threshold))
                transfer_f1_list[test_study].append(f1_score(test_y, test_y_proba>proba_threshold))
                try:
                    transfer_auroc_list[test_study].append(roc_auc_score(test_y, test_y_proba))
                except ValueError:
                    print("transfer_from_{}_to_{}_{}_{}_{}_{}_{} failed".format(train_study, test_study, stage, model, feat_type, fs_method, sampler))
                    pass
                #transfer_auroc_list[test_study].append(roc_auc_score(test_y, test_y_proba))
                transfer_tpr_list[test_study].append(interp_tpr)


    predict_matrix['mean_proba'] = predict_matrix.iloc[:, 0:repeats].mean(axis=1)
    # print(predict_matrix.columns.values.tolist())
    # print("predict_matrix columns: {}".format(predict_matrix.columns.values.tolist()))
    cv_prefix = "CV_{}_{}_{}_{}_{}_{}".format(train_study, stage, model, feat_type, fs_method, sampler)
    predict_matrix.to_csv(os.path.join(outdir, '{}_predict_metrix.csv'.format(cv_prefix)), index=False)
    ##reference: https://cmdlinetips.com/2018/01/how-to-create-pandas-dataframe-from-multiple-lists/
    classification_matrix = pd.DataFrame({'accuracy': accuracy_list,
                                          'precision': precision_list,
                                          'recall': recall_list,
                                          'f1': f1_list,
                                          'auroc': auroc_list})  ##True Positive Rate
    classification_matrix.to_csv(os.path.join(outdir, '{}_classification_metrix.csv'.format(cv_prefix)), index=False)
    np.savetxt(os.path.join(outdir, '{}_tpr.csv'.format(cv_prefix)), np.transpose(np.array(tpr_list)))
    ##
    feature_importance_df.to_csv(os.path.join(outdir, '{}_feature_importance.csv'.format(cv_prefix)), index=True,
                                 index_label="features")

    ###<=======================================================>###
    accuracy = np.array(accuracy_list).mean()
    precision = np.array(precision_list).mean()
    recall = np.array(recall_list).mean()
    f1 = np.array(f1_list).mean()
    auroc_mean = np.array(auroc_list).mean()
    auroc_std = np.array(auroc_list).std()
    mean_tpr = np.mean(tpr_list, axis=0)  ##np.mean(axis=0),压缩的是行。
    #mean_tpr[0] = 0
    mean_tpr[-1] = 1.0

    auroc_curve_dic = {}
    auroc_curve_dic[train_study+'(CV)'] = [auroc_mean, auroc_std, mean_tpr, mean_fpr]
    ###
    result_list = []
    result = [stage, train_study, train_study, model, feat_type, fs_method, sampler, accuracy, precision, recall, f1,  auroc_mean, auroc_std,
              n_features_in, n_features, selected_features]
    result_list.append(result)
    for test_study, value in test_study_dic.items():
        test_X = value['test_X']
        test_y = value['test_y']
        test_sample_list = value['test_sample_list']
        test_feature_list = value['test_feature_list']
        if len(Counter(test_y)) > 1:
            transfer_prefix = "transfer_from{}to{}_{}_{}_{}_{}_{}".format(train_study, test_study, stage, model, feat_type, fs_method, sampler)
            accuracy = np.array(transfer_accuracy_list[test_study]).mean()
            precision = np.array(transfer_precision_list[test_study]).mean()
            recall = np.array(transfer_recall_list[test_study]).mean()
            f1 = np.array(transfer_f1_list[test_study]).mean()
            auroc_mean = np.array(transfer_auroc_list[test_study]).mean()
            auroc_std = np.array(transfer_auroc_list[test_study]).std()
            #print("transfer_tpr_list[{}]: {}".format(test_study, transfer_tpr_list[test_study]))
            mean_tpr = np.mean(transfer_tpr_list[test_study], axis=0)  ##np.mean(axis=0),压缩的是行。
            #print("mean_tpr for test_study({}): {}".format(test_study, mean_tpr))
            #mean_tpr[0] = 0
            mean_tpr[-1] = 1.0
            np.savetxt(os.path.join(outdir, '{}_{}_tpr.csv'.format(transfer_prefix, test_study)), np.transpose(np.array(transfer_tpr_list[test_study])))
            predict_matrix = transfer_predict_matrix_dic[test_study]
            predict_matrix['mean_proba'] = predict_matrix.iloc[:, 0:repeats].mean(axis=1)
            predict_matrix.to_csv(os.path.join(outdir, '{}_predict_metrix.csv'.format(transfer_prefix)), index=False)
            auroc_curve_dic[test_study] = [auroc_mean, auroc_std, mean_tpr, mean_fpr]
            result = [stage, train_study, test_study, model, feat_type, fs_method, sampler, accuracy, precision, recall, f1, auroc_mean, auroc_std,
                      n_features_in, n_features, selected_features]
            result_list.append(result)

    plot_auroc_curve(train_study, stage, model, feat_type, fs_method, sampler, auroc_curve_dic, outdir)
    result_df = pd.DataFrame(result_list,
                             columns=['stage', 'train_study', 'test_study', 'model', 'feat_type', 'fs_method', 'sampler',
                                      'accuracy', 'precision', 'recall', 'f1', 'auroc_mean', 'auroc_std',
                                      'n_features_in', 'n_features', 'selected_features'])
    outfile = os.path.join(outdir, '01results_{}_{}_{}_{}_{}_{}.csv'.format(train_study, stage, model, feat_type, fs_method, sampler))
    result_df.to_csv(outfile, index=False)
    return result_df


