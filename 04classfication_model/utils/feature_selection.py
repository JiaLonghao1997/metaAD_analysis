import pandas as pd
import numpy as np
import scipy.stats as stats
import pymrmr
import os
from math import sqrt
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_classif
from boruta import BorutaPy

def transfer_feature_selection(taxon, stage, train_X, train_y, test_X, feature_list, fs_method='all'):
    if fs_method == 'all':
        print("select all features.")
        return train_X, test_X, feature_list
    if fs_method == 't-test':
        choose_feature = []
        for i in range(train_X.shape[1]):
            if stats.ttest_ind(train_X[:,i][train_y==1],train_X[:,i][train_y==0])[1]<0.05:
                choose_feature.append(i)
        print("{} selected by {}".format(len(choose_feature), fs_method))
        feature_list = [feature_list[index] for index in choose_feature]
        return train_X[:,choose_feature], test_X[:,choose_feature], feature_list
    if fs_method == 'Wilcox':
        choose_feature = []
        for i in range(train_X.shape[1]):
            var = train_X[:,i].tolist()
            if stats.ranksums(train_X[:,i][train_y==1],train_X[:,i][train_y==0])[1]<0.05:
                choose_feature.append(i)
        print("{} selected by {}".format(len(choose_feature), fs_method))
        feature_list = [feature_list[index] for index in choose_feature]
        return train_X[:,choose_feature], test_X[:,choose_feature], feature_list
    if fs_method == 'mRMR':
        train_X_df = pd.DataFrame(train_X, columns=feature_list)
        train_y_df = pd.DataFrame(train_y, columns=['group'])
        merge_df = pd.merge(train_y_df, train_X_df, left_index=True, right_index=True)
        #merge_df.to_csv("merge.df.csv", header=True, index_label='Sample', index=True)
        choose_feature = pymrmr.mRMR(merge_df, 'MIQ', round(sqrt(train_X_df.shape[1])))
        print("{} selected by {}".format(len(choose_feature), fs_method))
        test_X_df = pd.DataFrame(test_X, columns=feature_list)
        # print(train_X_df[:, choose_feature].values)
        # print(test_X_df[:, choose_feature].values)
        return train_X_df[choose_feature].values, test_X_df[choose_feature].values, choose_feature
    elif fs_method == "RandomForest":
        train_X_df = pd.DataFrame(train_X, columns=feature_list)
        train_y_df = pd.DataFrame(train_y, columns=['group'])
        test_X_df = pd.DataFrame(test_X, columns=feature_list)
        sel = SelectFromModel(RandomForestClassifier(n_estimators=500, class_weight="balanced"), max_features=round(sqrt(train_X_df.shape[1])))
        sel.fit(train_X_df, train_y_df)
        choose_feature = train_X_df.columns[(sel.get_support())]
        return train_X_df[choose_feature].values, test_X_df[choose_feature].values, choose_feature
    elif fs_method == 'all+RFECV':
        print("select all features.")
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, class_weight="balanced",  max_features="sqrt", n_jobs=-1)
        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)
        rfecv = RFECV(estimator=estimator, step=0.1, cv=cv, scoring='accuracy', n_jobs=-1)
        rfecv = rfecv.fit(train_X, train_y)
        print("Optimal number of features : %d" % rfecv.n_features_)
        train_X = rfecv.transform(train_X)
        test_X = rfecv.transform(test_X)
        feature_list = np.array(feature_list)[rfecv.support_]
        return train_X, test_X, feature_list
    elif fs_method == "RandomForest+RFECV":
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        train_X_df = pd.DataFrame(train_X, columns=feature_list)
        train_y_df = pd.DataFrame(train_y, columns=['group'])
        test_X_df = pd.DataFrame(test_X, columns=feature_list)
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, class_weight="balanced", max_features="sqrt", n_jobs=-1)
        sel = SelectFromModel(estimator=estimator, max_features=round(sqrt(train_X_df.shape[1])))
        sel.fit(train_X_df, train_y_df)
        feature_list = train_X_df.columns[(sel.get_support())]
        train_X = train_X_df.loc[:, feature_list].values
        test_X = test_X_df.loc[:, feature_list].values
        ###<-----------------------去除缺失值-------------------------------->###
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        ###<---------------RFECV-------------->###
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, class_weight="balanced", max_features="sqrt", n_jobs=-1)
        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)
        rfecv = RFECV(estimator=estimator, step=0.1, cv=cv, scoring='accuracy', n_jobs=-1)
        rfecv = rfecv.fit(train_X, train_y)
        print("Optimal number of features : %d" % rfecv.n_features_)
        train_X = rfecv.transform(train_X)
        test_X = rfecv.transform(test_X)
        feature_list = np.array(feature_list)[rfecv.support_]
        return train_X, test_X, feature_list
    elif fs_method == "RandomForest+RFECV_X2":
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        train_X_df = pd.DataFrame(train_X, columns=feature_list)
        train_y_df = pd.DataFrame(train_y, columns=['group'])
        test_X_df = pd.DataFrame(test_X, columns=feature_list)
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, class_weight="balanced", max_features="sqrt", n_jobs=-1)
        sel = SelectFromModel(estimator=estimator, max_features=round(2*sqrt(train_X_df.shape[1])))
        sel.fit(train_X_df, train_y_df)
        feature_list = train_X_df.columns[(sel.get_support())]
        train_X = train_X_df.loc[:, feature_list].values
        test_X = test_X_df.loc[:, feature_list].values
        ###<-----------------------去除缺失值-------------------------------->###
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        ###<---------------RFECV-------------->###
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, class_weight="balanced", max_features="sqrt", n_jobs=-1)
        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)
        rfecv = RFECV(estimator=estimator, step=0.1, cv=cv, scoring='accuracy', n_jobs=-1)
        rfecv = rfecv.fit(train_X, train_y)
        print("Optimal number of features : %d" % rfecv.n_features_)
        train_X = rfecv.transform(train_X)
        test_X = rfecv.transform(test_X)
        feature_list = np.array(feature_list)[rfecv.support_]
        return train_X, test_X, feature_list
    elif fs_method == "RandomForest+RFECV_X4":
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        train_X_df = pd.DataFrame(train_X, columns=feature_list)
        train_y_df = pd.DataFrame(train_y, columns=['group'])
        test_X_df = pd.DataFrame(test_X, columns=feature_list)
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, class_weight="balanced", max_features="sqrt", n_jobs=-1)
        sel = SelectFromModel(estimator=estimator, max_features=round(4*sqrt(train_X_df.shape[1])))
        sel.fit(train_X_df, train_y_df)
        feature_list = train_X_df.columns[(sel.get_support())]
        train_X = train_X_df.loc[:, feature_list].values
        test_X = test_X_df.loc[:, feature_list].values
        ###<-----------------------去除缺失值-------------------------------->###
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        ###<---------------RFECV-------------->###
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, class_weight="balanced", max_features="sqrt", n_jobs=-1)
        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)
        rfecv = RFECV(estimator=estimator, step=0.1, cv=cv, scoring='accuracy', n_jobs=-1)
        rfecv = rfecv.fit(train_X, train_y)
        print("Optimal number of features : %d" % rfecv.n_features_)
        train_X = rfecv.transform(train_X)
        test_X = rfecv.transform(test_X)
        feature_list = np.array(feature_list)[rfecv.support_]
        return train_X, test_X, feature_list
    elif fs_method == "RandomForest+RFECV20":
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        train_X_df = pd.DataFrame(train_X, columns=feature_list)
        train_y_df = pd.DataFrame(train_y, columns=['group'])
        test_X_df = pd.DataFrame(test_X, columns=feature_list)
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, class_weight="balanced", max_features="sqrt", n_jobs=-1)
        sel = SelectFromModel(estimator=estimator, max_features=20)
        sel.fit(train_X_df, train_y_df)
        feature_list = train_X_df.columns[(sel.get_support())]
        train_X = train_X_df.loc[:, feature_list].values
        test_X = test_X_df.loc[:, feature_list].values
        ###<-----------------------去除缺失值-------------------------------->###
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        ###<---------------RFECV-------------->###
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, class_weight="balanced", max_features="sqrt", n_jobs=-1)
        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)
        rfecv = RFECV(estimator=estimator, step=0.1, cv=cv, scoring='accuracy', n_jobs=-1)
        rfecv = rfecv.fit(train_X, train_y)
        print("Optimal number of features : %d" % rfecv.n_features_)
        train_X = rfecv.transform(train_X)
        test_X = rfecv.transform(test_X)
        feature_list = np.array(feature_list)[rfecv.support_]
        return train_X, test_X, feature_list
    elif fs_method == "RandomForest+RFECV80":
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        train_X_df = pd.DataFrame(train_X, columns=feature_list)
        train_y_df = pd.DataFrame(train_y, columns=['group'])
        test_X_df = pd.DataFrame(test_X, columns=feature_list)
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, class_weight="balanced", max_features="sqrt", n_jobs=-1)
        sel = SelectFromModel(estimator=estimator, max_features=80)
        sel.fit(train_X_df, train_y_df)
        feature_list = train_X_df.columns[(sel.get_support())]
        train_X = train_X_df.loc[:, feature_list].values
        test_X = test_X_df.loc[:, feature_list].values
        ###<-----------------------去除缺失值-------------------------------->###
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        ###<---------------RFECV-------------->###
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, class_weight="balanced", max_features="sqrt", n_jobs=-1)
        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)
        rfecv = RFECV(estimator=estimator, step=0.1, cv=cv, scoring='accuracy', n_jobs=-1)
        rfecv = rfecv.fit(train_X, train_y)
        print("Optimal number of features : %d" % rfecv.n_features_)
        train_X = rfecv.transform(train_X)
        test_X = rfecv.transform(test_X)
        feature_list = np.array(feature_list)[rfecv.support_]
        return train_X, test_X, feature_list
    elif fs_method == "MMUPHin+RFECV":
        taxon_prefix = taxon.strip().split('+')[0]
        MMUPHin_diff_dir = "/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff"
        infile = os.path.join(MMUPHin_diff_dir, taxon_prefix,  "MMUPHin_diff_merge.csv")
        MMUPHin_df = pd.read_csv(infile, header=0, index_col=0)
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        train_X_df = pd.DataFrame(train_X, columns=feature_list)
        train_y_df = pd.DataFrame(train_y, columns=['group'])
        test_X_df = pd.DataFrame(test_X, columns=feature_list)
        feature_list = MMUPHin_df.loc[MMUPHin_df.index.isin(train_X_df.columns.values.tolist()), ].index.values.tolist()
        print("feature_list from MMUPHin: {}".format(len(feature_list)))
        train_X = train_X_df.loc[:, feature_list].values
        test_X = test_X_df.loc[:,feature_list].values
        ###<-----------------------去除缺失值-------------------------------->###
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        ###<---------------RFECV-------------->###
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, class_weight="balanced", max_features="sqrt", n_jobs=-1)
        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)
        rfecv = RFECV(estimator=estimator, step=0.1, cv=cv, scoring='accuracy', n_jobs=-1)
        rfecv = rfecv.fit(train_X, train_y)
        print("Optimal number of features : %d" % rfecv.n_features_)
        train_X = rfecv.transform(train_X)
        test_X = rfecv.transform(test_X)
        feature_list = np.array(feature_list)[rfecv.support_]
        return train_X, test_X, feature_list
    elif fs_method == "MMUPHin_stage+RFECV":
        taxon_prefix = taxon.strip().split('+')[0]
        MMUPHin_diff_dir = "/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff"
        infile = os.path.join(MMUPHin_diff_dir, taxon_prefix,  "MMUPHin_diff_merge.csv")
        MMUPHin_df = pd.read_csv(infile, header=0, index_col=0)
        MMUPHin_df_stage = MMUPHin_df[MMUPHin_df[stage].notna()]
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        train_X_df = pd.DataFrame(train_X, columns=feature_list)
        train_y_df = pd.DataFrame(train_y, columns=['group'])
        test_X_df = pd.DataFrame(test_X, columns=feature_list)
        feature_list = MMUPHin_df_stage.loc[MMUPHin_df_stage.index.isin(train_X_df.columns.values.tolist()), ].index.values.tolist()
        print("feature_list from MMUPHin: {}".format(len(feature_list)))
        train_X = train_X_df.loc[:, feature_list].values
        test_X = test_X_df.loc[:,feature_list].values
        ###<-----------------------去除缺失值-------------------------------->###
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        ###<---------------RFECV-------------->###
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, class_weight="balanced", max_features="sqrt", n_jobs=-1)
        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)
        rfecv = RFECV(estimator=estimator, step=0.1, cv=cv, scoring='accuracy', n_jobs=-1)
        rfecv = rfecv.fit(train_X, train_y)
        print("Optimal number of features : %d" % rfecv.n_features_)
        train_X = rfecv.transform(train_X)
        test_X = rfecv.transform(test_X)
        feature_list = np.array(feature_list)[rfecv.support_]
        return train_X, test_X, feature_list
    elif fs_method == "MMUPHin2+RFECV":
        taxon_prefix = taxon.strip().split('+')[0]
        MMUPHin_diff_dir = "/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff"
        infile = os.path.join(MMUPHin_diff_dir, taxon_prefix, "MMUPHin_diff_merge.csv")
        MMUPHin_df = pd.read_csv(infile, header=0, index_col=0)
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        train_X_df = pd.DataFrame(train_X, columns=feature_list)
        train_y_df = pd.DataFrame(train_y, columns=['group'])
        test_X_df = pd.DataFrame(test_X, columns=feature_list)
        feature_list = MMUPHin_df.loc[(MMUPHin_df.index.isin(train_X_df.columns.values.tolist()) & (MMUPHin_df['Freq']>=2)),].index.values.tolist()
        print("feature_list from MMUPHin: {}".format(len(feature_list)))
        train_X = train_X_df.loc[:, feature_list].values
        test_X = test_X_df.loc[:, feature_list].values
        ###<-----------------------去除缺失值-------------------------------->###
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        ###<---------------RFECV-------------->###
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, class_weight="balanced", max_features="sqrt", n_jobs=-1)
        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)
        rfecv = RFECV(estimator=estimator, step=0.1, cv=cv, scoring='accuracy', n_jobs=-1)
        rfecv = rfecv.fit(train_X, train_y)
        print("Optimal number of features : %d" % rfecv.n_features_)
        train_X = rfecv.transform(train_X)
        test_X = rfecv.transform(test_X)
        feature_list = np.array(feature_list)[rfecv.support_]
        return train_X, test_X, feature_list
    elif fs_method == "SelectNonCollinear+RandomForest+RFECV":
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        train_X_df = pd.DataFrame(train_X, columns=feature_list)
        train_y_df = pd.DataFrame(train_y, columns=['group'])
        test_X_df = pd.DataFrame(test_X, columns=feature_list)
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, class_weight="balanced", max_features="sqrt", n_jobs=-1)
        sel = SelectFromModel(estimator=estimator, max_features=round(sqrt(train_X_df.shape[1])))
        sel.fit(train_X_df, train_y_df)
        feature_list = train_X_df.columns[(sel.get_support())]
        train_X = train_X_df.loc[:, feature_list].values
        test_X = test_X_df.loc[:, feature_list].values
        ####
        
        ###<-----------------------去除缺失值-------------------------------->###
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        ###<---------------RFECV-------------->###
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, class_weight="balanced", max_features="sqrt", n_jobs=-1)
        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)
        rfecv = RFECV(estimator=estimator, step=0.1, cv=cv, scoring='accuracy', n_jobs=-1)
        rfecv = rfecv.fit(train_X, train_y)
        print("Optimal number of features : %d" % rfecv.n_features_)
        train_X = rfecv.transform(train_X)
        test_X = rfecv.transform(test_X)
        feature_list = np.array(feature_list)[rfecv.support_]
        return train_X, test_X, feature_list
    elif fs_method == "all+Boruta":
        ##Feature Selection with Boruta in Python: https://towardsdatascience.com/feature-selection-with-boruta-in-python-676e3877e596
        ##package code: https://github.com/scikit-learn-contrib/boruta_py
        train_X = np.nan_to_num(train_X)
        test_X = np.nan_to_num(test_X)
        estimator = RandomForestClassifier(n_estimators=500, max_features="sqrt", class_weight="balanced",random_state=0, n_jobs=-1)
        sel = BorutaPy(estimator=estimator, n_estimators='auto', max_iter=500)  # number of iterations to perform
        sel.fit(train_X, train_y)
        feature_list = np.array(feature_list)[sel.support_].tolist()
        train_X = sel.transform(train_X)
        test_X = sel.transform(test_X)
        return train_X, test_X, feature_list
    elif fs_method == "MMUPHin+Boruta":
        taxon_prefix = taxon.strip().split('+')[0]
        MMUPHin_diff_dir = "/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff"
        infile = os.path.join(MMUPHin_diff_dir, taxon_prefix, "MMUPHin_diff_merge.csv")
        MMUPHin_df = pd.read_csv(infile, header=0, index_col=0)
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        train_X_df = pd.DataFrame(train_X, columns=feature_list)
        train_y_df = pd.DataFrame(train_y, columns=['group'])
        test_X_df = pd.DataFrame(test_X, columns=feature_list)
        feature_list = MMUPHin_df.loc[MMUPHin_df.index.isin(train_X_df.columns.values.tolist()),].index.values.tolist()
        print("feature_list from MMUPHin: {}".format(len(feature_list)))
        train_X = train_X_df.loc[:, feature_list].values
        test_X = test_X_df.loc[:, feature_list].values
        ###<-----------------------去除缺失值-------------------------------->###
        train_X = np.nan_to_num(train_X, nan=0.0)
        test_X = np.nan_to_num(test_X, nan=0.0)
        ###<---------------Boruta-------------->###
        estimator = RandomForestClassifier(n_estimators=500, max_features="sqrt", class_weight="balanced",random_state=0, n_jobs=-1)
        sel = BorutaPy(estimator=estimator, n_estimators='auto', max_iter=500)  # number of iterations to perform
        sel.fit(train_X, train_y)
        feature_list = np.array(feature_list)[sel.support_].tolist()
        train_X = sel.transform(train_X)
        test_X = sel.transform(test_X)
        return train_X, test_X, feature_list


def CVout_feature_selection(X, label, feature_list, fs_method):
    if fs_method == "all":
        X = X
        feature_list = feature_list
    elif fs_method == "Wilcox":
        choose_feature = []
        for i in range(X.shape[1]):
            if stats.ranksums(X[:, i][label == 1], X[:, i][label == 0])[1] < 0.05:
                choose_feature.append(i)
        X = X[:, choose_feature]
        feature_list = [feature_list[index] for index in choose_feature]
    elif fs_method == "t-test":
        choose_feature = []
        for i in range(X.shape[1]):
            if stats.ttest_ind(X[:, i][label == 1], X[:, i][label == 0])[1] < 0.05:
                choose_feature.append(i)
        X = X[:, choose_feature]
        feature_list = [feature_list[index] for index in choose_feature]
    elif fs_method == "mRMR":
        X_df = pd.DataFrame(X, columns=feature_list)
        label_df = pd.DataFrame(label, columns=['group'])
        merge_df = pd.merge(label_df, X_df, left_index=True, right_index=True)
        choose_feature = pymrmr.mRMR(merge_df, 'MIQ', round(sqrt(X_df.shape[1])))
        X = X_df.loc[:, choose_feature].values
        feature_list = choose_feature
    elif fs_method == "RandomForest":
        X_df = pd.DataFrame(X, columns=feature_list)
        label_df = pd.DataFrame(label, columns=['group'])
        sel = SelectFromModel(RandomForestClassifier(n_estimators=500), max_features=round(sqrt(X_df.shape[1])))
        sel.fit(X_df, label_df)
        feature_list = X_df.columns[(sel.get_support())]
        X = X_df.loc[:, feature_list].values
    elif fs_method == "RFECV+RandomForest":
        X_df = pd.DataFrame(X, columns=feature_list)
        label_df = pd.DataFrame(label, columns=['group'])
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, max_features="sqrt", n_jobs=-1)
        sel = SelectFromModel(estimator=estimator, max_features=round(sqrt(X_df.shape[1])))
        sel.fit(X_df, label_df)
        feature_list = X_df.columns[(sel.get_support())]
        X = X_df.loc[:, feature_list].values
        ###
        estimator = RandomForestClassifier(random_state=0, n_estimators=500, max_features="sqrt", n_jobs=-1)
        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)
        rfecv = RFECV(estimator=estimator, step=0.1, cv=cv, scoring='accuracy', n_jobs=-1)
        rfecv = rfecv.fit(X, label)
        print("Optimal number of features : %d" % rfecv.n_features_)
        X = rfecv.transform(X)
        feature_list = np.array(feature_list)[rfecv.support_]
    return X, feature_list
