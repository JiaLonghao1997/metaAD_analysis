from sklearn.feature_selection import RFECV
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import ElasticNet
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import GridSearchCV
from xgboost import XGBClassifier

from sklearn.feature_selection import RFE
from sklearn.metrics import roc_auc_score
import os
import numpy as np
import pandas as pd
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')


def plot_REFCV(clf, min_features_to_select, outdir, train_study, model, feat_type, fs_method, stage, n_features, auroc):
    n_scores = len(clf.cv_results_["mean_test_score"])
    #plt.figure()
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlabel("Number of features selected")
    ax.set_ylabel("Mean test AUROC")
    ax.errorbar(
        range(min_features_to_select, n_scores + min_features_to_select),
        clf.cv_results_["mean_test_score"],
        yerr=clf.cv_results_["std_test_score"],
    )
    ax.set(ylim=[-0.05, 1.05])
    ax.set_title("REF({}-{}-{}-{}-{},{} features, AUC={:.2f})".format(train_study, model, feat_type, fs_method, stage, n_features, auroc))
    plt.show()
    if not os.path.exists(os.path.join(outdir, train_study)):
        os.makedirs(os.path.join(outdir, train_study), exist_ok=True)
    fig.savefig(os.path.join(outdir, train_study, "{}_{}_{}_{}_{}_RFECV.pdf".format(stage, train_study, feat_type, fs_method, model)))
    fig.savefig(os.path.join(outdir, train_study, "{}_{}_{}_{}_{}_RFECV.png".format(stage, train_study, feat_type, fs_method, model)))

##pool.apply_async(func=RFECV_clf, args=(X, y, feature_list, outdir, stage, study, feat_type, model, repeats=1, fold=5))
##train_X, train_y, test_study_dic, feature_list, fs_method, outdir, stage, train_study, feat_type, model, repeats, fold
def RFECV_transfer(train_X, train_y, test_study_dic, feature_list, fs_method, outdir, stage, train_study, feat_type, model, repeats, fold):
    min_features_to_select = 1  # Minimum number of features to consider
    #cv = RepeatedStratifiedKFold(n_splits=fold, n_repeats=repeats)
    #cv = RepeatedStratifiedKFold(n_splits=fold, n_repeats=repeats)
    cv = RepeatedStratifiedKFold(n_splits=fold, n_repeats=repeats, random_state=0)
    estimator = ""
    param_grid = ""
    ##参数调整参考: https://www.kaggle.com/code/ldfreeman3/a-data-science-framework-to-achieve-99-accuracy/notebook#5.12-Tune-Model-with-Hyper-Parameters
    # if model == 'SVM':
    #     estimator = SVC(probability=True, random_state=0)
    #     param_grid = {'estimator__C': [1,2,3,4,5],
    #                   'estimator__gamma': ['auto', 'scale']}  ##SVM无法返回特征权重值，无法与RFECV结合。

    ##参数调整参考: https://www.kaggle.com/code/ldfreeman3/a-data-science-framework-to-achieve-99-accuracy/notebook#5.12-Tune-Model-with-Hyper-Parameters
    if model == 'SVM':
        estimator = SVC(probability=True, random_state=0)
        # param_grid = {'estimator__C': [1,2,3,4,5],
        #               'estimator__gamma': ['auto', 'scale']}  ##SVM无法返回特征权重值，无法与RFECV结合。

    if model == 'RandomForest':
        estimator = RandomForestClassifier(random_state=0, n_jobs=-1)
        # param_grid = {'estimator__n_estimators': [10, 50, 100, 300, 500],
        #               'estimator__max_depth': [1, 2, 4, 6, 8, 10, None]}

    if model == 'Lasso':
        estimator = LogisticRegression(penalty='l1', solver='liblinear', max_iter=1000, n_jobs=-1)
        # param_grid = {'estimator__C': [0.01, 0.1, 1, 10, 100, 1000]}

    if model == "ElasticNet":
        ##参考:https://stackoverflow.com/questions/66787845/how-to-perform-elastic-net-for-a-classification-problem
        estimator = LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5, random_state=0)
        # param_grid = {'estimator__C': [0.01, 0.1, 1, 10, 100, 1000],
        #               'estimator__l1_ratio': [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]}

    if model == "AdaBoost":  ##Adaptive Boosting, 自适应增强
        #Bagging(个体学习器间不存在强依赖关系、可同时生成的并行化 方法) =>随机森林。
        #Boosting(个体学习器间存在强依赖关系、必须串行生成的序列 化方法)。
        #以分类为例，Adaboost算法通过提高前一轮分类器分类错误的 样本的权值，而降低那 些被分类正确的样本的权值。
        estimator = AdaBoostClassifier(random_state=0)
        # param_grid = {'estimator__n_estimators': [10, 50, 100, 300, 500]}

    if model == "GradientBoosting": ##梯度提升决策树。
        # 和Adaboost不同，Gradient Boosting 在迭代的时候选择损失 函数在其梯度方向下降 的方向不停地改进模型。
        estimator = GradientBoostingClassifier(random_state=0)
        # param_grid = {'estimator__n_estimators': [10, 50, 100, 300, 500],
        #               'estimator__max_depth': [2, 4, 6, 8, 10, None]}

    if model == "XGBoost":
        # eXtreme Gradient Boosting，可译为极限梯度提升算法。
        estimator = XGBClassifier(n_jobs=8, seed=0)
        # param_grid = {'estimator__n_estimators': [10, 50, 100, 300, 500],
        #               'estimator__max_depth': [1, 2, 4, 6, 8, 10, None]}

    # selector = RFECV(
    #     estimator=estimator,
    #     step=1,
    #     cv=cv,
    #     scoring="roc_auc",
    #     min_features_to_select=min_features_to_select,
    #     n_jobs=-1,
    # )

    ##sklearn.model_selection.GridSearchCV
    #clf = GridSearchCV(selector, param_grid=param_grid, cv=cv, n_jobs=-1)

    clf = RFECV(
        estimator=estimator,
        step=1,
        cv=cv,
        scoring="roc_auc",
        min_features_to_select=min_features_to_select,
        n_jobs=8,
    )

    clf.fit(train_X, train_y)
    ###

    n_features_in = clf.n_features_in_
    n_features = clf.n_features_
    selected_features = ''
    try:
        selected_features = ",".join(np.array(feature_list)[clf.support_])
        #print("selected_features: ", selected_features)
    except:
        print("REF({}-{}-{}-{},{} features, feature_list: {})".format(model, feat_type, fs_method, stage, n_features, feature_list))
    cv_auroc_mean = clf.cv_results_['mean_test_score'][n_features-1]
    cv_auroc_std = clf.cv_results_['std_test_score'][n_features-1]

    print(f"Optimal number of features: {n_features}")  ##最佳特征数目是多少。
    print("REF({}-{}-{}-{}-{},{} features, AUC={:.2f}+/-{:.2f})".format(train_study, model, feat_type, fs_method, stage, n_features,
                                                               cv_auroc_mean, cv_auroc_std))
    plot_REFCV(clf, min_features_to_select, outdir, train_study, model, feat_type, fs_method, stage, n_features, cv_auroc_mean)

    result_list = []
    #selected_features = (X[:, clf.get_support()])
    for key, value in test_study_dic.items():
        test_study = key
        test_X = value['test_X']
        test_y = value['test_y']
        test_y_proba = clf.predict_proba(test_X)[:, 1]
        test_auroc = roc_auc_score(test_y, test_y_proba)
        print("REF(train_{}-test_{}-{}-{}-{}-{},{} features,transfer AUROC={:.2f}".format(train_study, test_study, model, feat_type, fs_method, stage, n_features, test_auroc))
        #results.append([stage, study, feat_type, model,n_features_in, n_features, selected_features, auroc])
        result = [stage, train_study, test_study, feat_type, fs_method, model,  cv_auroc_mean, cv_auroc_std, test_auroc, n_features_in, n_features, selected_features]
        result_list.append(result)
    result_df = pd.DataFrame(result_list,
                             columns=['stage', 'train_study', 'test_study', 'feat_type',
                                      'fs_method', 'model', 'cv_auroc_mean', 'cv_auroc_std', 'test_auroc',
                                      'n_features_in', 'n_features', 'selected_features'])
    outfile = os.path.join(outdir, '01RFECV_{}_{}_{}_{}_20230320_transfer.csv'.format(model, train_study, feat_type, fs_method))
    result_df.to_csv(outfile, index=False)
    return(result_df)
    #return(stage, study, feat_type, model,n_features_in,n_features,selected_features,auroc)
