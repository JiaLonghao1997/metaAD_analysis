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
import os
import numpy as np
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')


def plot_REFCV(clf, min_features_to_select, outdir, study, model, feat_type, fs_method, stage, n_features, auroc):
    n_scores = len(clf.best_estimator_.cv_results_["mean_test_score"])
    #plt.figure()
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlabel("Number of features selected")
    ax.set_ylabel("Mean test AUROC")
    ax.errorbar(
        range(min_features_to_select, n_scores + min_features_to_select),
        clf.best_estimator_.cv_results_["mean_test_score"],
        yerr=clf.best_estimator_.cv_results_["std_test_score"],
    )
    ax.set(ylim=[-0.05, 1.05])
    ax.set_title("REF({}-{}-{}-{}-{},{} features, AUC={:.2f})".format(study, model, feat_type, fs_method, stage, n_features, auroc))
    plt.show()
    if not os.path.exists(os.path.join(outdir, study)):
        os.makedirs(os.path.join(outdir, study), exist_ok=True)
    fig.savefig(os.path.join(outdir, study, "{}_{}_{}_{}_{}_RFECV.pdf".format(stage, study, feat_type, fs_method, model)))
    fig.savefig(os.path.join(outdir, study, "{}_{}_{}_{}_{}_RFECV.png".format(stage, study, feat_type, fs_method, model)))

##pool.apply_async(func=RFECV_clf, args=(X, y, feature_list, outdir, stage, study, feat_type, model, repeats=1, fold=5))
def RFECV_clf(X, y, feature_list, fs_method, outdir, stage, study, feat_type, model, repeats, fold):
    min_features_to_select = 1  # Minimum number of features to consider
    #cv = RepeatedStratifiedKFold(n_splits=fold, n_repeats=repeats)
    #cv = RepeatedStratifiedKFold(n_splits=fold, n_repeats=repeats)
    cv = StratifiedKFold(n_splits=fold)
    estimator = ""
    param_grid = ""
    ##参数调整参考: https://www.kaggle.com/code/ldfreeman3/a-data-science-framework-to-achieve-99-accuracy/notebook#5.12-Tune-Model-with-Hyper-Parameters
    # if model == 'SVM':
    #     estimator = SVC(probability=True, random_state=0)
    #     # param_grid = {'estimator__C': [1,2,3,4,5],
    #     #               'estimator__gamma': ['auto', 'scale']}  ##SVM无法返回特征权重值，无法与RFECV结合。

    if model == 'RandomForest':
        estimator = RandomForestClassifier(random_state=0, n_jobs=-1)
        param_grid = {'estimator__n_estimators': [10, 50, 100, 300, 500],
                      'estimator__max_depth': [1, 2, 4, 6, 8, 10, None]}

    if model == 'Lasso':
        estimator = LogisticRegression(penalty='l1', solver='liblinear', max_iter=1000, n_jobs=-1)
        param_grid = {'estimator__C': [0.01, 0.1, 1, 10, 100, 1000]}

    if model == "ElasticNet":
        ##参考:https://stackoverflow.com/questions/66787845/how-to-perform-elastic-net-for-a-classification-problem
        estimator = LogisticRegression(penalty='elasticnet', solver='saga', random_state=0)
        param_grid = {'estimator__C': [0.01, 0.1, 1, 10, 100, 1000],
                      'estimator__l1_ratio': [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]}

    if model == "AdaBoost":  ##Adaptive Boosting, 自适应增强
        #Bagging(个体学习器间不存在强依赖关系、可同时生成的并行化 方法) =>随机森林。
        #Boosting(个体学习器间存在强依赖关系、必须串行生成的序列 化方法)。
        #以分类为例，Adaboost算法通过提高前一轮分类器分类错误的 样本的权值，而降低那 些被分类正确的样本的权值。
        estimator = AdaBoostClassifier(random_state=0)
        param_grid = {'estimator__n_estimators': [10, 50, 100, 300, 500]}

    if model == "GradientBoosting": ##梯度提升决策树。
        # 和Adaboost不同，Gradient Boosting 在迭代的时候选择损失 函数在其梯度方向下降 的方向不停地改进模型。
        estimator = GradientBoostingClassifier(random_state=0)
        param_grid = {'estimator__n_estimators': [10, 50, 100, 300, 500],
                      'estimator__max_depth': [2, 4, 6, 8, 10, None]}

    if model == "XGBoost":
        # eXtreme Gradient Boosting，可译为极限梯度提升算法。
        estimator = XGBClassifier(seed=0)
        param_grid = {'estimator__n_estimators': [10, 50, 100, 300, 500],
                      'estimator__max_depth': [1, 2, 4, 6, 8, 10, None]}

    selector = RFECV(
        estimator=estimator,
        step=1,
        cv=cv,
        scoring="roc_auc",
        min_features_to_select=min_features_to_select,
        n_jobs=-1,
    )

    ##sklearn.model_selection.GridSearchCV
    clf = GridSearchCV(selector, param_grid=param_grid, cv=cv, n_jobs=-1)

    clf.fit(X, y)
    ###

    n_features_in = clf.best_estimator_.n_features_in_
    n_features = clf.best_estimator_.n_features_
    selected_features = ''
    try:
        selected_features = ",".join(np.array(feature_list)[clf.best_estimator_.get_support()])
        #print("selected_features: ", selected_features)
    except:
        print("REF({}-{}-{}-{},{} features, feature_list: {})".format(model, feat_type, fs_method, stage, n_features, feature_list))
    auroc = clf.best_score_

    #selected_features = (X[:, clf.best_estimator_.get_support()])

    print(f"Optimal number of features: {n_features}")   ##最佳特征数目是多少。
    print("REF({}-{}-{}-{}-{},{} features, AUC={:.2f})".format(study, model, feat_type, fs_method, stage, n_features, auroc))
    plot_REFCV(clf, min_features_to_select, outdir, study, model, feat_type, fs_method, stage, n_features, auroc)

    #results.append([stage, study, feat_type, model,n_features_in, n_features, selected_features, auroc])
    result = [stage, study, feat_type, fs_method, model, auroc,  n_features_in, n_features, selected_features]
    return(result)
    #return(stage, study, feat_type, model,n_features_in,n_features,selected_features,auroc)
