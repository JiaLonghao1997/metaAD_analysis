import pandas as pd
import numpy as np
import os
import sys
import warnings
warnings.filterwarnings('ignore')

sys.path.append("/home1/jialh/brain/01meta/multikingdom")
from scripts.utils.merge_feature import merge_features_transfer
from scripts.utils.feature_selection import transfer_feature_selection
#from scripts.utils.RFECV_transfer import RFECV_transfer
from scripts.utils.CV_transfer_after_RFECV import CV_transfer
from scripts.testwork import testwork
import multiprocessing
import scipy.stats as stats
import pymrmr
import matplotlib.pyplot as plt
import seaborn as sns
from scripts.utils.Solve_Imbalanced_Classes import imblearn
from scripts.utils.filter_samples_with_shap import filter_samples_with_shap, filter_samples_with_shap_5X
from scripts.utils.dataset_split import preprare_train_and_test
from collections import Counter

##preprare_train_and_test(workdir, meta_feats, stage, train_study, test_studies, taxon, model, feat_type, fs_method, sampler)
def run(workdir, taxon, meta_feats, stage, train_study,test_studies, feat_type, fs_method, sampler, outdir, model, repeats, kfold):
    outfile = os.path.join(outdir, '01results_{}_{}_{}_{}_{}_None.csv'.format(train_study, stage, model, feat_type, fs_method))
    #result_df = pd.DataFrame
    if os.path.exists(outfile):
        print("{} exists.".format(outfile))
        result_df = pd.read_csv(outfile, header=0, index_col=None)
    else:
        train_X, train_sample_list, train_y, test_study_dic, feature_list, n_features_in, n_features, selected_features = preprare_train_and_test(
            workdir, taxon, meta_feats, stage, train_study, test_studies,model, feat_type, fs_method, sampler)
        if train_X is not None:
            ##至少包括1个特征，3个样本，才进行分类分析。
            if train_X.shape[1]>=1 and train_X.shape[0]>=3:
                result_df = CV_transfer(train_X, train_y, test_study_dic, train_sample_list, feature_list, fs_method, sampler, outdir, stage,
                                    train_study, test_studies, feat_type, n_features_in, n_features, selected_features, model, repeats, kfold)
            else:
                result_df = None
        else:
            result_df = None
    # print("{}: {}".format(outfile, result_df))
    return result_df

def main():
    import warnings
    warnings.filterwarnings("ignore")
    workdir = "/home1/jialh/brain/01meta/multikingdom/06classification_20231115"
    inputdir = "/home1/jialh/brain/01meta/multikingdom/00profile"
    repeats = 20
    kfold = 5
    #model = sys.argv[1]
    stages = ['SCD', 'SCS', 'MCI', 'AD']
    # stages = ['AD', 'MCI']
    models = ['RandomForest']
    #feat_types = ['ABFV']
    train_studies = ['CHN+CHN2']
    # fs_methods = ['all+Boruta', 'MMUPHin+Boruta']
    # fs_methods = ['MMUPHin+RFECV','RandomForest+RFECV','RandomForest+RFECV_X2','RandomForest+RFECV_X4', 'all+Boruta', 'MMUPHin+Boruta']
    # fs_methods = ['RandomForest+RFECV','RandomForest+RFECV_X2','RandomForest+RFECV_X4']
    fs_methods = ['RandomForest+RFECV_X4']
    # fs_methods = ['all+RFECV', 'all+Boruta','MMUPHin+RFECV', 'MMUPHin+Boruta']
    samplers = ['None']
    #taxons = ['species', 'KOs', 'pathways', 'GBMs', 'GMMs']
    taxons = ['species+ConQuR', 'KOs+ConQuR', 'pathways+ConQuR', 'GMMs+ConQuR', 'GBMs+ConQuR']
    # taxons = ['KOs+ConQuR']
    # taxons = ['species+ConQuR', 'KOs+ConQuR', 'pathways+ConQuR', 'GMMs+ConQuR', 'GBMs+ConQuR']
    # taxons = ['species+ConQuR']
    outfile = os.path.join(workdir, '01merge_result_transfer.csv')
    result_list = []
    pool = multiprocessing.Pool(processes=16)
    ####********************************先选特征，然后再交叉验证******************************。
    for taxon in taxons:
        if "species" in taxon:
            feat_types = ['A', 'B', 'F', 'V', 'AB', 'AF', 'AV', 'BF', 'BV', 'FV',  'ABF', 'ABV', 'AFV', 'BFV', 'ABFV']
            # feat_types = ['ABFV']
        elif taxon == 'species+ConQuR':
            #feat_types = ['A', 'B', 'F', 'V', 'AB', 'AF', 'AV', 'BF', 'BV', 'FV', 'ABF', 'ABV', 'AFV', 'BFV', 'ABFV']
            feat_types = ['ABFV']
        elif taxon in ['pathways', 'KOs', 'pathways+ConQuR', 'KOs+ConQuR', 'GMMs+ConQuR', 'GBMs+ConQuR']:
            feat_types = [taxon.strip().split('+')[0]]
        elif taxon =='species+KOs+pathways':
            feat_types = ['ABFV+KOs+pathways']
        elif taxon == 'species+KOs+pathways+ConQuR':
            feat_types = ['ABFV+KOs+pathways+ConQuR']
        meta_feats = pd.read_csv(os.path.join(inputdir, "metadata_{}.csv".format(taxon)), header=0, index_col=0)
        #meta_feats = meta_feats.loc[~meta_feats.index.isin(['B1009', 'B1521']), ]

        #workdir = os.path.join(workdir, "06transfer_{}".format(taxon))

        print(meta_feats.columns.values.tolist()[0:15])
        print(meta_feats.groupby(['Batch', 'Group'])['Group'].count())
        ###<----------------------------->###
        ##https://www.analyticsvidhya.com/blog/2020/07/10-techniques-to-deal-with-class-imbalance-in-machine-learning/
        ###<--------------------------------------------------------###
        for model in models:
            for train_study in train_studies:
                # result_list = []
                # pool = multiprocessing.Pool(processes=16)
                test_studies = []
                if train_study == 'CHN':
                    test_studies = ['CHN2', 'GER', 'JPN']
                elif train_study == 'CHN2':
                    test_studies = ['CHN', 'GER', 'JPN']
                elif train_study == "CHN+CHN2":
                    test_studies = ['JPN', 'GER']
                elif train_study == 'CHN+CHN2+GER':
                    test_studies = ['JPN']
                outdir = os.path.join(workdir, "{}_{}".format(train_study, taxon), model)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
                for sampler in samplers:
                    for stage in stages:
                        for feat_type in feat_types:
                            for fs_method in fs_methods:
                                print("begin to deal with NC-{} by {}_{}".format(stage, model, train_study))
                                # result = run(workdir, taxon, meta_feats, stage,
                                # train_study, test_studies, feat_type, fs_method, sampler, outdir, model, repeats, kfold)
                                # print("result: {}".format(result))
                                result_list.append(pool.apply_async(func=run, args=(workdir, taxon, meta_feats, stage,
                                    train_study, test_studies, feat_type, fs_method, sampler, outdir, model, repeats, kfold)))
    pool.close()
    pool.join()
    merge_result_df = pd.DataFrame()
    for result in result_list:
        result_df = result.get()
        if merge_result_df is None:
            merge_result_df = result_df
        elif merge_result_df.shape[0]>1 and result is not None:
            merge_result_df = pd.concat([merge_result_df, result_df], axis=0)
    merge_result_df.to_csv(outfile, index=False)

if __name__ == '__main__':
    main()