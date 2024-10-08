import numpy as np
import pandas as pd
import os

##要绘制热图，首先需要通过/home1/jialh/brain/01meta/multikingdom/scripts/img/02performance_stat.py整理结果。
#os.system("python /home1/jialh/brain/01meta/multikingdom/scripts/img/02performance_stat.py")
workdir = "/home1/jialh/brain/01meta/multikingdom/06classification_20231115"
performance_stat_file = os.path.join(workdir, "performance_stat.csv")
data = pd.read_csv(performance_stat_file, header=0)
#data = data.loc[data['']]
###
test_studies = ['CHN+CHN2', 'JPN', 'GER']
fs_methods = ['MMUPHin+RFECV', 'MMUPHin_stage+RFECV', 'RandomForest+RFECV','RandomForest+RFECV_X2','RandomForest+RFECV_X4',
              'all+Boruta', "all+RFECV", 'MMUPHin+Boruta']

for fs_method in fs_methods:
    feats = ['A', 'B', 'F', 'V', 'AB', 'AF', 'AV', 'BF', 'BV', 'FV', 'ABF', 'ABV', 'AFV', 'BFV', 'ABFV', "KOs", 'GMMs',
             'GBMs', 'pathways']
    # stages = ['SCS', 'SCD', 'MCI', 'AD']
    stages = ['CHN+CHN2_SCS', 'CHN+CHN2_SCD', 'CHN+CHN2_MCI', 'CHN+CHN2_AD', 'JPN_MCI', 'JPN_AD', 'GER_AD']
    outdir = os.path.join(workdir, "heatmaps")
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    AUC_df = pd.DataFrame(np.zeros(shape=(len(feats), len(stages))), index=feats, columns=stages)
    feats_df = pd.DataFrame(np.zeros(shape=(len(feats), len(stages))), index=feats, columns=stages)

    for test_study in test_studies:
        data_filter = data.loc[(data['test_study']==test_study) & (data['fs_method']==fs_method), ]
        #c('A. Archaea', 'B. Bacteria', 'C. Fungi', 'D. Viruses'))
        #feat_type_dic = {'A':'Archaea', 'B':'Bacteria', 'F':'Fungi', 'V':'Viruses'}
        #data = data.replace({'feat_type': feat_type_dic})
        for index, row in data_filter.iterrows():
            feat = row['feat_type']
            stage = row['test_study']+ '_' + row['stage']
            AUC_df.loc[feat, stage] = row['auroc_mean']
            feats_df.loc[feat, stage] = row['n_features']

    ####
    zero_num = (AUC_df==0).astype(int).sum(axis=0)
    print("There are {} 0.".format(zero_num))
    AUC_df.to_csv(os.path.join(outdir, "{}_AUC.csv".format(fs_method)), index=True, index_label='feature')
    feats_df.to_csv(os.path.join(outdir, "{}_n_features.csv".format(fs_method)), index=True, index_label='feature')