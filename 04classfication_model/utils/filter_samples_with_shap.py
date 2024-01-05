import os
import pandas as pd

def filter_samples_with_shap(meta_feats, workdir, inputdir, taxon, studies, stages):
    #inputdir = "/home1/jialh/brain/01meta/multikingdom/06SHAP/AD_species0530"
    ##<-----------------AD vs NC--------------->##
    CHN2_AD_SHAP = pd.read_csv(
        os.path.join(inputdir, "JPN+GER+CHN", "kmeans", "JPN+GER+CHN_to_CHN2_ADvsNC_ABFV_AD_cluster.csv"), header=0,
        index_col=0)
    CHN2_AD_exclude = CHN2_AD_SHAP.loc[CHN2_AD_SHAP['cluster'].isin([2]),].index.values.tolist()
    CHN_AD_SHAP = pd.read_csv(
        os.path.join(inputdir, "JPN+GER+CHN2", "kmeans", "JPN+GER+CHN2_to_CHN_ADvsNC_ABFV_AD_cluster.csv"), header=0,
        index_col=0)
    CHN_AD_exclude = CHN_AD_SHAP.loc[CHN_AD_SHAP['cluster'].isin([1]),].index.values.tolist()

    ##<-----------------MCI vs NC--------------->##
    CHN2_MCI_SHAP = pd.read_csv(
        os.path.join(inputdir, "JPN+GER+CHN", "kmeans", "JPN+GER+CHN_to_CHN2_MCIvsNC_ABFV_MCI_cluster.csv"),
        header=0, index_col=0)
    CHN2_MCI_exclude = CHN2_MCI_SHAP.loc[CHN2_MCI_SHAP['cluster'].isin([4]),].index.values.tolist()
    CHN_MCI_SHAP = pd.read_csv(
        os.path.join(inputdir, "JPN+GER+CHN2", "kmeans", "JPN+GER+CHN2_to_CHN_MCIvsNC_ABFV_MCI_cluster.csv"),
        header=0, index_col=0)
    CHN_MCI_exclude = CHN_MCI_SHAP.loc[CHN_MCI_SHAP['cluster'].isin([3]),].index.values.tolist()

    ##<-----------------SCD vs NC--------------->##
    CHN2_SCD_SHAP = pd.read_csv(
        os.path.join(inputdir, "JPN+GER+CHN", "kmeans", "JPN+GER+CHN_to_CHN2_SCDvsNC_ABFV_SCD_cluster.csv"),
        header=0, index_col=0)
    CHN2_SCD_exclude = CHN2_SCD_SHAP.loc[CHN2_SCD_SHAP['cluster'].isin([3]),].index.values.tolist()
    CHN_SCD_SHAP = pd.read_csv(
        os.path.join(inputdir, "JPN+GER+CHN2", "kmeans", "JPN+GER+CHN2_to_CHN_SCDvsNC_ABFV_SCD_cluster.csv"),
        header=0, index_col=0)
    CHN_SCD_exclude = CHN_SCD_SHAP.loc[CHN_SCD_SHAP['cluster'].isin([3]),].index.values.tolist()

    ##<-----------------SCS vs NC--------------->##
    CHN2_SCS_SHAP = pd.read_csv(
        os.path.join(inputdir, "JPN+GER+CHN", "kmeans", "JPN+GER+CHN_to_CHN2_SCSvsNC_ABFV_SCS_cluster.csv"),
        header=0, index_col=0)
    CHN2_SCS_exclude = CHN2_SCS_SHAP.loc[CHN2_SCS_SHAP['cluster'].isin([3]),].index.values.tolist()
    CHN_SCS_SHAP = pd.read_csv(
        os.path.join(inputdir, "JPN+GER+CHN2", "kmeans", "JPN+GER+CHN2_to_CHN_SCSvsNC_ABFV_SCS_cluster.csv"),
        header=0, index_col=0)
    CHN_SCS_exclude = CHN_SCS_SHAP.loc[CHN_SCS_SHAP['cluster'].isin([4]),].index.values.tolist()

    with open(os.path.join(inputdir, "exclude_samples_stat.csv"), "w") as output, open(os.path.join(inputdir, "exclude_samples.csv"), "w") as exclude_samples:
        for study in studies:
            for stage in stages:
                count = len(eval('{}_{}_exclude'.format(study, stage)))
                samples = eval('{}_{}_exclude'.format(study, stage))
                print("exclude some samples: {}_{}: {}".format(study, stage, count))
                output.write("exclude some samples: {}_{}: {} \n".format(study, stage, count))
                for sample in samples:
                    exclude_samples.write("{},{},{}\n".format(study, stage, sample))
                #merge_exclude = merge_exclude + samples
    merge_exclude = CHN_AD_exclude + CHN2_AD_exclude + CHN_MCI_exclude + CHN2_MCI_exclude + CHN_SCD_exclude + CHN2_SCD_exclude + CHN_SCS_exclude + CHN2_SCS_exclude
    print("merge_exclude: {}".format(len(merge_exclude)))
    meta_feats2 = meta_feats.copy()
    meta_feats2['SHAP'] = 'include'
    meta_feats2.loc[meta_feats2.index.isin(merge_exclude), 'SHAP'] = 'exclude'
    meta_feats2.to_csv(os.path.join(workdir, "00profile", "metadata_{}.shap.csv".format(taxon)), index=True,
                       index_label='Sample')
    return merge_exclude

def filter_samples_with_shap_5X(meta_feats, workdir, inputdir, taxon, studies, stages):
    #inputdir = "/home1/jialh/brain/01meta/multikingdom/06SHAP/AD_species0530"
    ##<-----------------AD vs NC--------------->##
    CHN2_AD_SHAP = pd.read_csv(
        os.path.join(inputdir, "JPN+GER+CHN", "kmeans", "JPN+GER+CHN_to_CHN2_ADvsNC_ABFV_AD_cluster.csv"), header=0,
        index_col=0)
    CHN2_AD_exclude = CHN2_AD_SHAP.loc[CHN2_AD_SHAP['cluster'].isin([3]),].index.values.tolist()
    CHN_AD_SHAP = pd.read_csv(
        os.path.join(inputdir, "JPN+GER+CHN2", "kmeans", "JPN+GER+CHN2_to_CHN_ADvsNC_ABFV_AD_cluster.csv"), header=0,
        index_col=0)
    CHN_AD_exclude = CHN_AD_SHAP.loc[CHN_AD_SHAP['cluster'].isin([2]),].index.values.tolist()

    ##<-----------------MCI vs NC--------------->##
    CHN2_MCI_SHAP = pd.read_csv(
        os.path.join(inputdir, "JPN+GER+CHN", "kmeans", "JPN+GER+CHN_to_CHN2_MCIvsNC_ABFV_MCI_cluster.csv"),
        header=0, index_col=0)
    CHN2_MCI_exclude = CHN2_MCI_SHAP.loc[CHN2_MCI_SHAP['cluster'].isin([1]),].index.values.tolist()
    CHN_MCI_SHAP = pd.read_csv(
        os.path.join(inputdir, "JPN+GER+CHN2", "kmeans", "JPN+GER+CHN2_to_CHN_MCIvsNC_ABFV_MCI_cluster.csv"),
        header=0, index_col=0)
    CHN_MCI_exclude = CHN_MCI_SHAP.loc[CHN_MCI_SHAP['cluster'].isin([1]),].index.values.tolist()

    ##<-----------------SCD vs NC--------------->##
    CHN2_SCD_SHAP = pd.read_csv(
        os.path.join(inputdir, "JPN+GER+CHN", "kmeans", "JPN+GER+CHN_to_CHN2_SCDvsNC_ABFV_SCD_cluster.csv"),
        header=0, index_col=0)
    CHN2_SCD_exclude = CHN2_SCD_SHAP.loc[CHN2_SCD_SHAP['cluster'].isin([2]),].index.values.tolist()
    CHN_SCD_SHAP = pd.read_csv(
        os.path.join(inputdir, "JPN+GER+CHN2", "kmeans", "JPN+GER+CHN2_to_CHN_SCDvsNC_ABFV_SCD_cluster.csv"),
        header=0, index_col=0)
    CHN_SCD_exclude = CHN_SCD_SHAP.loc[CHN_SCD_SHAP['cluster'].isin([4]),].index.values.tolist()

    ##<-----------------SCS vs NC--------------->##
    CHN2_SCS_SHAP = pd.read_csv(
        os.path.join(inputdir, "JPN+GER+CHN", "kmeans", "JPN+GER+CHN_to_CHN2_SCSvsNC_ABFV_SCS_cluster.csv"),
        header=0, index_col=0)
    CHN2_SCS_exclude = CHN2_SCS_SHAP.loc[CHN2_SCS_SHAP['cluster'].isin([1]),].index.values.tolist()
    CHN_SCS_SHAP = pd.read_csv(
        os.path.join(inputdir, "JPN+GER+CHN2", "kmeans", "JPN+GER+CHN2_to_CHN_SCSvsNC_ABFV_SCS_cluster.csv"),
        header=0, index_col=0)
    CHN_SCS_exclude = CHN_SCS_SHAP.loc[CHN_SCS_SHAP['cluster'].isin([2]),].index.values.tolist()

    with open(os.path.join(inputdir, "exclude_samples_stat.csv"), "w") as output, open(os.path.join(inputdir, "exclude_samples.csv"), "w") as exclude_samples:
        for study in studies:
            for stage in stages:
                count = len(eval('{}_{}_exclude'.format(study, stage)))
                samples = eval('{}_{}_exclude'.format(study, stage))
                print("exclude some samples: {}_{}: {}".format(study, stage, count))
                output.write("exclude some samples: {}_{}: {} \n".format(study, stage, count))
                for sample in samples:
                    exclude_samples.write("{},{},{}\n".format(study, stage, sample))
                #merge_exclude = merge_exclude + samples
    merge_exclude = CHN_AD_exclude + CHN2_AD_exclude + CHN_MCI_exclude + CHN2_MCI_exclude + CHN_SCD_exclude + CHN2_SCD_exclude + CHN_SCS_exclude + CHN2_SCS_exclude
    print("merge_exclude: {}".format(len(merge_exclude)))
    meta_feats2 = meta_feats.copy()
    meta_feats2['SHAP'] = 'include'
    meta_feats2.loc[meta_feats2.index.isin(merge_exclude), 'SHAP'] = 'exclude'
    meta_feats2.to_csv(os.path.join(workdir, "00profile", "metadata_{}.shap.csv".format(taxon)), index=True,
                       index_label='Sample')
    return merge_exclude