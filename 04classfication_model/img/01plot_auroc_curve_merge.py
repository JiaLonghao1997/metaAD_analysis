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

##plot_auroc_curve(train_study, stage, model, fs_method, auroc_curve_dic, outdir)
def plot_auroc_curve(train_study, stage, model, fs_method, auroc_curve_dic, outdir):
    #mean_fpr = np.linspace(0, 1, 100)
    fig, ax = plt.subplots(figsize=(5.5, 5.5))

    ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="grey", label="Chance", alpha=0.8)
    # color_dic = {'Species': '#03ACDA', 'Species+ConQuR': '#03ACDA',
    #              'KOs': '#6AB42D', 'KOs+ConQuR': '#6AB42D',
    #              'Pathways': '#ff993f', 'Pathways+ConQuR': '#ff993f',
    #              'All':'#E03837', 'All+ConQuR':'#E03837'}

    # color_dic = {'Species': '#E03837',
    #              'KOs': '#ff993f',
    #              'GMMs': '#6AB42D',
    #              'GBMs':'#03ACDA',
    #              'Pathways': '#a374aa'}
    # color_dic = {'Species': '#E03837',
    #              'KOs': '#ff993f',
    #              'GMMs': '#6AB42D',
    #              'GBMs':'#03ACDA',
    #              'Pathways': '#a374aa'}
    color_dic = {'Species': '#E03837',
                 'KOs': '#ff993f',
                 'Pathways': '#03ACDA'}

    ##<---------------------如何保存绘制AUROC曲线所需数据。----------------------------->##
    ##auroc_mean, auroc_std, mean_tpr, mean_fpr
    for key, value in auroc_curve_dic.items():
        taxon = key
        color = color_dic[taxon]
        auroc_mean = value[0]
        auroc_std = value[1]
        mean_tpr = value[2]
        mean_fpr = value[3]
        ax.plot(
            mean_fpr,
            mean_tpr,
            color=color,
            label=r"%s (AUC = %0.2f)" % (taxon, auroc_mean),
            lw=2,
            alpha=0.8,
        )

    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05]
    )
    ax.legend(loc="lower right")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('False positive rate', fontsize=12)
    plt.ylabel('True positive rate', fontsize=12)
    plt.legend(fontsize=10, loc="lower right")
    plt.title(label="NC vs {} (transfer to JPN)".format(stage), fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "auroc_{}_{}_{}_merge.png".format(
        train_study, stage, model)), dpi=300)
    plt.savefig(os.path.join(outdir,"auroc_{}_{}_{}_merge.pdf".format(
        train_study, stage, model)), dpi=300)

####
workdir = "/home1/jialh/brain/01meta/multikingdom/06classification_20231115"
train_studies = ['CHN+CHN2']
test_studies = ['JPN']
#feat_types = ['species+ConQuR', 'KOs+ConQuR', 'species+KOs+ConQuR']
suffix = '+ConQuR'
#suffix = ''
taxons = ['species'+suffix, 'KOs'+suffix,  'pathways'+suffix]
#taxons = ['species'+suffix, 'KOs'+suffix]
#taxons = ['species+KOs+pathways']

fs_methods = ['RandomForest+RFECV', 'RandomForest+RFECV_X4']
models = ['RandomForest']
stages = ['AD', 'MCI', 'SCD', 'SCS']

# aucs = []
mean_fpr = np.linspace(0, 1, 100)
# outdir = os.path.join(workdir, "auroc_curves"+suffix)
# if not os.path.exists(outdir):
#     os.makedirs(outdir)

for train_study in train_studies:
    for test_study in test_studies:
        outdir = os.path.join(workdir, "auroc_curves_{}_{}".format(test_study, suffix))
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        for model in models:
            for stage in stages:
                # for fs_method in fs_methods:
                result_list = []
                auroc_curve_dic = {}
                for taxon in taxons:
                    if taxon == 'species+ConQuR':
                        feat_types = ['ABFV']
                        fs_method = 'RandomForest+RFECV_X4'
                    elif taxon == "KOs+ConQuR":
                        feat_types = [taxon.strip().split('+')[0]]
                        fs_method = 'RandomForest+RFECV'
                    elif taxon == "pathways+ConQuR":
                        feat_types = [taxon.strip().split('+')[0]]
                        fs_method = 'RandomForest+RFECV_X4'
                    for feat_type in feat_types:
                        inputdir = os.path.join(workdir, "{}_{}".format(train_study, taxon), model)
                        if test_study == train_study:
                            cv_prefix = "CV_{}_{}_{}_{}_{}_None".format(train_study, stage, model, feat_type, fs_method)
                            infile = os.path.join(inputdir, "{}_tpr.csv".format(cv_prefix))
                        else:
                            ##transfer_fromCHN+CHN2toGER_AD_RandomForest_ABFV_MMUPHin+RFECV_None_GER
                            transfer_prefix = "transfer_from{}to{}_{}_{}_{}_{}_None_{}".format(train_study, test_study, stage, model, feat_type, fs_method, test_study)
                            infile = os.path.join(inputdir, "{}_tpr.csv".format(transfer_prefix))
                        aucfile = os.path.join(inputdir, "01results_{}_{}_{}_{}_{}_None.csv".format(
                            train_study, stage, model, feat_type, fs_method))
                        print("infile: {}".format(infile))
                        print("aucfile: {}".format(aucfile))
                        if (os.path.exists(infile) and os.path.exists(aucfile)):
                            print("{} exists!".format(infile))
                            auc_df = pd.read_csv(aucfile, header=0)
                            auroc_mean = auc_df.loc[auc_df['test_study']==test_study, 'auroc_mean']
                            auroc_std = auc_df.loc[auc_df['test_study']==test_study, 'auroc_std']
                            tpr_list = np.transpose(np.loadtxt(infile))
                            mean_tpr = np.mean(tpr_list, axis=0)  ##np.mean(axis=0),压缩的是行。
                            # mean_tpr[0] = 0
                            mean_tpr[-1] = 1.0
                            if taxon in ["species+KOs+pathways+ConQuR", "species+KOs+pathways"]:
                                taxon = "All"
                            elif taxon in ["species+ConQuR", "species"]:
                                taxon = "Species"
                            elif taxon in ['pathways', 'pathways+ConQuR']:
                                taxon = "Pathways"
                            else:
                                taxon = taxon.rstrip('+ConQuR')
                            auroc_curve_dic[taxon] = [auroc_mean, auroc_std, mean_tpr, mean_fpr]
                    if len(auroc_curve_dic) > 0:
                        print(len(auroc_curve_dic))
                        plot_auroc_curve(train_study, stage, model, fs_method, auroc_curve_dic, outdir)
