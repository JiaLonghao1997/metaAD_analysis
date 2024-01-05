import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys
#sys.path.append("")
#from scripts.utils.merge_feature import merge_features
import seaborn as sns
import os

workdir="/home1/jialh/brain/01meta/multikingdom"
stages = ['AD', 'MCI', 'SCD', 'SCS']
models = ['RandomForest']
studies = ['CHN+CHN2']
fs_method = "RandomForest+RFECV"
# feat_type = "ABFV+KOs+pathways+ConQuR"
# taxon = "species+KOs+pathways+ConQuR"

feat_type = "KOs"
taxon = "KOs+ConQuR"
repeats = 20
kfold = 5

infile = os.path.join(workdir, "00profile", "metadata_KOs+ConQuR.csv")
meta_feats = pd.read_csv(infile, header=0, index_col=0)
outdir = os.path.join(workdir, "06classification", "feature_importance_boxplots")
if not os.path.exists(outdir):
    os.makedirs(outdir)

###
kegg_pathways = pd.read_table("/home1/jialh/brain/CAGs/KEGGdb/kegg_pathways.tsv", header=0, index_col=0)
#print(kegg_pathways)
kegg_pathways['Pathway_short_name'] = kegg_pathways['Pathway_name'].str.split(" [", n=1, expand=True, regex=False)[0]
#print(kegg_pathways)
kegg_pathways['pathway_desc'] = kegg_pathways['Pathway_short_name'].astype(str) + ' (' + kegg_pathways.index.values + ')'
###***********************只考虑代谢相关的通路******************************。
kegg_pathways_metabolism = kegg_pathways.loc[kegg_pathways['Pathway_category1']=="09100 Metabolism", ]
kegg_pathways_metabolism.to_csv(os.path.join(outdir, "kegg_pathways_metabolism.tsv"), sep="\t", header=True, index=True, index_label="pathways")


kegg_kos = pd.read_table("/home1/jialh/brain/CAGs/KEGGdb/kegg_kos.tsv", header=0, index_col=0)
kegg_kos['KO_short_name'] = kegg_kos['KO_name'].str.split(";", expand=True)[0].str.split(" ", n=1, expand=True)[1].str.strip()
kos_to_pathways = pd.read_table("/home1/jialh/brain/CAGs/KEGGdb/kegg_ko2pathways.tsv", header=None, index_col=0)
kos_to_pathways.columns = ['map']

merge_kos_to_pathways = pd.merge(kegg_kos, kos_to_pathways, left_index=True, right_index=True)
merge_kos_to_pathways_desc = pd.merge(merge_kos_to_pathways, kegg_pathways_metabolism, left_on="map", right_index=True)
merge_kos_to_pathways_desc.to_csv(os.path.join(outdir, "merge_kos_to_pathways_desc.tsv"), sep="\t",  header=True, index=True, index_label="KOs")
merge_kos_to_pathways_desc['KO_desc'] = merge_kos_to_pathways_desc.index.values + ' (' + merge_kos_to_pathways_desc['Pathway_short_name'] + ')'
merge_kos_to_pathways_desc.to_csv(os.path.join(outdir, "merge_kos_to_pathways_desc02.tsv"), sep="\t",  header=True, index=True, index_label="KOs")
##
##https://pandas.pydata.org/docs/reference/api/pandas.Index.duplicated.html#pandas.Index.duplicated
##pandas.Index.duplicated(keep='first')的意思是，在所有重复的index中, 除了第一个重复项为False外，其他重复项均为True，然后在取反，就是保留第一项。
merge_kos_to_pathways_desc_filter = merge_kos_to_pathways_desc[~merge_kos_to_pathways_desc.index.duplicated(keep='first')]
merge_kos_to_pathways_desc_filter['Pathway_category2_short'] = merge_kos_to_pathways_desc_filter['Pathway_category2'].str.split(pat=' ', n=1).str[1]
merge_kos_to_pathways_desc_filter.to_csv(os.path.join(workdir, "00profile", "KEGG_KOs_to_pathways_metabolism.csv"))
print("merge_kos_to_pathways_desc_filter.columns: {}".format(merge_kos_to_pathways_desc_filter.columns.values.tolist()))
#sys.exit()


# module_dir = "/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff/module_desc"
# for module_type in ['GMMs', 'GBMs']:
#     infile = os.path.join(module_dir, "{}_dataframe.txt".format(module_type))
#     outfile = os.path.join(module_dir, "KOs_to_{}.txt".format(module_type))
#     with open(infile) as input, open(outfile, "w") as output:
#         for line in input:
#             elements = line.strip().split("\t")
#             module_id = elements[0]
#             KO_list = elements[2]
#             module_desc = elements[3]
#             KOs = KO_list.split(",")
#             for KO in KOs:
#                 output.write("{}\t{}\t{}\n".format(KO, module_id, module_desc))
#
#     module_df = pd.read_table(outfile, header=0, index_col=0)
#     module_df.columns = ['Module', 'Description']
#     module_df['features'] = module_df.index.values + ' (' + module_df['Description'] + ')'
#     # print("kegg_pathways: ", kegg_pathways.iloc[0:5, ])
#     # print("kegg_kos: ", kegg_kos.iloc[0:5, ])


MMUPHin_df = pd.read_csv("/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff/KOs/MMUPHin_diff_merge.csv", header=0, index_col=0)
for study in studies:
    for stage in stages:
        meta_feats_filter = meta_feats.loc[(meta_feats['Group'].isin(['NC', stage])),]
        if study in ['CHN', 'CHN2', 'JPN', 'GER']:
            meta_feats_filter = meta_feats_filter.loc[meta_feats_filter['Batch'] == study,]
        elif study == "CHN+CHN2":
            meta_feats_filter = meta_feats_filter.loc[meta_feats_filter['Batch'].isin(['CHN', 'CHN2']),]

        meta_feats_filter = meta_feats_filter.loc[meta_feats_filter.iloc[:, 14:].sum(axis=1) > 0,]
        meta = meta_feats_filter.iloc[:, 0:14]
        features = meta_feats_filter.iloc[:, 14:]
        sample_list = meta_feats_filter.index.values.tolist()
        # X, feature_list = merge_features(features, feat_type)
        # X_df = pd.DataFrame(X, index=sample_list, columns=feature_list)
        # meta_feats_filter = pd.merge(meta, X_df, left_index=True, right_index=True)

        for model in models:
            outprefix = "CV_{}_{}_{}_{}_{}_None".format(study, stage, model, feat_type, fs_method)
            print("deal with {}".format(outprefix))
            infile = os.path.join(workdir, "06classification", '{}_{}'.format(study, taxon),
                                  model, '{}_feature_importance.csv'.format(outprefix))
            print("infile: {}".format(infile))
            if os.path.exists(infile):
                print("{} exists.".format(infile))
                data = pd.read_csv(infile, header=0, index_col=0)
                ##

                data['median'] = data.median(axis=1)
                data_sort = data.sort_values('median', ascending=False)
                #print(data['median'])
                data_sort = pd.merge(data_sort, merge_kos_to_pathways_desc_filter, left_index=True, right_index=True)
                # data_sort = data.sort_values('median', ascending=False)
                # data_sort = data_sort.iloc[0:20, ]

                #data_sort['features'] = data_sort.index.values.tolist()
                ####
                #data_sort['features'].update(kegg_pathways['features'])


                feature_list = data_sort.index.values.tolist()
                #print("feature_list: ", feature_list)
                for feature in feature_list:
                    if feature in MMUPHin_df.index.values.tolist():
                        coef = MMUPHin_df.loc[feature, stage]
                        if coef > 0:
                            data_sort.loc[feature, 'Significance'] = 'Elevation (P<0.05)'
                        elif coef < 0:
                            data_sort.loc[feature, 'Significance'] = "Depletion (P<0.05)"
                        else:
                            data_sort.loc[feature, 'Significance'] = "Not significant"
                    else:
                        data_sort.loc[feature, 'Significance'] = "Not significant"
                    # pvalue_greater = stats.ranksums(meta_feats_filter[feature][meta_feats_filter['Group']==stage],
                    #                         meta_feats_filter[feature][meta_feats_filter['Group']=='NC'],
                    #                         alternative='greater')[1]
                    # pvalue_less = stats.ranksums(meta_feats_filter[feature][meta_feats_filter['Group']==stage],
                    #                         meta_feats_filter[feature][meta_feats_filter['Group']=='NC'],
                    #                         alternative='less')[1]
                    # ##https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ranksums.html, 前者比后者less或者greater.
                    # if pvalue_greater < 0.05:
                    #     data_sort.loc[feature, 'Significance']  = 'Elevation (P<0.05)'
                    # elif pvalue_less < 0.05:
                    #     data_sort.loc[feature, 'Significance'] = "Depletion (P<0.05)"
                    # else:
                    #     data_sort.loc[feature, 'Significance'] = "Not significant"

                ####
                # data_sort.loc[data_sort.index.str.startswith("k__Archaea"), 'kingdom'] = "Species"
                # data_sort.loc[data_sort.index.str.startswith("k__Bacteria"), 'kingdom'] = "Species"
                # data_sort.loc[data_sort.index.str.startswith("k__Eukaryota|k__Fungi"), 'kingdom'] = "Species"
                # data_sort.loc[data_sort.index.str.startswith("k__Viruses"), 'kingdom'] = "Species"
                # ####
                # data_sort.loc[data_sort.index.str.startswith("k__Archaea"), 'features'] += " (A)"
                # data_sort.loc[data_sort.index.str.startswith("k__Bacteria"), 'features'] += "(B)"
                # data_sort.loc[data_sort.index.str.startswith("k__Eukaryota|k__Fungi"), 'features'] += " (F)"
                # data_sort.loc[data_sort.index.str.startswith("k__Viruses"), 'features'] += " (V)"
                #
                # ##
                # data_sort.loc[data_sort.index.str.startswith("K"), 'kingdom'] = "KOs"
                # data_sort.loc[data_sort.index.str.startswith("map"), 'kingdom'] = "Pathways"


                # print(data_sort.columns.values.tolist())
                #print(data_sort.iloc[0:5, ])  ###<------------------->###
                #data_sort['species'] = data_sort['KO_desc'].str.replace('.*s__', '', regex=True)
                data_sort['features'] = data_sort.index.values.tolist()
                data_sort.to_csv(os.path.join(outdir,
                                              "{}_{}_{}_feature_importance_sort.csv".format(study, stage, feat_type)))
                # data_sort.to_csv(os.path.join(outdir,
                #                               "{}_{}_{}_feature_importance_sort.csv".format(study, stage, feat_type)))
                print("data_sort.columns: {}".format(data_sort.columns.values.tolist()))
                data_sort_top20 = data_sort.iloc[0:20, ]
                data_melt = pd.melt(data_sort_top20, id_vars=['features', 'KO_name', 'KO_desc',  'Pathway_category2', 'Significance', 'median'],
                                    value_vars=[str(i) for i in range(0, repeats*kfold)], value_name='feature importance')
                #data_melt = data_melt.sort_values('median', ascending=False)
                data_melt['Contribution to model (%)'] = data_melt['feature importance'] * 100

                #print(data_melt.iloc[0:5, ])
                data_melt.to_csv(os.path.join(outdir, "{}_{}_{}_feature_importance.csv".format(study, stage, feat_type)))

                # fig, ax = plt.subplots(figsize=(4.2, 4.2))
                # # plt.xticks(fontsize=8)
                # # plt.yticks(fontsize=8)
                # plt.subplots_adjust(left=0.40, right=0.95, top=0.95, bottom=0.10)
                # palette = {"Elevation (P<0.05)": "#ff0000",
                #            "Depletion (P<0.05)": "#add8e6",
                #            "Not significant": "#999999"}
                # sns.boxplot(data=data_melt, hue='Significance', x='Contribution to model (%)', y='species',
                #             hue_order=["Elevation (P<0.05)", "Depletion (P<0.05)", "Not significant"],
                #             palette=palette, fliersize=1, dodge=False, ax=ax)
                # plt.xlim((0, 8))
                # plt.xticks(np.arange(0, 8, 2), fontsize=12)
                # plt.yticks(fontsize=10)
                # plt.xlabel('Contribution to model (%)', fontsize=12)
                # plt.ylabel('Features', fontsize=12)
                # plt.legend(fontsize=10)
                # #plt.tight_layout()
                # plt.title('{} vs NC (cross validation)'.format(stage), fontsize=14)
                # ax.set_yticklabels(labels=ax.get_yticklabels(), va='center')
                #
                # ##<--------------set xtick label colors based on label string---------------->##
                # ##https://stackoverflow.com/questions/71387431/set-xtick-label-colors-based-on-label-string
                # # for t in ax.get_yticklabels():
                # #     txt = t.get_text()
                # #     print(t.get_text())
                # #     if data_sort.loc[data_sort['species']==txt, 'kingdom'].values[0] == "Species":
                # #         t.set_color('#03ACDA')  ###蓝色
                # #     elif data_sort.loc[data_sort['species']==txt, 'kingdom'].values[0] == "KOs":
                # #         t.set_color('#6AB42D')  ##绿色
                # #     elif data_sort.loc[data_sort['species']==txt, 'kingdom'].values[0] == "Pathways":
                # #         t.set_color('#ff993f')  ##橙色
                #
                # ax.legend(loc="lower right", fontsize=8)
                # scatter_fig = fig.get_figure()
                # outfile = os.path.join(outdir, '{}_feature_importance.jpg'.format(outprefix))
                # scatter_fig.savefig(outfile, dpi=600)
                # outfile = os.path.join(outdir, '{}_feature_importance.pdf'.format(outprefix))
                # scatter_fig.savefig(outfile, dpi=600)

