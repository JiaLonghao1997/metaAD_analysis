import os
import sys

import pandas as pd
import numpy as np

def merge_diff(workdir, stages, taxon_type, feat_types, diff_file):
    diff_df = pd.DataFrame()
    for stage in stages:
        stage_df = pd.DataFrame()
        for feat_type in feat_types:
            infile = os.path.join(workdir, feat_type, "NC_{}".format(stage), "metaAD_MMUPHin_result_{}_CHN+CHN2_NC_{}.csv".format(feat_type, stage))
            MMUPHin_df = pd.read_csv(infile, header=0, index_col=0)
            MMUPHin_df_filter = MMUPHin_df.loc[MMUPHin_df['pval'] < 0.05, MMUPHin_df.columns.isin(['coef'])]
            MMUPHin_df_filter.rename(columns={"coef": stage}, inplace=True)
            if stage_df is None:
                stage_df = MMUPHin_df_filter
            else:
                stage_df = pd.concat([stage_df, MMUPHin_df_filter], axis=0)
        if diff_df is None:
            diff_df = stage_df
        else:
            diff_df = pd.merge(diff_df, stage_df, left_index=True, right_index=True, how="outer")

    diff_df.to_csv(diff_file, index=True, index_label=taxon_type)
    return diff_df

def merge_diff_pvalue(workdir, stages, taxon_type, feat_types, diff_file):
    diff_df_pvalue = pd.DataFrame()
    for stage in stages:
        stage_df = pd.DataFrame()
        for feat_type in feat_types:
            infile = os.path.join(workdir, feat_type, "NC_{}".format(stage), "metaAD_MMUPHin_result_{}_CHN+CHN2_NC_{}.csv".format(feat_type, stage))
            MMUPHin_df = pd.read_csv(infile, header=0, index_col=0)
            MMUPHin_df_filter = MMUPHin_df.loc[MMUPHin_df['pval'] < 0.05, MMUPHin_df.columns.isin(['pval'])]
            MMUPHin_df_filter.rename(columns={"pval": stage}, inplace=True)
            if stage_df is None:
                stage_df = MMUPHin_df_filter
            else:
                stage_df = pd.concat([stage_df, MMUPHin_df_filter], axis=0)
        if diff_df_pvalue is None:
            diff_df_pvalue = stage_df
        else:
            diff_df_pvalue = pd.merge(diff_df_pvalue, stage_df, left_index=True, right_index=True, how="outer")

    diff_df_pvalue.to_csv(diff_file, index=True, index_label=taxon_type)
    return diff_df_pvalue

def abundance_heatmap(infile, feat_types, diff_df, outfile):
    meta_feats = pd.read_csv(infile, header=0, index_col=0)
    meta_feats = meta_feats.loc[~meta_feats.index.isin(['B1009', 'B1521']),]
    meta_feats_filter = meta_feats.loc[meta_feats['Study'].isin(["CHN1", "CHN2"])]
    feats = meta_feats_filter.iloc[:, 14:]
    feats_rel_mean = pd.DataFrame()
    for feat_type in feat_types:
        feat_type_dic = {"A": "k__Archaea", "B": "k__Bacteria", "F": "k__Eukaryota|k__Fungi", "V": "k__Viruses"}
        feat_type_str = feat_type_dic[feat_type]
        feats_filter = feats.loc[:, feats.columns.str.contains(feat_type_str)]
        # convert to relative abundance
        feats_rel = feats_filter.div(feats_filter.sum(axis=1), axis=0)  ##row is samples, column is features.
        feats_rel = feats_rel.T
        feats_rel['mean'] = feats_rel.mean(axis=1)
        if feats_rel_mean is None:
            feats_rel_mean = feats_rel['mean']
        else:
            feats_rel_mean = pd.concat([feats_rel_mean, feats_rel['mean']], axis=0)

    feats_rel_mean.columns = ['mean']
    feats_rel_mean_filter = feats_rel_mean.loc[feats_rel_mean.index.isin(diff_df.index.values.tolist())]
    feats_rel_mean_filter.index = feats_rel_mean_filter.index.str.replace('sp.', 'sp', regex=False)
    feats_rel_mean_filter.index = feats_rel_mean_filter.index.str.replace('|', '.', regex=False)
    feats_rel_mean_filter.index = feats_rel_mean_filter.index.str.replace('k__Heunggongvirae.', '', regex=False)
    feats_rel_mean_filter.index = feats_rel_mean_filter.index.str.replace('k__Eukaryota.', '', regex=False)
    feats_rel_mean_filter.index = feats_rel_mean_filter.index.str.replace('.f__Siphoviridae.', '.o__.f__Siphoviridae.g__.', regex=False)
    feats_rel_mean_filter['mean'] = np.log10(feats_rel_mean_filter['mean'])
    feats_rel_mean_filter.to_csv(outfile, index=True, index_label="species")
    return feats_rel_mean_filter

def main():
    workdir = "/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff"
    #taxon_types = ['KOs', "GMMs", "GBMs"]
    taxon_types = ['species']
    for taxon_type in taxon_types:
        outdir = os.path.join(workdir, taxon_type)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        graphlan = "/home1/jialh/tools/anaconda3/envs/graphlan/bin/graphlan.py"
        graphlan_annotate  = "/home1/jialh/tools/anaconda3/envs/graphlan/bin/graphlan_annotate.py"
        #feat_types = ['A', 'B', 'F', 'V']
        if taxon_type in ['KOs', "GMMs", "GBMs", "pathways"]:
            feat_types = [taxon_type]
        stages = ['AD', 'MCI', 'SCD', 'SCS']
        diff_file = os.path.join(outdir, "MMUPHin_diff_merge.csv")
        #os.remove(diff_file)
        if not os.path.exists(diff_file):
            diff_df = merge_diff(workdir, stages, taxon_type, feat_types, diff_file)
        else:
            diff_df = pd.read_csv(diff_file, header=0, index_col=0)

        ##
        diff_pvalue_file = os.path.join(outdir, "MMUPHin_diff_merge_pvalue.csv")
        #os.remove(diff_pvalue_file)
        if not os.path.exists(diff_pvalue_file):
            diff_df_pvalue = merge_diff_pvalue(workdir, stages, taxon_type, feat_types, diff_pvalue_file)
        else:
            diff_df_pvalue = pd.read_csv(diff_pvalue_file, header=0, index_col=0)

        #print(diff_df_pvalue)

        if taxon_type in ['species', "KOs"]:
            # diff_df_filter = diff_df.loc[(diff_df<0.001).sum(axis=1) >= 1,]
            # print(diff_df_filter.columns.values.tolist())
            # diff_df_filter = diff_df_filter[['SCS', 'SCD', 'MCI', 'AD']]
            # diff_df_filter.to_csv(os.path.join(outdir, "MMUPHin_diff_merge_filter0.01.csv"))
            ###
            diff_df_pvalue_filter = diff_df_pvalue.loc[(diff_df_pvalue<0.001).sum(axis=1) >= 1,]
            print(diff_df_pvalue_filter.columns.values.tolist())
            diff_df_pvalue_filter = diff_df_pvalue_filter[['SCS', 'SCD', 'MCI', 'AD']]
            diff_df_pvalue_filter.to_csv(os.path.join(outdir, "MMUPHin_diff_merge_pvalue_filter0.001.csv"))
        elif taxon_type in ['GMMs', "pathways", "GBMs"]:
            diff_df_filter = diff_df #.loc[diff_df.isna().sum(axis=1) <= 2,]
            print(diff_df_filter.columns.values.tolist())
            diff_df_filter = diff_df_filter[['SCS', 'SCD', 'MCI', 'AD']]
            diff_df_filter.to_csv(os.path.join(outdir, "MMUPHin_diff_merge_filter.csv"))
            ###
            diff_df_pvalue_filter = diff_df_pvalue #.loc[diff_df_pvalue.isna().sum(axis=1) <= 2,]
            print(diff_df_pvalue_filter.columns.values.tolist())
            diff_df_pvalue_filter = diff_df_pvalue_filter[['SCS', 'SCD', 'MCI', 'AD']]
            diff_df_pvalue_filter.to_csv(os.path.join(outdir, "MMUPHin_diff_merge_pvalue_filter.csv"))


        sys.exit()
    #     diff_df2 = diff_df.copy()
    #     diff_df2.index = diff_df2.index.str.replace('|', '.', regex=False)
    #     diff_df2.index = diff_df2.index.str.replace('sp.', 'sp', regex=False)
    #     diff_df2.index = diff_df2.index.str.replace('k__Heunggongvirae.', '', regex=False)
    #     diff_df2.index = diff_df2.index.str.replace('k__Eukaryota.', '', regex=False)
    #     ##.f__Siphoviridae.
    #     diff_df2.index = diff_df2.index.str.replace('.f__Siphoviridae.', '.o__.f__Siphoviridae.g__.', regex=False)
    #
    #     ###abundance heatmap
    #     diff_heatmap_gradient = os.path.join(outdir, "diff_heatmap_gradient.txt")
    #     if not os.path.exists(diff_heatmap_gradient):
    #         inputdir = "/home1/jialh/brain/01meta/multikingdom/00profile"
    #         meta_feats = os.path.join(inputdir, "metadata_species.csv")
    #         diff_heatmap = os.path.join(outdir, "diff_heatmap.csv")
    #         feats_rel_mean_filter = abundance_heatmap(infile=meta_feats, feat_types=feat_types, diff_df=diff_df, outfile=diff_heatmap)
    #         print(feats_rel_mean_filter.describe())
    #         feats_rel_mean_filter['mean'] = 1 - (feats_rel_mean_filter['mean'] + 5) * 0.2
    #         ##如果相对丰度为1，则alpha为0
    #         ##如果相对丰度为0.1, 则alpha为1-(-1+5)*0.2=0.2.
    #         ##如果相对丰度为0.01， 则alpha为1-（-2+5)*0.2=0.4
    #         feats_rel_mean_filter.insert(loc=0, column='ring_label', value="ring_alpha")
    #         feats_rel_mean_filter.insert(loc=1, column='ring_order', value=1)
    #         feats_rel_mean_filter.to_csv(diff_heatmap_gradient, sep="\t", index=True, header=False)
    #
    #
    #
    # ##是否为差异丰度菌。
    # diff_ring = os.path.join(outdir, "diff_ring.txt")
    # if not os.path.exists(diff_ring):
    #     stage_dic = {'SCS':2, 'SCD':3, 'MCI':4, 'AD':5}
    #     with open(diff_ring, "w") as output:
    #         for index, row in diff_df2.iterrows():
    #             for stage in stages:
    #                 species = index
    #                 ring_order = stage_dic[stage]
    #                 cof = row[stage]
    #                 if cof > 0:
    #                     ring_color = "#ff6836"
    #                 elif cof < 0:
    #                     ring_color = "#36bccb"
    #                 else:
    #                     ring_color = "#ffffff"
    #                 output.write("{}\tring_color\t{}\t{}\n".format(species, ring_order, ring_color))
    #
    # diff_text = os.path.join(outdir, "MMUPHin_diff_species.csv")
    # if not os.path.exists(diff_text):
    #     with open(diff_text, "w") as output:
    #         for species in diff_df2.index.values.tolist():
    #             output.write(species + "\n")
    # ###
    # annotate_txt = os.path.join(outdir, "graphlan_annotation.txt")
    # annotate_template = os.path.join(outdir, "graphlan_annotation_template.txt")
    # os.system("cat {} > {}".format(annotate_template, annotate_txt))
    #
    # diff_df2_filter = diff_df2.loc[diff_df2.isna().sum(axis=1) <= 2, ]
    #
    # print("diff_df2_filter.shape: ", diff_df2_filter.shape)
    # with open(annotate_txt, "a+") as output:
    #     output.write("\n")
    #     for diff_in_two_stages in diff_df2_filter.index.values.tolist():
    #         output.write("{}\tclade_marker_size\t20\n".format(diff_in_two_stages))
    #         output.write("{}\tclade_marker_shape\t*\n".format(diff_in_two_stages))
    #     output.write("\n")
    # os.system("cat {} >> {}".format(diff_heatmap_gradient, annotate_txt))
    # os.system("cat {} >> {}".format(diff_ring, annotate_txt))
    #
    # tree_img = os.path.join(outdir, "tree.pdf")
    # tree_xml = os.path.join(outdir, "tree.xml")
    # os.system("{} --annot {} {} {}".format(graphlan_annotate, annotate_txt, diff_text, tree_xml))
    # os.system("{} {} {} --dpi 300 --size 4".format(graphlan, tree_xml, tree_img))


if __name__ == '__main__':
    main()