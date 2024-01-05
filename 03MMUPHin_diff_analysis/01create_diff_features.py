import pandas as pd
import os
taxon_types = ["species", "KOs", "GMMs", "GBMs", "pathways"]
for taxon_type in taxon_types:
    print("deal with {}".format(taxon_type))
    profiledir = "/home1/jialh/brain/01meta/multikingdom/00profile"
    meta_features_file = os.path.join(profiledir, "metadata_{}+ConQuR.csv".format(taxon_type))
    if os.path.exists(meta_features_file):
        meta_features = pd.read_csv(meta_features_file, header=0, index_col=0)
        meta_df = meta_features.iloc[:, 0:14]
        features_df = meta_features.iloc[:, 14:]

        MMUPHin_diff_dir = "/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff"
        MMUPHin_diff_file = os.path.join(MMUPHin_diff_dir, taxon_type, "MMUPHin_diff_merge.csv")
        MMUPHin_diff_df = pd.read_csv(MMUPHin_diff_file, header=0, index_col=0)
        diff_feature_list = MMUPHin_diff_df.index.values.tolist()

        features_df_diff = features_df.loc[:, features_df.columns.isin(diff_feature_list)]
        meta_features_diff = pd.merge(meta_df, features_df_diff, left_index=True, right_index=True)
        meta_features_diff.to_csv(os.path.join(profiledir, "metadata_{}+ConQuR_MMUPHin.csv".format(taxon_type)), index=True, index_label="Sample")
    else:
        print("{} does not exist.".format(meta_features_file))