import os
import pandas as pd
import sys
sys.path.append("/home1/jialh/brain/01meta/multikingdom/scripts/visualization")
from dimensionality_reduction import run_PCA, run_PCoA, run_tsne, run_umap

##参考:
##1.预处理：https://towardsdatascience.com/analyzing-microbiome-data-320728b56b8e
##2.高维数据可视化: https://towardsdatascience.com/visualizing-high-dimensional-microbiome-data-eacf02526c3a
def convert_tsv(df):
    """
    Helper function to take the output of the biom-format package from a .biom file
    Converts it to a feature table
    """
    length = df.shape[0]
    df = df.reset_index().T
    df.set_index(0)

    # We remove the DNA sequence for an encoding, taking less space
    new_header = ['Sample_id'] + list(range(length - 1))
    df = df[1:]
    df.columns = new_header
    df = df.reset_index().drop('index', axis=1)

    return df

def prepare_data(biom_file, feature_table_file, metadata_file):
    os.system("biom convert -i {} -o {} --to-tsv".format(biom_file, feature_table_file))

    feature_table_raw = pd.read_csv(feature_table_file, sep = '\t', )
    feature_table_raw.head()

    sample_table = pd.read_csv(metadata_file, sep = '\t')

    feature_table = convert_tsv(feature_table_raw)
    feature_table = feature_table.set_index('Sample_id').astype(float).astype(int)

    # Subset sample table
    # print(sample_table.iloc[0:5,])
    sample_table = sample_table.set_index('sample_name')
    sample_table = sample_table.loc[sample_table['host_subject_id'].isin(['M2', 'M3', 'M9'])]

    feature_table = feature_table.loc[sample_table.index]
    #print(unique(sample_table['host_subject_id']))
    print(sample_table.groupby(['host_subject_id', 'sample_type_qiita'])['host_subject_id'].count())

    # Make sure the sequence has more than 20,000 reads total
    feature_table = feature_table.loc[(feature_table.sum(axis='columns') > 500)]

    # Make sure the sequence has been seen in at least 3 different samples
    feature_table = feature_table[feature_table.columns[((feature_table > 0).sum() > 3)]]

    # Check dimensions of our data
    print("过滤后特征表大小: ", feature_table.shape)  ##(99, 598) 正好99个样本。

    return feature_table, sample_table

def main():
    workdir = "/home1/jialh/brain/01meta/multikingdom/scripts/visualization"
    outdir = os.path.join(workdir, "results")
    biom_file = os.path.join(workdir, "data/keyboard/table/46809/46809.table.biom")
    feature_table_file = os.path.join(workdir, "data/keyboard/feature_table.tsv")
    metadata_file = os.path.join(workdir, "data/keyboard/metadata/232/232.metadata.tsv")
    feature_table, sample_table = prepare_data(biom_file, feature_table_file, metadata_file)
    ##<-------------------------------------------------------->##
    run_PCA(feature_table, sample_table, outdir, color_column='host_subject_id', shape_column='sample_type_qiita',  n_components=2)
    run_PCoA(feature_table, sample_table, outdir, color_column='host_subject_id', shape_column='sample_type_qiita', metric='jaccard')
    run_tsne(feature_table, sample_table, outdir, color_column='host_subject_id', shape_column='sample_type_qiita',
             n_components=2, perplexity=30.0, metric='jaccard')
    run_umap(feature_table, sample_table, outdir, color_column='host_subject_id', shape_column='sample_type_qiita',
             n_components=2, n_neighbors=30, min_dist=0.1, metric='jaccard')

if __name__ == "__main__":
    main()


