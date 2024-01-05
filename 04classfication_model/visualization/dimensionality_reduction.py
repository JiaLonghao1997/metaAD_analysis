import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from skbio.stats.ordination import pcoa
from skbio.diversity import beta_diversity
import umap
import os

def run_PCA(feature_table, sample_table, outdir, title, color_column, shape_column,  n_components=2):
    pca = PCA(n_components=n_components)
    principalComponents = pca.fit_transform(np.array(feature_table))
    np.savetxt(os.path.join(outdir, "PCA_principalComponents.csv"), principalComponents,  delimiter=",")
    prop1 = 'PC1 ({:.2%})'.format(pca.explained_variance_ratio_[0])
    prop2 = 'PC2 ({:.2%})'.format(pca.explained_variance_ratio_[1])
    plot_df = pd.DataFrame(data=principalComponents, columns=[prop1, prop2], index=feature_table.index)
    # We can then plot this like so:
    plot_df = pd.concat([plot_df, sample_table[[color_column, shape_column]]], axis=1)
    PCA_img = sns.scatterplot(x=prop1, y=prop2, hue=color_column, style=shape_column, data=plot_df)
    PCA_img.set(title="PCA: {}".format(title))
    PCA_img = PCA_img.get_figure()
    PCA_img.savefig(os.path.join(outdir, "PCA_img_with_label.png"), dpi=400)
    plt.clf()
    #return PCA_img

##参考：https://readiab.org/machine-learning.html#principal-coordinates-analysis-pcoa
def run_PCoA(feature_table, sample_table, outdir, title, color_column, shape_column, metric='braycurtis'):
    #beta多样性: http://scikit-bio.org/docs/0.4.2/generated/skbio.diversity.beta_diversity.html
    #距离度量方法: https://scipy.github.io/devdocs/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist
    distance_matrix = beta_diversity(metric=metric, counts=feature_table, ids=feature_table.index)
    pcoa_ordination = pcoa(distance_matrix)
    # print("pcoa_ordination.proportion_explained:", pcoa_ordination.proportion_explained)
    # print("type(pcoa_ordination.proportion_explained):", type(pcoa_ordination.proportion_explained))
    # print("type(pcoa_ordination.samples):", type(pcoa_ordination.samples))
    pcoa_ordination.samples.to_csv(os.path.join(outdir, "PCoA_ordination.csv"))
    prop1 = 'PCoA_1 ({:.1%})'.format(pcoa_ordination.proportion_explained[0])
    prop2 = 'PCoA_2 ({:.1%})'.format(pcoa_ordination.proportion_explained[1])
    print("prop1: {}, prop2: {}".format(prop1, prop2))
    plot_df = pd.DataFrame(data={prop1: pcoa_ordination.samples['PC1'], prop2: pcoa_ordination.samples['PC2']},
                           columns=[prop1, prop2], index=feature_table.index)

    plot_df = pd.concat([plot_df, sample_table[[color_column, shape_column]]], axis=1)
    PCoA_img = sns.scatterplot(x=prop1, y=prop2, hue=color_column, style=shape_column, data=plot_df)
    PCoA_img.set(title="PCoA: {}".format(title))
    PCoA_img = PCoA_img.get_figure()
    PCoA_img.savefig(os.path.join(outdir, "PCoA_img_with_label.png"), dpi=400)
    plt.clf()
    #return PCoA_img

##https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html#sklearn.manifold.TSNE
def run_tsne(feature_table, sample_table, outdir, title, color_column, shape_column, n_components=2, perplexity=30.0, metric='braycurtis'):
    if perplexity < sample_table.shape[0]:
        perplexity = round(0.3*sample_table.shape[0])
    tsne = TSNE(n_components=n_components, metric=metric, perplexity=perplexity)
    embeddings = tsne.fit_transform(feature_table)
    np.savetxt(os.path.join(outdir, "tsne_embeddings.csv"), embeddings, delimiter=",")
    plot_df = pd.DataFrame(data=embeddings, columns=['dim1', 'dim2'], index=feature_table.index)

    plot_df = pd.concat([plot_df, sample_table[[color_column, shape_column]]], axis=1)
    tsne_img = sns.scatterplot(x='dim1', y='dim2', hue=color_column, style=shape_column, data=plot_df)
    tsne_img.set(title="t-SNE: {}".format(title))
    tsne_img = tsne_img.get_figure()
    tsne_img.savefig(os.path.join(outdir, "tsne_img_with_label.png"), dpi=400)
    plt.clf()
    #eturn tsne_img

##https://umap-learn.readthedocs.io/en/latest/
##UMAP包括4个主要的超参数：n_neighbors, min_dist, n_components, metric
def run_umap(feature_table, sample_table, outdir, title, color_column, shape_column, n_components=2, n_neighbors=15, min_dist=0.1, metric='braycurtis'):
    reducer = umap.UMAP(n_components=n_components, n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, random_state=0)
    embeddings = reducer.fit_transform(feature_table)
    np.savetxt(os.path.join(outdir, "umap_embeddings.csv"), embeddings, delimiter=",")
    plot_df = pd.DataFrame(data=embeddings, columns=['dim1', 'dim2'], index=feature_table.index)
    # Like before, we can plot our lower dimensional embedding with labels:
    plot_df = pd.concat([plot_df, sample_table[[color_column, shape_column]]], axis=1)
    umap_img = sns.scatterplot(x='dim1', y='dim2', hue=color_column, style=shape_column, data=plot_df)
    umap_img.set(title="UMAP: {}".format(title))
    umap_img = umap_img.get_figure()
    umap_img.savefig(os.path.join(outdir, "umap_img_with_label.png"), dpi=400)
    plt.clf()
    #return umap_img