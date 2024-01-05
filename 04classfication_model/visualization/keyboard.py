from functools import partial
from itertools import product
from scipy.spatial.distance import cdist
from scipy.sparse import csc_matrix
import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.metrics import silhouette_score

import os
##friendly_guacamole的定义：https://github.com/gwarmstrong/friendly-guacamole
import sys
sys.path.append("/home1/jialh/brain/01meta/multikingdom/scripts/visualization")
from friendly_guacamole.transforms import (
    FilterSamples,
    UniFrac,
    RarefactionBIOM,
    PCoA,
    AsDense,
    CLR,
)
from friendly_guacamole.datasets import KeyboardDataset
from skbio import DistanceMatrix
from umap import UMAP
from biom import Table
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

workdir="/home1/jialh/brain/01meta/multikingdom/scripts/visualization"
outdir=os.path.join(workdir, "results")
if not os.path.exists(outdir):
    os.makedirs(outdir)

keyboard_data = KeyboardDataset('data/keyboard')
tree = keyboard_data.apply('tree', 'path')

print("keyboard_data: ", keyboard_data)
print("tree: ", tree)
# Files already downloaded and verified
# keyboard_data:  <friendly_guacamole.datasets.KeyboardDataset object at 0x7ff155405e80>
# tree:  /tmp/tmp0zt3ljd8

RAREFACTION_DEPTH = 500

##
min_count_filter = FilterSamples(min_count=RAREFACTION_DEPTH)
kb_mf = keyboard_data['metadata']
kb_table = keyboard_data['table']
kb_mf = kb_mf.set_index('sample_name')

metadata = kb_mf[kb_mf.host_subject_id.isin(['M2','M3', 'M9'])]
kb_table = kb_table.filter(metadata.index)
table = min_count_filter.fit_transform(kb_table)  ##过滤掉低丰度的样本。
metadata = metadata.loc[table.ids('sample')]



##
print(type(table))
#table.to_csv(os.path.join(workdir, 'data', 'keyboard_table.csv'))
metadata.to_csv(os.path.join(workdir, 'data', 'keyboard_metadata.csv'))

##
rarefied_table = RarefactionBIOM(RAREFACTION_DEPTH).fit_transform(table)

def postprocess_umap(results):
    return pd.DataFrame(results, columns=[f'PC{i + 1}'
                                          for i in range(results.shape[1])])

aitchison_pipeline = Pipeline([
    ('asdense', AsDense()),
    ('clr', CLR()),
])

aitchison_prep_tables = [
    {
        'name': 'Aitchison',
        'metric': 'euclidean',
        'pipeline': aitchison_pipeline,
    },
]

prep_tables = [
    {
        'name': 'Aitchison',
        'metric': 'euclidean',
        'pipeline': aitchison_pipeline,
    },
]

def pcoa_amend_axes(transformer, axes_names):
    pe = transformer.ordination_.proportion_explained
    return [f'{axn} ({pexp:.1%})' for pexp, axn in zip(pe, axes_names)]

embedding_methods = [
    {
        'method': 'PCoA',
        'pipeline': PCoA,
        'axes': ['PCoA-1', 'PCoA-2'],
        'amend_axes': pcoa_amend_axes,
    },
    {
        'method': 'UMAP\nNeighbors=15',
        'pipeline': partial(UMAP,
                            min_dist=1,
                            random_state=724
                            ),
        'postprocess': postprocess_umap,
        'axes': ['UMAP-1', 'UMAP-2']
    },
    {
        'method': 'UMAP\nNeighbors=80',
        'pipeline': partial(UMAP,
                            min_dist=1,
                            n_neighbors=80,
                            random_state=825),
        'postprocess': postprocess_umap,
        'axes': ['UMAP-1', 'UMAP-2']
    },
    {
        'method': 'UMAP\nNeighbors=98',
        'pipeline': partial(UMAP,
                            min_dist=1,
                            n_neighbors=98,
                            random_state=901),
        'postprocess': postprocess_umap,
        'axes': ['UMAP-1', 'UMAP-2']
    },
]

results = dict()

for prep, emb in product(prep_tables, embedding_methods):
    metric = prep['metric']
    method = emb['method']
    name = prep['name']
    prepped_table = prep['pipeline'].fit_transform(table)
    transformer = emb['pipeline'](metric=metric)
    embedding = transformer.fit_transform(prepped_table)
    result = emb.get('postprocess', lambda x: x)(embedding)
    amend_axes = emb.get('amend_axes', lambda t, labels: labels)
    results[(name, method)] = {'ordination': result,
                               'axes': amend_axes(transformer, emb['axes']),
                               }

v_position_map = {x['name']: i for i, x in enumerate(prep_tables)}
h_position_map = {x['method']: i for i, x in enumerate(embedding_methods)}

metadata['host_surface'] = metadata['host_subject_id'] + \
                           metadata['sample_type']


fig, axs = plt.subplots(len(v_position_map), len(h_position_map),
                        figsize=(15, 5),
                        )
for (name, method), result in results.items():
        i = v_position_map[name]
        j = h_position_map[method]
        res = results[(name, method)]['ordination']
        res.index = metadata.index
        res = res.join(metadata)

        g = sns.scatterplot(
            x='PC1',
            y='PC2',
            hue='host_surface',
            hue_order=list(sorted(metadata['host_surface'].unique())),
            # style='sample_type',
            data=res,
            ax=axs[j],
            s=35,
            edgecolor='k',
            palette='Paired',
        )
        g.set_aspect('equal', 'datalim')
        g.legend().remove()

        g.set_title(f'{name}-{method}',
                    color='black',
                    fontsize=24)

        g.set_xlabel(result['axes'][0], color='black', fontsize=18)
        g.set_ylabel(result['axes'][1], color='black', fontsize=18)

plt.tight_layout()
plt.savefig(os.path.join(outdir, '2.0-real-data-keyboard-ordination-aitchison.png'), dpi=300)
plt.show()

