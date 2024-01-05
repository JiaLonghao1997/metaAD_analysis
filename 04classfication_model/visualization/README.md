# umap-microbiome-benchmarking

>   -   Download from：https://github.com/knightlab-analyses/umap-microbiome-benchmarking/tree/main
>   -   Paper: Armstrong G, Martino C, Rahman G, et al.Uniform manifold approximation and projection (UMAP) reveals composite patterns and resolves visualization artifacts in microbiome data[J]. Msystems, 2021, 6(5): e00691-21. 

## Installation

These notebooks' dependencies are listed in the `requirements.yml` file.

Additionally, using the notebooks will require Jupyter.

```bash
ENV_NAME=umap-benchmarking
conda create -n ${ENV_NAME}
conda activate ${ENV_NAME}
conda env update --file requirements.yml
jupyter notebook
```

## Figures

```bash
cd notebooks
mkdir results
jupyter notebook
```

Then these notebooks contain the analyses:
* `real-data-keyboard-benchmark.ipynb`
* `real-data-soil-benchmark.ipynb`
* `technical-effects-hmpv13v35.ipynb`
