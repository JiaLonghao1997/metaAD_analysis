# Metagenomic analysis reveals stage-specific roles of  the gut microbiota in Alzheimer's Disease 
### 0. Introduction

Here, we enrolled the 476 individuals across five stages of Alzheimer's Disease (AD) pathology (including 63 NC, 82 SCS, 90 SCD, 119 MCI, and 122 AD samples) and performed deep shotgun metagenomic sequencing. We then assessed stage-specific intestinal dysbiosis and identified differential microbial taxonomic and functional signatures between different AD stages and NC. Furthermore, we constructed taxonomic and functional classifiers with both cross-validation and external validation to evaluate the diagnostic potential of gut microbial signatures. We also interpreted the optimal classification model, identified the most discriminating features and discussed the potential stage-specific mechanisms of these features in AD pathology.

![Figure1-1225_1](https://jialh.oss-cn-shanghai.aliyuncs.com/img2/Figure1-1225_1.jpg)

### 1. Preprocessing, taxonomic and functional profiling

#### 1.1 Shotgun sequencing data preprocessing

-   The [hts_SuperDeduper](https://s4hts.github.io/HTStream/) (v.1.3.0) was employed to remove duplicated reads. 
-   The [Trim Galore](https://github.com/FelixKrueger/TrimGalore) (v.0.6.7) was used to remove adaptors, trim low-quality bases, and discarded short reads.
-   The remaining reads were then aligned to the mammalian genome, bacterial plasmids, complete plastomes and UniVec sequences to remove contaminant reads using [Burrows–Wheeler Aligner](https://bio-bwa.sourceforge.net/) (BWA) MEM v.0.7.17-r1188.

#### 1.2 Metagenomic taxonomic profiling

-   Preprocessed reads were assigned taxonomic classifications for bacteria, archaea, fungi, and viruses using [Kraken2](https://ccb.jhu.edu/software/kraken2/) v.2.1.2.
-   [Bracken](https://ccb.jhu.edu/software/bracken/) v.2.5.0 was utilized to accurately estimate taxonomic abundance, particularly at the species and genus levels.

#### 1.3 Metagenomic functional profiling

-   Preprocessed reads were assembled into contigs with [Megahit](https://github.com/voutcn/megahit) v.1.2.9 using ‘meta-sensitive’ parameters. 
-   Gene prediction was performed using [Prodigal](https://github.com/hyattpd/Prodigal) v.2.6.3 in the metagenome mode (-p meta). 
-   To construct a non-redundant microbial gene reference, [CD-HIT](http://cd-hit.org/) v. 4.8.1 was employed with a sequence identity cut-off of 0.95 and a minimum coverage cut-off of 0.9 for shorter sequences. 
-   The reference was annotated with [EggNOG mapper](http://eggnog-mapper.embl.de/) v.2.1.0 based on EggNOG orthology data. 
-   The abundance of microbial genes was estimated by mapping the high-quality reads to the reference sequences using CoverM v.0.4.0 (https://github.com/wwood/CoverM).

### 2. Diversity and confounder analysis

-   Alpha diversity metrics, such as Shannon and Simpson Indices of all kingdoms were calculated for each sample. 
-   In addition, beta diversity was assessed based on Bray–Curtis distance. 
-   Permutational multivariate analysis of variance (PERMANOVA) was performed by [Vegan](https://github.com/vegandevs/vegan) to investigate the microbial community differences between disease groups or batches with 999 permutations.

### 3. Differential analysis  by MMUPHin

MMUPHin was used to identify AD-related differential microbial species, which enables the normalization and combination of multiple microbial community studies.

### 4.  Microbial classification models for Alzheimer's Disease

To determine the diagnostic potential in differentiating AD stages from healthy controls using microbial signatures, we constructed classification models based on the taxonomic (86 Archaea, 343 Bacteria, 50 Fungi and 54 Viruses) and functional (6593 KOs, 103 GMMs, 46 GBMs and 407 pathways) profiles using random forest algorithm. An iterative feature elimination method was used for feature selection and performance was evaluated with 20 randomized 5-fold cross-validations using the Python package ‘[scikit-learn](https://scikit-learn.org/)’.  

To further test the generalizability of AD microbial markers, we validated the diagnostic model by cohort-to-cohort transfer validation as previously described. Briefly, the diagnostic model was trained on our dataset (CHN, 63 NC, 119 MCI, and 122 AD samples) and then validated on the Japanese dataset (JPN, 21 NC, 14 MCI and 7 AD samples). 

### 5. Mediation analysis between microbial species and gene families

For differential important KOs in the classification models, microbial species were associated with these KOs using Spearman correlation coefficients. To investigate the links between microbial species, KOs and AD stages, a mediation analysis was carried out using the mediate function in the R [mediation](https://cran.r-project.org/web/packages/mediation/) package, with batch, gender, age and BMI adjusted as covariates. 
