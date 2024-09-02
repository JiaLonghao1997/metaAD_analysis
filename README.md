# Metagenomic analysis reveals stage-specific roles of  the gut microbiota in Alzheimer's Disease 
### 0. Introduction

Here, we enrolled the 476 individuals across five stages of Alzheimer's Disease (AD) pathology (including 63 NC, 82 SCS, 90 SCD, 119 MCI, and 122 AD samples) and performed deep shotgun metagenomic sequencing. We then assessed stage-specific intestinal dysbiosis and identified differential microbial taxonomic and functional signatures between different AD stages and NC. Furthermore, we constructed taxonomic and functional classifiers with both cross-validation and external validation to evaluate the diagnostic potential of gut microbial signatures. We also interpreted the optimal classification model, identified the most discriminating features and discussed the potential stage-specific mechanisms of these features in AD pathology.

![Figure1_1](https://jialh.oss-cn-shanghai.aliyuncs.com/img2/Figure1_1.jpg)

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

### 3. Differential analysis  by MMUPHin and MaAsLin2

Batch effect correction was performed using the [MMUPHin](https://bioconductor.org/packages/release/bioc/html/MMUPHin.html) v.1.12.0, which employs a ComBat-like approach designed for zero-inflated microbial abundance data (MMUPHin_Correct). Subsequent analysis was performed using one of the most popular and powerful differential abundance tools [41], the [MaAsLin2](https://www.bioconductor.org/packages/release/bioc/html/Maaslin2.html) v.1.12.0, to identify microbial signatures associated with AD pathology. MaAsLin2 conducted multiple regression analyses using generalized linear and mixed models, incorporating default parameters and adjusting for covariates including sequencing batch, diet, age and BMI. The Benjamini–Hochberg procedure was used to adjust P values for multiple hypothesis testing. The default *q*-value threshold of 0.25 in MaAsLin2, which has been commonly accepted in microbial exploratory studies, was used to identify significant signatures.

### 4.  Microbial classification models for Alzheimer's Disease

To determine the diagnostic potential in differentiating AD stages from healthy controls using microbial signatures, we constructed classification models based on the taxonomic (86 Archaea, 343 Bacteria, 50 Fungi and 54 Viruses) and functional (6593 KOs, 103 GMMs, 46 GBMs and 407 pathways) profiles using random forest algorithm. An iterative feature elimination method was used for feature selection and performance was evaluated with 20 randomized 5-fold cross-validations using the Python package ‘[scikit-learn](https://scikit-learn.org/)’.  

To further test the generalizability of AD microbial markers, we validated the diagnostic model by cohort-to-cohort transfer validation as previously described. Briefly, the diagnostic model was trained on our dataset (CHN, 63 NC, 119 MCI, and 122 AD samples) and then validated on the Japanese dataset (JPN, 21 NC, 14 MCI and 7 AD samples). 

### 5. Validation by an in vivo gut simulator

To validate the robustness of the stage-specific taxonomic and functional alterations in AD pathology, we constructed an *in vitro* gut simulator. The gut simulator is an anaerobic single-chamber fermenter designed to mimic human digestion processes and has been widely used to explore the metabolic potential of gut microbiota in vitro (T&J-Minibox5 1.3L*4 Intelli-Ferm, see Section 1.7 of Supplementary Materials for details) [29, 30]. In our experiment, five fecal samples were selected from participants at different stages of AD pathology, including B2118 for NC, B2374 for SCS, B2238 for SCD, B1930 for MCI, and B2264 for AD. These fecal samples were first diluted to 4% liquid and then anaerobically cultured in a Coy chamber (5% H2, 10% CO2, and 85% N2) using 50 ml of brain-heart infusion (BHI) medium with tryptone (10.0 g/l), Na2HPO4 (2.5 g/l), beef heart powder (17.5 g/l), NaCl (5.0 g/l) and glucose (2.0 g/l). After 12 hours of incubation at 37°C, 5% of this pre-culture was transferred into the luminal compartment of the gut simulator. The cultures were sampled initially at baseline and then daily for one week. The samples were stored at −80°C and then subjected to shotgun sequencing and metagenomic analysis using the same methods as for the discovery cohort. If a differential signature identified in the discovery cohort was also identified as a differential signature with consistent alterations in the gut simulator (MaAsLin2, *q*-value < 0.25), it was considered to be validated.

