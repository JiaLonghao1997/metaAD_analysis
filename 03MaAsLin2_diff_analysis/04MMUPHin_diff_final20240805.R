#print("hello world!")
#install.packages("dplyr")
library(dplyr)
library(vegan)
library(ggplot2)
library(viridis)
library(cowplot)
library(plyr)
library("readxl")
library(ggpubr)
library(ggsignif)
#library(mia)
library(phyloseq)
library(MMUPHin)
library(Maaslin2)
#library(ConQuR)
rm(list = ls())
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

#data("CRC_abd", "CRC_meta")


MMUPHin_Batch <- function(metaAD_meta, metaAD_abd, outdir, group, feat_type, Batch){
    ####Step1: 批次矫正。
    ## feature_abd表的行为特征(species, KOs, pathways), 列为样本名称。
    fit_adjust_Batch <- adjust_batch(feature_abd = metaAD_abd,
                                     batch = "Batch",
                                     covariates = "Group", 
                                     data = metaAD_meta,
                                     control = list(verbose = FALSE))
    
    metaAD_abd_adj <- fit_adjust_Batch$feature_abd_adj
    
    ## ----permanova----------------------------------------------------------------
    library(vegan, quietly = TRUE)
    ##vegdist: 行为样本，类为特征。
    #any(is.na(metaAD_abd))
    any(rowSums(x=t(metaAD_abd), na.rm = TRUE) == 0)
    D_before <- vegdist(t(metaAD_abd), method="bray")
    D_after <- vegdist(t(metaAD_abd_adj), method="bray")
    
    set.seed(1)
    fit_adonis_before <- adonis2(D_before ~ Batch + Group, data = metaAD_meta, by="margin")
    fit_adonis_after <- adonis2(D_after ~ Batch + Group, data = metaAD_meta, by="margin")
    # print(fit_adonis_before)
    # print(fit_adonis_after)
    fit_adonis_before_df <- scores(fit_adonis_before)
    fit_adonis_after_df <- scores(fit_adonis_after)
    
    write.csv(fit_adonis_before_df, file=file.path(outdir, group, paste0("MMUPHin_before_adjust_Batch_", feat_type, '_', Batch, '_', group, ".csv")))
    write.csv(fit_adonis_after_df, file=file.path(outdir, group, paste0("MMUPHin_after_adjust_Batch_", feat_type, '_', Batch, '_', group, ".csv")))
    
    ###使用ConQuR来可视化批次矫正前后的结果。
    Batchid <- metaAD_meta[, "Batch"]
    pdf(file.path(outdir, group, paste0("MMUPHin_adjust_Batch_", feat_type, '_', Batch, '_', group, ".pdf")), width=6, height=6)
    par(mfrow=c(2, 2))
    
    ###Plot_PCoA的输入: TAX中行为样本,列为物种。
    ConQuR::Plot_PCoA(TAX=t(metaAD_abd), factor=Batchid, main="Before Correction, Bray-Curtis")
    ConQuR::Plot_PCoA(TAX=t(metaAD_abd_adj), factor=Batchid, main="MMUPHin, Bray-Curtis")
    
    ConQuR::Plot_PCoA(TAX=t(metaAD_abd), factor=Batchid, dissimilarity="Aitch", main="Before Correction, Aitchison")
    ConQuR::Plot_PCoA(TAX=t(metaAD_abd_adj), factor=Batchid, dissimilarity="Aitch", main="Plot_PCoA, Aitchison")
    dev.off()
    return(metaAD_abd_adj)
}

#BiocManager::install("MMUPHin")
##参考: https://microbiome.github.io/course_2022_turku/beta-diversity.html#pcoa-aggregated-to-phylum-level
MMUPHin_diff <- function(workdir, taxon_type, metadata_taxon, Batch, group, feat_type, outdir){
    # Batch = "CHN+CHN2"
    # taxon_type = "GBMs"
    # feat_type = "GBMs"
    # group = "NC_SCD"
    ###############################################################################################
    #outdir = file.path(workdir, "03MMUPHin_diff", taxon_type)
    dir.create(file.path(outdir, group), showWarnings = FALSE)
    
    #####
    case = strsplit(group, split="_", fixed=TRUE)[[1]][2]
    dir.create(file.path(outdir, group), showWarnings = FALSE)
    
    if(Batch=="all"){
        metadata_taxon_filter = metadata_taxon
    }else if(Batch=='CHN+CHN2'){
        metadata_taxon_filter = metadata_taxon[metadata_taxon$Batch %in% c('CHN', 'CHN2'), ]
        metadata_taxon_filter$Batch = factor(metadata_taxon_filter$Batch, levels=c('CHN', 'CHN2'))
    }
    
    if (group=="NC_AD"){
        metadata_taxon_filter = metadata_taxon_filter[metadata_taxon_filter$Group %in% c('NC', 'AD'), ]
        metadata_taxon_filter$Group = factor(metadata_taxon_filter$Group, levels=c('NC', 'AD'))
        case = "AD"
        #colors = c('#ea716d', '#b376b1')
    }else if (group=="NC_MCI"){
        metadata_taxon_filter = metadata_taxon_filter[metadata_taxon_filter$Group %in% c('NC', 'MCI'), ]
        metadata_taxon_filter$Group = factor(metadata_taxon_filter$Group, levels=c('NC', 'MCI'))
        case = "MCI"
        #colors = c('#ea716d', "#49a1d0")
    }else if (group=="NC_SCD"){
        metadata_taxon_filter = metadata_taxon_filter[metadata_taxon_filter$Group %in% c('NC', 'SCD'), ]
        metadata_taxon_filter$Group = factor(metadata_taxon_filter$Group, levels=c('NC', 'SCD'))
        case = "SCD"
        #colors = c('#ea716d', "#37b178")
    }else if (group=="NC_SCS"){
        metadata_taxon_filter = metadata_taxon_filter[metadata_taxon_filter$Group %in% c('NC', 'SCS'), ]
        metadata_taxon_filter$Group = factor(metadata_taxon_filter$Group, levels=c('NC', 'SCS'))
        case = "SCS"
        #colors = c('#ea716d', "#9f9c23")
    }else if (group=="all"){
        metadata_taxon_filter = metadata_taxon_filter[metadata_taxon_filter$Group %in% c('NC', 'SCS', 'SCD', 'MCI', 'AD'), ]
        metadata_taxon_filter$Group = factor(metadata_taxon_filter$Group, levels=c('NC', 'SCS', 'SCD', 'MCI', 'AD'))
    }
    
    metadata <- metadata_taxon_filter[, 1:15]
    taxon <- metadata_taxon_filter[, 16:ncol(metadata_taxon_filter)]
    
    source("D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\statistic\\01preprocess\\01feature_prepare.R")
    if(feat_type %in% c('A', 'B', 'F', 'V', 'ABFV')){
        taxon_rel = merge_feats(taxon, feat_type, min_abundance=1e-4, min_prevalence=0.1)
    }else if(feat_type %in% c('KOs', 'pathways', 'GMMs', 'GBMs')){
        taxon = round(taxon)
        taxon_rel = merge_feats(taxon, feat_type, min_abundance=1e-6, min_prevalence=0.1)
    }
    
    #taxon_rel = merge_feats(taxon, feat_type, min_abundance=1e-4, min_prevalence=0.1)
    rowSums(taxon_rel)
    taxon_filter = taxon[, colnames(taxon) %in% colnames(taxon_rel)]
    taxon_filter[is.na(taxon_filter)] <- 0
    taxon_filter <- taxon_filter[rowSums(taxon_filter)>0, ]
    #any(is.na(taxon_filter))
    metaAD_abd = t(taxon_filter)
    metaAD_meta = metadata[metadata$Sample %in% row.names(taxon_filter), ]
    print(dim(metaAD_abd))
    print(dim(metaAD_meta))
    ####Step1: 批次矫正。
    ## feature_abd表的行为特征(species, KOs, pathways), 列为样本名称。
    metaAD_abd_adj <- MMUPHin_Batch(metaAD_meta, metaAD_abd, outdir, group, feat_type, Batch)
    
    # ###中心对数化方法参考: https://malucalle.github.io/Microbiome-Variable-Selection/clr.html
    # ### margin=1表示对行操作，margin=2表示对列操作。
    # ### 报错: object ‘lvl_Batch’ not found ==> https://forum.biobakery.org/t/error-in-paste-lvl-Batch-ind-exposure-ind-exposure-cat-collapse-object-lvl-Batch-not-found/2561
    Maaslin2_dir = file.path(outdir, paste0("MaAsLin2_", feat_type, "_", group, "vsNC"))
    dir.create(Maaslin2_dir, showWarnings = FALSE, recursive = TRUE)
    outfile = file.path(outdir, paste0("MaAsLin2_", feat_type, "_", group, "vsNC.csv"))
    if(!file.exists(outfile)){
        fit_data <- Maaslin2(
            input_data=metaAD_abd_adj,
            input_metadata=metaAD_meta,
            output=Maaslin2_dir,
            #transform = "AST",
            fixed_effects = c("Group", "Batch", "Age", "BMI", "host_removed_count"),
            reference = c("Batch,CHN2", "Group,NC"),
            cores = 4,
            standardize = FALSE,
            plot_heatmap = FALSE,
            plot_scatter = FALSE,
            heatmap_first_n = 20)

        results = fit_data$results[fit_data$results$name %in% c("GroupSCS", "GroupSCD", "GroupMCI", "GroupAD"), ]
        write.csv(results, outfile, row.names=FALSE)
    }else{
        results = read.csv(outfile, header=TRUE)
    }

    ## https://forum.biobakery.org/t/q-value-threshold-for-significance/4520
    num_SCS = nrow(results[(results$value=="SCS") & (results$qval<0.25), ])
    num_SCD = nrow(results[(results$value=="SCD") & (results$qval<0.25), ])
    num_MCI = nrow(results[(results$value=="MCI") & (results$qval<0.25), ])
    num_AD = nrow(results[(results$value=="AD") & (results$qval<0.25), ])
    print(paste0("Diff ", feat_type, ": SCS=", num_SCS, ", SCD=", num_SCD, ", MCI=", num_MCI, ", AD=", num_AD))
    write.csv(results, outfile, row.names=FALSE)
}

theme_set(theme_bw())

#inputdir="D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\02profile"
#rm(list = ls())
workdir="D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom"
# source("D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\03MMUPHin_diff\\01MMUPHin_lm_meta.R")
# environment(MMUPHin_lm_meta) <- asNamespace('MMUPHin')
# assignInNamespace("lm_meta", MMUPHin_lm_meta, ns = "MMUPHin")
# options(warn=-1)


# taxon_types = c("KOs", "pathways", "GMMs", "GBMs")
# taxon_types = c("species", "KOs", "pathways", "GMMs", "GBMs")
taxon_types = c("species")
Batch = "CHN+CHN2"
group = "all"
for (taxon_type in taxon_types){
    #taxon_type = "species"
    
    metadata_taxon <- data.table::fread(file.path(workdir, "00profile", paste0('metadata_', taxon_type,'.csv')), sep=",")
    metadata_taxon <- as.data.frame(metadata_taxon)
    metadata_taxon <- metadata_taxon[metadata_taxon$Batch %in% c("CHN", "CHN2"), ]
    
    #metadata_taxon <- metadata_taxon[metadata_taxon$age>=60, ]
    print(dim(metadata_taxon))
    print(unique(metadata_taxon$Group))
    print(colnames(metadata_taxon)[1:15])
    metadata_taxon = metadata_taxon[metadata_taxon$Group %in% c('NC', 'SCS', 'SCD', 'MCI', 'AD'), ]
    metadata_taxon$Group <- factor(metadata_taxon$Group, levels = c('NC', 'SCS', 'SCD', 'MCI', 'AD'))
    print(unique(metadata_taxon$Batch))
    metadata_taxon$Batch <- factor(metadata_taxon$Batch, levels = c("CHN","CHN2"))
    
    rownames(metadata_taxon) <- metadata_taxon$Sample
    
    #str(metadata_taxon[, 1:15])
    
    #colnames(metadata_taxon)[which(names(metadata_taxon) == "genDEU")] <- "gender"
    metadata_taxon$BMI <- as.numeric(metadata_taxon$BMI)
    # metadata_taxon$BMI_num <- metadata_taxon$BMI
    # metadata_taxon$BMI_num <- as.numeric(metadata_taxon$BMI_num)
    # metadata_taxon$BMI[metadata_taxon$BMI_num<25.0 ] = "lean"
    # metadata_taxon$BMI[(metadata_taxon$BMI_num>=25.0) & (metadata_taxon$BMI_num<=30.0)] = "obese"
    # metadata_taxon$BMI[metadata_taxon$BMI_num>30.0 ] = "overweight"
    # metadata_taxon$BMI_num <- NULL
    # metadata_taxon$BMI = factor(metadata_taxon$BMI, levels=c("lean", "obese", "overweight"))
    # 
    # metadata_taxon$Age_num <- metadata_taxon$Age
    # metadata_taxon$Age_num <- as.numeric(metadata_taxon$Age_num)
    # metadata_taxon$Age[metadata_taxon$Age_num<65.0 ] = "<65"
    # metadata_taxon$Age[(metadata_taxon$Age_num>=65.0) & (metadata_taxon$Age_num<70.0)] = "65-70"
    # metadata_taxon$Age[(metadata_taxon$Age_num>=70.0) & (metadata_taxon$Age_num<75.0)] = "70-75"
    # metadata_taxon$Age[(metadata_taxon$Age_num>=75.0) & (metadata_taxon$Age_num<80.0)] = "75-80"
    # metadata_taxon$Age[metadata_taxon$Age_num>=80.0 ] = ">=80"
    # metadata_taxon$Age_num <- NULL
    # metadata_taxon$Age = factor(metadata_taxon$Age, levels=c("<65", "65-70", "70-75", "75-80", ">80"))
    
    
    metadata_taxon$Age <- as.numeric(metadata_taxon$Age)
    metadata_taxon$host_removed_count <- as.numeric(metadata_taxon$host_removed_count)
    metadata_taxon$ND_score <- as.numeric(metadata_taxon$ND_score)
    
    outdir = file.path(workdir, "03MMUPHin_diff", taxon_type)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
    imglist = list()
    if(taxon_type == "species"){
        all_A = MMUPHin_diff(workdir, taxon_type, metadata_taxon, Batch, group,  'A', outdir)
        all_B = MMUPHin_diff(workdir, taxon_type, metadata_taxon, Batch, group, 'B', outdir)
        all_F = MMUPHin_diff(workdir, taxon_type, metadata_taxon, Batch, group, 'F', outdir)
        all_V = MMUPHin_diff(workdir, taxon_type, metadata_taxon, Batch, group, 'V', outdir)
        #all_ABFV = MMUPHin_diff(workdir, taxon_type, metadata_taxon, Batch, group, 'ABFV', outdir)
    }else if(taxon_type %in% c("KOs", "pathways", "GMMs", "GBMs")){
        diff_features = MMUPHin_diff(workdir, taxon_type, metadata_taxon, Batch, group,  taxon_type, outdir)
        imglist[[group]] = diff_features
    }
}


