library(data.table)
library(ggplot2)

rm(list=ls())

cor_method = "spearman"
profile_dir = "/home1/jialh/brain/01meta/multikingdom/00profile"
species_abundance_df = as.data.frame(fread(file.path(profile_dir, "metadata_species.csv"), sep=","))
row.names(species_abundance_df) = species_abundance_df$Sample
species_abundance_df = species_abundance_df[species_abundance_df$Batch %in% c("CHN", "CHN2"), ]
print(species_abundance_df[1:5, 1:20])

KOs_abundance_df = as.data.frame(fread(file.path(profile_dir, "metadata_KOs.csv"), sep=","))
row.names(KOs_abundance_df) = KOs_abundance_df$Sample
KOs_abundance_df = KOs_abundance_df[KOs_abundance_df$Batch %in% c("CHN", "CHN2"), ]
print(KOs_abundance_df[1:5, 1:20])

#stage = "AD"
#species_abundance_df_filter = species_abundance_df[species_abundance_df$Group %in% c("NC", stage), ]
species_abundance_df_filter = species_abundance_df
species_abundance = species_abundance_df[, 16:ncol(species_abundance_df_filter)]
source("/home1/jialh/brain/01meta/multikingdom/statistic/01preprocess/01feature_prepare.R")
species_abundance_filter = merge_feats(taxon=species_abundance, feat_type="ABFV", min_abundance=1e-4, min_prevalence=0.1)

##
diff_species_df = read.csv("/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff/graphlan/MMUPHin_diff_merge.csv", header=1)
species_abundance_filter = species_abundance_filter[, colnames(species_abundance_filter) %in% c(diff_species_df$species)]




#KOs_abundance_df_filter = KOs_abundance_df[KOs_abundance_df$Group %in% c("NC", stage), ]
KOs_abundance_df_filter = KOs_abundance_df
KOs_abundance = KOs_abundance_df_filter[, 16:ncol(KOs_abundance_df_filter)]
source("/home1/jialh/brain/01meta/multikingdom/statistic/01preprocess/01feature_prepare.R")
KOs_abundance_filter = merge_feats(taxon=KOs_abundance, feat_type="KOs", min_abundance=1e-6, min_prevalence=0.1)

RF_important_KOs = NULL
for (stage in c("SCS", "SCD", "MCI", "AD")){
    #stage = "MCI"
    RFdir = "/home1/jialh/brain/01meta/multikingdom/06classification/CHN+CHN2_KOs+ConQuR/RandomForest"
    RF_KOs_stage_file = file.path(RFdir, paste0("CV_CHN+CHN2_",stage, "_RandomForest_KOs_RandomForest+RFECV_None_feature_importance.csv"))
    RF_KOs_stage_df = read.csv(RF_KOs_stage_file, header=1, row.names=1)
    RF_KOs_stage_df$median = apply(RF_KOs_stage_df, 1, median, na.rm=T)
    RF_KOs_stage_df_filter = data.frame(RF_KOs_stage_df$median)
    row.names(RF_KOs_stage_df_filter) = row.names(RF_KOs_stage_df)
    colnames(RF_KOs_stage_df_filter) = c(stage)
    if (is.null(RF_important_KOs)){
        RF_important_KOs = RF_KOs_stage_df_filter
    }else{
        RF_important_KOs = merge(RF_important_KOs, RF_KOs_stage_df_filter, by="row.names", all=TRUE)
        row.names(RF_important_KOs) = RF_important_KOs$Row.names
        RF_important_KOs = RF_important_KOs[, -1]
    }
    #print(RF_important_KOs)
}

###<----------------------------------------------->###
print(RF_important_KOs)
RF_important_KOs_filter = RF_important_KOs #[rowSums(is.na(RF_important_KOs))<=3, ]
#RF_important_KOs_filter = RF_important_KOs
KOs_abundance_important = KOs_abundance_filter[, colnames(KOs_abundance_filter) %in% row.names(RF_important_KOs_filter)]

##存储KOs和代谢物(包括极性代谢物和分子脂质)之间的关联。
##use="pairwise.complete.obs"通过一种不太合理的方式来处理缺失值。https://bwlewis.github.io/covar/missing.html
##当指定时，R使用通过在成对的基础上省略缺失值的行形成的向量来计算每对列的相关性。因此，每个列向量可能会根据其配对而变化，从而导致相关性值甚至无法比较。
##cor(x, y = NULL, use = "everything", method = c("pearson", "kendall", "spearman"))
##cor{stats} 用于计算两个向量的相关性(correlation)或者两个矩阵各个列的相关矩阵(correlations)。
KOs_species_cor = matrix (NA, nrow = ncol (KOs_abundance_important), ncol = ncol (species_abundance_filter))
rownames (KOs_species_cor) = colnames (KOs_abundance_important)
colnames (KOs_species_cor) = colnames (species_abundance_filter)

# species_abundance_filter = as.numeric(species_abundance_filter)
# KOs_abundance_important = as.numeric(KOs_abundance_important)

for (species in colnames (species_abundance_filter)) {
    KOs_species_cor [ , species] = apply (KOs_abundance_important, MARGIN = 2, FUN = function (x)
        cor (x, species_abundance_filter [, species],
             method = cor_method, use = "pairwise.complete.obs"))
}


####<---------------------------->####
###
diff_KOs = read.csv("/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff/KOs/MMUPHin_diff_merge.csv", header=1, row.names=1)
diff_species = read.csv("/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff/graphlan/MMUPHin_diff_merge.csv", header=1, row.names=1)
diff_KOs_filter = diff_KOs[rowSums(is.na(diff_KOs))<=2, ]
diff_species_filter = diff_species  #[rowSums(is.na(diff_species))<=2, ]
print(KOs_species_cor[1:5, 1:3])
KOs_species_cor = KOs_species_cor[rownames(KOs_species_cor) %in% rownames(diff_KOs_filter), ]
print(colnames(KOs_species_cor)[1:5])
print(rownames(diff_species_filter)[1:5])
KOs_species_cor = KOs_species_cor[, colnames(KOs_species_cor) %in% rownames(diff_species_filter)]

KOs_species_cor_filter =  KOs_species_cor[, colMeans(abs(KOs_species_cor)>0.4)>=0.1]
colnames(KOs_species_cor_filter) = gsub("^.*s__", "", colnames(KOs_species_cor_filter))


###
KOs_to_pathways = read.csv(file.path(profile_dir, "KEGG_KOs_to_pathways_metabolism.csv"), row.name=1, header=1)
KOs_to_pathways$KO_desc[KOs_to_pathways$KO_desc=="K00260 (Alanine, aspartate and glutamate metabolism)"] <- "K00260 (Ala, Asp and Glu metabolism)"
KOs_to_pathways$KO_desc[KOs_to_pathways$KO_desc=="K06208 (Phenylalanine, tyrosine and tryptophan biosynthesis)"] <- "K06208 (Phe, Tyr and  Trp biosynthesis)"
# KOs_to_pathways$KO_desc[KOs_to_pathways$KO_desc=="K00260 (Alanine, aspartate and glutamate metabolism)"] <- "K00260 (Ala, Asp and Glu metabolism)"
# KOs_to_pathways$KO_desc[KOs_to_pathways$KO_desc=="K00260 (Alanine, aspartate and glutamate metabolism)"] <- "K00260 (Ala, Asp and Glu metabolism)"


print(colnames(KOs_to_pathways))
KOs_to_pathways_filter = data.frame(KOs_to_pathways[, "KO_desc"])
rownames(KOs_to_pathways_filter) = rownames(KOs_to_pathways)
colnames(KOs_to_pathways_filter) = c("KO_desc")

KOs_species_cor_heatmap = merge(KOs_species_cor_filter, KOs_to_pathways_filter, by="row.names")
rownames(KOs_species_cor_heatmap) = KOs_species_cor_heatmap$KO_desc
KOs_species_cor_heatmap = KOs_species_cor_heatmap[, !(colnames(KOs_species_cor_heatmap) %in% c("Row.names", "KO_desc"))]

#KOs_species_cor_heatmap = KOs_species_cor_heatmap[, colSums(KOs_species_cor_heatmap)]

###
library(ComplexHeatmap)
library(circlize)
marker_df = data.frame(data=array(NA, dim(KOs_species_cor_filter)))
rownames(marker_df) = rownames(KOs_species_cor_filter)
colnames(marker_df) = colnames(KOs_species_cor_filter)

for(i in rownames(KOs_species_cor_filter)){
    for (j in colnames(KOs_species_cor_filter)){
        if(KOs_species_cor_filter[i,j]>0.3){
            marker_df[i, j] = "+"
        }else if(KOs_species_cor_filter[i,j]< -0.3){
            marker_df[i, j] = "-"
        }else{
            marker_df[i, j] = ""
        }
    }
}

print(marker_df[1:5, 1:5])

marker_df = merge(marker_df, KOs_to_pathways_filter, by="row.names")
rownames(marker_df) = marker_df$KO_desc
marker_df = marker_df[, !(colnames(marker_df) %in% c("Row.names", "KO_desc"))]



# color mapping for -log10(pvalue)
###<-------------------行名的相关注释---------------->###
diff_KOs = read.csv("/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff/KOs/MMUPHin_diff_merge.csv", header=1, row.names=1)
diff_KOs_annotation = merge(KOs_species_cor_filter, diff_KOs, by="row.names", all.x=TRUE)
rownames(diff_KOs_annotation) = diff_KOs_annotation$Row.names
diff_KOs_annotation = diff_KOs_annotation[, c("SCS", "SCD", "MCI", "AD")]
diff_KOs_annotation[diff_KOs_annotation>0] = 1
diff_KOs_annotation[diff_KOs_annotation<0] = -1
diff_KOs_annotation[is.na(diff_KOs_annotation)] = 0
diff_KOs_annotation = merge(diff_KOs_annotation, KOs_to_pathways_filter, by="row.names")
rownames(diff_KOs_annotation) = diff_KOs_annotation$features
diff_KOs_annotation = diff_KOs_annotation[, !(colnames(diff_KOs_annotation) %in% c("Row.names", "features"))]

KOs_col_fun = colorRamp2(c(-1, 0, 1), c("#36bccb", "white", "#ff6836"))
ha = rowAnnotation(
    SCS = anno_simple(diff_KOs_annotation$SCS, col = KOs_col_fun, gp = gpar(col = "black", lwd = 0.5)),
    SCD = anno_simple(diff_KOs_annotation$SCD, col = KOs_col_fun, gp = gpar(col = "black", lwd = 0.5)),
    MCI = anno_simple(diff_KOs_annotation$MCI, col = KOs_col_fun, gp = gpar(col = "black", lwd = 0.5)),
    AD = anno_simple(diff_KOs_annotation$AD, col = KOs_col_fun, gp = gpar(col = "black", lwd = 0.5)),
    annotation_name_side = "top", gp = gpar(col = "black", lwd = 0.5))
draw(ha)


###<-----------------------列名的相关注释------------------------>###
diff_species = read.csv("/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff/graphlan/MMUPHin_diff_merge.csv", header=1, row.names=1)
diff_species$species = gsub("^.*s__", "", rownames(diff_species))
diff_species$phylum = gsub("^.*p__", "", rownames(diff_species))
diff_species$phylum = gsub("c__.*$", "", diff_species$phylum)
diff_species$phylum = gsub("|", "", diff_species$phylum, fixed=TRUE)
diff_species_annotation = merge(t(KOs_species_cor_filter), diff_species, by.x="row.names", by.y="species", all.x=TRUE)
#####
rownames(diff_species_annotation) = diff_species_annotation$Row.names
diff_species_annotation = diff_species_annotation[, c("SCS", "SCD", "MCI", "AD")]
diff_species_annotation[diff_species_annotation>0] = 1
diff_species_annotation[diff_species_annotation<0] = -1
diff_species_annotation[is.na(diff_species_annotation)] = 0
# diff_species_annotation = merge(diff_species_annotation, KOs_to_pathways_filter, by="row.names")
# rownames(diff_species_annotation) = diff_species_annotation$features
# diff_species_annotation = diff_species_annotation[, !(colnames(diff_species_annotation) %in% c("Row.names", "features"))]

species_col_fun = colorRamp2(c(-1, 0, 1), c("#36bccb", "white", "#ff6836"))
species_ha = HeatmapAnnotation(
    SCS = anno_simple(diff_species_annotation$SCS, col = species_col_fun, gp = gpar(col = "black", lwd = 0.5)),
    SCD = anno_simple(diff_species_annotation$SCD, col = species_col_fun, gp = gpar(col = "black", lwd = 0.5)),
    MCI = anno_simple(diff_species_annotation$MCI, col = species_col_fun, gp = gpar(col = "black", lwd = 0.5)),
    AD = anno_simple(diff_species_annotation$AD, col = species_col_fun, gp = gpar(col = "black", lwd = 0.5)),
    annotation_name_side = "right", gp = gpar(col = "black", lwd = 0.5))
draw(species_ha)

###16x24的矩阵。
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("#4575b4", "#e9f6e6", "#d73027"))
print(rownames(KOs_species_cor_heatmap)[1:5])
print(colnames(KOs_species_cor_heatmap)[1:5])

print(rownames(marker_df)[1:5])
print(colnames(marker_df)[1:5])
#KOs_species_cor_heatmap[abs(KOs_species_cor_heatmap) < 0.1] = 0
marker_df = marker_df[rownames(KOs_species_cor_heatmap), ]
markder_df = marker_df[, colnames(KOs_species_cor_heatmap)]


img <- Heatmap(KOs_species_cor_heatmap, col = col_fun, 
        cluster_rows = TRUE, cluster_columns = TRUE, 
        rect_gp = gpar(col = "white", lwd = 1),
        left_annotation = ha, bottom_annotation = species_ha, 
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(marker_df[i, j], x, y, gp = gpar(fontsize = 12))},
        #column_split = c('A', 'A', 'A', 'B'),
        row_names_side = "left", column_names_side = "bottom",
        column_title_side = "bottom", 
        column_title = "species", row_title = "KOs",
        column_title_gp = gpar(fontsize=12),
        row_title_gp = gpar(fontsize=12),
        row_dend_gp = gpar(), 
        column_dend_gp = gpar(),
        row_names_gp = gpar(fontsize = 12),
        column_names_rot = 30,
        column_names_gp = gpar(fontsize = 12),
        row_names_max_width = unit(12, "cm"),
        heatmap_legend_param = list(title = "SCC", legend_width = unit(6, "cm"), 
                                    position="top", direction = "horizontal", 
                                    title_position = "lefttop",
                                    title_gp = gpar(fontsize = 12),
                                    labels_gp = gpar(fontsize = 12)))

####<-------------------------------------------------->####
#getwd()
outdir = "/home1/jialh/brain/01meta/multikingdom/03multi_omics"
write.csv(KOs_species_cor_heatmap, file.path(outdir, "KOs_species_cor_heatmap.csv"))
png(file.path(outdir, "KOs_and_species_associations.png"), width=1200*3, height=800*3, res=300)
draw(img, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side="top", show_heatmap_legend=TRUE)
dev.off()


pdf(file.path(outdir, "KOs_and_species_associations.pdf"), width=13, height=8)
draw(img, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side="top", show_heatmap_legend=TRUE)
dev.off()

