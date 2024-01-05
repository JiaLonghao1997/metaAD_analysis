# ##############################################################################
#
##  Alpha and Beta Diversity
#
# ##############################################################################

# Packages
library("labdsv")
library("coin")
library("vegan")
library("yaml")
library("ggpubr")
library("cowplot")
library("tidyverse")
library("data.table")

rm(list=ls())
pipelinedir = "D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\statistic\\Rscripts"
parameters <- yaml.load_file(file.path(pipelinedir, 'parameters.yaml'))
ref.studies <- parameters$ref.studies
Batch.cols <- unlist(parameters$plotting$study.cols)
group.cols <- unlist(parameters$plotting$group.cols)
Batch.shapes <- unlist(parameters$plotting$study.shapes)
start.time <- proc.time()[1]

# args = commandArgs(trailingOnly=TRUE)
# if (length(args)==0) {
#     stop("The analysis tag needs to be provided! Exiting...\n")
# }
# tag <- args[1]
# stage <- args[2]
# feat_type <- args[3]


tag = "species"
stage = "stages"
feat_type = 'ABFV'

workdir = "D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\02EMBLmeta"
#datadir = file.path(workdir, "data", tag, paste0("NCvs", stage))
#inputdir = file.path(workdir, "files", tag, paste0("NCvs", stage))
outdir = file.path(workdir, "02diversity")
#dir.create(inputdir, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# ##############################################################################
# Get Data
# fn.feat <- file.path(datadir, paste0('feat_rar_', feat_type, '.tsv'))
# feat.all <- as.matrix(read.table(fn.feat, sep='\t', header=TRUE, row.names=1, 
#                                  stringsAsFactors = FALSE, 
#                                  check.names = FALSE, quote=''))
# feat.all <- prop.table(feat.all)
# 
# meta <- read_tsv(file.path(datadir, 'meta.tsv'))
inputdir = "D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\00profile"
setwd(workdir)

metadata_features = fread(file.path(inputdir, "metadata_species.csv"), sep=",")
metadata_features = as.data.frame(metadata_features)

rownames(metadata_features) = metadata_features$Sample
print(unique(metadata_features$Batch))
metadata_features_filter = metadata_features[metadata_features$Batch %in% c('CHN', 'CHN2'), ]


###
metadata_features_filter$Batch[metadata_features_filter$Batch=="CHN"] = "CHN1"
metadata = metadata_features_filter[, 1:15]
features = metadata_features_filter[, 16:ncol(metadata_features_filter)]
print(rownames(features)[0:10])
print(unique(metadata_features_filter$Batch))

source("D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\statistic\\01preprocess\\01feature_prepare.R")
features_rel = merge_feats(taxon=features, feat_type="ABFV", min_abundance=1e-4, min_prevalence=0.1)
features_filter = features[, colnames(features) %in% colnames(features_rel)]

###
# data(BCI)  ##A data frame with 50 plots (rows) of 1 hectare with counts of trees on each plot with total of 225 species (columns). 
# S <- specnumber(BCI) # observed number of species
# (raremax <- min(rowSums(BCI)))
# Srare <- rarefy(BCI, raremax)
# plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
# abline(0, 1)
# rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)
(raremax <- min(rowSums(features_filter)))

features_rar_file = file.path(workdir, paste0("rrarefy_", feat_type, ".csv"))
if(file.exists(features_rar_file)){
    features_rar = read.csv(features_rar_file, header=1, row.names=1, stringsAsFactors = FALSE)
}else{
    features_rar = rrarefy(floor(features_filter), sample=raremax)
    
    write.csv(features_rar, file.path(workdir, paste0("rrarefy_", feat_type, ".csv")))
}

###
feat.rar = as.matrix(features_rar)
feat.rel <- prop.table(feat.rar, margin=1)
meta <- metadata

print(rowSums(feat.rel)[1:10])
# ##############################################################################
# Compute Alpha diversity
print(str(meta))

stage = "stages"
##https://github.com/zellerlab/crc_meta/blob/master/src/figure_a_b_diversity.R
plot_alpha_diversity <- function(meta, feat.rar, feat.rel, alpha_method="shannon"){
    #alpha_method = "fisher"
    if (alpha_method %in% c("shannon", "simpson", "invsimpson")){
        ##diversity是vegan中的，行为样本，列为物种名称。
        df.div <- meta %>% 
            mutate(`All species` = vegan::diversity(feat.rel, index = alpha_method)) %>%
            mutate(`Archaea` = vegan::diversity(feat.rel[, grep(pattern = 'k__Archaea', x=colnames(feat.rel), value = TRUE)], index = alpha_method)) %>%
            mutate(`Bacteria` = vegan::diversity(feat.rel[, grep(pattern = 'k__Bacteria',
                                                          x=colnames(feat.rel), value = TRUE)], index = alpha_method)) %>%
            mutate(`Fungi` = vegan::diversity(feat.rel[, grep(pattern = 'k__Eukaryota.k__Fungi',
                                                       x=colnames(feat.rel), value = TRUE)], index = alpha_method)) %>%
            mutate(`Viruses` = vegan::diversity(feat.rel[, grep(pattern = 'k__Viruses', x=colnames(feat.rel), value = TRUE)], index = alpha_method)) %>%
            
            select(Group, Batch, `All species`, `Archaea`, `Bacteria`, `Fungi`, `Viruses`) %>% 
            gather(key=type, value=diversity, -c(Batch, Group)) %>% 
            mutate(Batch=str_remove(Batch, pattern=paste0('-', stage))) %>%
            mutate(Batch=factor(Batch, levels = str_remove(ref.studies, pattern=paste0('-', stage)))) %>%
            mutate(Group=factor(Group, levels=c('NC', 'SCS', 'SCD', 'MCI', 'AD'))) %>% 
            mutate(type=factor(type, levels=c('All species', 
                                              'Archaea', 'Bacteria', 'Fungi', 'Viruses')))
        
         
    }else if(alpha_method=="fisher"){
        #feat.rar <- floor(feat.rar)
        df.div <- meta %>% 
            mutate(`All species` = fisher.alpha(feat.rar)) %>%
            mutate(`Archaea` = fisher.alpha(feat.rar[, grep(pattern = 'k__Archaea', x=colnames(feat.rar), value = TRUE)])) %>%
            mutate(`Bacteria` = fisher.alpha(feat.rar[, grep(pattern = 'k__Bacteria',
                                                          x=colnames(feat.rar), value = TRUE)])) %>%
            mutate(`Fungi` = fisher.alpha(feat.rar[, grep(pattern = 'k__Eukaryota.k__Fungi',
                                                       x=colnames(feat.rar), value = TRUE)])) %>%
            mutate(`Viruses` = fisher.alpha(feat.rar[, grep(pattern = 'k__Viruses', x=colnames(feat.rar), value = TRUE)])) %>%
            
            select(Group, Batch, `All species`, `Archaea`, `Bacteria`, `Fungi`, `Viruses`) %>% 
            gather(key=type, value=diversity, -c(Batch, Group)) %>% 
            mutate(Batch=str_remove(Batch, pattern=paste0('-', stage))) %>%
            mutate(Batch=factor(Batch, levels = str_remove(ref.studies, pattern=paste0('-', stage)))) %>%
            mutate(Group=factor(Group, levels=c('NC', 'SCS', 'SCD', 'MCI', 'AD'))) %>% 
            mutate(type=factor(type, levels=c('All species', 
                                              'Archaea', 'Bacteria', 'Fungi', 'Viruses')))
    }
    
    #blocked wilcoxon
    # df.div %>%
    #   filter(type=='All species') %>%
    #   wilcox_test(diversity~Group|Batch, data=.)
    # df.div %>%
    #   filter(type=='Archaea') %>%
    #   wilcox_test(diversity~Group|Batch, data=.)
    # df.div %>%
    #   filter(type=='Bacteria') %>%
    #   wilcox_test(diversity~Group|Batch, data=.)
    # df.div %>%
    #     filter(type=='Fungi') %>%
    #     wilcox_test(diversity~Group|Batch, data=.)
    # df.div %>%
    #     filter(type=='Viruses') %>%
    #     wilcox_test(diversity~Group|Batch, data=.)
    
    
    # anova
    summary(aov(rank~Group+Batch, data=df.div %>% 
                    filter(type=='All species') %>% 
                    mutate(rank=rank(diversity))))
    summary(aov(rank~Group+Batch, data=df.div %>% 
                    filter(type=='Archaea') %>% 
                    mutate(rank=rank(diversity))))
    summary(aov(rank~Group+Batch, data=df.div %>%  
                    filter(type=='Bacteria') %>% 
                    mutate(rank=rank(diversity))))
    summary(aov(rank~Group+Batch, data=df.div %>%  
                    filter(type=='Fungi') %>%  
                    mutate(rank=rank(diversity))))
    summary(aov(rank~Group+Batch, data=df.div %>% 
                    filter(type=='Viruses') %>% 
                    mutate(rank=rank(diversity))))
    # ##############################################################################
    # plot
    #print(group.cols)
    
    if(stage=="stages"){
        my_comparisons <- list( c("NC", "SCS"), c("NC", "SCD"), c("NC", "MCI"), c("NC", "AD"))
    }else{
        my_comparisons <- list( c("NC", stage))
    }
    
    #print(df.div$Group)
    library(stringr)
    max_diversity <- max(df.div$diversity)
    g <- df.div %>% 
      ggplot(aes(x=Group, fill=Group, y=diversity)) +
        geom_boxplot() +
        facet_wrap(~type, ncol=5, scales = 'fixed') + 
        theme_bw() + 
        ylab(paste0(str_to_title(alpha_method), " diversity")) +
        xlab('') +
        theme(
            axis.text.x = element_text(size=12, hjust = 1, angle =30),
            axis.text.y = element_text(size=12, hjust = 0.5),
            axis.title.x = element_text(size=12, hjust = 0.5),
            axis.title.y = element_text(size=12, hjust = 0.5),
            plot.title = element_text(size=14, hjust=0.5),
            strip.text = element_text(size=12, hjust = 0.5)) + 
        scale_fill_manual(values=group.cols) + 
        stat_compare_means(label = "p.format", method = "wilcox.test", 
                           comparisons = my_comparisons, size=2)
    
    print(g)
    return(g)
}

####<------------------------------------->###
Shannon_img <- plot_alpha_diversity(metadata, feat.rar, feat.rel, alpha_method="shannon")
Simpson_img <- plot_alpha_diversity(metadata, feat.rar, feat.rel, alpha_method="simpson")
Fisher_img  <- plot_alpha_diversity(metadata, feat.rar, feat.rel, alpha_method="fisher")

pdf(file.path(outdir, paste0('alpha_diversity_', feat_type, '.pdf')), width = 12, height = 8)
merge_img = plot_grid(Shannon_img, Simpson_img, Fisher_img, ncol = 1, nrow=3, align="hv", labels=letters[1:3])
print(merge_img)
dev.off()


cat('Successfully plotted alpha diversity in',
    proc.time()[1]-start.time, 'second...\n')



# ##############################################################################
# Beta Diversity
# species
#feat.all <- rrarefy(features_filter, sample=min(rowSums(features_filter)))
feat.all <- feature_filter

# ##############################################################################
# compute PCoA
dist = vegdist(feat.all, method = 'bray')
pco.results = pco(dist, k=2)

##通过置换多因素方差分析，检查Batch和group分别在beta多样性中p-value差异。
##https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permanova/
##置换因素方差分析顺序：https://blog.csdn.net/qazplm12_3/article/details/120520772
print(dim(feat.all))
print(dim(metadata))
print(colnames(metadata))
colnames(metadata)[colnames(metadata)=="host_removed_count"] <- "ReadCount"
inputdir = "D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\02EMBLmeta\\02diversity"
print(inputdir)
permanova_file = file.path(inputdir, paste0("permanova_beta_diversity_", feat_type, ".Rdata"))
# if(file.exists(permanova_file)){
#     load(permanova_file)
# }else{
#     results = adonis2(formula = feat.all ~ Batch + Group + Gender + ReadCount, 
#                       data = metadata, method = 'bray', permutations = 999, by = "margin")
#     save(results, file=permanova_file)
# }

results = adonis2(formula = feat.all ~ Batch + Group + Gender + ReadCount,
                  data = metadata, method = 'bray', permutations = 999, by = "margin")

print(results)
print(results[rownames(results)=="Batch"]$`Pr(>F)`)

Batch.pvalue = signif(results$`Pr(>F)`[1], 2)
group.pvalue = signif(results$`Pr(>F)`[2], 2)
gender.pvale = signif(results$`Pr(>F)`[3], 2)

axis.1.title <- paste('PCoA1 [', 
                      round((pco.results$eig[1]/sum(pco.results$eig))*100,1),
                      '%]', sep='')
axis.2.title <- paste('PCoA2 [', 
                      round((pco.results$eig[2]/sum(pco.results$eig))*100,1),
                      '%]', sep='')

df.plot <- tibble(Axis1 = -1*pco.results$points[,1],
                  Axis2 = pco.results$points[,2],
                  Sample_ID = rownames(pco.results$points),
                  Group=metadata$Group,
                  Batch=metadata$Batch)

##
aixs1.group.pvalue <- signif(kruskal.test(Axis1 ~ Group, data = df.plot)$p.value, 2)
aixs1.Batch.pvalue <- signif(kruskal.test(Axis1 ~ Batch, data = df.plot)$p.value, 2)
aixs2.group.pvalue <- signif(kruskal.test(Axis2 ~ Group, data = df.plot)$p.value, 2)
aixs2.Batch.pvalue <- signif(kruskal.test(Axis2 ~ Batch, data = df.plot)$p.value, 2)

print(min(df.plot$Axis1))
print(min(df.plot$Axis2))
# ##############################################################################
# subplots
print(df.plot$Group)
df.plot$Group <- factor(df.plot$Group, levels=c("NC", "SCS", "SCD", "MCI", "AD"))

# main plot
g.main <- df.plot %>% 
  ggplot(aes(x=Axis1, y=Axis2, shape=Batch, col=Group)) +
  geom_point() + 
  scale_colour_manual(values=group.cols, guide="none") + 
  scale_shape_manual(values=Batch.shapes) + 
  scale_x_continuous(position='top') +
  xlab(axis.1.title) + ylab(axis.2.title) +
  annotate("text", x = min(df.plot$Axis1), y = min(df.plot$Axis1), hjust = 0, vjust = 1, 
           label= paste0("Batch: ", 'P', "=", Batch.pvalue, 
           "\nGroup: ", 'P', "=", group.pvalue)) + 
  theme(panel.background = element_rect(fill='white', color = 'black'),
        legend.background = element_blank(),
        axis.text = element_blank(), axis.ticks=element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank(),
        panel.grid = element_blank())
print(g.main)
# dev.off()

# Batch boxplot axis 1
g.s.1 <- df.plot %>% 
    mutate(Batch=factor(Batch, levels=names(Batch.cols))) %>% 
    ggplot(aes(y=Axis1, x=Batch, fill=Batch)) + 
    geom_boxplot() + 
    scale_fill_manual(values=Batch.cols) +
    xlab(paste0("Batch\n", 'P', "=",  aixs1.Batch.pvalue)) +
    theme(axis.ticks = element_blank(),
        panel.background = element_rect(fill='white', color = 'black'),
        axis.text = element_blank(), axis.title.x = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank(),
        panel.grid = element_blank()) + 
    coord_flip()

# Batch boxplot axis 2
g.s.2 <- df.plot %>% 
  mutate(Batch=factor(Batch, levels=names(Batch.cols))) %>% 
  ggplot(aes(y=Axis2, x=Batch, fill=Batch)) + 
  geom_boxplot() + 
  scale_fill_manual(values=Batch.cols, guide="none") +
  xlab(paste0("Batch\n", 'P', "=",  aixs2.Batch.pvalue)) +
  scale_x_discrete(position='top') +
  theme(axis.ticks=element_blank(), 
        panel.background = element_rect(fill='white', color = 'black'),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank())

# group plot axis1
g.g.1 <- df.plot %>% 
  ggplot(aes(x=Group, y=Axis1, fill=Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.cols) + 
  xlab(paste0("Group\n", 'P', "=",  aixs1.group.pvalue)) +
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank(),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank()) + 
  coord_flip()
# group plot axis2
g.g.2 <- df.plot %>% 
  ggplot(aes(x=Group, y=Axis2, fill=Group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.cols, guide="none") + 
  scale_x_discrete(position='top') + 
  scale_y_continuous(position = 'right') +
  xlab(paste0("Group\n", 'P', "=",  aixs2.group.pvalue)) +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank())

##Extract Legends from a ggplot object: https://rpkgs.datanovia.com/ggpubr/reference/get_legend.html
#g.main.legend <- get_legend(g.main)
pdf(file.path(outdir, paste0('beta_diversity_legend_', feat_type, '.pdf')), width=6, height=3)
Batch.shape.legend <- as_ggplot(get_legend(g.main))
Batch.color.legend <- as_ggplot(get_legend(g.s.1))
group.color.legend <- as_ggplot(get_legend(g.g.1))
merge_legend1 <- plot_grid(Batch.color.legend, Batch.shape.legend, nrow=2, rel_heights = c(0.2, 0.3))
merge_legend2 <- plot_grid(merge_legend1, group.color.legend, ncol=2, rel_widths = c(0.2, 0.3)) 
print(merge_legend2)
dev.off()

# ##############################################################################
# Plot everyting together
pdf(file.path(outdir, paste0('beta_diversity_', feat_type, '.pdf')), width=6.5, height=6.5)
g.main <- g.main + theme(legend.position="none")
g.s.1 <- g.s.1 + theme(legend.position="none")
g.g.1 <- g.g.1 + theme(legend.position="none")
# merge_img1 <- plot_grid(g.main, g.s.2, g.g.2, g.s.1, NULL, NULL, g.g.1, NULL, NULL,
#           nrow=3,
#           rel_widths = c(0.8, 0.2, 0.3), rel_heights = c(0.8, 0.2, 0.3))
merge1 <- plot_grid(g.main, g.s.2, g.g.2, nrow=1,  rel_widths = c(1.6, 0.2, 0.5))
merge2 <- plot_grid(g.s.1, g.g.1, ncol=1,  rel_heights=c(0.2, 0.5))
merge3 <- plot_grid(merge2, merge_legend2, ncol=2,  rel_widths=c(1.6, 0.7))
merge4 <- plot_grid(merge1, merge3, ncol=1,  rel_heights=c(1.6, 0.7))
print(merge4)
dev.off()

cat('Successfully plotted beta diversity in',
    proc.time()[1]-start.time, 'second...\n')

# ######################################################################################################
# End of script
# ######################################################################################################
# 比较:
inputdir = "D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\00profile"
setwd(workdir)
lineage="phylum"
metadata_features = fread(file.path(inputdir, paste0("metadata_", lineage, ".csv")), sep=",")
metadata_features = as.data.frame(metadata_features)
rownames(metadata_features) = metadata_features$Sample
print(unique(metadata_features$Batch))
metadata_features_filter = metadata_features[metadata_features$Batch %in% c('CHN', 'CHN2'), ]
###
metadata_features_filter$Batch[metadata_features_filter$Batch=="CHN"] = "CHN1"
metadata = metadata_features_filter[, 1:15]
features = metadata_features_filter[, 16:ncol(metadata_features_filter)]
print(rownames(features)[0:10])
print(unique(metadata_features_filter$Batch))

source("D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\statistic\\01preprocess\\01feature_prepare.R")
features_rel = merge_feats(taxon=features, feat_type="ABFV", min_abundance=1e-4, min_prevalence=0.1)

metadata_features_rel = merge(metadata, features_rel, by.x="row.names", by.y="row.names")
##The Firmicutes/Bacteroidetes ratio
metadata_features_rel$FB_ratio = metadata_features_rel$`k__Bacteria|p__Firmicutes` / metadata_features_rel$`k__Bacteria|p__Bacteroidota`
print(unique(metadata_features_rel$Group))
metadata_features_rel$Group = factor(metadata_features_rel$Group, levels=c("NC", "SCS", "SCD", "MCI", "AD"))
anova_result <- aov(FB_ratio ~ Group, data = metadata_features_rel)
print(summary(anova_result))
Pvalue <- summary(anova_result)[[1]][["Pr(>F)"]][1]

my_comparisons <- list( c("NC", "SCS"), c("NC", "SCD"), c("NC", "MCI"), c("NC", "AD"))

FB_image <- ggplot(metadata_features_rel, aes(x=Group, y=FB_ratio)) + 
    geom_violin(stat='ydensity',aes(fill=Group), trim = TRUE) +
    geom_boxplot(width=0.3,alpha=0.3)+
    geom_point(size = 0.5, position=position_jitter(width=0.2))+
    stat_compare_means(label = "p.format", method = "wilcox.test", 
                       label.x = c(0.6, 0.7, 0.8, 0.9, 1.0), label.y = c(0.6, 0.7, 0.8, 0.9, 1.0),
                       comparisons = my_comparisons, size=2) + 
    theme_bw() +
    labs(title=paste0('Firmicutes/Bacteroidetes ratio (ANOVA P-value=', signif(Pvalue,2), ")"),  y = 'Ratio') +
    guides(fill=guide_legend(title="Group")) +
    scale_y_continuous(limits = c(0, 6)) + 
    scale_fill_manual(values=group.cols) + 
    geom_hline(aes(yintercept=1), size=0.5, colour="red", linetype="dashed") + 
    theme(axis.text.x = element_text(size=12,angle = 30, hjust = 1),
          axis.text.y = element_text(size=12),
          axis.title.x = element_text(size=14),
          strip.text.x = element_text(size=12),
          axis.title.y = element_text(size=14),
          legend.title = element_text(size=12),
          legend.text = element_text(size=10),
          legend.position = "right", 
          plot.title = element_text(size = 14, hjust=0.5)
    )

print(FB_image)
pdf(file.path(workdir, "FB_image.pdf"), width=8, height=4)
print(FB_image)
dev.off()

# ######################################################################################################
# End of script
# ######################################################################################################

#绘制genus水平的柱状图:
inputdir = "D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\00profile"
setwd(workdir)
lineage="genus"
metadata_features = fread(file.path(inputdir, paste0("metadata_", lineage, ".csv")), sep=",")
metadata_features = as.data.frame(metadata_features)

rownames(metadata_features) = metadata_features$Sample
print(unique(metadata_features$Batch))
metadata_features_filter = metadata_features[metadata_features$Batch %in% c('CHN', 'CHN2'), ]
###
metadata_features_filter$Batch[metadata_features_filter$Batch=="CHN"] = "CHN1"
metadata = metadata_features_filter[, 1:15]
features = metadata_features_filter[, 16:ncol(metadata_features_filter)]
print(rownames(features)[0:10])
print(unique(metadata_features_filter$Batch))

source("D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\statistic\\01preprocess\\01feature_prepare.R")
features_rel = merge_feats(taxon=features, feat_type="ABFV", min_abundance=1e-4, min_prevalence=0.1)
features_filter = features[, colnames(features) %in% colnames(features_rel)]

###
#taxon_types = c("Archaea", "Bacteria", "Fungi", "Viruses")
plot_top_genus <- function(metadata, features_filter, taxon_type, lineage, plot_type="sample"){
    # taxon_type = "Archaea"
    # lineage = "genus"
    # plot_type = "sample"
    if (taxon_type == "Archaea"){
        features_single_kingdom = features_filter[, grep(pattern = 'k__Archaea', x=colnames(features_filter), value = TRUE)]
        viridis_option = "D"
    }else if(taxon_type == "Bacteria"){
        features_single_kingdom = features_filter[, grep(pattern = 'k__Bacteria', x=colnames(features_filter), value = TRUE)]
        viridis_option = "E"
    }else if(taxon_type == "Fungi"){
        features_single_kingdom = features_filter[, grep(pattern = 'k__Eukaryota.k__Fungi', x=colnames(features_filter), value = TRUE)]
        viridis_option = "F"
    }else if (taxon_type=="Viruses"){
        features_single_kingdom = features_filter[, grep(pattern = 'k__Viruses', x=colnames(features_filter), value = TRUE)]
        viridis_option = "G"
    }
    #features_single_kingdom = features_filter[, grep(pattern = 'k__Archaea', x=colnames(features_filter), value = TRUE)]
    features_single_kingdom = as.matrix(features_single_kingdom)
    features_single_kingdom_rel = prop.table(features_single_kingdom, margin=1)
    features_single_kingdom_rel[is.na(features_single_kingdom_rel)] <- 0

    ###
    metadata_feature_rel = merge(metadata, features_single_kingdom_rel, by.x="row.names", by.y="row.names")
    row.names(metadata_feature_rel) <- metadata_feature_rel$Row.names
    metadata_feature_rel <- metadata_feature_rel[, !colnames(metadata_feature_rel) %in% "Row.names"]
    if(plot_type == "group"){
        metadata_feature_mean = aggregate(metadata_feature_rel[, 16:ncol(metadata_feature_rel)], list(metadata_feature_rel$Group), FUN=mean) 
        row.names(metadata_feature_mean) = metadata_feature_mean$Group.1
        metadata_feature_mean = metadata_feature_mean[, !colnames(metadata_feature_mean) %in% c("Group.1")]
        taxon_rel = as.data.frame(t(metadata_feature_mean))
    }else if(plot_type == "sample"){
        taxon_rel = as.data.frame(t(metadata_feature_rel[, 16:ncol(metadata_feature_rel)]))
    }
    
    taxon_rel$sum <- rowSums(taxon_rel)
    taxon_rel <- taxon_rel[order(taxon_rel$sum, decreasing = TRUE), ]
    seleted_num <- min(10, nrow(taxon_rel))
    taxon_top10 <- taxon_rel[1:seleted_num, -ncol(taxon_rel)]
    taxon_top10['others', ] <- 1 - colSums(taxon_top10)
    
    taxon_top10$Genus <- factor(rownames(taxon_top10), levels = rev(rownames(taxon_top10)))
    print(taxon_top10$Genus)
    taxon_top10 <- reshape2::melt(taxon_top10, id = 'Genus')
    colnames(taxon_top10) <- c("Genus", "Sample", "abundance") 
    #taxon_top10$Genus = factor(taxon_top10$Genus)
    if(lineage=="genus"){
        taxon_top10$Genus <- gsub("^.*g__", "", taxon_top10$Genus, perl=TRUE)
    }else if (lineage=="phylum"){
        taxon_top10$Genus <- gsub("^.*p__", "", taxon_top10$Genus, perl=TRUE)
    }
    #taxon_top10$Genus <- gsub("^.*g__", "", taxon_top10$Genus, perl=TRUE)
    Genus_list <- unique(taxon_top10$Genus[order(taxon_top10$abundance, decreasing = T)])
    Genus_list <- c("others", setdiff(Genus_list, 'others'))
    taxon_top10$Genus = factor(taxon_top10$Genus, levels = rev(Genus_list))
    print(unique(taxon_top10$Genus))
    taxon_top10_filter = taxon_top10[taxon_top10$Genus==Genus_list[1], ]
    ###
    sample_list = unique(taxon_top10_filter$Sample[order(taxon_top10_filter$abundance, decreasing = T)])
    #b = c(-1, 0, 1)
    #taxon_list = c("other", setdiff(taxon_list, 'other'))
    taxon_top10$Sample = factor(taxon_top10$Sample, levels = sample_list)
    
    meta_filter = metadata[, colnames(metadata) %in% c("Sample", "Group")]
    meta_taxon_top10 = merge(taxon_top10, meta_filter, by.x="Sample", by.y="Sample")
    
    meta_taxon_top10$Group <- factor(meta_taxon_top10$Group, levels=c("NC", "SCS", "SCD", "MCI", "AD"))
    #taxon_top10$Group <- factor(taxon_top10$Group, levels=c("NC", "SCS", "SCD", "MCI", "AD"))
    #print(taxon_top10)
    
    # taxon_top10$Genus = factor(taxon_top10$Genus)
    # taxon_top10$Genus <- gsub("^.*g__", "", taxon_top10$Genus, perl=TRUE)
    # taxon_top10$Genus = factor(taxon_top10$Genus, levels = rev(unique(taxon_top10$Genus[order(taxon_top10$abundance, decreasing = T)])))
    
    #getPalette = colorRampPalette(brewer.pal(10, "Set3"))
    #colors = getPalette(21)
    
    library("viridis") 
    #print(outdir)
    pdf(file = file.path(outdir, paste0("taxon_barplot_", lineage, "_", taxon_type, ".pdf")),width=6.5, height=4)
    p <- ggplot(meta_taxon_top10, aes(x=Sample, y=abundance, fill = Genus)) +
        geom_col(position = 'stack', width = 1.0) + 
        facet_wrap(~Group, ncol=5, scales = "free_x") + 
        scale_fill_viridis(discrete = TRUE, begin = 0.2, end = 1, option = viridis_option) +
        labs(x = 'Groups', y = 'Relative abundance', title = paste0("Top 10 ", lineage, " in ", taxon_type)) +
        theme_bw() + 
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=12, hjust = 0.5),
            axis.title.x = element_text(size=12, hjust = 0.5),
            axis.title.y = element_text(size=12, hjust = 0.5),
            plot.title = element_text(size = 12, hjust=0.5),
            strip.text.x = element_text(size = 12, hjust=0.5),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) 
    
    print(p)
    dev.off()
    return(p)
}

###
taxon_types = c("Archaea", "Bacteria", "Fungi", "Viruses") 
A = plot_top_genus(metadata, features_filter, taxon_type="Archaea", lineage=lineage)
B = plot_top_genus(metadata, features_filter, taxon_type="Bacteria", lineage=lineage)
F = plot_top_genus(metadata, features_filter, taxon_type="Fungi", lineage=lineage)
V = plot_top_genus(metadata, features_filter, taxon_type="Viruses", lineage=lineage)

pdf(file.path(outdir, paste0('merge_', lineage, '_boxplot.pdf')), width = 13, height = 8)
merge_img = plot_grid(A, B, F, V, ncol = 2, nrow=2, align="hv", labels=letters[1:4])
print(merge_img)
dev.off()

