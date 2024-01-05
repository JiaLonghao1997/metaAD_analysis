library("migest")
library(data.table)
library(dplyr)
library(tidyverse)


rm(list=ls())
workdir = "/home1/jialh/brain/01meta/multikingdom/04mediation"
data = as.data.frame(fread(file.path(workdir, "merge__Group_KOs_mediation.forward.csv"), sep=",", header=TRUE))
print(colnames(data))
data = data[data$Treat_Mediator_Y!="Treat_Mediator_Y", ]

data_filter = data[data$ACME.p<0.05, ]

data_filter[c('Treat', 'Mediator', 'Y')] <- stringr::str_split_fixed(data_filter$Treat_Mediator_Y, '\\|', 3)

###
profile_dir = "/home1/jialh/brain/01meta/multikingdom/00profile"
KOs_to_pathways = read.csv(file.path(profile_dir, "KEGG_KOs_to_pathways_metabolism.csv"), row.name=1, header=1)
KOs_to_pathways$KO_desc[KOs_to_pathways$KO_desc=="K00260 (Alanine, aspartate and glutamate metabolism)"] <- "K00260 (Ala, Asp and Glu metabolism)"
KOs_to_pathways$KO_desc[KOs_to_pathways$KO_desc=="K06208 (Phenylalanine, tyrosine and tryptophan biosynthesis)"] <- "K06208 (Phe, Tyr and  Trp biosynthesis)"
# KOs_to_pathways$KO_desc[KOs_to_pathways$KO_desc=="K00260 (Alanine, aspartate and glutamate metabolism)"] <- "K00260 (Ala, Asp and Glu metabolism)"
# KOs_to_pathways$KO_desc[KOs_to_pathways$KO_desc=="K00260 (Alanine, aspartate and glutamate metabolism)"] <- "K00260 (Ala, Asp and Glu metabolism)"


print(colnames(KOs_to_pathways))
KOs_to_pathways_filter = data.frame(KOs_to_pathways[, "KO_desc"])
rownames(KOs_to_pathways_filter) = rownames(KOs_to_pathways)
colnames(KOs_to_pathways_filter) = c("KO_desc")

data_filter_desc = merge(data_filter, KOs_to_pathways_filter, by.x="Mediator", by.y="row.names")

####
heatmap_df = read.csv("/home1/jialh/brain/01meta/multikingdom/03multi_omics/KOs_species_cor_heatmap.csv", header=1, row.names=1)

diff_species = read.csv("/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff/graphlan/MMUPHin_diff_merge.csv", header=1, row.names=1)
diff_species_filter = diff_species #[rowSums(is.na(diff_species))<=2, ]
print(rownames(diff_species_filter)[1:5])
diff_species_filter$genus = gsub("^.*g__(.*)\\|s__.*", "\\1", rownames(diff_species_filter))

rownames(diff_species_filter) = gsub("^.*s__", "", rownames(diff_species_filter))
diff_species_filter = diff_species_filter["genus"]

data_filter_desc_filter = merge(data_filter_desc, diff_species_filter, by.x="Treat", by.y="row.names")

data_filter_desc_filter = data_filter_desc_filter[data_filter_desc_filter$Treat %in% colnames(heatmap_df), ]
data_filter_desc_filter = data_filter_desc_filter[data_filter_desc_filter$KO_desc %in% rownames(heatmap_df), ]


data_filter_desc_filter$ACME.p = as.numeric(data_filter_desc_filter$ACME.p)
data_filter_desc_filter$ACME.logp = -log10(data_filter_desc_filter$ACME.p)

#log10(0.064)
write.csv(data_filter_desc_filter, file.path(workdir, "metaAD_circos.csv"))

###################################################################
dat = data_filter_desc_filter
print(colnames(dat))
colnames(dat)[colnames(dat)=="ACME.p"] = "Pval"
colnames(dat)[colnames(dat)=="ACME.logp"] = "log10P"
colnames(dat)[colnames(dat)=="Treat"] = "microbe"
#colnames(dat)[colnames(dat)=="genus"] = "MicroType"


# adjust the log10P value to a range 
library(scales)
dat$flow <-  rescale(dat$log10P, to = c(1.5, 5))


# dat <- dat %>%
#     mutate(genus_KO_desc = paste(genus, KO_desc,sep = "_"))

dat.d <- dat %>% 
    mutate(orig = microbe,
           #   flow=log10P,
           dest = KO_desc ) %>%
    dplyr::select(orig, dest, flow, Pval, genus)


#mig_chord(x = dat.d)
#参考: https://cloud.tencent.com/developer/article/2206509
library(ggalluvial)

print(unique(dat.d$genus))
dat.d$genus <- factor(dat.d$genus, levels=c("Bacteroides", "Barnesiella", "Faecalibacterium", "Parabacteroides", "Phocaeicola"))
colors = c("#cb1414", "#fd9a9a", "#16a3a7", "#7beeff", "#ffc1f6")
img <- ggplot(data=dat.d, aes(axis1=orig,axis2=dest, y=flow))+
        geom_alluvium(aes(fill=genus), width = 0.1, aes.bind = "flows") +
        scale_fill_manual(values=colors) + 
        geom_stratum(width=0.1, fill="#eeeeee")+
        geom_label(stat = "stratum", aes(label = after_stat(stratum)), size=3)+
        scale_x_continuous(breaks = 1:2, labels = c("microbes", "KOs"), expand = c(0.1,0.2)) +
        ggtitle("Mediation analysis") + 
        theme_void() +
        theme(axis.text.x = element_text(size=12),
              legend.position = "top", 
              plot.title = element_text(size=14, hjust = 0.5))
print(img)

print(workdir)
pdf(file.path(workdir, "01sankey_plot.pdf"), width=7, height= 4)
print(img)
dev.off()

