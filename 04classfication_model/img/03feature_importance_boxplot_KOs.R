library(ggplot2)
library(ggpubr)

stages = c("SCS", "SCD", "MCI", "AD")
inputdir = "/home1/jialh/brain/01meta/multikingdom/06classification_20231115/feature_importance_boxplots_KOs"

imglist = list()
for (stage in stages){
    # stage = "AD"
    infile = file.path(inputdir, paste0("CHN+CHN2_", stage, "_KOs_feature_importance.csv"))
    data = read.csv(infile, header=1, row.names=1)
    print(colnames(data))
    print(unique(data$Significance))
    data$Significance <- factor(data$Significance, levels=c("Depletion (P<0.05)", "Elevation (P<0.05)", "Not significant"))
    colors = c("#36bccb", "#ff6836", "#dddddd")
    colnames(data)[colnames(data)=="Contribution.to.model...."] = "contribution_to_model"
    print(unique(data$KO_desc))
    data$KO_desc = factor(data$KO_desc, levels=rev(unique(data$KO_desc)))
    data$pathway_class = "other"
    print(colnames(data))
    #print(unique(data$Pathway_category2))
    print(table(data$KO_desc))
    # ###氨基酸代谢。
    # data$pathway_class[data$Pathway_category2 %in% c("09105 Amino acid metabolism", "09106 Metabolism of other amino acids")]  = "Amino acid metabolism"
    # ###碳水化合物代谢。
    # data$pathway_class[data$Pathway_category2 %in% c("09101 Carbohydrate metabolism", "09107 Glycan biosynthesis and metabolism")]  = "Carbohydrate/Glycan metabolism"
    # ###核酸代谢
    # data$pathway_class[data$Pathway_category2 == "09104 Nucleotide metabolism"]  = "Nucleotide metabolism"
    # ###能量代谢
    # data$pathway_class[data$Pathway_category2 == "09102 Energy metabolism"]  = "Energy metabolism"
    # 
    # data$pathway_class <- factor(data$pathway_class, levels=c("Amino acid metabolism", "Carbohydrate/Glycan metabolism", "Nucleotide metabolism", "Energy metabolism", "other"))
    #####<---------------------------->####
    p = ggplot(data, aes(x=KO_desc, y=contribution_to_model)) + 
        geom_boxplot(aes(fill=Significance)) + 
        scale_fill_manual(values=colors) +
        coord_flip() +
        labs(title=paste0("NC vs ", stage), y ="Contribution to model (%)", x="") +
        theme_bw() + 
        theme(
            axis.text.x = element_text(size=12),
            axis.title.x = element_text(size = 12),
            axis.text.y = element_text(size = 10),
            axis.title.y = element_text(size = 12),
            plot.title = element_text(hjust=0.5, size=16),
            plot.subtitle = element_text(hjust=0.5, size=12),
            strip.text.x = element_text(size=12),
            legend.position = "none",
            legend.text=element_text(size=10),
            legend.title=element_text(size=12)
        )
    print(p)
    pdf(file.path(inputdir, paste0(stage, "_important_boxplot.pdf")), width=6, height=6)
    print(p)
    dev.off()
    
    imglist[[stage]] = p
}

merge_img = ggarrange(plotlist=imglist, nrow=2, ncol=2, align="hv", labels=letters[1:4])
pdf(file.path(inputdir, "merge_important_KOs.pdf"), width=13, height=8)
print(merge_img)
dev.off()
