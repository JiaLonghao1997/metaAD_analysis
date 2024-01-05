rm(list=ls())
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggpubr)

stages = c("SCS", "SCD", "MCI", "AD")
inputdir = "/home1/jialh/brain/01meta/multikingdom/06classification/feature_importance_boxplots"

imglist = list()
for (stage in stages){
    stage = "AD"
    infile = file.path(inputdir, paste0("CHN+CHN2_", stage, "_KOs_feature_importance_sort.csv"))
    data = read.csv(infile, header=1, row.names=1)
    print(colnames(data))
    print(unique(data$Significance))
    data$Significance <- factor(data$Significance, levels=c("Depletion (P<0.05)", "Elevation (P<0.05)", "Not significant"))
    colors = c("#36bccb", "#ff6836", "#dddddd")
    colnames(data)[colnames(data)=="Contribution.to.model...."] = "contribution_to_model"
    # print(unique(data$KO_desc))
    data$KO_desc = factor(data$KO_desc, levels=rev(unique(data$KO_desc)))
    data$pathway_class = "other"
    print(colnames(data))
    print(unique(data$Pathway_short_name))
    data_stat = as.data.frame(table(data$Pathway_short_name))
    print(data_stat)
    ##绘制pie图: https://r-graph-gallery.com/piechart-ggplot2.html
    colnames(data_stat) = c("pathways", "count")
    data_stat <- data_stat[order(data_stat$count, decreasing = TRUE),]
    print(data_stat)
    
    
    #print(data_stat$pathways)
    # data_stat$pathways = factor(data_stat$pathways, levels=c(
    #     "Metabolism of cofactors and vitamins",  
    #     "Carbohydrate metabolism", 
    #     "Energy metabolism", 
    #     "Amino acid metabolism",
    #     "Metabolism of other amino acids", 
    #     "Lipid metabolism",
    #     "Biosynthesis of other secondary metabolites",
    #     "Nucleotide metabolism", 
    #     "Xenobiotics biodegradation and metabolism",
    #     "Glycan biosynthesis and metabolism",
    #     "Metabolism of terpenoids and polyketides"))
    # colors = c(
    #     "#cb1414",        ##"Metabolism of cofactors and vitamins",  
    #     "#fd9a9a", ##"Carbohydrate metabolism"
    #     "#fdb103",  ##"Energy metabolism" 
    #     "#16a3a7", #"Amino acid metabolism"
    #     "#fdb103", # Metabolism of other amino acids
    #     "#c3c00e", # "Lipid metabolism",
    #     "#7a4b03", # "Biosynthesis of other secondary metabolites"
    #     "#6992d8",  ##"Nucleotide metabolism"
    #     "#c9e4ff",  ###"Xenobiotics biodegradation and metabolism"
    #     "#7beeff", ##"Glycan biosynthesis and metabolism"
    #     "#ffc1f6"  ##"Metabolism of terpenoids and polyketides"
    # )
    colors = unlist(list(
        "Metabolism of cofactors and vitamins"="#cb1414",
        "Carbohydrate metabolism"="#fd9a9a",
        "Energy metabolism"="#cb7014",
        "Amino acid metabolism"="#16a3a7",
        "Metabolism of other amino acids"="#fdb103",
        "Lipid metabolism"="#c3c00e",
        "Biosynthesis of other secondary metabolites"="#7a4b03",
        "Nucleotide metabolism"="#6992d8",
        "Xenobiotics biodegradation and metabolism"="#c9e4ff",
        "Glycan biosynthesis and metabolism"="#7beeff",
        "Metabolism of terpenoids and polyketides"="#ffc1f6"
    ))
    
    
    # Compute the position of labels
    data_stat <- data_stat %>% 
        arrange(desc(pathways)) %>%
        mutate(prop = count / sum(data_stat$count) *100) %>%
        mutate(ypos = cumsum(prop)- 0.5*prop )
    
    
    
    # Basic piechart
    img <- ggplot(data_stat, aes(x="", y=prop, fill=pathways)) +
        geom_bar(stat="identity", width=1, color="white") +
        #geom_text(aes(x=1.2, label = paste0(round(prop, 2), " %")), position = position_stack(vjust = 0.5), color = "white", size=4) +
        coord_polar("y", start=0) +
        geom_label_repel(data = data_stat,
                         aes(y = ypos, label = paste0(round(prop,2), "%")),
                         size = 4, nudge_x = 0.5, show.legend = FALSE) +
        theme_void() + 
        guides(fill=guide_legend(ncol=3)) + 
        theme(
              plot.title = element_text(size=14, hjust=0.5),
              axis.text = element_blank(),
              axis.title = element_blank(),
              legend.position="bottom",
              legend.text = element_text(size=12),
              legend.title = element_text(size=12)) +
        labs(title=paste0("NC vs ", stage)) + 
        scale_fill_manual(values=colors)
    
  
    print(img)
    dev.off()
    
    pdf(file.path(inputdir, paste0(stage, "_important_pie.pdf")), width=4, height=4)
    print(img)
    dev.off()
    
    imglist[[stage]] = img
}

setwd("/home1/jialh/brain/01meta/multikingdom")
merge_img = ggarrange(plotlist=imglist, nrow=1, ncol=4, align="hv", labels=letters[1:4], legend="bottom", common.legend = TRUE)
pdf(file.path(inputdir, "merge_important_pie.pdf"), width=13, height=5)
print(merge_img)
dev.off()
