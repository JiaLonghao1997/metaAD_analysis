####
library(ggplot2)
library(data.table)

rm(list=ls())

plot_Acc_F1_AUC <- function(outdir, data, stage, fs_method, test_study){
    # stage = 'AD'
    # fs_method = 'RandomForest+RFECV'
    # test_study = 'JPN'
    data_filter = data[(data$stage==stage) & (data$fs_method==fs_method) & (data$test_study==test_study), ]
    
    data_filter = data_filter[data_filter$feat_type %in% c("ABFV", "KOs+ConQuR", "pathways+ConQuR", "ABFV+KOs+pathways+ConQuR"), ]
    

    data_filter$feat_type=gsub("+ConQuR", "", data_filter$feat_type, fixed = TRUE)
    data_filter = data_filter[data_filter$feat_type %in% c("ABFV", "KOs"), ]
    
    data_filter$feat_type[data_filter$feat_type=='ABFV'] = 'Species'
    data_filter$feat_type[data_filter$feat_type=='pathways'] = 'Pathways'
    #data_filter$feat_type[data_filter$feat_type=='ABFV+KOs+pathways'] = 'All'
    
    #data_filter$measure[data_filter$measure=='auroc_mean'] = 'auroc'
    colnames(data_filter)[which(names(data_filter) == "accuracy")] <- "Accuracy"
    colnames(data_filter)[which(names(data_filter) == "precision")] <- "Precision"
    colnames(data_filter)[which(names(data_filter) == "recall")] <- "Recall"
    colnames(data_filter)[which(names(data_filter) == "f1")] <- "F1"
    colnames(data_filter)[which(names(data_filter) == "auroc_mean")] <- "AUC"
    
    
    print(colnames(data_filter))
    
    ####
    if (dim(data_filter)[1]>0){
        data_filter = melt(as.data.table(data_filter), id.vars=c('stage', 'train_study', 'test_study', 'model', 'feat_type', 'fs_method', 'sampler'),
                           measure.vars=c('Accuracy', "Precision", "Recall", 'F1', 'AUC'), 
                           variable.name = "measure", value.name = "value")
        
        print(unique(data_filter$measure))
        data_filter$measure = factor(data_filter$measure, levels=c('Accuracy', "Precision", "Recall", 'F1', 'AUC'))  
        
        print(unique(data_filter$feat_type))
        #data_filter$feat_type = factor(data_filter$feat_type, levels=c("Species", "KOs", "Pathways", "All"))
        data_filter$feat_type = factor(data_filter$feat_type, levels=c("Species", "KOs"))
        #colors = c('#03ACDA', '#6AB42D', '#ff993f', '#E03837')
        colors = c('#03ACDA', '#E03837')
        
        ##
        img <- ggplot(data=data_filter, aes(x=feat_type, y=value)) + 
            geom_col(aes(fill=feat_type)) +
            theme_bw() +
            facet_wrap(.~measure, ncol=5) +
            scale_fill_manual(values=colors) + 
            labs(x = 'Features', y = 'Count', title= paste0("NCvs ", stage, '( tranfer to ', test_study, ')')) +
            scale_y_continuous(breaks=seq(0,1.0,0.2), labels=seq(0,1.0,0.2), limits=c(0, 1.0)) + 
            theme(
                panel.grid = element_blank(),
                axis.text.x = element_text(size=12, angle=30, hjust=1),
                axis.text.y = element_text(size=12),
                axis.title.x = element_text(size=12),
                axis.title.y = element_text(size=12),
                strip.text.x = element_text(size = 12),
                legend.position = "none",  
                plot.title = element_text(size = 14, hjust=0.5))
        print(img)
        pdf(file.path(outdir, paste0("Acc_F1_AUC_", "NCvs ", stage, '_', test_study, '+', fs_method, '.pdf')), width=7, height=4)
        print(img)
        dev.off()
        
        png(file.path(outdir, paste0("Acc_F1_AUC_", "NCvs ", stage, '_', test_study, '+', fs_method, '.png')), width=1600, height=1600, res=300)
        print(img)
        dev.off()
    }
}

workdir = "/home1/jialh/brain/01meta/multikingdom/06transfer0710DCU"
data = read.csv(file.path(workdir, "performance_stat_clean.csv"))
outdir = file.path(workdir, "Acc_F1_AUC")
if(! file.exists(outdir)){
    dir.create(outdir)
}

###
#group[group[, .I[which.max(pt)], by=Subject]$V1]

print(unique(data$stage))
print(unique(data$fs_method))
print(unique(data$test_study))
for(stage in c("AD", "MCI")){
    for(fs_method in c("RandomForest+RFECV")){
        for(test_study in c("CHN+CHN2", "JPN", "GER")){
            plot_Acc_F1_AUC(outdir, data, stage, fs_method, test_study)
        }
    }
}
    
