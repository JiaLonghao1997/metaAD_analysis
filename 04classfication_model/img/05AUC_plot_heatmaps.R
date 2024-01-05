# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(cowplot)
library(grid)
library(circlize)


rm(list=ls())

plot_heatmap <- function(infile, metric, model, outdir){
    data = read.csv(infile, header=1, row.names = 1, stringsAsFactors = FALSE) 
    #data = data[row.names(data) %in% binners,  ]
    #print(typeof(data))
    #print(class(data))
    #data = as.matrix(data)
    #print(data)
    data$Means <-apply(data,1,0.5*sum)
    #print(data)
    if(metric=="AUC"){
        data$Means <- round(data$Means, 2)
        col_fun = colorRamp2(c(0, 0.5,1), c("#4575b4", "#e9f6e6", "#d73027"))
    }else if(metric=="n_features"){
        data$Means <- round(data$Means)
        col_fun = colorRamp2(c(0, 0.5*max(data), max(data)), c("#4575b4", "#e9f6e6", "#d73027"))
    }
    #data$Means <- round(data$Means, 2)
    data = as.matrix(data)
    #print(data)
    #print(max(data))
    #col_fun <- colorRamp2(c(0, 0.5*max(data), max(data)), c("blue", "white", "red"))
    #lgd <- Legend(col_fun = col_fun, title = "foo", direction = "horizontal")
    outfile = file.path(outdir, paste0(model, '_', metric,  "_multi_kingdom.png"))
    print(paste0("heatmap: ", outfile))
    png(outfile, width=800, height=600, res=300)
    Heatmap(data, name = paste0(model, '_', metric), col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, 
                 rect_gp = gpar(col = "white", lwd = 1),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%0.2f", data[i, j]), x, y, gp = gpar(fontsize = 7))},
                 #column_split = c('A', 'A', 'A', 'B'),
                 row_names_side = "left", column_names_side = "bottom",
                 column_title_side = "bottom", 
                 column_title = "studies", row_title = "features",
                 column_title_gp = gpar(fontsize=7),
                 row_title_gp = gpar(fontsize=7),
                 row_dend_gp = gpar(), 
                 column_dend_gp = gpar(),
                 row_names_gp = gpar(fontsize = 6),
                 column_names_rot = 30,
                 column_names_gp = gpar(fontsize = 6),
                 heatmap_legend_param = list(title = "AUC", legend_width = unit(1, "cm"), 
                                             position="top", #direction = "horizontal", 
                                             #title_position = "lefttop",
                                             title_gp = gpar(fontsize = 6),
                                             labels_gp = gpar(fontsize = 6))
    )
    dev.off()
    #return(p)
}

workdir = "/home1/jialh/brain/01meta/multikingdom/06classification_20231115/heatmaps"
outdir="/home1/jialh/brain/01meta/multikingdom/06classification_20231115/heatmaps"
fs_methods = c('MMUPHin+RFECV', 'MMUPHin_stage+RFECV', 'RandomForest+RFECV','RandomForest+RFECV_X2','RandomForest+RFECV_X4', 'all+Boruta', 'all+RFECV', 'MMUPHin+Boruta')
imglist = list()

for (fs_method in fs_methods){
    for(metric in c('AUC')){
        # fs_method = "RandomForest+RFECV"
        # metric = "AUC"
        #model = "randomforest"
        #metric = "n_features"
        infile <- file.path(workdir, paste0(fs_method, '_', metric,  ".csv"))
        if (file.exists(infile)){
            #plot_heatmap(infile, metric, model, outdir)
            data = read.csv(infile, header=1, row.names = 1, stringsAsFactors = FALSE) 
            #data = data[row.names(data) %in% binners,  ]
            #print(typeof(data))
            #print(class(data))
            #data = as.matrix(data)
            #print(data)
            data$Means <-apply(data, 1, mean)
            #data$Means <- 1.5 * data$Means
            #print(data)
            if(metric=="AUC"){
                data$Means <- round(data$Means, 2)
                col_fun = colorRamp2(c(0, 0.5,1), c("#4575b4", "#e9f6e6", "#d73027"))
                #col_fun = colorRamp2(c(0.5,1), c("#e9f6e6", "#d73027"))
                text_format = "%0.2f"
                
            }else if(metric=="n_features"){
                data$Means <- round(data$Means)
                col_fun = colorRamp2(c(0, 0.5*max(data), max(data)), c("#4575b4", "#e9f6e6", "#d73027"))
                text_format = "%d"
            }
            #data$Means <- round(data$Means, 2)
            data = as.matrix(data)
            #print(data)
            #print(max(data))
            #col_fun <- colorRamp2(c(0, 0.5*max(data), max(data)), c("blue", "white", "red"))
            #lgd <- Legend(col_fun = col_fun, title = "foo", direction = "horizontal")
            ##行按照指定顺序而不是自动排序: https://github.com/jokergoo/ComplexHeatmap/issues/959
            row_split <- factor(c(rep("A", 4), rep("B", 11), rep("C", 4)), levels=c(1,2,3))
            print(row_split)
            img <- Heatmap(data, name = metric, col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, 
                    rect_gp = gpar(col = "white", lwd = 1),
                    cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(sprintf(text_format, data[i, j]), x, y, gp = gpar(fontsize = 12))},
                    column_split = c('A', 'A', 'A', 'A', 'B', 'B', 'C', 'D'),
                    row_split = c(rep("A", 4), rep("B", 11), rep("C", 4)),
                    column_gap=unit(.02, "npc"),
                    row_gap=unit(.02, "npc"),
                    #row_split = factor(c(rep("single_kingdom", 4), rep("multi_kingdom", 11), rep("functional", 4))), 
                    
                    cluster_row_slices = FALSE,
                    row_names_side = "left", column_names_side = "bottom",
                    column_title_side = "top", 
                    column_title = fs_method, row_title = "Features",
                    column_title_gp = gpar(fontsize=14),
                    row_title_gp = gpar(fontsize=14),
                    row_dend_gp = gpar(fontsize=12), 
                    column_dend_gp = gpar(fontsize=12),
                    row_names_gp = gpar(fontsize = 12),
                    column_names_rot = 30,
                    column_names_gp = gpar(fontsize = 12),
                    heatmap_legend_param = list(title = "AUC", legend_width = unit(3, "cm"), 
                                                position="top", direction = "horizontal", 
                                                title_position = "lefttop",
                                                title_gp = gpar(fontsize = 12),
                                                labels_gp = gpar(fontsize = 12))
            )
            
            imglist[[fs_method]] = img
            draw(img, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side="top", show_heatmap_legend=TRUE)
            outfile = file.path(outdir, paste0(fs_method, '_', metric,  "_heatmap.pdf"))
            print(paste0("heatmap: ", outfile))
            pdf(outfile, width=7, height=12)
            draw(img, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side="top", show_heatmap_legend=TRUE)
            dev.off()
        }else{
            print(paste0(infile, " does not exist."))
        }
    }
}

print(length(names(imglist)))
pdf(file.path(outdir, "merge_heatmap.pdf"), width=24, height=8)
draw(imglist[[1]]+imglist[[2]]+imglist[[3]]+imglist[[4]]+imglist[[5]]+imglist[[6]]+imglist[[7]]+imglist[[8]], padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side="top", show_heatmap_legend=TRUE)
dev.off()
