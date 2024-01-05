#Install package to get relevant information from KEGG database

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("KEGGREST")

#Note: Uncomment the above code once the installation is successful

#Load package
library(KEGGREST)
library(utils)

workdir = "/home1/jialh/brain/01meta/multikingdom/03multi_omics"
module_dir = file.path(workdir, "modules")
if(! file.exists(module_dir)){
    dir.create(module_dir)
}

#Get list of modules in KEGG
print(module_list)
listDatabases()
# pathway_list = keggList("pathway")
# pathway_info <-keggGet(names(pathway_list)[1])

module_list <- keggList("module")


download_pathways_or_modules
## 第一个位置：新建一个其实进度条
pb <- txtProgressBar(style=3)


i = 1
merge_mod_df = NULL
for(module in names(module_list)){
    i = i + 1
    #module = "M00814"
    module_file = file.path(module_dir, paste0(module, "csv"))
    if(file.exists(module_file)){
        mod_df = read.csv(module_file, header=1)
    }else{
        print(paste0("deal with: ", module))
        #Search for corresponding ortholog
        module_info <-keggGet(module)
        #Save list of orthologs as a string separated by ","
        module <- module_info[[1]]$ENTRY
        module_name <- module_info[[1]]$NAME
        module_class <- module_info[[1]]$CLASS
        pathways <- paste(names(module_info[[1]]$PATHWAY), collapse = ",")
        KOs<-paste(names(module_info[[1]]$ORTHOLOGY),collapse = ",")
        mod_df = data.frame(module,module_name, module_class, pathways, KOs)
        colnames(mod_df)  <- c("module","module_name", "module_class", "pathways", "KOs")
        write.csv(mod_df, module_file, row.names=FALSE)
    }
    
    if(is.null(merge_mod_df)){
        merge_mod_df = mod_df
    }else{
        merge_mod_df = rbind(merge_mod_df, mod_df)
    }
    
    ###
    ## 第二个位置：实时反映进度
    setTxtProgressBar(pb, i/length(module_list))
}

## 第三个位置关闭进度条
close(pb)

# #Name columns
# colnames(mod_df)<-c("module","module_name", "module_class", "pathways", "KOs")

#Display first few entries in the table
head(merge_mod_df)

#Save table to csv file
write.csv(merge_mod_df,"Module_KOs.csv")

