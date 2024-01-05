####<--------------------------->####
workdir = "/home1/jialh/brain/01meta/multikingdom/03multi_omics"
setwd(workdir)
### - KEGG_modules.tab
### 本文件包括3列，分别是ModuleID, Module功能注释，Module所包含的KOs。
prepare_KEGG_modules <- function(KEGG_modules_infile){
    tmp = read.table (KEGG_modules_infile, sep = "\t")
    koann = strsplit (tmp[,3], split = ";")
    names (koann) = tmp[,1]
    
    module_mapping = tmp[,2] ### description of Kegg modules
    names (module_mapping) = tmp[,1] ; rm (tmp)
    
    ### remove KEGG references in square brackets for more clean names for plotting
    module_mapping_clean = sapply(module_mapping, function(x) strsplit(x, " \\[")[[1]][1])
    ## module_mapping_clean_df = as.data.frame(module_mapping_clean) #输出可以转化为dataframe,行名为moduleID, 列为module注释。
    
    ### Filtering out KEGG modules that are eukaryotic only
    ### 过滤掉只在真核生物中出现的modules。
    euk = unique (c ("M00352", "M00355", "M00354", "M00285", "M00295", "M00341", 
                     "M00177", "M00160", "M00359", "M00391", "M00182", "M00340", 
                     "M00359", "M00182", "M00160", "M00177", "M00391", "M00180", 
                     ### "eukaryotes" in name
                     "M00351", "M00352", "M00355", "M00354", "M00353", 
                     "M00427")) ### spliceosome or nuclear functions
    incl = setdiff (names (koann), euk)
    koann = koann [incl] ; rm (incl)
    ### Eukaryote filtering done
    return(list(koann=koann, module_mapping_clean=module_mapping_clean))
}

KEGG_modules_infile = file.path(workdir, "data", "KEGG_modules.tab")
modules = prepare_KEGG_modules(KEGG_modules_infile)
koann = modules$koann
module_mapping_clean = modules$module_mapping_clean
koann_df = as.data.frame(koann) ##记录了每个模块所包含的KOs信息。

tmpMat = array (NA, c (length (koann), 4, 2))
dimnames (tmpMat) [[1]] = names (koann)
dimnames (tmpMat) [[2]] = c ("NC_vs_SCS", "NC_vs_SCD", "NC_vs_MCI", "NC_vs_AD")
dimnames (tmpMat) [[3]] = c ("estimate", "p.value")


# ko_abundance_df = read.csv("/home1/jialh/brain/01meta/multikingdom/00profile/metadata_KOs.csv", header=1, row.names=1)
# ####<------------------------------------->####
# stage = "AD"
#ko_abundance_stage_df = ko_abundance_df[ko_abundance_df$Gender]

### 通过MUPPHin获得KOs与AD发病各个阶段的关联。  
### First, obtain the correlation coefficients for Spearman correlation between KOs and HOMA-IR
for(stage in c("SCS", "SCD", "MCI", "AD")){
    inputdir = paste0("/home1/jialh/brain/01meta/multikingdom/03MMUPHin_diff/KOs/NC_", stage)
    KO_cor_df = read.csv(file.path(inputdir, paste0("metaAD_MMUPHin_result_KOs_CHN+CHN2_NC_", stage, ".csv")), header=1, row.names = 1)
    KO_cor =  KO_cor_df['coef']
    #KO_cor_HOMA_df = as.data.frame(KO_cor_HOMA)  ##记录了总共6785个丰度KOs与HOMA的相关性信息。
    
    ### Then, test for difference in correlation coefficients between KOs in the
    ### KEGG module and all other KOs.
    ### 首先计算位于KEGG模块中的KOs和其他的KOs的相关系数是否存在差异。               
    for (k in names (koann)) {
        #k = names (koann)[1]
        incat =    na.omit (KO_cor_df[rownames (KO_cor) %in% koann [[k]],  ]$coef) ### select all correlations between HOMA-IR and KOs in the KEGG module
        notincat = na.omit (KO_cor_df[!(rownames (KO_cor) %in% koann [[k]]), ]$coef) ### select all correlations between HOMA-IR and KOs NOT in the KEGG module
        #print(paste0("incat: ", incat))
        if (length (incat) > 0 & length (notincat) > 0) {
            x = wilcox.test (incat, notincat)
            tmpMat [k, paste0("NC_vs_", stage), "p.value"] = x$p.value
            tmpMat [k, paste0("NC_vs_", stage), "estimate"] = (median (incat, na.rm = T) - median (notincat, na.rm = T))
        }
    }
}

tmpMat_df = as.data.frame(tmpMat)
cor_HOMA.IR [["keggmodules"]] <- tmpMat
#rm (tmpMat, ko.subset, ctrl.no.na.4KOs)
###保留keggmodules的分析结果，如下所示:
write.csv(tmpMat, "metaAD_keggmodules_association.csv")     

