library(data.table)
library(dplyr)
library(mediation)

rm(list=ls())
workdir = "/home1/jialh/brain/01meta/multikingdom/04mediation"
setwd(workdir)

outdir = file.path(workdir, "metaAD_results")
if(! file.exists(outdir)){
    dir.create(outdir, showWarnings = FALSE)
}

# load data ============================================
# load(file.path(workdir, "data/microbialFeatures_standardized.RData"))
# load(file.path(workdir, "data/meta_phenotypes.RData"))

## 用力肺活量（forced vital capacity, FVC），是指尽力最大吸气后，尽力尽快呼气所能呼出的最大气量。
## 略小于没有时间限制条件下测得的肺活量。该指标是指将测定肺活量的气体用最快速呼出的能力。
## 其中，开始呼气第一秒内的呼出气量为一秒钟用力呼气容积（forced expiratory volume in one second, FEV1.0），其临床应用较广，常以FEV1.0/FVC%表示。
##meta.cat包括34列, 主要是性别、是否吸烟、是否慢性阻塞性肺疾病(COPD)、生物燃料暴露、职业污染、是否患有其他疾病等分类型(category)变量。
##meta.num包括7列，主要是年龄、BMI、FEV1.0/FVC%、预测的一秒钟用力呼气容积、慢阻肺评估测试（COPD assessment test，CAT）、年平均PM2.5浓度等。
# print(summary(meta.cat))  
# print(summary(meta.num))
inputdir = "/home1/jialh/brain/01meta/multikingdom/00profile"
meta_taxon = as.data.frame(fread(file.path(inputdir, "metadata_species.csv"), sep=","))
meta_taxon = meta_taxon[meta_taxon$Batch %in% c("CHN", "CHN2"), ]
rownames(meta_taxon) = meta_taxon$Sample
print(rownames(meta_taxon)[1:5])


meta = meta_taxon[, 1:15]
species_df = meta_taxon[, 16:ncol(meta_taxon)]
source("/home1/jialh/brain/01meta/multikingdom/statistic/01preprocess/01feature_prepare.R")
print(colnames(species_df)[1:5])
print(rownames(species_df)[1:5])
species_df_filter = merge_feats(taxon=species_df, feat_type="ABFV", min_abundance=1e-4, min_prevalence=0.1)

species_df_filter <- tibble::rownames_to_column(species_df_filter, "Sample")
colnames(species_df_filter) = gsub("^.*s__", "", colnames(species_df_filter))

# exposures are the treators: 
Exposures <- colnames(species_df_filter)[2:ncol(species_df_filter)]

# lung functions:
print(colnames(meta))
Outcomes <- c("Group","ECOG_score","MMSE", "MoCA_B", "ACEIII_score")
#print(unique(meta$Group))
meta$Group <- stringr::str_replace_all(meta$Group, c("NC" = "0", "SCS" = "1", "SCD"= "2", "MCI" = "3", "AD" = "4"))
meta$Group <- as.integer(meta$Group)

# Covariates:
Covars <- c("Batch", "Age",  "BMI", "Gender", "ND_score")
covar_df <- meta %>% dplyr::select(Sample, all_of(Covars))
print(unique(covar_df$Batch))
covar_df$Batch <- as.integer(stringr::str_replace_all(covar_df$Batch, c("CHN" = "1", "CHN2" = "2")))
print(unique(covar_df$Gender))
covar_df$Gender <- as.integer(stringr::str_replace_all(covar_df$Gender, c("M" = "0", "F" = "1")))
covar_df$BMI <- as.numeric(covar_df$BMI)
covar_df$ND_score <- as.numeric(covar_df$ND_score)
#covar_df$Medication <- as.integer(sub( "Y", 1, sub("N",0, covar_df$Medication)))


# mediation analysis =================
merge.Mediation.forward <- NULL
merge.Mediation.reverse <- NULL

for(exp in Exposures){
  #exp=Exposures[1]
  species.expo <- species_df_filter %>% dplyr::select(Sample, all_of(exp))
  #species.expo$Biofuel_exposure <- as.integer(sub("N",0,sub("Y",1,species.expo$Biofuel_exposure)))
  
  for(outcome in Outcomes){
     #outcome = Outcomes[1]
     ##all_of() is for strict selection. If any of the variables in the character vector is missing, an error is thrown.
     meta.outcome = meta %>% dplyr::select(Sample, all_of(outcome))
     #print(meta.outcome)
     for(function_type in c("pathways", "KOs")){
       #function_type = "KOs"
       profile_dir = "/home1/jialh/brain/01meta/multikingdom/00profile"
       infile = file.path(profile_dir, paste0("metadata_", function_type, ".csv"))
       meta_KO <- as.data.frame(fread(infile, sep=","))  ##对于细菌，总共有1651个样本，100个变量(98个genus以及alpha diversity, PCoA1)。
       meta_KO = meta_KO[meta_KO$Batch %in% c("CHN", "CHN2"), ]
       rownames(meta_KO) = meta_KO$Sample
       #print(rownames(meta_KO)[1:5])
       
       
       meta = meta_KO[, 1:15]
       KO_df = meta_KO[, 16:ncol(meta_KO)]
       source("/home1/jialh/brain/01meta/multikingdom/statistic/01preprocess/01feature_prepare.R")
       # print(colnames(KO_df)[1:5])
       # print(rownames(KO_df)[1:5])
       KO.data = merge_feats(taxon=KO_df, feat_type=function_type, min_abundance=1e-6, min_prevalence=0.1)
       #meta_KO_filter <- tibble::rownames_to_column(meta_KO_filter, "Sample")
       if(function_type == "KOs"){
           important_KOs = read.csv("/home1/jialh/brain/01meta/multikingdom/06classification/important_features.csv", header=1)
           KO.data = KO.data[, colnames(KO.data) %in% important_KOs$features]
       }
       
       colnames(KO.data) <- gsub("-", "_", colnames(KO.data), fixed=TRUE)
       # print(paste0("Exposure: ", exp, "; Outcome: ", outcome, "; Feature: ", function_type))
       ########################################################################################
       Mediation.forward <- NULL
       Mediation.reverse <- NULL
       Mediation.forward.file = file.path(outdir, paste(exp, outcome, function_type, "mediation.forward.csv", sep="_"))
       Mediation.reverse.file = file.path(outdir, paste(exp, outcome, function_type, "mediation.reverse.csv", sep="_"))
       ########################################################################################
       if(file.exists(Mediation.forward.file) & file.exists(Mediation.reverse.file)){
           Mediation.forward = read.csv(Mediation.forward.file)
           Mediation.reverse = read.csv(Mediation.reverse.file)
       }else{
           pb <- utils::txtProgressBar(min = 0, max = 100, style = 3)
           for(i in 1:ncol(KO.data)){
             #i = 1
             #print(paste0("deal with ", colnames(KO.data)[i]))
             utils::setTxtProgressBar(pb, i/ncol(KO.data)*100)     # 更新进度条的值
             KO.dat.tmp <- KO.data[,i, drop=F]
             KO = colnames(KO.data)[i]
             
             dat <- merge(merge(merge(species.expo, meta.outcome, by = "Sample"),
                                covar_df, by="Sample"),
                          KO.dat.tmp, by.x="Sample", by.y = 0)
             
             dat <- dat[complete.cases(dat),]
             
             # forward mediation --------------------
             id = paste(exp, KO, outcome, sep = "|")
             #print(id)
             # 基于中介的线性回归模型，"bact.alpha.st ~ Biofuel_exposure + Age + BMI + Gender + Medication"。
             # "K00001 ~ Sample + Batch + Age + BMI + Gender + ND_score"
             model.m = glm(as.formula( paste(KO, " ~ ", exp, " + ", 
                                            paste(colnames(dat)[!(colnames(dat) %in% c("Sample",outcome, KO, exp))], collapse = " + "),
                                            sep = "") ), data = dat)
             # 基于结果的线性回归分析, "A_FEV1_FVC_Post ~ Biofuel_exposure + bact.alpha.st + Age + BMI + Gender + Medication"
             # Group ~ Candidatus_Methanomethylophilus_alvus + K00001 + Batch + Age + BMI + Gender + ND_score
             # print(dat$Group)
             model.y = glm(as.formula( paste(outcome, " ~ ", exp," + ", KO," + ", 
                                            paste(colnames(dat)[!(colnames(dat) %in% c("Sample",outcome, KO, exp))], collapse = " + "),
                                            sep = "") ), data = dat)
             
             # print(model.m$coefficients)
             # print(model.y$coefficients)
             # model.m: a fitted model object for mediator.
             # model.y: a fitted model object for outcome.
             if(anyNA(model.m$coefficients) | anyNA(model.y$coefficients)){
                 print(paste0("Exposure: ", exp, "; Outcome: ", outcome, "; Feature_type: ", function_type, "; feature: ", KO, " get NA coefficients in forward mediation."))
             }else{
                 summary = summary(mediate(model.m,model.y,treat=exp, mediator=KO, boot=F,sims=1000))
                 res <- capture.output(summary,append=FALSE)
                 
                 tmp <-  base::strsplit(res[grep("ACME",res)],"\\s")[[1]]
                 tmp <- tmp[tmp != "" & tmp!="."]
                 tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
                 ACME.p <- tmp[length(tmp)]
                 
                 tmp <- base::strsplit(res[grep("ADE",res)],"\\s")[[1]]
                 tmp <- tmp[tmp != "" & tmp!="."]
                 tmp <- tmp[!grepl("*",tmp,fixed = T) ]
                 ADE.p <- tmp[length(tmp)]
                 
                 tmp <- base::strsplit(res[grep("Prop. Mediated", res)],"\\s")[[1]]
                 tmp <- tmp[tmp != "" ]
                 i_str = which(grepl("Mediated", tmp))
                 prop.mediated <- tmp[(i_str + 1)]
                 
                 forw_vec = c(id, ACME.p, ADE.p, prop.mediated)
                 names(forw_vec) <- c("Treat_Mediator_Y", "ACME.p", "ADE.p", "prop.mediated")
                 
                 Mediation.forward <- bind_rows(Mediation.forward, forw_vec)
             }
             
             
             # reverse mediation --------------------
             id = paste(exp, outcome, KO, sep = "|")
             
             model.m = glm(as.formula( paste(outcome, " ~ ", exp, " + ", 
                                            paste(colnames(dat)[!(colnames(dat) %in% c("Sample",outcome, KO, exp))], collapse = " + "),
                                            sep = "") ), data = dat)
             
             model.y = glm(as.formula( paste(KO, " ~ ", exp," + ", outcome," + ", 
                                            paste(colnames(dat)[!(colnames(dat) %in% c("Sample",outcome, KO, exp))], collapse = " + "),
                                            sep = "") ), data = dat)
             if(anyNA(model.m$coefficients) | anyNA(model.y$coefficients)){
                 print(paste0("Exposure: ", exp, "; Outcome: ", outcome, "; Feature_type: ", function_type, "; feature: ", KO, " get NA coefficients in reverse mediation."))
             }else{
                 summary = summary(mediate(model.m,model.y,treat=exp,mediator=outcome,boot=F,sims=1000))
                 res <- capture.output(summary,append=FALSE)
                 
                 #sub( "^()\\s", "\\1", res[7])
                 tmp <-  base::strsplit(res[grep("ACME",res)],"\\s")[[1]]
                 tmp <- tmp[tmp != "" & tmp!="."]
                 tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
                 ACME.p <- tmp[length(tmp)]
                 
                 tmp <- base::strsplit(res[grep("ADE",res)],"\\s")[[1]]
                 tmp <- tmp[tmp != "" & tmp!="."]
                 tmp <- tmp[!grepl("*",tmp,fixed = T) ]
                 ADE.p <- tmp[length(tmp)]
                 
                 tmp <- base::strsplit(res[grep("Prop. Mediated", res)],"\\s")[[1]]
                 tmp <- tmp[tmp != "" ]
                 i_str = which(grepl("Mediated", tmp))
                 prop.mediated <- tmp[(i_str + 1)]
                 
                 revers_vec = c(id, ACME.p, ADE.p, prop.mediated)
                 names(revers_vec) <- c("Treat_Mediator_Y", "ACME.p", "ADE.p", "prop.mediated")
                 
                 Mediation.reverse <- bind_rows(Mediation.reverse, revers_vec)
             }
           }# loop through individual microbial features
           ##################################################################
           write.csv(Mediation.forward, Mediation.forward.file)
           write.csv(Mediation.reverse, Mediation.reverse.file)
           ##################################################################
       }
       merge.Mediation.forward <- bind_rows(merge.Mediation.forward, Mediation.forward)
       merge.Mediation.reverse <- bind_rows(merge.Mediation.reverse, Mediation.reverse)
     }# loop through microbial data frames
  }# loop through lung functions
}# loop through exposures

write.csv(merge.Mediation.forward, file.path(outdir, "merge.Mediation.forward.csv"))
write.csv(merge.Mediation.reverse, file.path(outdir, "merge.Mediation.reverse.csv"))