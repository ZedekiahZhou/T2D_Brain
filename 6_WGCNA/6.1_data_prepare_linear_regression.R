# ====================================================================================================
# Author:   Aaron
# Function:  regress out covariates and SVs, filter out low expressed genes
# Version:  1.0
# Date:     Oct 15, 2018
# ====================================================================================================

rm(list=ls())
library(sva)
library(DESeq2)
setwd("/data/MyProgram/Final_diabrain/7.scripts/6_WGCNA/")
dir.create("covar_prepare/")

#-------------------------------------------------------------------------------------------------
# I. Prepare the datExpr
#-------------------------------------------------------------------------------------------------

#mytissue <- "Brain_Caudate"
#nsv <- 15
fname <- dir("/data/MyProgram/Final_diabrain/2.results/Cov_optimal/trans_RData", full.names = TRUE)
#sname <- dir("/data/MyProgram/Final_diabrain/1.clean/matched_sam_optimal", full.names = TRUE)
tname <- sub("[[:print:]]+trans_RData/", "", sub("_trans.RData", "", fname))

for (mytissue in tname) {
    print(paste0(mytissue, " begin at ", Sys.time()))
    linput <- load(fname[grep(mytissue, fname, fixed = TRUE)])
    #tsaminfo <- read.table(sname[grep(mytissue, sname, fix=TRUE)], header = T, sep="\t", quote="", stringsAsFactors = F) 
    #colnames(tsaminfo) <- c(colnames(tsaminfo)[-8], "DIABETES")
    
    tres05 <- results(tdds, alpha = 0.05)
    idx <- !is.na(tres05$padj)
    print(table(idx))
    tres05 <- tres05[idx, ]
    rld <- rld[idx, ]
    tdds <- tdds[idx, ]
    mycol <- data.frame(colData(tdds))
    
    # 2. covariate + SV
    dds_residual <- assay(rld)
    if (length(table(mycol$RACE))>1) {
        if (length(table(mycol$SEX))>1) {
            myformula <- "geneexp~SEX+AGE+RACE+BMI+SMRIN"
        } else {
            myformula <- "geneexp~AGE+RACE+BMI+SMRIN"
        }
    } else if (length(table(mycol$SEX))>1) {
        myformula <- "geneexp~SEX+AGE+BMI+SMRIN"
    } else {
        myformula <- "geneexp~AGE+BMI+SMRIN"
    }
    print(myformula)
    sv_factor <- names(mycol)[grep("SV", names(mycol))]
    myformula <- paste0(myformula, "+", paste0(sv_factor, collapse = "+"))
    print(myformula)
    for (i in 1:dim(rld)[1]) {
        geneexp <- dds_residual[i, ]
        fit1 <- lm(myformula, data = mycol)
        dds_residual[i, ] <- fit1$residual
    }
    save(tdds, rld, vsd, dds_residual, mycol, file = paste0("covar_prepare/", mytissue, ".RData"))
    print(paste0(mytissue, " over at ", Sys.time()))
}





# if (0) {
#     # 1. covariate + 15SVss
#     dds_residual <- assay(rld)
#     for (i in 1:dim(rld)[1]) {
#         geneexp <- dds_residual[i, ]
#         myformula <- paste0("geneexp~SEX+AGE+RACE+BMI+SMRIN+", paste0("SV", 1:nsv, collapse = "+"))
#         fit1 <- lm(myformula, data = mycol)
#         dds_residual[i, ] <- fit1$residual
#     }
#     dir.create("./covar_15SV")
#     save.image(file = paste0("covar_15SV/", mytissue, ".RData"))
# }
# if(0) {
# # 3.rld
# dds_residual <- assay(rld)
# 
# dir.create("./rld")
# save.image(file = paste0("rld/", mytissue, ".RData"))
# }
