# ====================================================================================================
# Author:   Aaron
# Function: transform raw counts to vsd or rlog
# Version:  1.0
# Date:     Oct 15, 2018
# ====================================================================================================

rm(list=ls())
library(DESeq2)

setwd("/data/MyProgram/Final_diabrain/2.results/Cov_optimal/RData/")
dir.create("../trans_RData")
fname <- dir("./")
tname <- sub("\\.[[:alnum:]]+", "", fname)
for (i in 1:length(fname)) {
    print(paste("No.", i, fname[i], "is begin at", Sys.time()))
    
    load(fname[i])
    vsd <- vst(tdds, blind = FALSE)
    rld <- rlog(tdds, blind = FALSE)
    save(tdds, KeepGene, vsd, rld, file = paste0("../trans_RData/", tname[i], "_trans.RData"))
    
    print(paste("No.", i, fname[i], "is over at", Sys.time()))
}