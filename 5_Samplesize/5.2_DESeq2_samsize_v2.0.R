# ====================================================================================================
# Author:   Aaron
# Function: randomly choose subset of T2D samples with matched control samples, use DESeq2 to identify
#           diferentially expressed genes. 
# Version:  2.0
# Date:     Mar 5, 2019
# ====================================================================================================


rm(list=ls())
library(DESeq2)
library(sva)
#library(foreach)
#library(doParallel)
#cl <- makeCluster(16)
#registerDoParallel(cl)
library(BiocParallel)
register(MulticoreParam(6))
options(stringsAsFactors = FALSE)
posnum <- commandArgs(T)

# ----------------------------------------------------------------------------------------------
# Initiation. 
# ---------------------------------------------------------------------------------------------
setwd("/home/zhouzhe/Final_diabrain/1.clean/data/")
cmethod <- "optimal"
fname <- dir("./")
tname <- sub("_[0-9]+.[[:alnum:]]+", "", fname)
PTIMES <- 100
MAX_SAMSIZE <- 35
CUT_OFF <- 0.05

# Initial an array to store the result of differentially expression analysis
degenes <- array(NA, dim = c(length(fname), MAX_SAMSIZE/5-1, PTIMES), 
                 dimnames = list(sub("Brain_", "", tname), 
                                 paste0("s", seq(from = 10, to = MAX_SAMSIZE, by = 5)), 
                                 c(1:PTIMES)))


# ----------------------------------------------------------------------------------------------
# Select sample according to input, and perform DESeq2 to identify diferentially expressed genes. 
# ----------------------------------------------------------------------------------------------
DEfun <- function(select_sam) {
    #library(DESeq2)
    #library(sva)
    # Select the samples -----------------------------------------------------------------------
    sams <- matchMatrix[select_sam, ]
    sams <- c(sams$T2D, sams$Ctrl1, sams$Ctrl2)
    
    tsaminfo <- hsaminfo[match(sams, rownames(hsaminfo)), ]
    tdata <- hdata[, rownames(tsaminfo)]
    
    diab <- tsaminfo$DIABETES
    tsaminfo$DIABETES <- as.factor(tsaminfo$DIABETES)
    tsaminfo$RACE <- as.factor(tsaminfo$RACE)
    tsaminfo$SEX <- as.factor(tsaminfo$SEX)
    
    # Body -------------------------------------------------------------------------------------
    ## Construct DESeqDataSet and filter out lower expression
    if (length(table(tsaminfo$RACE))>1) {
        if (length(table(tsaminfo$SEX))>1) {
            dds_design <- formula(~ AGE+RACE+SEX+BMI+SMRIN+DIABETES)
            dds_design0 <- formula(~ AGE+RACE+SEX+BMI+SMRIN)
        } else {
            dds_design <- formula(~ AGE+RACE+BMI+SMRIN+DIABETES)
            dds_design0 <- formula(~ AGE+RACE+BMI+SMRIN)
        }
    } else if (length(table(tsaminfo$SEX))>1) {
        dds_design <- formula(~ AGE+SEX+BMI+SMRIN+DIABETES)
        dds_design0 <- formula(~ AGE+SEX+BMI+SMRIN)
    } else {
        dds_design <- formula(~ AGE+BMI+SMRIN+DIABETES)
        dds_design0 <- formula(~ AGE+BMI+SMRIN)
    }
    
    tdds <- DESeqDataSetFromMatrix(countData = tdata,
                                   colData = tsaminfo,
                                   design = dds_design)
    keep <- rowMeans(counts(tdds)) >= 2
    tdds <- tdds[keep, ]
    tdds$DIABETES <- factor(tdds$DIABETES, levels = c("0", "1"))
    
    ## svaseq 
    tdds <- estimateSizeFactors(tdds)
    dat  <- counts(tdds, normalized = TRUE)
    mod  <- model.matrix(dds_design, colData(tdds))
    mod0 <- model.matrix(dds_design0, colData(tdds))
    svseq <- svaseq(dat, mod, mod0, n.sv = 3)
    
    ## remove sva factor correlate with diab
    svmatrix <- data.frame(svseq$sv, diab)
    svcomp <- split(svmatrix, svmatrix$diab)
    pcor <- numeric(dim(svseq$sv)[2])
    for (j in 1:dim(svseq$sv)[2]) {
        pcor[j] <- t.test(svcomp[[1]][, j], svcomp[[2]][, j])$p.value
    }
    idx <- pcor > 0.05

    # if only 0 or 1 surrogate variable kept, there will cause error, so..
    if (all(!idx)) {
	nsv <- 0
    } else {
	tmp <- svseq$sv[, idx, drop = FALSE]
    	nsv <- dim(tmp)[2]
    
	colnames(tmp) <- paste0("SV", 1:nsv)
        colData(tdds) <- cbind(colData(tdds), tmp)
        design(tdds) <- as.formula(paste("~", dds_design0[2], "+", 
                                     paste(colnames(tmp), collapse = "+"), "+DIABETES"))
    } 
    print(paste("There are", nsv, "surrogate variables identified which is not correlated with DIABETES"))
    print(design(tdds))
    
    #Defferential Analysis
    tdds <- DESeq(tdds, parallel = TRUE)
    tres05 <- results(tdds, alpha=0.05)
    
    #Summary--------------------------------------------------------------------
    sum(tres05$padj < CUT_OFF, na.rm = TRUE)
}


system.time(
for (i in posnum[1]:posnum[2]) {
#for (i in 1:length(fname)) {
    print(paste("No.", i, fname[i], "is begin at", Sys.time()))
    
    # -----------------------------------------------------------------------------------------------
    # Prepare the expression matrix, match matrix and sample info
    # -----------------------------------------------------------------------------------------------
    ## expression matrix
    hdata <- read.table(fname[i], header=T, sep=" ", quote="")
    row.names(hdata) <- hdata$Name
    hdata <- hdata[,c(-1,-2)]
    
    ## sample info and match matrix
    tmp <- dir(paste0("/home/zhouzhe/Final_diabrain/1.clean/math_matrix/"), full.names = T)
    matchMatrix <- read.table(tmp[grep(tname[i], tmp)], header = TRUE)
    tmp <- dir(paste0("/home/zhouzhe/Final_diabrain/1.clean/matched_sam_", cmethod), full.names = T)
    hsaminfo <- read.table(tmp[grep(tname[i], tmp, fix=TRUE)], header = T, sep = "\t") 
    
    ## Check if the match matrix is right
    rownames(hsaminfo) <- make.names(hsaminfo$SAMPID)
    hdata <- hdata[, rownames(hsaminfo)]
    colnames(hsaminfo) <- c(colnames(hsaminfo)[-8], "DIABETES")
    idx <- hsaminfo$DIABETES == 0
    print(all(sort(rownames(hsaminfo)[!idx]) == sort(matchMatrix$T2D)))
    print(all(sort(rownames(hsaminfo)[idx]) == sort(c(matchMatrix$Ctrl1, matchMatrix$Ctrl2))))
    
    # -----------------------------------------------------------------------------------------------
    # Sample size 
    # -----------------------------------------------------------------------------------------------
    nT2D <- nrow(matchMatrix)
    for (sample_size in seq(from = 10, to = MAX_SAMSIZE, by = 5)) {
        if (sample_size > nT2D) next
        if (choose(nT2D, sample_size) <= PTIMES) {
            all_comb <- combn(nT2D, sample_size)
            maxn <- choose(nT2D, sample_size)
        } else {
            all_comb <- replicate(PTIMES, sample(c(1:nT2D), sample_size))
            maxn <- PTIMES
        }
        for (j in 1:maxn) {
	    degenes[i, sample_size/5-1, j] <- DEfun(all_comb[, j])
	}

	#tmp <- foreach(j=1:maxn, .combine = "c") %dopar% DEfun(all_comb[, j]) 
        #degenes[i, sample_size/5-1, 1:length(tmp)] <- tmp 
        write.table(degenes[i, sample_size/5-1, ], paste0("/home/zhouzhe/Final_diabrain/DESeq2data/", tname[i], "_s", sample_size), 
                    quote = F, row.names = F, col.names = F, na = "NA")
    }
    print(paste(tname[i], "is over at", Sys.time()))
}    
)

# ----------------------------------------------------------------------------------------------
# Save the result. 
# ----------------------------------------------------------------------------------------------
#demeans <- apply(degenes, MARGIN = c(1, 2), mean, na.rm=T, simplify = T)
#demedians <- apply(degenes, MARGIN = c(1, 2), median, na.rm=T, simplify = T)
save.image(paste0("/home/zhouzhe/Final_diabrain/DESeq2data/", posnum[1],"_", posnum, ".RData"))
