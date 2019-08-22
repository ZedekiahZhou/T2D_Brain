# ====================================================================================================
# Author:   Aaron
# Function: Identified deferentially expressed genes using DESeq2 in all tissues with maxium number of 
#           samples.
# Version:  1.0
# Date:     Sep 17, 2018
# ====================================================================================================


rm(list=ls())
sink(file = "/data/MyProgram/Final_diabrain/7.interdata/4.1_Find_DG_opt.log")
library(DESeq2)
library(sva)
library(BiocParallel)
register(MulticoreParam(10))

fpath <- "/data/MyProgram/Final_diabrain/2.results/Cov_optimal/"
dir.create(fpath, recursive = T)
cmethod <- "optimal"
setwd(fpath)

dir.create("./RData");dir.create("./result1");dir.create("./result05")
setwd("./result1");dir.create("./Total_Gene");dir.create("./Up_Gene");dir.create("./Down_Gene");dir.create("./res")
setwd("../result05");dir.create("./Total_Gene");dir.create("./Up_Gene");dir.create("./Down_Gene");dir.create("./res")

setwd("../../../1.clean/data/")
fname <- dir("./")
tname <- sub("_[0-9]+.[[:alnum:]]+", "", fname)
#codinggene <- read.table("../codinggene_20034.list", stringsAsFactors = F)[[1]]
#seqgene <- read.table("../RNAseq_gene_symbol_56202.list", header=T, stringsAsFactors = F)[[1]]
#coding_num <- seqgene %in% codinggene
id2symbol <- read.table("../id2symbol.tab", header = T, sep = "\t", 
                        quote = "", stringsAsFactors = F)


genesummary <- matrix(nrow=length(fname), ncol=9)
colnames(genesummary) <- c("total1", "up1", "down1", 
                           "total05", "up05", "down05",
                           "total01", "up01", "down01")
genesummary <- data.frame(genesummary, Tissue=tname)
tsam_size <- matrix(nrow=length(fname), ncol=3)
colnames(tsam_size) <- c("Ctrl", "T2D", "Total")
tsam_size <- data.frame(tsam_size, Tissue=tname)

for (i in 1:length(fname)) {
    print(paste("No.", i, fname[i], "is begin at", Sys.time()))
    
    #Prepare--------------------------------------------------------------------
    tdata <- read.table(fname[i], header=T, sep=" ", quote="")
    row.names(tdata) <- tdata$Name
    tdata <- tdata[,c(-1,-2)]
    tmp <- dir(paste0("../matched_sam_", cmethod), full.names = T)
    tsaminfo <- read.table(tmp[grep(tname[i], tmp, fix=TRUE)], 
                           header = T, sep="\t", quote="", stringsAsFactors = F) 
    #tsaminfo <- tsaminfo[, c("SAMPID", "AGE", "RACE", "SEX", "BMI",  "DIABETES")]
    
    #Check if the rownames of tsaminfo and colnames of tdata are the same
    rownames(tsaminfo) <- make.names(tsaminfo$SAMPID)
    tdata <- tdata[, rownames(tsaminfo)]
    colnames(tsaminfo) <- c(colnames(tsaminfo)[-8], "DIABETES")
    
    #filter out type I Diabetes
    diab <- tsaminfo$DIABETES
    tsaminfo$DIABETES <- as.factor(tsaminfo$DIABETES)
    tsaminfo$RACE <- as.factor(tsaminfo$RACE)
    tsaminfo$SEX <- as.factor(tsaminfo$SEX)
    
    #Body ----------------------------------------------------------------------
    #Construct DESeqDataSet and filter out lower expression
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
    print(design(tdds))
    
    #tdds
    keep <- rowMeans(counts(tdds)) >= 5
    sum(keep)
    tdds <- tdds[keep, ]
    KeepGene <- row.names(counts(tdds))
    tdds$DIABETES <- factor(tdds$DIABETES, levels = c("0", "1"))
    tsam_size[i, 1:2] <- table(tdds$DIABETES)
    tsam_size[i, 3] <- dim(tsaminfo)[1]
    

    #sva 
    tdds <- estimateSizeFactors(tdds)
    dat  <- counts(tdds, normalized = TRUE)
    mod  <- model.matrix(dds_design, colData(tdds))
    mod0 <- model.matrix(dds_design0, colData(tdds))
    #n.sv = num.sv(dat,mod,method="leek")
    svseq <- svaseq(dat, mod, mod0, n.sv = 3)
    #plot(svseq, pch=19, col="blue")
    
    #pcor <- apply(svseq$sv, 2, function(x) {cor.test(x, diab)$p.value})
    svmatrix <- data.frame(svseq$sv, diab)
    svcomp <- split(svmatrix, svmatrix$diab)
    pcor <- numeric(dim(svseq$sv)[2])
    for (j in 1:dim(svseq$sv)[2]) {
      pcor[j] <- t.test(svcomp[[1]][, j], svcomp[[2]][, j])$p.value
    }
    
    idx <- pcor > 0.05
    tmp <- svseq$sv[, idx]
    nsv <- dim(tmp)[2]
    colnames(tmp) <- paste0("SV", 1:nsv)
    colData(tdds) <- cbind(colData(tdds), tmp)
    design(tdds) <- as.formula(paste("~", dds_design0[2], "+", 
                                   paste(colnames(tmp), collapse = "+"), "+DIABETES"))
    print(paste("There are", nsv, "surrogate variables identified which is not correlated with DIABETES"))
    print(design(tdds))
    
    #Defferential Analysis
    system.time(tdds <- DESeq(tdds, parallel = TRUE))
    tres1 <- results(tdds) #tres
    tres05 <- results(tdds, alpha=0.05)
    tres01 <- results(tdds, alpha=0.01)
    summary(tres1)
    summary(tres05)
    
    resultsNames(tdds)
    #tresLFC05 <- lfcShrink(dds=tdds, coef=length(resultsNames(tdds)), res=tres05, parallel=TRUE) #tresLFC
    
    # ------------------------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------------------------
    genesummary[i, 1] <- sum(tres1$padj < 0.1, na.rm = T)
    genesummary[i, 2] <- sum((tres1$padj < 0.1)&(tres1$log2FoldChange>0), na.rm = T)
    genesummary[i, 3] <- sum((tres1$padj < 0.1)&(tres1$log2FoldChange<0), na.rm = T)
    genesummary[i, 4] <- sum(tres05$padj < 0.05, na.rm = T)
    genesummary[i, 5] <- sum((tres05$padj < 0.05)&(tres05$log2FoldChange>0), na.rm = T)
    genesummary[i, 6] <- sum((tres05$padj < 0.05)&(tres05$log2FoldChange<0), na.rm = T)
    genesummary[i, 7] <- sum(tres01$padj < 0.01, na.rm = T)
    genesummary[i, 8] <- sum((tres01$padj < 0.01)&(tres01$log2FoldChange>0), na.rm = T)
    genesummary[i, 9] <- sum((tres01$padj < 0.01)&(tres01$log2FoldChange<0), na.rm = T)
    
    save(tdds, KeepGene,  
         file=paste0(fpath, "RData/", tname[i], ".RData"))
    
    #Extract genelist
    DEix <- (!is.na(tres1$padj))&(tres1$padj<0.1)
    UpDEix <- (!is.na(tres1$padj))&(tres1$padj<0.1)&(tres1$log2FoldChange>0)
    DownDEix <- (!is.na(tres1$padj))&(tres1$padj<0.1)&(tres1$log2FoldChange<0)
    DEres1 <- data.frame(signif(data.frame(tres1[DEix, ]), 3), Gene=KeepGene[DEix])
    ToGene <- id2symbol$Description[match(rownames(DEres1),id2symbol$Name)]
    print(paste("Are the Genes and rows match for p=0.1? :", all(ToGene==DEres1$Gene)))
    write(rownames(tdds)[DEix], file=paste0(fpath, "result1/Total_Gene/", tname[i], ".DE"))
    write(rownames(tdds)[UpDEix], file=paste0(fpath, "result1/Up_Gene/", tname[i], ".up"))
    write(rownames(tdds)[DownDEix], file=paste0(fpath, "result1/Down_Gene/", tname[i], ".down"))
    write.table(DEres1, quote = F, row.names = F, sep="\t", 
                file=paste0(fpath, "result1/res/", tname[i], ".res"))
    
    #Extract genelist
    DEix <- (!is.na(tres05$padj))&(tres05$padj<0.05)
    UpDEix <- (!is.na(tres05$padj))&(tres05$padj<0.05)&(tres05$log2FoldChange>0)
    DownDEix <- (!is.na(tres05$padj))&(tres05$padj<0.05)&(tres05$log2FoldChange<0)
    DEres05 <- data.frame(signif(data.frame(tres05[DEix, ]), 3), Gene=KeepGene[DEix])
    ToGene <- id2symbol$Description[match(rownames(DEres05),id2symbol$Name)]
    print(paste("Are the Genes and rows match for p=0.05? :", all(ToGene==DEres05$Gene)))
    write(rownames(tdds)[DEix], file=paste0(fpath, "result05/Total_Gene/", tname[i], ".DE"))
    write(rownames(tdds)[UpDEix], file=paste0(fpath, "result05/Up_Gene/", tname[i], ".up"))
    write(rownames(tdds)[DownDEix], file=paste0(fpath, "result05/Down_Gene/", tname[i], ".down"))
    write.table(DEres05, quote = F, row.names = F, sep="\t", 
                file=paste0(fpath, "result05/res/", tname[i], ".res"))

    print(paste("No.", i, fname[i], "is over at", Sys.time()))
}

write.table(genesummary, file=paste0(fpath, "DEsummary.tab"), 
            quote=F, row.names = F, sep="\t")
write.table(tsam_size, file=paste0(fpath, "sample_size.tab"), 
            quote=F, row.names = F, sep="\t")
detach("package:DESeq2")
detach("package:sva")
sink()
