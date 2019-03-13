# ====================================================================================================
# Author:   Aaron
# Function: remove samples with low quality
#           remove samples which is not control or type 2 diabetes
#           split sample information into tissues
# Version:  1.0
# Date:     unknown
# ====================================================================================================

rm(list=ls())
setwd("/data/MyProgram/Final_diabrain/0.raw/")

# load data--------------------------------------------------------------------
sam_attr <- read.table("Sample_Attributes.txt", header = TRUE, sep = "\t",
                          quote = "", fill = 1, stringsAsFactors = FALSE)
sub_pheno <- read.table("Subject_Phenotypes.txt", header = TRUE, sep = "\t",
                        quote = "", fill = 1, stringsAsFactors = FALSE)
sam_in_data <- read.table("GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_metrics.tsv", 
                          header = TRUE, sep = "\t", quote = "", fill = 1, 
                          stringsAsFactors = FALSE)

# filter-------------------------------------------------------------------------
sam_attr_seq <- sam_attr[(sam_attr$SAMPID %in% sam_in_data$Sample), ]
## Confirm if the samples in data matrixs are of high quality
table(sam_attr_seq$SMAFRZE)
## exclude low quality samples
sam_attr_seq <- sam_attr_seq[sam_attr_seq$SMAFRZE == "RNASEQ", ]
table(sam_attr_seq$ANALYTE_TYPE)
table(sam_attr_seq$SMTORMVE)
table(sam_attr_seq$SMFLGRMRK)

# merge sample and subject data--------------------------------------------------
samid <- sam_attr_seq$SAMPID
subid <- sapply(strsplit(samid, "-"), function(x) {paste(c(x[1], x[2]), collapse="-")})
sam_attr_seq <- data.frame(sam_attr_seq, SUBJID=subid)
sam_attr_seq <- merge(sam_attr_seq, sub_pheno, by="SUBJID")

# filter samples with low RIN----------------------------------------------------
sam_attr_seq <- sam_attr_seq[sam_attr_seq$SMRIN >6, ]

# assign samples which is unknown or NA or both in type I&II diabetesto 99-------
lsam <- nrow(sam_attr_seq)
t1d <- sam_attr_seq$MHT1D
t2d <- sam_attr_seq$MHT2D
DIABETES <- rep(0, lsam)
for (i in 1:lsam) {
    if (is.na(t1d[i]+t2d[i])) {
        DIABETES[i] <- 99
    } else if ((t1d[i]+t2d[i])>=2) {
        DIABETES[i] <- 99
    } else if (t1d[i]==1) {
        DIABETES[i] <- 1
    } else if (t2d[i]==1) {
        DIABETES[i] <-2
    }
}
sam_attr_seq <- data.frame(sam_attr_seq, DIABETES)
ix <- ((DIABETES==0)|(DIABETES==2))&((sam_attr_seq$RACE==2)|(sam_attr_seq$RACE==3))
sam_attr_seq <- sam_attr_seq[ix,]

# split sample info into list of tissue--------------------------------------------
tissuelist <- split(sam_attr_seq, sam_attr_seq$SMTSD)
tname <- names(tissuelist)
tissuelist <- tissuelist[grep("Brain", tname)]
names(tissuelist)

l <- length(tissuelist)
tissuename <- names(tissuelist)
tissuename
tissuename <- gsub(" ", "_", gsub(" - ", "_", gsub("\\ \\(.*\\)", "", tissuename)))
tissuename
tissue_sam_count <- matrix(nrow=l, ncol=3)
colnames(tissue_sam_count) <- c("Healthy", "T2D", "Total")
rownames(tissue_sam_count) <- tissuename

dir.create("../1.clean/sam_info")
dir.create("../1.clean/samid")

#save tissue samples (total)---------------------------------------------------- 
for (i in 1:l) {
    file_name <- paste0(tissuename[i], "_", nrow(tissuelist[[i]]), ".tab")
    write.table(tissuelist[[i]], file=paste0("../1.clean/sam_info/", file_name), 
                sep="\t", quote=FALSE, row.names=FALSE)
    write(tissuelist[[i]]$SAMPID, file=paste0("../1.clean/samid/", file_name))
    diab <- tissuelist[[i]]$DIABETES
    tissue_sam_count[i, 1] <- sum(diab==0)
    tissue_sam_count[i, 2] <- sum(diab==2)
    tissue_sam_count[i, 3] <- nrow(tissuelist[[i]])
}

write.table(data.frame(tissue_sam_count, Tissue=tissuename), 
            "../1.clean/tissue_sam_count.txt", 
            row.names = F, quote=FALSE, sep="\t")

