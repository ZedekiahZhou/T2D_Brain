# ====================================================================================================
# Author:   Aaron
# Function: match the sample to control co-variables
# Version:  1.0
# Date:     Jun 28, 2018
# ====================================================================================================

rm(list=ls())
library(MatchIt)
library(extrafont)
library(RColorBrewer)

fpath <- "/data/MyProgram/Final_diabrain/1.clean/sam_info/"
setwd(fpath)

sink("../../7.interdata/3_Matchit_optimal.log")
# dir.create("../matched_sam_optimal/")
# dir.create("../matched_summary/")
# dir.create("../math_matrix/")
# dir.create("../../4.plots/match_plot_optimal/", recursive = T)
# dir.create("../../4.plots/check_outliers_optimal/")

# Specify the name of each tissue
fname <- dir("./")
tname <- sub("_[0-9]+.[A-Za-z]+", "", fname)
short_name <- sub("Brain_", "", tname)
short_name <- gsub("_", " ", short_name)

loadfonts()
pdf(paste0("../../4.plots/match_plot_optimal/All_tissues.pdf"), family = "Arial", height = 15, width = 9)
par(mfcol = c(5, 3))

# -------------------------------------------------------------------------------------
# I. The function to summary co-cavariate info for matched samples
# ------------------------------------------------------------------------------------
sam_summary <- function(x) {
    res <- numeric(7)
    names(res) <- c("Male", "Female", "AGE", "White", "Black", "BMI", "RIN")
    res[1] <- table(x$SEX)["1"]
    res[2] <- table(x$SEX)["2"]
    res[3] <- mean(x$AGE)
    res[4] <- table(x$RACE)["3"]
    res[5] <- table(x$RACE)["2"]
    res[6] <- mean(x$BMI)
    res[7] <- mean(x$SMRIN)
    round(res, digits = 2)
}

#----------------------------------------------------------------------------------------------------------------
# II. match the sample for each tissue
#----------------------------------------------------------------------------------------------------------------

for (i in 1:length(fname)) {
    print(paste("No.", i, fname[i], "is begin at", Sys.time()))
    
    tsam_info <- read.table(fname[i], header = T, sep = "\t", 
                            quote = "", stringsAsFactors = F)
    tsam_info <- subset(tsam_info, 
                        select=c(SUBJID,SAMPID,SEX,AGE,RACE,BMI,SMRIN,MHT2D))
    rownames(tsam_info) <- make.names(tsam_info$SAMPID)
    
    # Use optimal method to match the sample with Ctrl/T2D = 2:1
    tout <- matchit(MHT2D ~ SEX + AGE + RACE + BMI + SMRIN, data = tsam_info, method = "optimal",  ratio = 2)
    # extract matched samples with all covariates
    tsam_match <- match.data(tout)
    
    # Get the matrix which sample matched to which sample
    matchMatrix <- cbind(row.names(tout$match.matrix), tout$match.matrix)
    colnames(matchMatrix) <- c("T2D", "Ctrl1", "Ctrl2")
    write.table(matchMatrix, file = paste0("../math_matrix/", tname[i], "_matchsam_", dim(tsam_match)[1], ".tab"), 
                quote = F, row.names = F, sep = "\t")
    write.table(tsam_match[, 1:8], file = paste0("../matched_sam_optimal/", tname[i], "_", dim(tsam_match)[1], ".tab"), 
                quote = F, row.names = F, sep="\t")
    
    print(summary(tout))
   
    ## plot the match result of each tissue
    plot(tout,type = "jitter", interactive = FALSE, sub = short_name[i], col = brewer.pal(9, "PuBu")[7])
    
    # generate the summary and statistics for covariates of matched samples
    tdataT2D <- match.data(tout, group = "treat")
    tdataCtrl <- match.data(tout, group = "control")
    match_summary <- data.frame(Control = sam_summary(tdataCtrl), T2D = sam_summary(tdataT2D), p.value = NA)
    match_summary$p.value[3] <- t.test(tdataCtrl$AGE, tdataT2D$AGE)$p.value
    match_summary$p.value[6] <- t.test(tdataCtrl$BMI, tdataT2D$BMI)$p.value 
    match_summary$p.value[7] <- t.test(tdataCtrl$SMRIN, tdataT2D$SMRIN)$p.value
    match_summary <- data.frame(Variables = row.names(match_summary), match_summary, stringsAsFactors = F)
    write.table(match_summary, file = paste0("../matched_summary/", tname[i], "_", dim(tsam_match)[1], ".tab"), 
                quote = F, row.names = F, sep = "\t", na = "")
    
    ## check outliers
    # tdata <- read.table(paste0("../data/", fname[i]), header=T, sep=" ", quote="")
    # row.names(tdata) <- tdata$Name
    # tdata <- tdata[,c(-1,-2)]
    # rownames(tsam_match) <- make.names(tsam_match$SAMPID)
    # tdata <- tdata[, rownames(tsam_match)]
    # 
    # hc <- hclust(dist(t(log2(tdata+1))), method = "ward.D2")
    # sam_label <- strsplit(colnames(tdata), "\\.")
    # sam_label <- sapply(sam_label, function(x) {paste(x[1], x[2], sep="-")})
    # pdf(paste0("../../4.plots/check_outliers_optimal/", tname[i], ".pdf"), 
    #     height = 8, width = 12)
    # plot(hc, labels = sam_label, hang=-1)
    # dev.off()
    
    print(paste("No.", i, fname[i], "is over at", Sys.time()))
}

dev.off()
# embed_fonts("/data/MyProgram/Final_diabrain/4.plots/match_plot_optimal/All_tissues.pdf")

sink()

