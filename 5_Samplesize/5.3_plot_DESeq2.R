# ====================================================================================================
# Author:   Aaron
# Function: plot the sample size use DESeq2
# Version:  1.0
# Date:     May 9, 2019
# ====================================================================================================

rm(list=ls())
library(ggplot2)
library(extrafont)
loadfonts()

setwd("/data/MyProgram/Final_diabrain/7.interdata/5_Samplesize_from_server/DESeq2data/")
options(stringsAsFactors = FALSE)
# load("./run_tsamsize.RData")
fname <- dir("/data/MyProgram/Final_diabrain/1.clean/data/")
tname <- sub("Brain_", "", sub("_[0-9]+\\.tab", "", fname))
sname <- c("Amygdala", "Anterior Cingulate Cortex", "Caudate", "Cerebellar Hemisphere", "Cerebellum", "Cortex", "Frontal Cortex", 
           "Hippocampus", "Hypothalamus", "Nucleus Accumbens", "Putamen", "Spinal Cord", "Substantia Nigra")
PTIMES <- 1000

# ----------------------------------------------------------------------------------------------------
# I. sample size increase from 15-35, for plot
# ----------------------------------------------------------------------------------------------------
# Function to calculate se
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)
    # 计算长度
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }
    # 以 groupvars 为组,计算每组的长度,均值,以及标准差
    # ddply 就是 dplyr 中的 group_by + summarise
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                       c(N    = length2(xx[[col]], na.rm=na.rm),
                         mean = mean   (xx[[col]], na.rm=na.rm),
                         sd   = sd     (xx[[col]], na.rm=na.rm)
                       )
                   },
                   measurevar
    )
    # 重命名  
    datac <- plyr::rename(datac, c("mean" = measurevar))
    
    # 计算标准偏差
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    # 计算置信区间
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
}

# construct degenes array from seprated files
degenes <- array(NA, dim = c(13, 6, PTIMES), 
                 dimnames = list(sname, paste0("s", seq(from = 10, to = 35, by = 5)), c(1:PTIMES)))
flist <- dir("./")
for (i in 1:dim(degenes)[1]) {
    for (j in 1:6) {
        tmpfname <- paste("Brain", tname[i], dimnames(degenes)[[2]][j], sep = "_")
        if (tmpfname %in% flist) {
            tmp <- read.table(tmpfname)[[1]]
            degenes[i, j, 1:length(tmp)] <- tmp
        }
    }
}

# only bootstrapped 100 times accross sample size
degenes_raw <- degenes
degenes <- degenes[,,1:100]

# transform degenes array to ggplot data.frame
# copy the permutations less than 100 but not zero to 100
pre_data <- data.frame(DEs = numeric(), Tissue = character(), Sample_size = character())
l <- dim(degenes)[3]
for (i in 1:dim(degenes)[1]) {
    for (j in 1:dim(degenes)[2]) {
        tmp <- degenes[i, j, ]
        tmp <- tmp[!is.na(tmp)]
        if (length(tmp)==0) {
            next
        } else if (length(tmp)<dim(degenes)[3]) {
            tmp1 <- rep_len(tmp, length.out = l)
            tmp <- tmp1
        }
        sz <- as.integer(sub("s", "", dimnames(degenes)[[2]][j]))
        tmpdf <- data.frame(DEs = tmp, Tissue = rep(dimnames(degenes)[[1]][i], l), Sample_size = rep(paste("Ctrl", 2*sz, "\nT2D", sz), l))
        pre_data <- rbind(pre_data, tmpdf)
    }
}
plot_data <- summarySE(pre_data, measurevar = "DEs", groupvars = c("Tissue", "Sample_size"))

# plot the degenes vs sample size
pdf("../../../4.plots/Sample_size/samplesize.pdf", family = "Arial", height = 6, width = 6)
ggplot(plot_data, aes(x=Sample_size, y=DEs, colour=Tissue, group = Tissue)) + 
    geom_errorbar(aes(ymin=DEs-ci, ymax=DEs+ci), width=.1) +
    geom_line() + geom_point() +
    xlab("Sample Size") + ylab("Numbers of T2D-associated Genes") +
    theme(axis.title.x=element_text(size=8,face="bold"),
          axis.title.y=element_text(size=8,face="bold"),
          axis.text=element_text(size=8),
          axis.text.y = element_text(face = "bold"),
          legend.title=element_text(size=8,face="bold"),
          legend.text=element_text(size=8,face="bold"),
          #panel.grid.major  = element_blank(),            # remove backgroud 
          #panel.background = element_blank(),             # remove backgroud 
          # axis.line=element_line(size = .3, colour="black"),  # add axis
          strip.text.x=element_text(size=8,face="bold"))
dev.off()
# embed_fonts("/data/MyProgram/Final_diabrain/4.plots/Sample_size/samplesize.pdf")


# ----------------------------------------------------------------------------------------------------
# II. calculate means in sample size T2D 10 vs Ctrl 20
# ----------------------------------------------------------------------------------------------------
# construct another degenes array from seprated files contain more info
degenes <- array(NA, dim = c(13, 9, PTIMES), 
                 dimnames = list(sname, c("total1", "up1", "down1", "total05", "up05", "down05", "total01", "up01", "down01"), c(1:PTIMES)))
demeans <- matrix(NA, nrow = 13, ncol = 9)
rownames(demeans) <- sname

for (i in 1:dim(degenes)[1]) {
    tmpfname <- paste("Brain", tname[i], "moreinfo_s10", sep = "_")
    if (tmpfname %in% flist) {
        tmp <- as.matrix(read.table(tmpfname, row.names = 1))
        degenes[i,,] <- tmp
        demeans[i, ] <- rowMeans(tmp, na.rm = T)
    }
}
colnames(demeans) <- c("total1", "up1", "down1", "total05", "up05", "down05", "total01", "up01", "down01")
write.table(demeans, "../../../4.plots/Sample_size/sample_size.txt", sep = "\t", quote = F)


# pd <- position_dodge(0.1)
# ggplot(plot_data, aes(x=Sample_size, y=DEs, colour=Tissue, group=Tissue)) + 
#     geom_errorbar(aes(ymin=DEs-se, ymax=DEs+se), colour="black", width=.1, position=pd) +
#     geom_line(position=pd) +
#     geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
#     xlab("Sample size") +
#     ylab("Means of number of DEgenes") +
#     scale_colour_hue(name="Tissue",    # Legend label, use darker colors
#                      #breaks=c("OJ", "VC"),
#                      #labels=c("Orange juice", "Ascorbic acid"),
#                      l=40) +                    # Use darker colors, lightness=40
#     ggtitle("Ttest Sample size") +
#     expand_limits(y=0) +                        # Expand y range
#     scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
#     theme_bw() +
#     theme(legend.justification=c(1,0),# 这一项很关键,如果没有这个参数,图例会偏移,读者可以试一试
#           legend.position=c(1,0))     
