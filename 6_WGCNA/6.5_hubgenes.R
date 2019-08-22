# ====================================================================================================
# Author:   Aaron
# Function: Hub genes
# Version:  1.0
# Date:     May 13, 2019
# ====================================================================================================

library(WGCNA)
library(igraph)
library(extrafont)
loadfonts()
rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/")
tmp <- load("/data/MyProgram/Final_diabrain/7.interdata/6_WGCNA/Brain_Caudate_covar_WGCNAres_v1.0.RData")
id2symbol <- read.table("/data/MyProgram/Final_diabrain/1.clean/id2symbol.tab", header = T)

#-------------------------------------------------------------------------------------------------
# I. for cytoscape
#-------------------------------------------------------------------------------------------------
# prepare cytoscape input edge files and node files
TOM <- TOMsimilarityFromExpr(datExpr, networkType = "signed", power = 5, nThreads = 10)
dissTOM <- 1-TOM
# save(TOM, file = "./cytoscape/TOM_caudate.Rdata")
modules <- bw_label2color$Colors
modules <- modules[modules != "grey"]

for (mymod in modules) {
    inModule <- bwnetColors == mymod
    modTOM <- TOM[inModule, inModule]
    modID <- colnames(datExpr)[inModule]
    modGene <- id2symbol$Description[match(modID, id2symbol$Name)]
    dimnames(modTOM) <- list(modID, modID)
    
    cyt = exportNetworkToCytoscape(modTOM,
                                   edgeFile = paste("./cytoscape/CytoscapeInput-edges-", mymod, ".txt", sep=""),
                                   nodeFile = paste("./cytoscape/CytoscapeInput-nodes-", mymod, ".txt", sep=""),
                                   weighted = TRUE,
                                   threshold = 0,
                                   nodeNames = modID,
                                   altNodeNames = modGene,
                                   nodeAttr = bwnetColors[inModule])
}

#-------------------------------------------------------------------------------------------------
# II. top connections plot
#-------------------------------------------------------------------------------------------------

# ...................................................................................................
# function to calculate in layout in concentric circles
# Arguments: g is an igraph object
#            group is an factor used to split node into various circles, radius increase with levels
# ...................................................................................................
layout_in_circles <- function(g, group) {
    layout <- lapply(split(V(g), group), function(x) {
        layout_in_circle(induced_subgraph(g,x))
    })
    layout <- Map(`*`, layout, seq_along(layout))
    x <- matrix(0, nrow=vcount(g), ncol=2)
    split(x, group) <- layout
    x
}

# ...................................................................................................
# function to plot the top connections
# Arguments: xmod id the mod to plot
#            ntop is the number of top connections to plot
#            prop_hub spectify the top proportion in number of connections to be considered as hub genes
# ...................................................................................................
plot_topConnections <- function(xmod, ntops = 500, prop_hub = 0.05, label.dist.cex = 1.1, vertex.color = "black", edge.color = "darkgrey") {
    modEdge <- read.table(paste0("./cytoscape/CytoscapeInput-edges-", xmod, ".txt"), header = TRUE)
    modEdge <- modEdge[order(modEdge$weight, decreasing = TRUE), ]
    # extract ntops connections for plot
    topConnec <- modEdge[1:ntops, c("fromAltName", "toAltName", "weight")]
    modnet <- graph_from_data_frame(d=topConnec, directed = F)
    
    # Identify hub genes with 5% top connectivity
    modNodesConnec <- sort(table(c(topConnec$fromAltName, topConnec$toAltName)), decreasing = T)
    nhubs <- floor(length(modNodesConnec)*prop_hub)
    hubgenes <- names(modNodesConnec)[1:nhubs]
    V(modnet)$is_hub <- factor(V(modnet)$name %in% hubgenes,  levels = c("TRUE", "FALSE"))
    
    # other attributes of vertex
    V(modnet)$Connectivity <- modNodesConnec[match(V(modnet)$name, names(modNodesConnec))]
    V(modnet)$size <- log2(V(modnet)$Connectivity)
    # V(modnet)$label <- ifelse(V(modnet)$is_hub == TRUE, V(modnet)$name, NA)
    
    # layout and label degree
    modlayout <- layout_in_circles(modnet, V(modnet)$is_hub)
    # modlayout <- layout_in_circle(modnet)
    tmp <- atan(-modlayout[, 1]/modlayout[, 2])
    modlabel.degree = ifelse(tmp < 0,  90 + tmp*180/pi, 270 + tmp*180/pi)
    
    plot(modnet, layout=modlayout, vertex.label="", vertex.color = vertex.color, vertex.frame.color = vertex.color, edge.color = edge.color)
    x = modlayout[,1]/2*label.dist.cex
    y = modlayout[,2]/2*label.dist.cex
    for (i in 1:length(x)) {
        text(x=x[i], y=y[i], labels=V(modnet)$name[i], adj=NULL, pos=NULL, cex=0.7, col="black", srt=modlabel.degree[i], xpd=T)
    }
    
    # deprecated, this way cann't rotate label text itself 
    # tmp <- atan(modlayout[, 2]/modlayout[, 1])
    # modlabel.degree <- ifelse(modlayout[, 1] > 0, -tmp, pi-tmp)
    # plot(modnet,  layout = modlayout, 
    #      vertex.label.dist = 1, vertex.label.color = "black", vertex.label.degree = modlabel.degree,
    #      vertex.color = "black", vertex.frame.color = "black", edge.color = "grey")
    
}

library(RColorBrewer)
mycolor <- brewer.pal(11, "Spectral")
pdf("./hubgenes.pdf", family = "Arial", width = 7, height = 7)
# par(ps = 7)
plot_topConnections("turquoise", ntops = 500, prop_hub = 0.05, label.dist.cex = 1.25)
plot_topConnections("turquoise", ntops = 500, prop_hub = 0.05, label.dist.cex = 1.25, vertex.color = blues9[8], edge.color = blues9[4])
plot_topConnections("turquoise", ntops = 500, prop_hub = 0.05, label.dist.cex = 1.25, vertex.color = blues9[7], edge.color = "grey")
plot_topConnections("turquoise", ntops = 500, prop_hub = 0.05, label.dist.cex = 1.25, vertex.color = mycolor[2], edge.color = "grey")
dev.off()
# embed_fonts("/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/hubgenes.pdf")



#-------------------------------------------------------------------------------------------------
# III. mini visualization plot for fig1
#-------------------------------------------------------------------------------------------------
nSelect <- 800       # random select 400 genes for visualization
nGenes <- dim(TOM)[1]
set.seed(10)
select <- sample(nGenes, size = nSelect)
selectTOM <- dissTOM[select, select]
selectTree <- hclust(as.dist(selectTOM), method = "average")
selectColors <- bwnetColors[select]
plotDiss <- selectTOM^3
diag(plotDiss) <- NA
mycol <- brewer.pal(9, "YlOrRd")[9:1] 
pdf("./mini_visualization.pdf")
TOMplot(plotDiss, selectTree, selectColors, col = mycol)
dev.off()
















# MEs <- moduleEigengenes(datExpr, bwnetColors)$eigengenes
# MEs <- MEs[, colnames(MEs)!="MEgrey"]
# MEs <- orderMEs(MEs)
# modNames <- substring(names(MEs), 3)
# 
# geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
# names(geneModuleMembership) <- paste("MM", modNames, sep = "")
# 
# turquoiseMM <- geneModuleMembership[bwnetColors == "turquoise",]
# turquoiseMM <- turquoiseMM[order(turquoiseMM$MMturquoise, decreasing = T), ]
# top1 <- rownames(turquoiseMM)[1:100]
# top1 <- id2symbol$Description[match(top1, id2symbol$Name)]
# top1
# 
# totalkIM <- intramodularConnectivity.fromExpr(datExpr, bwnetColors, networkType = "signed", scaleByMax = T)
# rownames(totalkIM) <- rownames(geneModuleMembership)
# turquoisekIM <- totalkIM[bwnetColors=="turquoise", ]
# turquoisekIM <- turquoisekIM[order(turquoisekIM$kWithin, decreasing = T), ]
# top2 <- rownames(turquoisekIM)[1:100]
# top2 <- id2symbol$Description[match(top2, id2symbol$Name)]
# top2
# 
# x <- merge(data.frame(turquoisekIM, Name = rownames(turquoisekIM)), data.frame(turquoiseMM, Name = rownames(turquoiseMM)), by = "Name")
# plot(x$kWithin, x$MMturquoise)
# cor.test(x$kWithin, x$MMturquoise)