library(dplyr)
library(edgeR)
library(ggplot2)
library(RColorBrewer) 
library(gplots)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(pheatmap)



## Preparing Count Data
mypath <- "/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/STAR+HTSEQ+DESEQ2/input/Q5/"

filenames=list.files(path=mypath)
names <- substr(filenames, 1,8) 
lista <- list(names)

for (i in names){
  filepath <-file.path(mypath, paste(i, ".Q5.txt", sep=""))
  assign(i,read.csv(filepath, header = FALSE, sep ="", stringsAsFactors = FALSE))
}

DF <- cbind(`18NR0293`, `18NR0295`$V2, `18NR0297`$V2, `18NR0299`$V2,`18NR0301`$V2, `18NR0303`$V2, 
            `18NR0305`$V2, `18NR0307`$V2, `18NR0309`$V2, `18NR0311`$V2, `18NR0313`$V2, `18NR0315`$V2)
DF <- DF[-c(58885:58889),]
setwd("/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/STAR+HTSEQ+DESEQ2/input/")
targets <- read.csv(file = "targets2.csv", row.names = 1)
rownames(DF) <- DF[,1]
DF <- dplyr::select(DF, -V1)
colnames(DF) <- rownames(targets)

#GSEA
write.csv(DF, file = "star+htseq+desq2_gsea.csv")

# Create coundata and coldata
# targets2.csv, file with rows=samples, and columns=sample information
countdata <- as.matrix(read.csv("star+htseq+desq2_gsea.csv", row.names = 1, check.names = F))
coldata <- read.csv(file = "targets2.csv", row.names = 1)
# countdata2 annotated with gene symbols in order to annotate the heatmap
countdata2 <- countdata
countdata2 <- cbind(ensembl=rownames(countdata2), countdata2)
countdata2 <- as.data.frame(countdata2)
countdata2$genes <- mapIds(org.Hs.eg.db, keys = as.vector(countdata2$ensembl), 
                                 keytype = "ENSEMBL", column="SYMBOL")
countdata2$ID <- mapIds(org.Hs.eg.db, keys = as.vector(countdata2$ensembl), 
                           keytype = "ENSEMBL", column="ENTREZID")
# how many genes are annotated
sum(!is.na(countdata2$ID))

# Check all samples IDs are the same in colData and countData and in the same order
all(rownames(coldata) %in% colnames(countdata)) # must be TRUE
all(rownames(coldata) == colnames(countdata)) # must be TRUE
# Create data set
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, 
                              design = ~ clon + tratamiento + clon:tratamiento)
dds$group <- factor(paste0(dds$clon, dds$tratamiento))
design(dds) <- ~ group
dds <- DESeq(dds)

# Results by group contrasts
resultsNames(dds)

results1 <- results(dds, contrast = c("group", "sen0h", "sen6h"))
# LFC shrinkage it looks at the largest fold changes that are not due to low counts 
# and uses these to inform a prior distribution. So the large fold changes from genes 
# with lots of statistical information are not shrunk, while the imprecise fold changes are shrunk.
#This allows you to compare all estimated LFC across experiments, for example, which is not really feasible without the use of a prior.
results1 <- lfcShrink(dds, contrast = c("group", "sen0h", "sen6h"), res = results1)
results1ordered <- results1[order(results1$padj),]
mcols(results1ordered)$description
write.csv(as.data.frame(results1ordered), file="star+htseq+desq2_sen0hvssen6h.csv")
## Venn's diagram
# topadj1 <- read.csv(file="star+htseq+desq2_sen0hvssen6h.csv")
# topadj1 <- topadj1 %>% filter(padj<0.05)
# topadj1 <- topadj1 %>% dplyr::select(X, log2FoldChange, padj)
# colnames(topadj1) <- c("ID", "logFC", "adj.P.Val")
# topadj1$ID <- mapIds(org.Hs.eg.db, keys = as.vector(topadj1$ID), 
#                            keytype = "ENSEMBL", column="SYMBOL")
# write.csv(topadj1, file = "star_sen0sen6_venn.csv")

results2 <- results(dds, contrast = c("group", "res0h", "res6h"))
results2 <- lfcShrink(dds, contrast = c("group", "res0h", "res6h"), res = results2)
results2ordered <- results2[order(results2$padj),]
mcols(results2ordered)$description
write.csv(as.data.frame(results2ordered), file="star+htseq+desq2_res0hvsres6h.csv")

results3 <- results(dds, contrast = c("group", "sen0h", "res0h"))
results3 <- lfcShrink(dds, contrast = c("group", "sen0h", "res0h"), res = results3)
results3ordered <- results3[order(results3$padj),]
mcols(results3ordered)$description
write.csv(as.data.frame(results3ordered), file="star+htseq+desq2_sen0hvsres0h.csv")
## Venn's diagram
# topadj3 <- read.csv(file="star+htseq+desq2_sen0hvsres0h.csv")
# topadj3 <- topadj3 %>% filter(padj<0.05)
# topadj3 <- topadj3 %>% dplyr::select(X, log2FoldChange, padj)
# colnames(topadj3) <- c("ID", "logFC", "adj.P.Val")
# topadj3$ID <- mapIds(org.Hs.eg.db, keys = as.vector(topadj3$ID), 
#                      keytype = "ENSEMBL", column="SYMBOL")
# write.csv(topadj3, file = "star_sen0res0_venn.csv")

results4 <- results(dds, contrast = c("group", "sen6h", "res6h"))
results4 <- lfcShrink(dds, contrast = c("group", "sen6h", "res6h"), res = results4)
results4ordered <- results4[order(results4$padj),]
mcols(results4ordered)$description
write.csv(as.data.frame(results4ordered), file="star+htseq+desq2_sen6hvsres6h.csv")

# plotMA
dev.off()
dev.new()
par(mfrow=c(2,2))
plotMA(results1ordered, main = "sen0h_vs_sen6h", ylim = c(-2, 2))
plotMA(results2ordered, main = "res0h_vs_res6h", ylim = c(-2, 2))
plotMA(results3ordered, main = "sen0h_vs_res0h", ylim = c(-2, 2))
plotMA(results4ordered, main = "sen6h_vs_res6h", ylim = c(-2, 2))


# Visualization and clustering
# Previous transformation of count data (Variance stabilizing transformation)
vsd <- vst(dds, blind=FALSE)
# Heatmap
dev.off()
dev.new()
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
rv <- order(apply(assay(vsd), MARGIN=1, FUN=var, na.rm=TRUE),decreasing = T)[1:20]
df <- as.data.frame(coldata)
pheatmap(assay(vsd)[rv,], cluster_rows=T, show_rownames=F, main = "STAR + HTSeq + DESeq2",
         cluster_cols=T, annotation_col=df, cutree_cols = 2, annotation_names_col = F,
         annotation_colors = mycolors)

# Sample-to-sample distances
dev.off()
dev.new()
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$tratamiento, vsd$clon, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, main = "STAR + HTSeq + DESeq2",
         col=colors)

# MDS (multidimensional scaling plot)
dev.off()
dev.new()
x <- plotPCA(vsd, intgroup=c("clon", "tratamiento"), returnData = T)
percentVar <- round(100 * attr(x, "percentVar"))
rownames(sampleDistMatrix) <- paste(vsd$clon, vsd$tratamiento, sep="-")
colnames(sampleDistMatrix) <- NULL
mds <- data.frame(cmdscale(sampleDistMatrix))
rld <- rlog(dds)
mds <- cbind(mds, as.data.frame(colData(rld)))
ggplot(mds, aes(X1, X2, color=clon, shape=tratamiento, label = rownames(df))) +  geom_point(size=3) + geom_text(size = 2.8, vjust=-1, nudge_x = -0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("STAR + HTSeq + DESeq2")
coord_fixed()

## Heatmap top 100 genes sen0-sen6 
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
results1ordered_heatmap <- results1ordered[1:100,]
pheatmap(assay(vsd)[rownames(results1ordered_heatmap),c("18NR0293","18NR0301","18NR0309","18NR0297","18NR0305","18NR0313")],
         scale = "row", labels_row = countdata2$genes, annotation_col = as.data.frame(coldata),
         treeheight_col = 20,treeheight_row = 20, annotation_names_col = F, fontsize = 8,
         fontsize_row = 2.5, cellwidth = 15,cutree_cols = 2,
         annotation_colors = mycolors,
         main = "STAR + HTSeq + DESeq2\nSEN0 vs SEN6\ntop 100 genes\n")

# dev.off()
# dev.new()
# mypalette <- brewer.pal(11,"RdYlBu")
# morecols <- colorRampPalette(mypalette)
# results1ordered_heatmap <- results1ordered[1:100,]
# design <- model.matrix(~0 + dds$group)
# rownames(design) <- rownames(coldata)
# design_df <- as.data.frame(design)[c("18NR0293","18NR0301","18NR0309","18NR0297","18NR0305","18NR0313"),
#                                    3:4]
# col1 <- factor(as.factor(design_df$`dds$groupsen0h`), labels = c(NA, "#FF99FF"))
# col2 <- factor(as.factor(design_df$`dds$groupsen6h`), labels = c(NA, "#990099"))
# col3 <- coalesce(as.vector(col1), as.vector(col2))
# lmat = rbind(c(5,4), c(0,1), c(3,2))
# lhei <- c(1,0.2,5)
# lwid <- c(0.5, 1)
# heatmap.2(assay(vsd)[rownames(results1ordered_heatmap),c("18NR0293","18NR0301","18NR0309","18NR0297","18NR0305","18NR0313")], 
#           col=rev(morecols(100)), ColSideColors = col3, trace='none',scale='row',
#           main = "STAR+HTSeq+DESeq2",
#           dendrogram = "column", labRow = countdata2$genes,density.info = "histogram",key.title = NA,keysize = 0.5,
#           cexRow = 0.6, cexCol = 1, revC = T, na.rm = T,
#           margins = c(4,4), srtCol = 0, adjCol = c(0.5,0.5),
#           lmat = lmat, lhei = lhei, lwid = lwid,
#           reorderfun=function(d, w) reorder(d, w, agglo.FUN = var))
# legend("left", legend = c("sensible-0h", "sensible-6h"), fill = c("#FF99FF","#990099"), 
#        border = c("#FF99FF","#990099"), box.col = "#CCCCCC",
#        cex = 0.8, inset = .01, xjust = 0.5)
# 
# # Heatmap top 100 genes sen0-res0
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
results3ordered_heatmap <- results3ordered[1:100,]
pheatmap(assay(vsd)[rownames(results3ordered_heatmap),c("18NR0293","18NR0301","18NR0309","18NR0295","18NR0303","18NR0311")],
         scale = "row", labels_row = countdata2$genes, annotation_col = as.data.frame(coldata),
         treeheight_col = 20,treeheight_row = 20, annotation_names_col = F, fontsize = 8,
         fontsize_row = 2.5, cellwidth = 15,cutree_cols = 2,
         annotation_colors = mycolors,
         main = "STAR + HTSeq + DESeq2\nSEN0 vs RES0\ntop 100 genes\n")

# dev.off()
# dev.new()
# mypalette <- brewer.pal(11,"RdYlBu")
# morecols <- colorRampPalette(mypalette)
# results3ordered_heatmap <- results3ordered[1:100,]
# design <- model.matrix(~0 + dds$group)
# rownames(design) <- rownames(coldata)
# design_df2 <- as.data.frame(design)[c("18NR0293","18NR0301","18NR0309","18NR0295","18NR0303","18NR0311"),
#                                     c(1,3)]
# lmat = rbind(c(5,4), c(0,1), c(3,2))
# lhei <- c(1,0.2,5)
# lwid <- c(0.5, 1)
# col4 <- factor(as.factor(design_df2$`dds$groupsen0h`), labels = c(NA, "#99FF33"))
# col5 <- factor(as.factor(design_df2$`dds$groupres0h`), labels = c(NA, "#336600"))
# col6 <- coalesce(as.vector(col4), as.vector(col5))
# heatmap.2(assay(vsd)[rownames(results3ordered_heatmap),c("18NR0293","18NR0301","18NR0309","18NR0295","18NR0303","18NR0311")], 
#           col=rev(morecols(100)), ColSideColors = col6, trace='none',scale='row', 
#           main = "STAR+HTSeq+DESeq2",
#           dendrogram = "column", labRow = countdata2$genes,density.info = "histogram",key.title = NA,keysize = 0.5,
#           cexRow = 0.6, cexCol = 1, revC = T, na.rm = T,
#           margins = c(4,4), srtCol = 0, adjCol = c(0.5,0.5),
#           lmat = lmat, lhei = lhei, lwid = lwid,
#           reorderfun=function(d, w) reorder(d, w, agglo.FUN = var))
# legend("left", legend = c("sensible-0h", "resistente-0h"), fill = c("#99FF33","#336600"),
#        border = c("#99FF33","#336600"), box.col = "#CCCCCC",
#        cex = 0.8, inset = .01, xjust = 0.5)
# 
# 
