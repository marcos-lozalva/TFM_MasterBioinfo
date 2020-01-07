## Installation

library("dplyr")

# BiocManager::install("DESeq2")
library("DESeq2")

# BiocManager::install("ReportingTools")
library("ReportingTools")

# BiocManager::install("ggplot2")
library("ggplot2")

# BiocManager::install("tximport")
library("tximport")

# BiocManager::install("tximportData")
library("tximportData")

# BiocManager::install("readr")
library("readr")

# BiocManager::install("biomaRt")
library("biomaRt")

# BiocManager::install("AnnotationDbi")
library("AnnotationDbi")

library(pheatmap)
library(RColorBrewer) 



## Set working directory
setwd("/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/HISAT2+STRINGTIE+DESEQ2/output/limpios/")

## Import otuput from Stringtie to Deseq2
#1. Using tximport()
# targets2.csv, file with rows=samples, and columns=sample information
dir <- system.file("extdata", package = "tximportData")
list.files(dir)
files <- list.files(pattern = "t_data*", recursive = T)
tmp <- read.delim(files[1])
tx2gene <- tmp[, c("t_name", "gene_id")]
txi <- tximport(files=files, type = "stringtie", tx2gene = tx2gene)
coldata <- read.csv(file = "targets2.csv", row.names = 1)
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ clon + tratamiento + clon:tratamiento)
dds$group <- factor(paste0(dds$clon, dds$tratamiento))
design(dds) <- ~ group
dds <- DESeq(dds)

## Export count matrix to GSEA (no necesario, se puede hacer el .gct a partir del gene_count_matrix.csv)
# With txi$counts
matrixgsea <- txi$counts
matrixgsea_log <- log2(txi$counts)
colnames(matrixgsea) <- c("18NR0293", "18NR0295", "18NR0297", "18NR0299", "18NR0301", "18NR0303", "18NR0305",
                          "18NR0307", "18NR0309", "18NR0311", "18NR0313", "18NR0315")
write.csv(matrixgsea, file = "txicounts_gsea.csv")

# With txi$abundance
matrixgsea2 <- txi$abundance
colnames(matrixgsea2) <- c("18NR0293", "18NR0295", "18NR0297", "18NR0299", "18NR0301", "18NR0303", "18NR0305",
                          "18NR0307", "18NR0309", "18NR0311", "18NR0313", "18NR0315")
write.csv(matrixgsea2, file = "txiabundance_gsea.csv")

# With dds 
matrixgsea3 <- assays(dds)$counts
write.csv(matrixgsea3, file = "ddscounts_gsea.csv")

#2. Using python prepDE.py
# Create coundata and coldata
# targets2.csv, file with rows=samples, and columns=sample information
countdata <- as.matrix(read.csv("gene_count_matrix.csv", row.names = "gene_id", check.names = F))
coldata <- read.csv(file = "targets2.csv", row.names = 1)
# Check all samples IDs are the same in colData and countData and in the same order
all(rownames(coldata) %in% colnames(countdata)) # must be TRUE
all(rownames(coldata) == colnames(countdata)) # must be TRUE
# Create data set
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ clon + tratamiento + clon:tratamiento)
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
write.csv(as.data.frame(results1ordered), file="deseq2_sen0hvssen6h.csv")


# Preranked list of genes to GSEA
# library("org.Hs.eg.db")
# results1ordered <- as.data.frame(results1ordered)
# results1prerank <- cbind(rownames(results1ordered), results1ordered$pvalue)
# results1prerank <- as.data.frame(results1prerank)
# results1filtered <- dplyr::filter(results1ordered, !grepl("MSTRG*",V1))
# colnames(results1prerank) <- c("gene_id", "pvalue")
# results1prerank$symbol <- mapIds(org.Hs.eg.db, as.vector(results1prerank$gene_id), 
#                                  keytype = "ENSEMBL", column="SYMBOL")

prerank1 <- select(results1filtered, symbol, pvalue)
colnames(prerankedlist1) <- NULL
rownames(prerankedlist1) <- NULL
write.csv(prerankedlist1,file = "prerank_sen0sen6.csv", row.names = F)

results2 <- results(dds, contrast = c("group", "res0h", "res6h"))
results2 <- lfcShrink(dds, contrast = c("group", "res0h", "res6h"), res = results2)
results2ordered <- results2[order(results2$padj),]
mcols(results2ordered)$description
write.csv(as.data.frame(results2ordered), file="deseq2_res0hvsres6h.csv")

results3 <- results(dds, contrast = c("group", "sen0h", "res0h"))
results3 <- lfcShrink(dds, contrast = c("group", "sen0h", "res0h"), res = results3)
results3ordered <- results3[order(results3$padj),]
mcols(results3ordered)$description
write.csv(as.data.frame(results3ordered), file="deseq2_sen0hvsres0h.csv")

results4 <- results(dds, contrast = c("group", "sen6h", "res6h"))
results4 <- lfcShrink(dds, contrast = c("group", "sen6h", "res6h"), res = results4)
results4ordered <- results4[order(results4$padj),]
mcols(results4ordered)$description
write.csv(as.data.frame(results4ordered), file="deseq2_sen6hvsres6h.csv")

# plotMA
dev.off()
dev.new()
par(mfrow=c(2,2))
plotMA(results1ordered, main = "sen0h_vs_sen6h", ylim = c(-2, 2))
plotMA(results2ordered, main = "res0h_vs_res6h", ylim = c(-2, 2))
plotMA(results3ordered, main = "sen0h_vs_res0h", ylim = c(-2, 2))
plotMA(results4ordered, main = "sen6h_vs_res6h", ylim = c(-2, 2))


## Visualization and clustering
# Previous transformation of count data (Variance stabilizing transformation)
vsd <- vst(dds, blind=FALSE)

## Heatmap
dev.off()
dev.new()
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
rv <- order(apply(assay(vsd), MARGIN=1, FUN=var, na.rm=TRUE),decreasing = T)[1:20]
df <- as.data.frame(coldata)
colnames(df) <- c("clon", "tratamiento")
pheatmap(assay(vsd)[rv,], cluster_rows=T, show_rownames=F, main = "HISAT2 + Stringtie + DESeq2",
         cluster_cols=T, annotation_col=df, cutree_cols = 2, annotation_names_col = F,
         annotation_colors = mycolors)

## Sample-to-sample distances
dev.off()
dev.new()
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$clon, vsd$tratamiento, sep="-")
colnames(sampleDistMatrix) <- rownames(coldata)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,  main = "HISAT2 + Stringtie + DESeq2",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

## MDS (multidimensional scaling plot)
dev.off()
dev.new()
x <- plotPCA(vsd, intgroup=c("clon", "tratamiento"), returnData = T)
percentVar <- round(100 * attr(x, "percentVar"))
rownames(sampleDistMatrix) <- paste(vsd$clon, vsd$tratamiento, sep="-")
colnames(sampleDistMatrix) <- NULL
mds <- data.frame(cmdscale(sampleDistMatrix))
rld <- rlog(dds)
mds <- cbind(mds, as.data.frame(colData(rld)))
ggplot(mds, aes(X1, X2, color=clon, shape=tratamiento, label = rownames(df))) + geom_point(size=3) + geom_text(size = 2.8, vjust=-1, nudge_x = -0.5) + ggtitle("HISAT2 + Stringtie + DESeq2") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
coord_fixed()

