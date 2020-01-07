library(Rsubread)
library(edgeR)
library(ggplot2)
library(RColorBrewer) 
library(gplots)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(pheatmap)


## Summarizing ouput from HISAT2
setwd("/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/HISAT2+FC+DESEQ2/output/BAM/")
alignments = list.files(path = "/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/HISAT2+FC+DESEQ2/output/BAM/",
                        pattern = ".bam$", recursive = T)
fc <- featureCounts(files = alignments , annot.inbuilt="hg38")

## DGEList
dge <- DGEList(counts = fc$counts, genes = fc$annotation[c("GeneID", "Length")])

## Filtering (keeps genes that have at least min.count reads in a worthwhile number samples)
keep <- filterByExpr(dge)
dge_filt <- dge[keep,]

## Normalization
dge <- calcNormFactors(dge_filt)
colnames(dge) <- c("18NR0293", "18NR0295", "18NR0297", "18NR0299", "18NR0301", "18NR0303", 
                   "18NR0305", "18NR0307", "18NR0309", "18NR0311", "18NR0313", "18NR0315")

## Annotation with Gene Symbol
dge$genes$GeneSymbol <- mapIds(org.Hs.eg.db, keys = as.vector(as.factor(dge$genes$GeneID)), 
                               keytype = "ENTREZID", column="SYMBOL")

## GSEA
write.csv(dge$counts, file = "../dge_hisat2fc_gsea.csv")

 ## Create coundata and coldata
# targets2.csv, file with rows=samples, and columns=sample information
write.csv(cbind(dge$genes$GeneSymbol, dge$counts), file = "../dge_hisat2fc.csv", row.names = F)
countdata <- as.matrix(read.csv("../dge_hisat2fc_gsea.csv", row.names = 1, check.names = F))
coldata <- read.csv(file = "/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/HISAT2+FC+DESEQ2/input/targets2.csv", 
                    row.names = 1)
# countdata2 annotated with gene symbols in order to annotate the heatmap
countdata2 <- countdata
countdata2 <- cbind(geneid=rownames(countdata2), countdata2)
countdata2 <- as.data.frame(countdata2)
countdata2$genes <- mapIds(org.Hs.eg.db, keys = as.vector(countdata2$geneid), 
                           keytype = "ENTREZID", column="SYMBOL")

## Check all samples IDs are the same in colData and countData and in the same order
all(rownames(coldata) %in% colnames(countdata)) # must be TRUE
all(rownames(coldata) == colnames(countdata)) # must be TRUE

## Create data set
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, 
                              design = ~ clon + tratamiento + clon:tratamiento)
dds$group <- factor(paste0(dds$clon, dds$tratamiento))
design(dds) <- ~ group
dds <- DESeq(dds)

## Results by group contrasts
resultsNames(dds)

results1 <- results(dds, contrast = c("group", "sen0h", "sen6h"))
# LFC shrinkage it looks at the largest fold changes that are not due to low counts 
# and uses these to inform a prior distribution. So the large fold changes from genes 
# with lots of statistical information are not shrunk, while the imprecise fold changes are shrunk.
#This allows you to compare all estimated LFC across experiments, for example, 
#which is not really feasible without the use of a prior.
results1 <- lfcShrink(dds, contrast = c("group", "sen0h", "sen6h"), res = results1)
results1ordered <- results1[order(results1$padj),]
mcols(results1ordered)$description
write.csv(as.data.frame(results1ordered), file="hisat2_fc_deseq2_sen0hvssen6h.csv")
# # Venn's diagram
# topadj1 <- read.csv(file="hisat2_fc_deseq2_sen0hvssen6h.csv")
# topadj1 <- topadj1 %>% filter(padj<0.05)
# topadj1 <- topadj1 %>% dplyr::select(X, log2FoldChange, padj)
# colnames(topadj1) <- c("ID", "logFC", "adj.P.Val")
# topadj1$ID <- mapIds(org.Hs.eg.db, keys = as.vector(as.character(topadj1$ID)), 
#                      keytype = "ENTREZID", column="SYMBOL")
# write.csv(topadj1, file = "hisat2_sen0sen6_venn.csv")

#############################################################################
# Preranked list of genes to GSEA
#results1ordered <- as.data.frame(results1ordered)
#results1prerank <- cbind(rownames(results1ordered), results1ordered$pvalue)
#results1prerank <- as.data.frame(results1prerank)
#results1filtered <- dplyr::filter(results1ordered, !grepl("MSTRG*",V1))
#colnames(results1prerank) <- c("gene_id", "pvalue")
#results1prerank$symbol <- mapIds(org.Hs.eg.db, as.vector(results1prerank$gene_id), 
                                 #keytype = "ENSEMBL", column="SYMBOL")
#prerank1 <- select(results1filtered, symbol, pvalue)
#colnames(prerankedlist1) <- NULL
#rownames(prerankedlist1) <- NULL
#write.csv(prerankedlist1,file = "prerank_sen0sen6.csv", row.names = F)
#############################################################################

results2 <- results(dds, contrast = c("group", "res0h", "res6h"))
results2 <- lfcShrink(dds, contrast = c("group", "res0h", "res6h"), res = results2)
results2ordered <- results2[order(results2$padj),]
mcols(results2ordered)$description
write.csv(as.data.frame(results2ordered), file="hisat2_fc_deseq2_res0hvsres6h.csv")

results3 <- results(dds, contrast = c("group", "sen0h", "res0h"))
results3 <- lfcShrink(dds, contrast = c("group", "sen0h", "res0h"), res = results3)
results3ordered <- results3[order(results3$padj),]
mcols(results3ordered)$description
write.csv(as.data.frame(results3ordered), file="hisat2_fc_deseq2_sen0hvsres0h.csv")
# # Venn's diagram
# topadj3 <- read.csv(file="hisat2_fc_deseq2_sen0hvsres0h.csv")
# topadj3 <- topadj3 %>% filter(padj<0.05)
# topadj3 <- topadj3 %>% dplyr::select(X, log2FoldChange, padj)
# colnames(topadj3) <- c("ID", "logFC", "adj.P.Val")
# topadj3$ID <- mapIds(org.Hs.eg.db, keys = as.vector(as.character(topadj3$ID)), 
#                      keytype = "ENTREZID", column="SYMBOL")
# write.csv(topadj3, file = "hisat2_sen0res0_venn.csv")

results4 <- results(dds, contrast = c("group", "sen6h", "res6h"))
results4 <- lfcShrink(dds, contrast = c("group", "sen6h", "res6h"), res = results4)
results4ordered <- results4[order(results4$padj),]
mcols(results4ordered)$description
write.csv(as.data.frame(results4ordered), file="hisat2_fc_deseq2_sen6hvsres6h.csv")

## plotMA
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
pheatmap(assay(vsd)[rv,], cluster_rows=T, show_rownames=F, main = "HISAT2 + featureCounts + DESeq2",
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
pheatmap(sampleDistMatrix,  main = "HISAT2 + featureCounts + DESeq2",
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
ggplot(mds, aes(X1, X2, color=clon, shape=tratamiento, label = rownames(df))) + geom_point(size=3) + geom_text(size = 2.8, vjust=-1, nudge_x = -0.5) + ggtitle("HISAT2 + featureCounts + DESeq2") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
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
         main = "HISAT2 + featureCounts + DESeq2\nSEN0 vs SEN6\ntop 100 genes\n")



# ## Heatmap top 100 genes sen0-res0
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
results3ordered_heatmap <- results3ordered[1:100,]
pheatmap(assay(vsd)[rownames(results3ordered_heatmap),c("18NR0293","18NR0301","18NR0309","18NR0295","18NR0303","18NR0311")],
         scale = "row", labels_row = countdata2$genes, annotation_col = as.data.frame(coldata),
         treeheight_col = 20,treeheight_row = 20, annotation_names_col = F, fontsize = 8,
         fontsize_row = 2.5, cellwidth = 15,cutree_cols = 2,
         annotation_colors = mycolors,
         main = "HISAT2 + featureCounts + DESeq2\nSEN0 vs RES0\ntop 100 genes\n")

## Heatmap TOLLPATHWAY sen6-sen0
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
toll <- read.csv("/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/HISAT2+FC+DESEQ2/output/gsea/SEN-6_versus_SEN-0_hisat2.Gsea.1572283209475/TOLLPATHWAY.csv")
toll <- toll %>% dplyr::select(PROBE)
toll$ID <-  mapIds(org.Hs.eg.db, keys = as.vector(as.factor(toll$PROBE)), 
                   keytype = "SYMBOL", column="ENTREZID")
pheatmap(assay(vsd)[rownames(vsd) %in% toll$ID,c("18NR0293","18NR0301","18NR0309","18NR0297","18NR0305","18NR0313")],
         scale = "row", labels_row = toll$PROBE, color = rev(morecols(100)),  annotation_col = as.data.frame(coldata),
         annotation_colors = mycolors, border_color = NA,
         treeheight_col = 20,treeheight_row = 20, annotation_names_col = F, fontsize = 8,
         fontsize_row = 6, cellwidth = 15, cutree_rows = 2, cutree_cols = 2,
         main = "TOLLPATHWAY\nSEN6 vs SEN0")

## Heatmap NFKBPATHWAYCANONIC sen6-sen0
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
nfkb <- read.csv("/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/HISAT2+FC+DESEQ2/output/gsea/SEN-6_versus_SEN-0_hisat2.Gsea.1572283209475/NFKBPATHWAYCANONIC.csv")
nfkb <- nfkb %>% dplyr::select(PROBE)
nfkb$ID <-  mapIds(org.Hs.eg.db, keys = as.vector(as.factor(nfkb$PROBE)), 
                   keytype = "SYMBOL", column="ENTREZID")
pheatmap(assay(vsd)[rownames(vsd) %in% nfkb$ID,c("18NR0293","18NR0301","18NR0309","18NR0297","18NR0305","18NR0313")],
         scale = "row", labels_row = nfkb$PROBE, annotation_col = as.data.frame(coldata),
         annotation_colors = mycolors, border_color = NA,
         treeheight_col = 20,treeheight_row = 20, annotation_names_col = F, fontsize = 8,
         fontsize_row = 6, cellwidth = 15, cutree_rows = 2, cutree_cols = 2,
         main = "NFKBPATHWAYCANONIC\nSEN6 vs SEN0")

## Heatmap IL1RPATHWAY sen6-sen0
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
il1r <- read.csv("/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/HISAT2+FC+DESEQ2/output/gsea/SEN-6_versus_SEN-0_hisat2.Gsea.1572283209475/IL1RPATHWAY.csv")
il1r <- il1r %>% dplyr::select(PROBE)
il1r$ID <-  mapIds(org.Hs.eg.db, keys = as.vector(as.factor(il1r$PROBE)), 
                   keytype = "SYMBOL", column="ENTREZID")
pheatmap(assay(vsd)[rownames(vsd) %in% il1r$ID,c("18NR0293","18NR0301","18NR0309","18NR0297","18NR0305","18NR0313")],
         scale = "row", labels_row = il1r$PROBE, annotation_col = as.data.frame(coldata),
         annotation_colors = mycolors, border_color = NA,
         treeheight_col = 20,treeheight_row = 20, annotation_names_col = F, fontsize = 8,
         fontsize_row = 6, cellwidth = 15, cutree_rows = 2, cutree_cols = 2,
         main = "IL1RPATHWAY\nSEN6 vs SEN0")

## Heatmap GCBCELL sen0-sen6
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
gcbcell <- read.csv("/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/HISAT2+FC+DESEQ2/output/gsea/SEN-0_versus_SEN-6_hisat2.Gsea.1572283028419/GC BCELL.csv")
gcbcell <- gcbcell %>% dplyr::select(PROBE)
gcbcell$ID <-  mapIds(org.Hs.eg.db, keys = as.vector(as.factor(gcbcell$PROBE)), 
                   keytype = "SYMBOL", column="ENTREZID")
pheatmap(assay(vsd)[rownames(vsd) %in% gcbcell$ID,c("18NR0293","18NR0301","18NR0309","18NR0297","18NR0305","18NR0313")],
         scale = "row", labels_row = gcbcell$PROBE, annotation_col = as.data.frame(coldata),
         annotation_colors = mycolors, border_color = NA,
         treeheight_col = 20,treeheight_row = 20, annotation_names_col = F, fontsize = 8,
         fontsize_row = 1.5, cellwidth = 15, cutree_rows = 2, cutree_cols = 2,
         main = "GCBCELL\nSEN0 vs SEN6")

## Heatmap BLOOD PAN-BCELL sen0-sen6
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
blood <- read.csv("/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/HISAT2+FC+DESEQ2/output/gsea/SEN-0_versus_SEN-6_hisat2.Gsea.1572283028419/GC BCELL.csv")
blood <- blood %>% dplyr::select(PROBE)
blood$ID <-  mapIds(org.Hs.eg.db, keys = as.vector(as.factor(blood$PROBE)), 
                      keytype = "SYMBOL", column="ENTREZID")
pheatmap(assay(vsd)[rownames(vsd) %in% blood$ID,c("18NR0293","18NR0301","18NR0309","18NR0297","18NR0305","18NR0313")],
         scale = "row", labels_row = blood$PROBE, annotation_col = as.data.frame(coldata),
         annotation_colors = mycolors, border_color = NA,
         treeheight_col = 20,treeheight_row = 20, annotation_names_col = F, fontsize = 8,
         fontsize_row = 1.4, cellwidth = 15, cutree_rows = 2, cutree_cols = 2,
         main = "BLOOD PAN-BCELL\nSEN0 vs SEN6")

## Heatmap MHC-II sen0-res0
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
mhc2 <- read.csv("/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/HISAT2+FC+DESEQ2/output/gsea/SEN-0_versus_RES-0_hisat2.Gsea.1572283379002/MHC-II.csv")
mhc2 <- mhc2 %>% dplyr::select(PROBE)
mhc2$ID <-  mapIds(org.Hs.eg.db, keys = as.vector(as.factor(mhc2$PROBE)), 
                    keytype = "SYMBOL", column="ENTREZID")
pheatmap(assay(vsd)[rownames(vsd) %in% mhc2$ID,c("18NR0293","18NR0301","18NR0309","18NR0295","18NR0303","18NR0311")],
         scale = "row", labels_row = mhc2$PROBE, annotation_col = as.data.frame(coldata),
         annotation_colors = mycolors, border_color = NA,
         treeheight_col = 20,treeheight_row = 20, annotation_names_col = F, fontsize = 8,
         fontsize_row = 6, cellwidth = 15, cellheight = 15, cutree_rows = 2, cutree_cols = 2,
         main = "MHC-II\nSEN0 vs RES0")

## Heatmap MHC-I sen0-res0
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
mhc1 <- read.csv("/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/HISAT2+FC+DESEQ2/output/gsea/SEN-0_versus_RES-0_hisat2.Gsea.1572283379002/MHC-I.csv")
mhc1 <- mhc1 %>% dplyr::select(PROBE)
mhc1$ID <-  mapIds(org.Hs.eg.db, keys = as.vector(as.factor(mhc1$PROBE)), 
                   keytype = "SYMBOL", column="ENTREZID")
pheatmap(assay(vsd)[rownames(vsd) %in% mhc1$ID,c("18NR0293","18NR0301","18NR0309","18NR0295","18NR0303","18NR0311")],
         scale = "row", labels_row = mhc1$PROBE, annotation_col = as.data.frame(coldata),
         annotation_colors = mycolors, border_color = NA,
         treeheight_col = 20,treeheight_row = 20, annotation_names_col = F, fontsize = 8,
         fontsize_row = 6, cellwidth = 15, cellheight = 15, cutree_rows = 2, cutree_cols = 2,
         main = "MHC-I\nSEN0 vs RES0")

## Heatmap STROMALCCG sen0-res0
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
strom <- read.csv("/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/HISAT2+FC+DESEQ2/output/gsea/SEN-0_versus_RES-0_hisat2.Gsea.1572283379002/STROMALCCG.csv")
strom <- strom %>% dplyr::select(PROBE)
strom$ID <-  mapIds(org.Hs.eg.db, keys = as.vector(as.factor(strom$PROBE)), 
                   keytype = "SYMBOL", column="ENTREZID")
pheatmap(assay(vsd)[rownames(vsd) %in% strom$ID,c("18NR0293","18NR0301","18NR0309","18NR0295","18NR0303","18NR0311")],
         scale = "row", labels_row = strom$PROBE, annotation_col = as.data.frame(coldata),
         annotation_colors = mycolors, border_color = NA,
         treeheight_col = 20,treeheight_row = 20, annotation_names_col = F, fontsize = 8,
         fontsize_row = 6, cellwidth = 15, cellheight = 15, cutree_rows = 2, cutree_cols = 2,
         main = "STROMALCCG\nSEN0 vs RES0")

## Heatmap BCL-6 TARGETS sen0-res0
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
bcl6 <- read.csv("/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/HISAT2+FC+DESEQ2/output/gsea/SEN-0_versus_RES-0_hisat2.Gsea.1572283379002/BCL-6 TARGETS.csv")
bcl6 <- bcl6 %>% dplyr::select(PROBE)
bcl6$ID <-  mapIds(org.Hs.eg.db, keys = as.vector(as.factor(bcl6$PROBE)), 
                    keytype = "SYMBOL", column="ENTREZID")
pheatmap(assay(vsd)[rownames(vsd) %in% bcl6$ID,c("18NR0293","18NR0301","18NR0309","18NR0295","18NR0303","18NR0311")],
         scale = "row", labels_row = bcl6$PROBE, annotation_col = as.data.frame(coldata),
         annotation_colors = mycolors, border_color = NA,
         treeheight_col = 20,treeheight_row = 20, annotation_names_col = F, fontsize = 8,
         fontsize_row = 6, cellwidth = 15, cutree_rows = 2,
         main = "BCL-6 TARGETS\nSEN0 vs RES0")
