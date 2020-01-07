library(Rsubread)
library(limma)
library(edgeR)
library(ggplot2)
library(RColorBrewer) 
library(gplots)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(genefilter)
library(pheatmap)
library(ggfortify)
library(tidyverse)


# setwd("/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/limma+voom/input/")

## Index
# gappedIndex para poder funcionar con 16GB de RAM en mac
# buildindex(basename = "GRCh38", reference = "hg38.fa.gz", gappedIndex = T)
# buildindex(basename = "hg19", reference = "hg19.fa")

setwd("/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/fastp")

## Mapping
reads1 = list.files(path = "/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/fastp",
                   pattern = "*_R1_001.fastq.gz", recursive = T)
reads2 = list.files(path = "/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/fastp",
                    pattern = "*_R2_001.fastq.gz*", recursive = T)

align(index="GRCh38", readfile1 = reads1, readfile2 = reads2,
      type = "rna", input_format="gzFASTQ",
      output_format="BAM", nthreads = 2)


## Summarize featurecounts
setwd("/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/limma+voom/output/")
alignments = list.files(path = "/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/limma+voom/output/",
                        pattern = ".BAM$", recursive = T)
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
write.csv(dge$counts, file = "dge_limma_gsea.csv")

## Differential expresssion
setwd("/Volumes/blue_kinex/TFM/GEP_CellLines/RNAseq_marcos/limma+voom/input")
targets <- readTargets("targets2.txt")
rownames(targets) <- targets$X
targets <- dplyr::select(targets, clon, tratamiento)
group <- factor(paste0(targets$clon, targets$tratamiento))
design <- model.matrix(~0 + group)
rownames(design) <- rownames(targets)
voom_dge <- voom(dge,design,plot=TRUE)
boxplot(voom_dge$E, cex.axis = 0.7)


cont_matrix1 <- makeContrasts(sen0hvssen6h = groupsen0h-groupsen6h, levels = design)
cont_matrix2 <- makeContrasts(res0hvsres6h = groupres0h-groupres6h, levels = design)
cont_matrix3 <- makeContrasts(sen0hvsres0h = groupsen0h-groupres0h, levels = design)
cont_matrix4 <- makeContrasts(sen6hvsres6h = groupsen6h-groupres6h, levels = design)

fit <- lmFit(voom_dge, design)
fit1 <- contrasts.fit(fit, cont_matrix1)
fit1 <- eBayes(fit1)
fit2 <- contrasts.fit(fit, cont_matrix2)
fit2 <- eBayes(fit2)
fit3 <- contrasts.fit(fit, cont_matrix3)
fit3 <- eBayes(fit3)
fit4 <- contrasts.fit(fit, cont_matrix4)
fit4 <- eBayes(fit4)

toptable1 <- topTable(fit1, number = nrow(fit$genes), genelist = fit$genes$GeneSymbol, adjust.method = "BH", sort.by = "p")
sum(toptable1$adj.P.Val<0.05)
# Venn's diagram
topadj1 <- toptable1 %>% filter(adj.P.Val<0.05)
topadj1 <- topadj1 %>% dplyr::select(ID, logFC, adj.P.Val)
write.csv(topadj1, file = "limma_sen0sen6_venn.csv")

toptable2 <- topTable(fit2, number = nrow(fit$genes), genelist = fit$genes$GeneSymbol, adjust.method = "BH", sort.by = "p")
sum(toptable2$adj.P.Val<0.05)

toptable3 <- topTable(fit3, number = nrow(fit$genes), genelist = fit$genes$GeneSymbol, adjust.method = "BH", sort.by = "p")
sum(toptable3$adj.P.Val<0.05)
# Venn's diagram
topadj3 <- toptable3 %>% filter(adj.P.Val<0.05)
topadj3 <- topadj3 %>% dplyr::select(ID, logFC, adj.P.Val)
write.csv(topadj3, file = "limma_sen0res0_venn.csv")

toptable4 <- topTable(fit4, number = nrow(fit$genes), genelist = fit$genes$GeneSymbol, adjust.method = "BH", sort.by = "p")
sum(toptable4$adj.P.Val<0.05)

## plotMA
dev.off()
dev.new()
par(mfrow=c(2,2))
summa_fit1 <- decideTests(fit1, p.value = 0.05) 
summa_fit2 <- decideTests(fit2, p.value = 0.05) 
summa_fit3 <- decideTests(fit3, p.value = 0.05) 
summa_fit4 <- decideTests(fit4, p.value = 0.05) 
plotMD(fit1, status = summa_fit1, main = "sen0h_vs_sen6h", ylim = c(-2, 2), legend = F)
plotMD(fit2, status = summa_fit2, main = "res0h_vs_res6h", ylim = c(-2, 2), legend = F)
plotMD(fit3, status = summa_fit3, main = "sen0h_vs_res0h", ylim = c(-2, 2), legend = F)
plotMD(fit4, status = summa_fit4, main = "sen6h_vs_res6h", ylim = c(-2, 2), legend = F)


## MDS (multidimensional scaling plot)
dev.off()
dev.new()
col <- as.numeric(group)
plotMDS(voom_dge$E, col = col, cex = 0.8, xlab = "PC1", ylab = "PC2", main = "limma + voom")
legend("right", legend = c("sensible-0h", "sensible-6h", "resistente-0h", "resistente-6h"),
       fill = c("green","blue","black","red"), 
       border = c("green","blue","black","red"),
       cex = 0.7)

## Sample to sample distances heatmap
dev.off()
dev.new()
sampleDistsvoom <- dist(t(voom_dge$E))
sampleDistMatrixvoom <- as.matrix(sampleDistsvoom)
rownames(sampleDistMatrixvoom) <- paste(targets$clon, targets$tratamiento, sep="-")
colnames(sampleDistMatrixvoom) <- rownames(targets)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrixvoom, main = "limma + voom",
         clustering_distance_rows=sampleDistsvoom,
         clustering_distance_cols=sampleDistsvoom,
         col=colors)

## Heatmap
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
rv <- order(apply(voom_dge, MARGIN=1, FUN=var, na.rm=TRUE),decreasing = T)[1:20]
df <- as.data.frame(targets)
pheatmap(voom_dge$E[rv,], cluster_rows=T, show_rownames=F, main = "limma + voom",
         cluster_cols=T, annotation_col=df, cutree_cols = 2, annotation_names_col = F,
         annotation_colors = mycolors)

## Heatmap top 100 genes sen0-sen6
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
toptable1_hm <- toptable1[1:100,]
pheatmap(voom_dge$E[rownames(toptable1_hm),c("18NR0293","18NR0297","18NR0301","18NR0305","18NR0309","18NR0313")],
         scale = "row", labels_row = toptable1_hm$ID, annotation_col = as.data.frame(targets),
         treeheight_col = 20,treeheight_row = 20, annotation_names_col = F, fontsize = 8,
         fontsize_row = 2.5, cellwidth = 15,cutree_cols = 2,
         annotation_colors = mycolors,
         main = "limma + voom\nSEN0 vs SEN6\ntop 100 genes\n")

# dev.off()
# dev.new()
# toptable1_hm <- toptable1[1:100,]
# mypalette <- brewer.pal(11,"RdYlBu")
# morecols <- colorRampPalette(mypalette)
# design_df <- as.data.frame(design)[c("18NR0293","18NR0297","18NR0301","18NR0305","18NR0309","18NR0313"),c("groupsen0h","groupsen6h")]
# col1 <- factor(as.factor(design_df$groupsen0h), labels = c(NA, "#FF99FF"))
# col2 <- factor(as.factor(design_df$groupsen6h), labels = c(NA, "#990099"))
# col3 <- coalesce(as.vector(col1), as.vector(col2))
# lmat = rbind(c(5,4), c(0,1), c(3,2))
# lhei <- c(1,0.2,5)
# lwid <- c(0.5, 1)
# heatmap.2(voom_dge$E[rownames(toptable1_hm),c("18NR0293","18NR0297","18NR0301","18NR0305","18NR0309","18NR0313")], 
#           col=rev(morecols(100)), ColSideColors = col3, trace='none',scale='row', 
#           main = "limma + voom\ntop 100 genes",
#           dendrogram = "column", labRow = toptable1_hm$ID, density.info = "histogram",key.title = NA,keysize = 0.5,
#           cexRow = 0.6, cexCol = 1, revC = T, na.rm = T,
#           margins = c(6,6), srtCol = 0, adjCol = c(0.5,0.5),
#           lmat = lmat, lhei = lhei, lwid = lwid,
#           reorderfun=function(d, w) reorder(d, w, agglo.FUN = var))
# legend("left", legend = c("sensible-0h", "sensible-6h"), fill = c("#FF99FF","#990099"), 
#        border = c("#FF99FF","#990099"), box.col = "#CCCCCC",
#        cex = 0.8, inset = .001, xjust = 0.5)
#                       
# ## Heatmap top 100 genes sen0-res0
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycolors = list(tratamiento=c("0h"="#FF99FF","6h"="#9900CC"),clon=c(sen="#99FF66",res="#006633"))
toptable3_hm <- toptable3[1:100,]
pheatmap(voom_dge$E[rownames(toptable3_hm),c("18NR0293","18NR0301","18NR0309","18NR0295","18NR0303","18NR0311")],
         scale = "row", labels_row = toptable3_hm$ID, annotation_col = as.data.frame(targets),
         treeheight_col = 20,treeheight_row = 20, annotation_names_col = F, fontsize = 8,
         fontsize_row = 2.5, cellwidth = 15,cutree_cols = 2,
         annotation_colors = mycolors,
         main = "limma + voom\nSEN0 vs RES0\ntop 100 genes\n")

# dev.off()
# dev.new()
# toptable3_hm <- toptable3[1:100,]
# mypalette <- brewer.pal(11,"RdYlBu")
# morecols <- colorRampPalette(mypalette)
# design_df2 <- as.data.frame(design)[c("18NR0293","18NR0301","18NR0309","18NR0295","18NR0303","18NR0311"),c("groupsen0h","groupres0h")]
# lmat = rbind(c(5,4), c(0,1), c(3,2))
# lhei <- c(1,0.2,5)
# lwid <- c(0.5, 1)
# col4 <- factor(as.factor(design_df2$groupsen0h), labels = c(NA, "#99FF33"))
# col5 <- factor(as.factor(design_df2$groupres0h), labels = c(NA, "#336600"))
# col6 <- coalesce(as.vector(col4), as.vector(col5))
# heatmap.2(voom_dge$E[rownames(toptable3_hm),c("18NR0293","18NR0301","18NR0309","18NR0295","18NR0303","18NR0311")], 
#           col=rev(morecols(100)), ColSideColors = col6, trace='none',scale='row', 
#           main = "limma + voom\ntop 100 genes",
#           dendrogram = "column", labRow = toptable3_hm$ID, density.info = "histogram",key.title = NA,keysize = 0.5,
#           cexRow = 0.6, cexCol = 1, revC = T, na.rm = T,
#           margins = c(6,6), srtCol = 0, adjCol = c(0.5,0.5),
#           lmat = lmat, lhei = lhei, lwid = lwid,
#           reorderfun=function(d, w) reorder(d, w, agglo.FUN = var))
# legend("left", legend = c("sensible-0h", "resistente-0h"), fill = c("#99FF33","#336600"),
#        border = c("#99FF33","#336600"), box.col = "#CCCCCC",
#        cex = 0.8, inset = .001,xjust = 0.5)
