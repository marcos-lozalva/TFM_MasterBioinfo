## Install limma
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("limma")
library(limma)
library(genefilter)

## Install annotation package
# BiocManager::install("hgug4112a.db") # Human whole genome
# library("hgug4112a.db")
# install.packages("~/Documents/GEP_CellLines/arrays_marcos/", repos=NULL, type="source")
library(org.Hs.eg.db)
library(hga039494.db) # annotation db from desing ID 03949

# Set working directory where raw data is located
setwd("/Volumes/blue_kinex/TFM/GEP_CellLines/arrays_marcos/raw_data/")

# Reading data with file targets.txt 
targets <- readTargets("/Volumes/blue_kinex/TFM/GEP_CellLines/arrays_marcos/input/targets.txt")
data <- read.maimages(targets, source="agilent", green.only=TRUE)

## Checking control probes
spottypes <- readSpotTypes()
data$genes$Status <- controlStatus(spottypes, data)
plotMD(data)

# Background correct
bgdata <- backgroundCorrect(data, method = "normexp", offset = 1)

## Normalization ("For single-channel arrays, within array normalization is not usually relevant and so 
# normalizeBetweenArrays is the sole normalization step")
normdata <- normalizeBetweenArrays(bgdata,method = "quantile")

# Create expset for running it in GSEA 
expset <- cbind(normdata$genes$ProbeName, normdata$E)
write.csv(expset, "/Volumes/blue_kinex/TFM/GEP_CellLines/arrays_marcos/input/expset.csv", row.names = F)

# Gene annotation
normdata$genes$EntrezID <- mapIds(hga039494.db, normdata$genes$ProbeName, keytype="PROBEID", column="ENTREZID")
normdata$genes$Symbol <- mapIds(hga039494.db, normdata$genes$ProbeName, keytype="PROBEID", column="SYMBOL")
normdata$genes$ENS <- mapIds(hga039494.db, normdata$genes$ProbeName, keytype="PROBEID", column="ENSEMBL")

expset_symbol <- cbind(normdata$genes$Symbol, normdata$E)
write.csv(expset_symbol, "/Volumes/blue_kinex/TFM/GEP_CellLines/arrays_marcos/input/expset_symbol.csv", row.names = F)

# QC assesment
par(mfrow=c(1,2))
par(cex.axis = 0.3)
boxplot(data$E, main="Boxplot Before Normalization", col = "lightgrey")
boxplot(normdata$E, main="Boxplot After Normalization", col = "lightgrey", )

# Design matrix
design1 <- cbind(res0h=c(1,1,1,0,0,0,0,0,0,0,0,0), res6h=c(0,0,0,1,1,1,0,0,0,0,0,0))
rownames(design1) <- targets$FileName
design2 <- cbind(sen0h=c(0,0,0,0,0,0,1,1,1,0,0,0), sen6h=c(0,0,0,0,0,0,0,0,0,1,1,1))
rownames(design2) <- targets$FileName
design3 <- cbind(res0h=c(1,1,1,0,0,0,0,0,0,0,0,0), sen0h=c(0,0,0,0,0,0,1,1,1,0,0,0))
rownames(design3) <- targets$FileName
design4 <- cbind(res6h=c(0,0,0,1,1,1,0,0,0,0,0,0), sen6h=c(0,0,0,0,0,0,0,0,0,1,1,1))
rownames(design4) <- targets$FileName


# Contrast matrix
cont.matrix1 <- makeContrasts(res0h_vs_res6h = res0h-res6h, levels = design1)
cont.matrix2 <- makeContrasts(sen0h_vs_sen6h = sen0h-sen6h, levels = design2)
cont.matrix3 <- makeContrasts(res0h_vs_sen0h = res0h-sen0h, levels = design3)
cont.matrix4 <- makeContrasts(res6h_vs_sen6h = res6h-sen6h, levels = design4)

# Lineal model and ebayes
fit_design1 <- lmFit(normdata, design1)
fit_design2 <- lmFit(normdata, design2)
fit_design3 <- lmFit(normdata, design3)
fit_design4 <- lmFit(normdata, design4)
fit2_design1 <- contrasts.fit(fit_design1, cont.matrix1)
fit2_design1 <- eBayes(fit2_design1)
fit2_design2 <- contrasts.fit(fit_design2, cont.matrix2)
fit2_design2 <- eBayes(fit2_design2)
fit2_design3 <- contrasts.fit(fit_design3, cont.matrix3)
fit2_design3 <- eBayes(fit2_design3)
fit2_design4 <- contrasts.fit(fit_design4, cont.matrix4)
fit2_design4 <- eBayes(fit2_design4)

# Toptable
toptable1 <- topTable(fit2_design1, number = dim(normdata), adjust.method = "BH", sort.by = "p")
toptable1_NA <- na.omit(toptable1)
results <- decideTests(fit2_design1)
toptable2 <- topTable(fit2_design2, number = dim(normdata), adjust.method = "BH", sort.by = "p")
toptable2_NA <- na.omit(toptable2)
toptable3 <- topTable(fit2_design3, number = dim(normdata), adjust.method = "BH", sort.by = "p")
toptable3_NA <- na.omit(toptable3)
toptable4 <- topTable(fit2_design4, number = dim(normdata), adjust.method = "BH", sort.by = "p")
toptable4_NA <- na.omit(toptable4)

# Save results
save(toptable1,file="RES_basalvsRES_tto.RData")
save(toptable2,file="SEN_basalvsSEN_tto")
save(toptable3,file="RES_basalvsSEN_basal")
save(toptable4,file="RES_ttovsSEN_tto.RData")



