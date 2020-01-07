library("AnnotationDbi")
library("AnnotationForge")
available.db0pkgs()
library(BiocManager)
install("human.db0")
install("org.Hs.eg.db")
library(org.Hs.eg.db)
available.dbschemas()

makeDBPackage("HUMANCHIP_DB",
              affy=FALSE,
              prefix="hga039494",
              fileName="~/Documents/GEP_CellLines/arrays_marcos/raw_data/genelist_array_accesion_agilent039494.txt",
              baseMapType="gb",
              outputDir = ".",
              version="1.0.0",
              manufacturer = "Agilent",
              chipName = "hga039494")

install.packages("~/Documents/GEP_CellLines/hga039494.db/", repos=NULL, type="source")
library(hga039494.db)
