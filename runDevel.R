# library(WEIN)

.libPaths()

# .libPaths("c:/Rlibs/")

# BiocManager::install("DEFormats",ask = F, update = T)
# install.packages(c("ashr","readr"), dependencies = T )
# BiocManager::install(c('DESeq2', 'pcaExplorer', 'IHW', 'goseq', 'GOstats', 'GO.db', 'rentrez', 'rintrojs'))
# BiocManager::install(c('airway', 'apeglm', 'BiocStyle', 'org.Hs.eg.db','TxDb.Hsapiens.UCSC.hg38.knownGene'))
library(ashr)
 library(airway)
 library(shinyjqui)
 library(apeglm)
library(DESeq2) 
library(SummarizedExperiment) 
library(GenomicRanges) 
library(reshape2)
library(IRanges)
library(S4Vectors) 
library(ggplot2) 
library(heatmaply) 
library(plotly)
library(pcaExplorer) 
library(IHW) 
library(gplots) 
library(UpSetR) 
library(goseq) 
library(stringr) 
library(dplyr)
library(limma) 
library(GOstats) 
library(GO.db) 
library(AnnotationDbi) 
library(shiny)
library(shinydashboard) 
library(shinyBS) 
library(DT) 
library(rentrez) 
library(rintrojs) 
library(ggrepel) 
library(knitr)
library(rmarkdown) 
library(shinyAce) 
library(BiocParallel) 
library(grDevices) 
library(base64enc)
library(methods)
library(testthat) 
 library(BiocStyle) 
 library(airway) 
 library(org.Hs.eg.db)
 library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
library(DEFormats) 
library(edgeR)
library(RColorBrewer)
library(topGO)


cp = load(file = "idealState_20231118_123859.RData")
# WEIN(dds_obj = dds_airway, annotation_obj = annotation_airway, res_obj = res_airway)
WEIN(state_file = "idealState_20231118_123859.RData")
  # WEIN()
# cp = load(file = "test.Rdata")
# dds_obj = ideal_env$ideal_values_20210526_152858$dds_obj 
# res_obj = ideal_env$ideal_values_20210526_152858$res_obj
# annotation_obj = ideal_env$ideal_values_20210526_152858$annotation_obj
# countmatrix = ideal_env$ideal_values_20210526_152858$countmatrix
# expdesign = ideal_env$ideal_values_20210526_152858$expdesign
# gene_signatures = ideal_env$ideal_values_20210526_152858$gene_signatures
# 
# 
# options(browser = "/Applications/Google\ Chrome.app/Contents/MacOS/Google\ Chrome")
# 
# source("R/ideal.R", local = T);
# WEIN(dds_obj = dds_obj,
#       res_obj = res_obj,
#       annotation_obj = annotation_obj,
#       countmatrix = countmatrix,
#       expdesign = expdesign, 
#       gene_signatures = gene_signatures
#     )
# 
