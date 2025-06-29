#paquetes necesarios
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma", "Biobase", "pheatmap", "EnhancedVolcano", "dplyr", "VennDiagram", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "enrichR", "edgeR", "ggfortify", "uwot2", "caret", "glmnet", "e1071", "pROC"))

#librerías
library(GEOquery)
library(limma)
library(Biobase)
library(pheatmap)
library(EnhancedVolcano)
library(dplyr)
library(tidyr)
library(matrixStats)
library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(enrichR)
library(edgeR)
library(ggplot2)
library(ggfortify)
library(uwot)
library(caret)
library(glmnet)
library(e1071)
library(pROC)

#directorio de trabajo 
setwd("~/Documents/BIOINF_UNIR/TFM/")

#cargamos datos descargados
counts <- read.csv("GSE138852_counts.csv", row.names = 1)
covariates <- read.csv("GSE138852_covariates.csv")

#ajustamos nombres de columnas para que coincidan
colnames(counts) <- gsub("\\.", "-", colnames(counts))
covariates$X <- gsub("\\.", "-", covariates$X)

#verificamos tipos de pacientes
unique(covariates$oupSample.subclustCond)

#ignoramos los undetermined
covariates <- covariates %>%
  filter(oupSample.subclustCond %in% c("AD", "ct"))

#establecemos IDs para cada paciente
covariates <- covariates %>%
  mutate(
    Condition = oupSample.subclustCond,
    PatientID = paste0(Condition, "_", oupSample.subclustID)
  )

#verificamos cada paciente
unique(covariates$PatientID)

#vemos qué tipos celulares tenemos y cuántas células de cada tipo
table(covariates$oupSample.cellType)

#vemos el número de células en los subgrupos de pacientes y el número en cada
table(covariates$oupSample.batchCond)

#vemos tipos celulares para cada tipo de paciente
table(covariates$oupSample.cellType_batchCond)

#estudiamos la distribución del dataset 
table(covariates$oupSample.subclustCond,
      covariates$oupSample.subclustID)

#filtramos el tipo de células y eliminamos 'OPC', 'endo',doublet' y 'unID'
cell_types_valid <- c("astro", "oligo", "mg", "neuron")
covariates <- covariates[covariates$oupSample.cellType %in% cell_types_valid, ]

table(covariates$oupSample.subclustCond)
unique(covariates$oupSample.subclustCond)
unique(covariates$PatientID)
