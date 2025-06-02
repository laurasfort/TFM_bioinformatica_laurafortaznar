#Instalación y carga de paquetes
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "Biobase", "pheatmap", "EnhancedVolcano"))

library(GEOquery)
library(limma)
library(Biobase)
library(pheatmap)
library(EnhancedVolcano)
library(tidyverse)
library(dplyr)

#Directorio de trabajo
setwd("~/Documents/BIOINF_UNIR/TFM/")

#Carga de datos
counts <- read.csv("GSE138852_counts.csv", row.names = 1)
covariates <- read.csv("GSE138852_covariates.csv")

#Ajuste de nombres de columnas
colnames(counts) <- gsub("\\.", "-", colnames(counts))
covariates$X <- gsub("\\.", "-", covariates$X)

#Filtrado de células válidas y condición
covariates <- covariates %>% 
  filter(oupSample.subclustCond %in% c("AD", "ct"))

cell_types_valid <- c("astro", "oligo", "mg", "endo", "neuron", "OPC")
covariates <- covariates[covariates$oupSample.cellType %in% cell_types_valid, ]
