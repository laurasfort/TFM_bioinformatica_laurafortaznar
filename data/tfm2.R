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

### ---- OLIGOS ----### 

#filtramos solo un tipo celular, aquí oligodendrocitos
cov_oligos <- covariates %>% 
  filter(oupSample.cellType_batchCond %in% c("oligo_AD", "oligo_Control"))

#añadimos el identificador de grupo por paciente
cov_oligos$group <- cov_oligos$PatientID

#subset de la matriz de expresión
counts_oligos <- counts[, cov_oligos$X]

#alineamos las columnas
counts_oligos <- counts_oligos[, match(cov_oligos$X, colnames(counts_oligos))]

#agrupar por paciente (suma de UMI por gen)
group_vector_oligo <- cov_oligos$group
counts_oligo_pseudobulk <- rowsum(t(counts_oligos), group = group_vector_oligo)
counts_oligo_pseudobulk <- t(counts_oligo_pseudobulk)  # Volver a genes x muestras

#matriz de diseño
group_labels_oligo <- ifelse(grepl("AD", colnames(counts_oligo_pseudobulk)), "EA", "Control")
group_oligo <- factor(group_labels_oligo)
design_oligo <- model.matrix(~ 0 + group_oligo)
colnames(design_oligo) <- levels(group_oligo)

#normalizar
log_counts_oligo <- log2(counts_oligo_pseudobulk + 1)

#análisis con limma
fit_oligo <- lmFit(log_counts_oligo, design_oligo)
contrast_oligo <- makeContrasts(EA - Control, levels = design_oligo)
fit2_oligo <- contrasts.fit(fit_oligo, contrast_oligo)
fit2_oligo <- eBayes(fit2_oligo)

# QQ plot para evaluación de ajuste (limma)
png("QQplot_oligo_EA_vs_Control.png", width = 800, height = 800)
qqt(fit2_oligo$t, df = fit2_oligo$df.total, main = "QQ plot - Oligodendrocitos")
abline(0, 1, col = "red")
dev.off()

#resultados
results_oligo <- topTable(fit2_oligo, adjust = "fdr", number = Inf)
head(results_oligo)
write.csv(results_oligo, "DEG_oligos_EA_vs_Control.csv")

#resumen de DEGs sobreexpresados
sum(results_oligo$adj.P.Val < 0.05 & abs(results_oligo$logFC) > 1)
up_oligo <- results_oligo %>% filter(adj.P.Val < 0.05 & logFC > 1)
down_oligo <- results_oligo %>% filter(adj.P.Val < 0.05 & logFC < -1)

nrow(up_oligo) #sobreexpresados
nrow(down_oligo) #infraexpresados

pdf("Volcano_EA_vs_ct_Oligos.pdf", width = 10, height = 8)
EnhancedVolcano(results_oligo,
                lab = rownames(results_oligo),
                x = "logFC",
                y = "P.Value",
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Oligodendrocitos: EA vs Control")
dev.off()

results_oligo <- as.data.frame(results_oligo)

# Top 20 genes con menor FDR
top20_genes_oligo <- results_oligo[order(results_oligo$adj.P.Val), ][1:20, ]

colnames(log_counts_oligo) <- gsub("^AD", "EA", colnames(log_counts_oligo))

#anotación por paciente y condición
annotation_col_oligo <- data.frame(
  Condicion = group_labels_oligo,
  Paciente = colnames(log_counts_oligo)
)
rownames(annotation_col_oligo) <- colnames(log_counts_oligo)

# Subset expresión para los top 20 genes y escalado

genes_heatmap_oligo <- rownames(top20_genes_oligo)
heat_data_oligo <- log_counts_oligo[genes_heatmap_oligo, ]
heat_data_scaled_oligo <- t(scale(t(heat_data_oligo)))

orden_oligos <- order(annotation_col_oligo$Condicion, decreasing = FALSE)
heat_data_scaled_oligo <- heat_data_scaled_oligo[, orden_oligos]
annotation_col_oligo <- annotation_col_oligo[orden_oligos, ]

pdf("Heatmap_EA_vs_ct_Oligos.pdf", width = 10, height = 8)
pheatmap(heat_data_scaled_oligo,
         annotation_col = annotation_col_oligo,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 5,
         main = "Heatmap - Top 20 genes Oligos")
dev.off()

# Filtrar genes significativos con logFC alto
biomarcadores_oligo <- results_oligo %>%
  filter(adj.P.Val < 0.01 & abs(logFC) > 2) %>%
  arrange(desc(abs(logFC))) 

biomarcadores_oligo <- biomarcadores_oligo[1:20, ]

rownames(biomarcadores_oligo)

head(biomarcadores_oligo, 20)

---------------------------
------ # MICROGLÍA # ------
  
#filtramos por microglia
cov_micro <- covariates %>%
  filter(oupSample.cellType_batchCond %in% c("mg_AD", "mg_Control"))

cov_micro$group <- cov_micro$PatientID

counts_micro <- counts[, cov_micro$X]

counts_micro <- counts_micro[, match(cov_micro$X, colnames(counts_micro))]

group_vector_micro <- cov_micro$group
counts_micro_pseudobulk <- rowsum(t(counts_micro), group = group_vector_micro)
counts_micro_pseudobulk <- t(counts_micro_pseudobulk)

group_labels_micro <- ifelse(grepl("AD", colnames(counts_micro_pseudobulk)), "EA", "Control")
group_micro <- factor(group_labels_micro)
design_micro <- model.matrix(~ 0 + group_micro)
colnames(design_micro) <- levels(group_micro)

log_counts_micro <- log2(counts_micro_pseudobulk + 1)

fit_micro <- lmFit(log_counts_micro, design_micro)
contrast_micro <- makeContrasts(EA - Control, levels = design_micro)
fit2_micro <- contrasts.fit(fit_micro, contrast_micro)
fit2_micro <- eBayes(fit2_micro)

png("QQplot_micro_EA_vs_Control.png", width = 800, height = 800)
qqt(fit2_micro$t, df = fit2_micro$df.total, main = "QQ plot - Microglia")
abline(0, 1, col = "red")
dev.off()

results_micro <- topTable(fit2_micro, adjust = "fdr", number = Inf)
head(results_micro)
write.csv(results_micro, "DEG_micro_EA_vs_Control.csv")

#resumen de DEGs sobreexpresados
sum(results_micro$adj.P.Val < 0.05 & abs(results_micro$logFC) > 1)
up_micro <- results_micro %>% filter(adj.P.Val < 0.05 & logFC > 1)
down_micro <- results_micro %>% filter(adj.P.Val < 0.05 & logFC < -1)

nrow(up_micro) #sobreexpresados
nrow(down_micro) #infraexpresados

pdf("Volcano_micro_EA_vs_Control.pdf", width = 10, height = 8)
EnhancedVolcano(results_micro,
                lab = rownames(results_micro),
                x = "logFC",
                y = "P.Value",
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Microglía: EA vs Control")
dev.off()

# Top 20 genes con menor FDR

variable_genes_micro <- apply(log_counts_micro, 1, var)
top20_genes_micro <- names(sort(variable_genes_micro, decreasing = TRUE))[1:20]

colnames(log_counts_micro) <- gsub("^AD", "EA", colnames(log_counts_micro))

# Subset y escalado
heat_data_micro <- log_counts_micro[top20_genes_micro, ]
heat_data_scaled_micro <- t(scale(t(heat_data_micro)))

colnames(log_counts_micro) <- gsub("^AD", "EA", colnames(log_counts_micro))

annotation_col_micro <- data.frame(
  Condicion = group_labels_micro,
  Paciente = colnames(log_counts_micro)
)
rownames(annotation_col_micro) <- colnames(log_counts_micro)

orden_micro <- order(annotation_col_micro$Condicion, decreasing = FALSE)
heat_data_scaled_micro <- heat_data_scaled_micro[, orden_micro]
annotation_col_micro <- annotation_col_micro[orden_micro, ]

pdf("Heatmap_EA_vs_ct_Microglia.pdf", width = 10, height = 8)
pheatmap(heat_data_scaled_micro,
         annotation_col = annotation_col_micro,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 6,
         main = "Heatmap – Top 20 genes Microglía")
dev.off()

biomarcadores_micro <- results_micro %>%
  filter(adj.P.Val < 0.01 & abs(logFC) > 2) %>%
  arrange(desc(abs(logFC))) 

biomarcadores_micro <- biomarcadores_micro[1:20, ]

rownames(biomarcadores_micro)

head(biomarcadores_micro, 20)

---------------------------
------ # ASTROCITOS # -----

#filtramos por astrocitos
cov_astro <- covariates %>%
  filter(oupSample.cellType_batchCond %in% c("astro_AD", "astro_Control"))

cov_astro$group <- cov_astro$PatientID

counts_astro <- counts[, cov_astro$X]
counts_astro <- counts_astro[, match(cov_astro$X, colnames(counts_astro))]

group_vector_astro <- cov_astro$group
counts_astro_pseudobulk <- rowsum(t(counts_astro), group = group_vector_astro)
counts_astro_pseudobulk <- t(counts_astro_pseudobulk)

group_labels_astro <- ifelse(grepl("AD", colnames(counts_astro_pseudobulk)), "EA", "Control")
group_astro <- factor(group_labels_astro)
design_astro <- model.matrix(~ 0 + group_astro)
colnames(design_astro) <- levels(group_astro)

log_counts_astro <- log2(counts_astro_pseudobulk + 1)

fit_astro <- lmFit(log_counts_astro, design_astro)
contrast_astro <- makeContrasts(EA - Control, levels = design_astro)
fit2_astro <- contrasts.fit(fit_astro, contrast_astro)
fit2_astro <- eBayes(fit2_astro)

png("QQplot_astro_EA_vs_Control.png", width = 800, height = 800)
qqt(fit2_astro$t, df = fit2_astro$df.total, main = "QQ plot - Astrocitos")
abline(0, 1, col = "red")
dev.off()

results_astro <- topTable(fit2_astro, adjust = "fdr", number = Inf)
write.csv(results_astro, "DEG_astro_EA_vs_Control.csv")

#resumen de DEGs sobreexpresados
sum(results_astro$adj.P.Val < 0.05 & abs(results_astro$logFC) > 1)
up_astro <- results_astro %>% filter(adj.P.Val < 0.05 & logFC > 1)
down_astro <- results_astro %>% filter(adj.P.Val < 0.05 & logFC < -1)

nrow(up_astro) #sobreexpresados
nrow(down_astro) #infraexpresados

top20_genes_astro <- results_astro[order(results_astro$adj.P.Val), ][1:20, ]
top20_genes_astro

pdf("Volcano_astro_EA_vs_Control.pdf", width = 10, height = 8)
EnhancedVolcano(results_astro,
                lab = rownames(results_astro),
                x = "logFC",
                y = "P.Value",
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Astrocitos: EA vs Control")
dev.off()

# Subset expresión para los top 20 genes
variable_genes_astro <- apply(log_counts_astro, 1, var)
top_20_genes_astro <- names(sort(variable_genes_astro, decreasing = TRUE))[1:20]

colnames(log_counts_astro) <- gsub("^AD", "EA", colnames(log_counts_astro))

#anotación por paciente y condición
annotation_col_astro <- data.frame(
  Condicion = group_labels_astro,
  Paciente = colnames(log_counts_astro)
)
rownames(annotation_col_astro) <- colnames(log_counts_astro)

# Subset expresión para los top 20 genes y escalado
genes_heatmap_astro <- top_20_genes_astro
heat_data_astro <- log_counts_astro[genes_heatmap_astro, ]
heat_data_scaled_astro <- t(scale(t(heat_data_astro)))

colnames(log_counts_astro) <- gsub("^AD", "EA", colnames(log_counts_astro))

orden_astro <- order(annotation_col_astro$Condicion, decreasing = FALSE)
heat_data_scaled_astro <- heat_data_scaled_astro[, orden_astro]
annotation_col_astro <- annotation_col_astro[orden_astro, ]

pdf("Heatmap_EA_vs_ct_Astro.pdf", width = 10, height = 8)
pheatmap(heat_data_scaled_astro,
         annotation_col = annotation_col_astro,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 5,
         main = "Heatmap - Top 20 genes Astro")
dev.off()

biomarcadores_astro <- results_astro %>%
  filter(adj.P.Val < 0.01 & abs(logFC) > 2) %>%
  arrange(desc(abs(logFC))) 

biomarcadores_astro <- biomarcadores_astro[1:20, ]

rownames(biomarcadores_astro)

head(biomarcadores_astro, 20)

----------------------------
  #TOP GENES + FUNCIÓN

  get_top10_with_function <- function(result_df, tipo) {
    top_genes <- result_df %>%
      filter(adj.P.Val < 0.05) %>%
      arrange(desc(abs(logFC))) %>%
      head(10) %>%
      mutate(Gene = rownames(.))
    
    gene_symbols <- top_genes$Gene
    
    # Anotar funciones usando org.Hs.eg.db
    annotated <- bitr(gene_symbols,
                      fromType = "SYMBOL",
                      toType = c("ENTREZID", "GENENAME"),
                      OrgDb = org.Hs.eg.db)
    
    top_genes_annotated <- left_join(top_genes, annotated, by = c("Gene" = "SYMBOL"))
    top_genes_annotated$CellType <- tipo
    
    # Redondear logFC y adj.P.Val a 3 cifras decimales
    top_genes_annotated <- top_genes_annotated %>%
      mutate(
        logFC = round(logFC, 3),
        adj.P.Val = signif(adj.P.Val, 3)  # mantiene formato científico si es necesario
      )
    
    return(top_genes_annotated[, c("Gene", "logFC", "adj.P.Val", "GENENAME", "CellType")])
  }

# Aplicar por tipo celular
top10_oligo <- get_top10_with_function(results_oligo, "Oligodendrocitos")
top10_micro <- get_top10_with_function(results_micro, "Microglía")
top10_astro <- get_top10_with_function(results_astro, "Astrocitos")

# Combinar en una sola tabla
top10_combined <- bind_rows(top10_oligo, top10_micro, top10_astro)

# Guardar
write.csv(top10_combined, "Top10_DEGs_funciones_por_celula.csv", row.names = FALSE)


----------------------------
----- # VENN DIAGRAM # -----

# Cargar resultados (usa los CSVs generados por topTable)
deg_oligo <- read.csv("DEG_oligos_EA_vs_Control.csv", row.names = 1)
deg_micro <- read.csv("DEG_micro_EA_vs_Control.csv", row.names = 1)
deg_astro <- read.csv("DEG_astro_EA_vs_Control.csv", row.names = 1)

# Filtrar por FDR y logFC
genes_oligo <- rownames(deg_oligo[deg_oligo$adj.P.Val < 0.01 & abs(deg_oligo$logFC) > 1, ])
genes_micro <- rownames(deg_micro[deg_micro$adj.P.Val < 0.01 & abs(deg_micro$logFC) > 1, ])
genes_astro <- rownames(deg_astro[deg_astro$adj.P.Val < 0.01 & abs(deg_astro$logFC) > 1, ])

# Crear lista de conjuntos
gene_lists <- list(
  Oligodendrocitos = genes_oligo,
  Microglia = genes_micro,
  Astrocitos = genes_astro
)

# Generar diagrama de Venn
png("venndiagram.png", width = 800, height = 800)
venn.plot <- venn.diagram(
  x = gene_lists,
  filename = NULL,
  fill = c("cornflowerblue", "darkolivegreen3", "goldenrod1"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "Genes diferencialmente expresados (FDR<0.01, logFC>1)"
)

# Mostrar el diagrama
grid::grid.newpage()
grid::grid.draw(venn.plot)
dev.off()

common_oligo_micro <- intersect(genes_oligo, genes_micro)
common_oligo_astro <- intersect(genes_oligo, genes_astro)
common_micro_astro <- intersect(genes_micro, genes_astro)

common_all_three <- Reduce(intersect, list(genes_oligo, genes_micro, genes_astro))

length(common_oligo_micro)     # Oligo ∩ Microglía
length(common_oligo_astro)     # Oligo ∩ Astrocitos
length(common_micro_astro)     # Microglía ∩ Astrocitos
length(common_all_three)       # Comunes a los tres

common_oligo_micro
common_oligo_astro
common_micro_astro
common_all_three

write.csv(common_all_three, "genes_comunes_tres_tipos.csv", row.names = FALSE)

#reducimos número de genes
intersect_all <- Reduce(intersect, list(genes_oligo, genes_micro, genes_astro))
length(intersect_all)

common_genes <- intersect_all

tabla_comun <- data.frame(
  Gene = common_genes,
  logFC_Oligo = deg_oligo[common_genes, "logFC"],
  FDR_Oligo = deg_oligo[common_genes, "adj.P.Val"],
  logFC_Micro = deg_micro[common_genes, "logFC"],
  FDR_Micro = deg_micro[common_genes, "adj.P.Val"],
  logFC_Astro = deg_astro[common_genes, "logFC"],
  FDR_Astro = deg_astro[common_genes, "adj.P.Val"]
)

tabla_comun_top <- tabla_comun[order(tabla_comun$FDR_Oligo), ][1:20, ]

tabla_comun_redondeada <- tabla_comun_top %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

write.csv(tabla_comun_redondeada, "tabla_comun_20.csv", row.names = FALSE)

# Convertir a ENTREZ IDs 
entrez_ids <- bitr(intersect_all,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

entrez_vector <- entrez_ids$ENTREZID

go_enrichment <- enrichGO(gene = entrez_vector,
                          OrgDb = org.Hs.eg.db,
                          ont = "BP",  # "MF" o "CC" también posibles
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          readable = TRUE)

kegg_enrichment <- enrichKEGG(gene = entrez_vector,
                              organism = "hsa",
                              pvalueCutoff = 0.05)

# GO dotplot
dotplot(go_enrichment, showCategory = 15) + ggtitle("Análisis de enriquecimiento GO (biological process")

---------------------------- ----------------------------
  ----- # METODOS NO SUPERVISADOS Y SUPERVISADOS # -----

#PCA por tipo celular
#oligo
log_counts_t_oligo <- t(log_counts_oligo)
varianzas_gen_oligo <- apply(log_counts_t_oligo, 2, var)
log_counts_t_filtered_oligo <- log_counts_t_oligo[, varianzas_gen_oligo > 0]

pca_res <- prcomp(log_counts_t_filtered_oligo, scale. = TRUE)

autoplot(pca_res, data = data.frame(Grupo = group_labels_oligo), colour = 'Grupo') + 
  ggtitle("PCA - Oligodendrocitos") +
  theme_minimal()

#microglía
log_counts_t_micro <- t(log_counts_micro)
varianzas_gen_micro <- apply(log_counts_t_micro, 2, var)
log_counts_t_filtered_micro <- log_counts_t_micro[, varianzas_gen_micro > 0]

pca_res_micro <- prcomp(log_counts_t_filtered_micro, scale. = TRUE)

autoplot(pca_res_micro, data = data.frame(Grupo = group_labels_micro), colour = 'Grupo') + 
  ggtitle("PCA - Microglía") +
  theme_minimal()

#astrocito
log_counts_t_astro <- t(log_counts_astro)
varianzas_gen_astro <- apply(log_counts_t_astro, 2, var)
log_counts_t_filtered_astro <- log_counts_t_astro[, varianzas_gen_astro > 0]

pca_res_astro <- prcomp(log_counts_t_filtered_astro, scale. = TRUE)

autoplot(pca_res_astro, data = data.frame(Grupo = group_labels_astro), colour = 'Grupo') + 
  ggtitle("PCA - Astrocito") +
  theme_minimal()

#PCA all
#todos los tipos celulares
common_genes <- Reduce(intersect, list(
  rownames(counts_oligo_pseudobulk),
  rownames(counts_micro_pseudobulk),
  rownames(counts_astro_pseudobulk)
))

# Subset con genes comunes
oligo_mat <- counts_oligo_pseudobulk[common_genes, ]
micro_mat <- counts_micro_pseudobulk[common_genes, ]
astro_mat <- counts_astro_pseudobulk[common_genes, ]

#Combinar columnas
expr_matrix_all <- cbind(oligo_mat, micro_mat, astro_mat)

#Log-transformar
log_expr_matrix_all <- log2(expr_matrix_all + 1)

#transponemos para PCA
log_expr_matrix_t <- t(log_expr_matrix_all)

varianzas_gen_all <- apply(log_expr_matrix_t, 2, var)
log_expr_matrix_t_filtered <- log_expr_matrix_t[, varianzas_gen_all > 0]

pca_all <- prcomp(log_expr_matrix_t_filtered, scale. = TRUE)

#extraemos grupos y etiquetas
grupos <- c(
  rep("Oligo", ncol(oligo_mat)),
  rep("Microglia", ncol(micro_mat)),
  rep("Astro", ncol(astro_mat))
)

condiciones <- c(
  ifelse(grepl("AD", colnames(oligo_mat)), "EA", "Control"),
  ifelse(grepl("AD", colnames(micro_mat)), "EA", "Control"),
  ifelse(grepl("AD", colnames(astro_mat)), "EA", "Control")
)

pca_df <- as.data.frame(pca_all$x)
pca_df$TipoCelular <- grupos
pca_df$Condicion <- condiciones

png("PCA_all.png", width = 800, height = 800)
ggplot(pca_df, aes(x = PC1, y = PC2, color = TipoCelular, shape = Condicion)) +
  geom_point(size = 3) +
  labs(title = "PCA conjunto - Tipos celulares (pseudobulk)",
       x = paste0("PC1 (", round(summary(pca_all)$importance[2, 1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca_all)$importance[2, 2] * 100, 2), "%)")) +
  theme_minimal() +
  scale_color_manual(values = c("Oligo" = "darkorchid1", "Microglia" = "darkolivegreen3", "Astro" = "coral1"))
dev.off()


#UMAP
# Ejecutar UMAP (necesita samples como filas)
# Usamos la misma expr_matrix que en PCA

set.seed(42)

umap_res <- umap(log_expr_matrix_t_filtered)

umap_df <- data.frame(
  UMAP1 = umap_res[, 1],
  UMAP2 = umap_res[, 2],
  TipoCelular = grupos,
  Condicion = condiciones
)

# Gráfico
png("UMAP_all.png", width = 800, height = 800)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = TipoCelular, shape = Condicion)) +
  geom_point(size = 3) +
  labs(title = "UMAP conjunto - Tipos celulares (pseudobulk)",
       x = "UMAP1", y = "UMAP2") +
  theme_minimal() +
  scale_color_manual(values = c("Oligo" = "#B03060", "Microglia" = "#00BFFF", "Astro" = "#FFD700")) +
  scale_shape_manual(values = c("Control" = 16, "EA" = 17))  # círculos y triángulos
dev.off()

#datos para los modelos
expr_matrix_all <- cbind(oligo_mat, micro_mat, astro_mat)
expr_matrix_all <- log2(expr_matrix_all + 1)

grupos <- c(
  rep("Oligo", ncol(oligo_mat)),
  rep("Microglia", ncol(micro_mat)),
  rep("Astro", ncol(astro_mat))
)

condiciones <- c(
  ifelse(grepl("AD", colnames(oligo_mat)), "EA", "Control"),
  ifelse(grepl("AD", colnames(micro_mat)), "EA", "Control"),
  ifelse(grepl("AD", colnames(astro_mat)), "EA", "Control")
)

gene_var <- apply(expr_matrix_all, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:200]

expr_subset <- expr_matrix_all[top_genes, ]
expr_df <- as.data.frame(t(expr_subset))

expr_df$Condicion <- condiciones

colnames(expr_df) <- make.names(colnames(expr_df))

dim(expr_df)
table(expr_df$Condicion)

#'Condicion' es factor
expr_df$Condicion <- as.factor(expr_df$Condicion)

#validación cruzada
set.seed(123)

ctrl <- trainControl(
method = "repeatedcv", 
number = 5,
repeats = 20,
classProbs = TRUE,
savePredictions = "final",
summaryFunction = twoClassSummary
)

#rf
model_rf <- train(Condicion ~ ., data = expr_df, method = "rf", trControl = ctrl, metric = "ROC")

#(SVM) lineal
model_svm <- train(Condicion ~ ., data = expr_df, method = "svmLinear", trControl = ctrl, metric = "ROC")

# Regresión logística penalizada (Lasso/GLM)
model_glm <- train(Condicion ~ ., data = expr_df, method = "glmnet", trControl = ctrl, metric = "ROC")

pred_rf <- model_rf$pred
pred_svm <- model_svm$pred
pred_glm <- model_glm$pred

pred_rf$obs <- factor(pred_rf$obs, levels = c("Control", "EA"))
pred_svm$obs <- factor(pred_svm$obs, levels = c("Control", "EA"))
pred_glm$obs <- factor(pred_glm$obs, levels = c("Control", "EA"))

roc_rf <- roc(pred_rf$obs, pred_rf$AD)
roc_svm <- roc(pred_svm$obs, pred_svm$AD)
roc_glm <- roc(pred_glm$obs, pred_glm$AD)

print(auc(roc_rf))
print(auc(roc_svm))
print(auc(roc_glm))

plot(roc_rf, col = "#AB82FF", main = "Curvas ROC - Modelos Comparados", legacy.axes = TRUE)
lines(roc_svm, col = "#B3EE3A")
lines(roc_glm, col = "palevioletred2")
legend("bottomright",
       legend = c(
         paste("Random Forest (AUC =", round(auc(roc_rf), 2), ")"),
         paste("SVM (AUC =", round(auc(roc_svm), 2), ")"),
         paste("GLMnet (AUC =", round(auc(roc_glm), 2), ")")
       ),
       col = c("#AB82FF", "#B3EE3A", "palevioletred2"), lwd = 2)

varImp(model_rf)
varImp(model_svm)
varImp(model_glm)

confusionMatrix(pred_rf$pred, pred_rf$obs)
confusionMatrix(pred_svm$pred, pred_svm$obs)
confusionMatrix(pred_glm$pred, pred_glm$obs)

#genes clave
rf_importance <- varImp(model_rf)
plot(rf_importance, top = 20, main = "Importancia de genes - Random Forest")

svm_importance <- varImp(model_svm)
plot(svm_importance, top = 20, main = "Importancia de genes - SVM")

coef(model_glm$finalModel, model_glm$bestTune$lambda)

imp_rf <- varImp(model_rf)$importance
imp_svm <- varImp(model_svm)$importance

genes_rf <- rownames(imp_rf)
genes_svm <- rownames(imp_svm)

common_rf_genes <- intersect(genes_rf, common_all_three)
common_svm_genes <- intersect(genes_svm, common_all_three)

common_rf_genes
common_svm_genes
