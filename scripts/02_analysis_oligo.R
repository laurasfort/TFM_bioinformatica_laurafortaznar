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
group_labels_oligo <- ifelse(grepl("AD", colnames(counts_oligo_pseudobulk)), "AD", "Control")
group_oligo <- factor(group_labels_oligo)
design_oligo <- model.matrix(~ 0 + group_oligo)
colnames(design_oligo) <- levels(group_oligo)

#normalizar
log_counts_oligo <- log2(counts_oligo_pseudobulk + 1)

#análisis con limma
fit_oligo <- lmFit(log_counts_oligo, design_oligo)
contrast_oligo <- makeContrasts(AD - Control, levels = design_oligo)
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

#anotación por paciente y condición
annotation_col_oligo <- data.frame(
  Condicion = group_labels_oligo,
  Paciente = colnames(log_counts_oligo)
)
rownames(annotation_col_oligo) <- colnames(log_counts_oligo)

# Subset expresión para los top 50 genes y escalado
genes_heatmap_oligo <- rownames(top20_genes_oligo)
heat_data_oligo <- log_counts_oligo[genes_heatmap_oligo, ]
heat_data_scaled_oligo <- t(scale(t(heat_data_oligo)))

pdf("Heatmap_EA_vs_ct_Oligos.pdf", width = 10, height = 8)
pheatmap(heat_data_scaled_oligo,
         annotation_col = annotation_col_oligo,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE,
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
