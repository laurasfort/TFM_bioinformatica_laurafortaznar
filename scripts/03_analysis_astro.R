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

group_labels_astro <- ifelse(grepl("AD", colnames(counts_astro_pseudobulk)), "AD", "Control")
group_astro <- factor(group_labels_astro)
design_astro <- model.matrix(~ 0 + group_astro)
colnames(design_astro) <- levels(group_astro)

log_counts_astro <- log2(counts_astro_pseudobulk + 1)

fit_astro <- lmFit(log_counts_astro, design_astro)
contrast_astro <- makeContrasts(AD - Control, levels = design_astro)
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

# Subset expresión para los top 50 genes
variable_genes_astro <- apply(log_counts_astro, 1, var)
top_20_genes_astro <- names(sort(variable_genes_astro, decreasing = TRUE))[1:20]

#anotación por paciente y condición
annotation_col_astro <- data.frame(
  Condicion = group_labels_astro,
  Paciente = colnames(log_counts_astro)
)
rownames(annotation_col_astro) <- colnames(log_counts_astro)

# Subset expresión para los top 50 genes y escalado
genes_heatmap_astro <- top_20_genes_astro
heat_data_astro <- log_counts_astro[genes_heatmap_astro, ]
heat_data_scaled_astro <- t(scale(t(heat_data_astro)))

pdf("Heatmap_EA_vs_ct_Astro.pdf", width = 10, height = 8)
pheatmap(heat_data_scaled_astro,
         annotation_col = annotation_col_astro,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize = 5,
         main = "Heatmap - Top 20 genes Astro")
dev.off()

biomarcadores_astro <- results_astro %>%
  filter(adj.P.Val < 0.01 & abs(logFC) > 2) %>%
  arrange(desc(abs(logFC))) 

biomarcadores_astro <- biomarcadores_astro[1:20, ]

rownames(biomarcadores_astro)

head(biomarcadores_astro, 20)
