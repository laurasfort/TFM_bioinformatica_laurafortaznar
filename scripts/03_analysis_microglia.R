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

group_labels_micro <- ifelse(grepl("AD", colnames(counts_micro_pseudobulk)), "AD", "Control")
group_micro <- factor(group_labels_micro)
design_micro <- model.matrix(~ 0 + group_micro)
colnames(design_micro) <- levels(group_micro)

log_counts_micro <- log2(counts_micro_pseudobulk + 1)

fit_micro <- lmFit(log_counts_micro, design_micro)
contrast_micro <- makeContrasts(AD - Control, levels = design_micro)
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

# Subset y escalado
heat_data_micro <- log_counts_micro[top20_genes_micro, ]
heat_data_scaled_micro <- t(scale(t(heat_data_micro)))

annotation_col_micro <- data.frame(
  Condicion = group_labels_micro,
  Paciente = colnames(log_counts_micro)
)
rownames(annotation_col_micro) <- colnames(log_counts_micro)

pdf("Heatmap_EA_vs_ct_Microglia.pdf", width = 10, height = 8)
pheatmap(heat_data_scaled_micro,
         annotation_col = annotation_col_micro,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize = 6,
         main = "Heatmap – Top 20 genes Microglía")
dev.off()

biomarcadores_micro <- results_micro %>%
  filter(adj.P.Val < 0.01 & abs(logFC) > 2) %>%
  arrange(desc(abs(logFC))) 

biomarcadores_micro <- biomarcadores_micro[1:20, ]

rownames(biomarcadores_micro)

head(biomarcadores_micro, 20)
