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
  ifelse(grepl("AD", colnames(oligo_mat)), "AD", "Control"),
  ifelse(grepl("AD", colnames(micro_mat)), "AD", "Control"),
  ifelse(grepl("AD", colnames(astro_mat)), "AD", "Control")
)

pca_df <- as.data.frame(pca_all$x)
pca_df$TipoCelular <- grupos
pca_df$Condicion <- condiciones

ggplot(pca_df, aes(x = PC1, y = PC2, color = TipoCelular, shape = Condicion)) +
  geom_point(size = 3) +
  labs(title = "PCA conjunto - Tipos celulares (pseudobulk)",
       x = paste0("PC1 (", round(summary(pca_all)$importance[2, 1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca_all)$importance[2, 2] * 100, 2), "%)")) +
  theme_minimal() +
  scale_color_manual(values = c("Oligo" = "darkorchid1", "Microglia" = "darkolivegreen3", "Astro" = "coral1"))


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
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = TipoCelular, shape = Condicion)) +
  geom_point(size = 3) +
  labs(title = "UMAP conjunto - Tipos celulares (pseudobulk)",
       x = "UMAP1", y = "UMAP2") +
  theme_minimal() +
  scale_color_manual(values = c("Oligo" = "#B03060", "Microglia" = "#00BFFF", "Astro" = "#FFD700")) +
  scale_shape_manual(values = c("Control" = 16, "AD" = 17))  # círculos y triángulos
