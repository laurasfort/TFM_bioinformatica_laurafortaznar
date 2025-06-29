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
