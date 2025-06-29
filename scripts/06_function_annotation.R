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
