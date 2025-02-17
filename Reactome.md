# Reactome analysis for intestinal organoids under different treatment conditions (R script)
## SPI1 vs Untreated
```
# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(reactome.db)

# Convert Log2 fold change and FDR p-value to numeric
Updated_SPI1vsUninfectedsheet1$`Log2 fold change` <- as.numeric(as.character(Updated_SPI1vsUninfectedsheet1$`Log2 fold change`))
Updated_SPI1vsUninfectedsheet1$`FDR p-value` <- as.numeric(as.character(Updated_SPI1vsUninfectedsheet1$`FDR p-value`))

# Remove rows with NA values
data_clean_SPI1 <- Updated_SPI1vsUninfectedsheet1[!is.na(Updated_SPI1vsUninfectedsheet1$`Log2 fold change`) & !is.na(Updated_SPI1vsUninfectedsheet1$`FDR p-value`), ]

# Create a named vector of log2 fold changes
gene_list_SPI1 <- setNames(data_clean$`Log2 fold change`, data_clean$Name)

# Convert gene symbols to Entrez IDs
gene_ids_SPI1 <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Keep only the genes that were successfully mapped to Entrez IDs
gene_list_SPI1 <- gene_list_SPI1[names(gene_list_SPI1) %in% gene_ids_SPI1$SYMBOL]
gene_ids_SPI1 <- gene_ids_SPI1[match(names(gene_list_SPI1), gene_ids_SPI1$SYMBOL),]

# Create a named vector of log2 fold changes with Entrez IDs as names
gene_list_entrez_SPI1 <- setNames(gene_list_SPI1, gene_ids_SPI1$ENTREZID)

# Perform GO enrichment analysis
go_result_SPI1 <- enrichGO(gene = names(gene_list_entrez_SPI1),
                      universe = names(gene_list_entrez_SPI1),
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)

# Convert the result to a data frame
print(class(go_result_SPI1))

go_df_SPI1 <- as.data.frame(go_result_SPI1)
if (!is.data.frame(go_result_SPI1)) {
  go_df_alt_SPI1 <- go_result@result
  print(class(go_df_alt_SPI1))
  print(str(go_df_alt_SPI1))
}

# Display the top 20 enriched pathways
print(head(go_df_alt_SPI1, 20))

# Create a dot plot of enriched pathways
dotplot(go_df_alt_SPI1, showCategory = 20) +
  ggtitle("Top 20 Enriched GO Biological Processes SPI1 VS Uninfected")

ggsave("go_enrichment_dotplot.png", width = 12, height = 10)

# Create an enrichment map
emapplot(ggo_df_alt, showCategory = 20) +
  ggtitle("Enrichment Map of GO Biological Processes")

ggsave("go_enrichment_emapplot.png", width = 12, height = 10)

# Separate upregulated and downregulated genes
up_genes_SPI1 <- names(gene_list[gene_list_SPI1 > 0])
down_genes_SPI1 <- names(gene_list[gene_list_SPI1 < 0])


# Convert gene symbols to Entrez IDs using a more robust method
entrez_ids_GOUP_SPI1 <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = up_genes_SPI1,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL")
# Convert gene symbols to Entrez IDs using a more robust method
entrez_ids_GODOWN_SPI1 <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = down_genes_SPI1,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL")
if (nrow(entrez_ids_GODOWN_SPI1) == 0) {
  stop("No gene symbols could be converted to Entrez IDs. Please check the gene symbols.")
}
# Perform GO enrichment analysis for upregulated genes
go_up_SPI1 <- enrichPathway(gene =  entrez_ids_GOUP_SPI1$ENTREZID, 
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

# Perform GO enrichment analysis for downregulated genes
go_down_SPI1 <- enrichPathway(gene = entrez_ids_GODOWN_SPI1$ENTREZID,
                  pAdjustMethod = "BH",
                    readable = TRUE)


SPI1DOTPLOT20 <- dotplot(go_up_SPI1, showCategory = 20, font.size = 20) + ggtitle("Top 20 Upregulated GO Biological Processes SPI1 vs Uninfected") + theme(plot.title = element_text(size = 25, face = "bold"), axis.text = element_text(size = 16),  # Increase axis text size
           legend.text = element_text(size = 16),  # Increase legend text size
           legend.title = element_text(size = 16))
ggsave("Top 20 Upregulated GO Biological Processes SPI1 vs Uninfected.png", SPI1DOTPLOT20 , width = 16, height = 16)
 #dotplot(go_down_SPI1, showCategory = 10, title = "Top 10 Downregulated GO Biological Processes")
 

# Combine the plots
#library(patchwork)
#combined_plot_SPI1 <- p1 / p2

# Save the combined plot
ggsave("up_down_regulated_pathways_SPI1.png", combined_plot, width = 12, height = 16)

#print("Analysis complete. Check the generated PNG files for visualizations.")

```
## SP2 vs Untreated

```
# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)


# Convert Log2 fold change and FDR p-value to numeric
Updated_SPI2_vs_UninfectedSheet1$`Log2 fold change` <- as.numeric(as.character(Updated_SPI2_vs_UninfectedSheet1$`Log2 fold change`))
Updated_SPI2_vs_UninfectedSheet1$`FDR p-value` <- as.numeric(as.character(Updated_SPI2_vs_UninfectedSheet1$`FDR p-value`))

# Remove rows with NA values
data_clean_SPI2 <- Updated_SPI2_vs_UninfectedSheet1[!is.na(Updated_SPI2_vs_UninfectedSheet1$`Log2 fold change`) & !is.na(Updated_SPI2_vs_UninfectedSheet1$`FDR p-value`), ]

# Create a named vector of log2 fold changes
gene_list_SPI2 <- setNames(data_clean_SPI2$`Log2 fold change`, data_clean_SPI2$Name)

# Convert gene symbols to Entrez IDs
gene_ids_SPI2 <- bitr(names(gene_list_SPI2), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Keep only the genes that were successfully mapped to Entrez IDs
gene_list_SPI2 <- gene_list_SPI2[names(gene_list_SPI2) %in% gene_ids_SPI2$SYMBOL]
gene_ids_SPI2 <- gene_ids_SPI2[match(names(gene_list_SPI2), gene_ids_SPI2$SYMBOL),]

# Create a named vector of log2 fold changes with Entrez IDs as names
gene_list_entrez_SPI2 <- setNames(gene_list_SPI2, gene_ids_SPI2$ENTREZID)

# Perform GO enrichment analysis
go_result_SPI2 <- enrichGO(gene = names(gene_list_entrez_SPI2),
                      universe = names(gene_list_entrez_SPI2),
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)

# Convert the result to a data frame
print(class(go_result_SPI1))

go_df_SPI2 <- as.data.frame(go_result_SPI2)
if (!is.data.frame(go_result_SPI2)) {
  go_df_alt_SPI2 <- go_result@result
  print(class(go_df_alt_SPI2))
  print(str(go_df_alt_SPI2))
}

# Display the top 20 enriched pathways
print(head(go_df_alt_SPI2, 20))

# Create a dot plot of enriched pathways
dotplot(go_df_alt_SPI2, showCategory = 20) +
  ggtitle("Top 20 Enriched GO Biological Processes SPI1 VS Uninfected")

ggsave("go_enrichment_dotplot.png", width = 12, height = 10)

# Create an enrichment map
#emapplot(ggo_df_alt, showCategory = 20) +
  #ggtitle("Enrichment Map of GO Biological Processes")

#ggsave("go_enrichment_emapplot.png", width = 12, height = 10)

# Separate upregulated and downregulated genes
up_genes_SPI2 <- names(gene_list[gene_list_SPI2 > 0])
down_genes_SPI2 <- names(gene_list[gene_list_SPI2 < 0])


# Convert gene symbols to Entrez IDs using a more robust method
entrez_ids_GOUP_SPI2 <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = up_genes_SPI2,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL")
# Convert gene symbols to Entrez IDs using a more robust method
entrez_ids_GODOWN_SPI2 <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = down_genes_SPI2,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL")
if (nrow(entrez_ids_GODOWN_SPI2) == 0) {
  stop("No gene symbols could be converted to Entrez IDs. Please check the gene symbols.")
}
# Perform GO enrichment analysis for upregulated genes
go_up_SPI2 <- enrichPathway(gene =  entrez_ids_GOUP_SPI2$ENTREZID, 
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

# Perform GO enrichment analysis for downregulated genes
go_down_SPI2 <- enrichPathway(gene = entrez_ids_GODOWN_SPI2$ENTREZID,
                    pAdjustMethod = "BH",
                    readable = TRUE)



SPI2DOTPLOT20 <- dotplot(go_up_SPI2, showCategory = 20, font.size = 20) + ggtitle("Top 20 Upregulated GO Biological Processes SPI2 vs Uninfected") + theme(plot.title = element_text(size = 25, face = "bold"), axis.text = element_text(size = 16),  # Increase axis text size
           legend.text = element_text(size = 16),  # Increase legend text size
           legend.title = element_text(size = 16))
ggsave("Top 20 Upregulated GO Biological Processes SPI2 vs Uninfected.png", SPI2DOTPLOT20 , width = 16, height = 16)
 #dotplot(go_down_SPI1, showCategory = 10, title = "Top 10 Downregulated GO Biological Processes")
 

# Combine the plots
#library(patchwork)
#combined_plot_SPI1 <- p1 / p2

# Save the combined plot
#ggsave("up_down_regulated_pathways_SPI1.png", combined_plot, width = 12, height = 16)

#print("Analysis complete. Check the generated PNG files for visualizations.")

```
# DM_vs_EM

```
# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)


# Convert Log2 fold change and FDR p-value to numeric
Updated_LynnEnteroidP25DM_vs_EMsheet1$`Log2 fold change` <- as.numeric(as.character(Updated_LynnEnteroidP25DM_vs_EMsheet1$`Log2 fold change`))
Updated_LynnEnteroidP25DM_vs_EMsheet1$`FDR p-value` <- as.numeric(as.character(Updated_LynnEnteroidP25DM_vs_EMsheet1$`FDR p-value`))

# Remove rows with NA values
data_clean_DM_vs_EM <- Updated_LynnEnteroidP25DM_vs_EMsheet1[!is.na(Updated_LynnEnteroidP25DM_vs_EMsheet1$`Log2 fold change`) & !is.na(Updated_LynnEnteroidP25DM_vs_EMsheet1$`FDR p-value`), ]

# Create a named vector of log2 fold changes
gene_list_DM_vs_EM <- setNames(data_clean_DM_vs_EM$`Log2 fold change`, data_clean_DM_vs_EM$Name)

# Convert gene symbols to Entrez IDs
gene_ids_DM_vs_EM <- bitr(names(gene_list_DM_vs_EM), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Keep only the genes that were successfully mapped to Entrez IDs
gene_list_DM_vs_EM <- gene_list_DM_vs_EM[names(gene_list_DM_vs_EM) %in% gene_ids_DM_vs_EM$SYMBOL]
gene_ids_DM_vs_EM <- gene_ids_DM_vs_EM[match(names(gene_list_DM_vs_EM), gene_ids_DM_vs_EM$SYMBOL),]

# Create a named vector of log2 fold changes with Entrez IDs as names
gene_list_entrez_DM_vs_EM <- setNames(gene_list_DM_vs_EM, gene_ids_DM_vs_EM$ENTREZID)

# Perform GO enrichment analysis
go_result_DM_vs_EM <- enrichpathway(gene = names(gene_list_entrez_DM_vs_EM),
                      universe = names(gene_list_entrez_DM_vs_EM),
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)

# Convert the result to a data frame
print(class(go_result_DM_vs_EM))

go_df_DM_vs_EM <- as.data.frame(go_result_DM_vs_EM)
if (!is.data.frame(go_result_DM_vs_EM)) {
  go_df_alt_DM_vs_EM <- go_result@result
  print(class(go_df_alt_DM_vs_EM))
  print(str(go_df_alt_DM_vs_EM))
}

# Display the top 20 enriched pathways
print(head(go_df_alt_DM_vs_EM, 20))

# Create a dot plot of enriched pathways
dotplot(go_df_alt_DM_vs_EM, showCategory = 20) +
  ggtitle("Top 20 Enriched GO Biological Processes SPI1 VS Uninfected")

ggsave("go_enrichment_dotplot.png", width = 12, height = 10)

# Create an enrichment map
#emapplot(ggo_df_alt, showCategory = 20) +
  #ggtitle("Enrichment Map of GO Biological Processes")

#ggsave("go_enrichment_emapplot.png", width = 12, height = 10)

# Separate upregulated and downregulated genes
up_genes_DM_vs_EM <- names(gene_list_DM_vs_EM[gene_list_DM_vs_EM > 0])
down_genes_DM_vs_EM <- names(gene_list_DM_vs_EM[gene_list_DM_vs_EM < 0])


# Convert gene symbols to Entrez IDs using a more robust method
entrez_ids_GOUP_DM_vs_EM <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = up_genes_DM_vs_EM,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL")
#Convert gene symbols to Entrez IDs using a more robust method
entrez_ids_GODOWN_DM_vs_EM <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = down_genes_DM_vs_EM,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL")

# Perform GO enrichment analysis for upregulated genes
go_up_DM_vs_EM <- enrichPathway(gene =  entrez_ids_GOUP_DM_vs_EM$ENTREZID, 
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

 

#Perform GO enrichment analysis for downregulated genes
go_down_DM_vs_EM <- enrichPathway(gene = entrez_ids_GODOWN_DM_vs_EM$ENTREZID,
                    pAdjustMethod = "BH",
                    readable = TRUE)



DM_vs_EMDOTPLOT20 <- dotplot(go_up_DM_vs_EM, showCategory = 20, font.size = 20) + ggtitle("Top 20 Upregulated GO Biological Processes DM vs EM ") + theme(plot.title = element_text(size = 25, face = "bold"), axis.text = element_text(size = 16),  # Increase axis text size
           legend.text = element_text(size = 16),  # Increase legend text size
           legend.title = element_text(size = 16))
ggsave("Top 20 Upregulated GO Biological Processes DM_vs_EM.png", DM_vs_EMDOTPLOT20 , width = 16, height = 16) 

DM_vs_EMDOTPLOT20down <- dotplot(go_down_DM_vs_EM, showCategory = 20, font.size = 17) + ggtitle("Top 20 Down-regulated GO Biological Processes DM vs EM ") + theme(plot.title = element_text(size = 25, face = "bold"), axis.text = element_text(size = 16),  # Increase axis text size
           legend.text = element_text(size = 16),  # Increase legend text size
           legend.title = element_text(size = 16))
ggsave("Top 20 Down-regulated GO Biological Processes DM vs EM.png", DM_vs_EMDOTPLOT20down , width = 16, height = 16) 


```
# Dome vs Transwell

```
library(readxl)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)


# Convert Log2 fold change and FDR p-value to numeric
Updated_Lynn_DM_Dome_V_TranswellCLCsheet1$`Log2 fold change` <- as.numeric(as.character(Updated_Lynn_DM_Dome_V_TranswellCLCsheet1$`Log2 fold change`))
Updated_Lynn_DM_Dome_V_TranswellCLCsheet1$`FDR p-value` <- as.numeric(as.character(Updated_Lynn_DM_Dome_V_TranswellCLCsheet1$`FDR p-value`))

# Remove rows with NA values
data_clean_DM_Dome_V_TranswellCLC <- Updated_Lynn_DM_Dome_V_TranswellCLCsheet1[!is.na(Updated_Lynn_DM_Dome_V_TranswellCLCsheet1$`Log2 fold change`) & !is.na(Updated_Lynn_DM_Dome_V_TranswellCLCsheet1$`FDR p-value`), ]

# Create a named vector of log2 fold changes
gene_list_DM_Dome_V_TranswellCLC <- setNames(data_clean_DM_Dome_V_TranswellCLC$`Log2 fold change`, data_clean_DM_Dome_V_TranswellCLC$Name)

# Convert gene symbols to Entrez IDs
gene_ids_DM_Dome_V_TranswellCLC <- bitr(names(gene_list_DM_Dome_V_TranswellCLC), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Keep only the genes that were successfully mapped to Entrez IDs
gene_list_DM_Dome_V_TranswellCLC <- gene_list_DM_Dome_V_TranswellCLC[names(gene_list_DM_Dome_V_TranswellCLC) %in% gene_ids_DM_Dome_V_TranswellCLC$SYMBOL]
gene_ids_DM_Dome_V_TranswellCLC <- gene_ids_DM_Dome_V_TranswellCLC[match(names(gene_list_DM_Dome_V_TranswellCLC), gene_ids_DM_Dome_V_TranswellCLC$SYMBOL),]

# Create a named vector of log2 fold changes with Entrez IDs as names
gene_list_entrez_DM_Dome_V_TranswellCLC <- setNames(gene_list_DM_Dome_V_TranswellCLC, gene_ids_DM_Dome_V_TranswellCLC$ENTREZID)

# Perform GO enrichment analysis
go_result_DM_Dome_V_TranswellCLC <- enrichpathway(gene = names(gene_list_entrez_DM_Dome_V_TranswellCLC),
                      universe = names(gene_list_entrez_DM_Dome_V_TranswellCLC),
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)

# Convert the result to a data frame
print(class(go_result_DM_Dome_V_TranswellCLC))

go_df__DM_Dome_V_TranswellCLC <- as.data.frame(go_result_DM_Dome_V_TranswellCLC)
if (!is.data.frame(go_result_DM_Dome_V_TranswellCLC)) {
  go_df_alt__DM_Dome_V_TranswellCLC <- go_result@result
  print(class(go_df_alt_DM_Dome_V_TranswellCLC))
  print(str(go_df_alt_DM_Dome_V_TranswellCLC))
}

# Display the top 20 enriched pathways
print(head(go_df_alt_DM_Dome_V_TranswellCLC, 20))

# Create a dot plot of enriched pathways
dotplot(go_df_alt_DM_Dome_V_TranswellCLC, showCategory = 20) +
  ggtitle("Top 20 Enriched GO Biological Processes SPI1 VS Uninfected")

ggsave("go_enrichment_dotplot.png", width = 12, height = 10)

# Create an enrichment map
#emapplot(ggo_df_alt, showCategory = 20) +
  #ggtitle("Enrichment Map of GO Biological Processes")

#ggsave("go_enrichment_emapplot.png", width = 12, height = 10)

# Separate upregulated and downregulated genes
up_genes_DM_Dome_V_TranswellCLC <- names(gene_list_DM_Dome_V_TranswellCLC[gene_list_DM_Dome_V_TranswellCLC > 0])
down_genes_DM_Dome_V_TranswellCLC <- names(gene_list_DM_Dome_V_TranswellCLC[gene_list_DM_Dome_V_TranswellCLC < 0])


# Convert gene symbols to Entrez IDs using a more robust method
entrez_ids_GOUP_DM_Dome_V_TranswellCLC <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = up_genes_DM_Dome_V_TranswellCLC,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL")
#Convert gene symbols to Entrez IDs using a more robust method
entrez_ids_GODOWN_DM_Dome_V_TranswellCLC <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = down_genes_DM_Dome_V_TranswellCLC,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL")

# Perform GO enrichment analysis for upregulated genes
go_up_DM_Dome_V_TranswellCLC <- enrichPathway(gene =  entrez_ids_GOUP_DM_Dome_V_TranswellCLC$ENTREZID, 
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

 

#Perform GO enrichment analysis for downregulated genes
go_down_DM_Dome_V_TranswellCLC<- enrichPathway(gene = entrez_ids_GODOWN_DM_Dome_V_TranswellCLC$ENTREZID,
                    pAdjustMethod = "BH",
                    readable = TRUE)



DM_Dome_V_TranswellCLCDOTPLOT20 <- dotplot(go_up_DM_Dome_V_TranswellCLC, showCategory = 20, font.size = 20) + ggtitle("Top 20 Upregulated GO Biological Processes DM Dome Vs Transwell") + theme(plot.title = element_text(size = 25, face = "bold"), axis.text = element_text(size = 16),  # Increase axis text size
           legend.text = element_text(size = 16),  # Increase legend text size
           legend.title = element_text(size = 16))
ggsave("Top 20 Upregulated GO Biological Processes DM Dome Vs Transwell.png", DM_Dome_V_TranswellCLCDOTPLOT20 , width = 16, height = 16) 

DM_Dome_V_TranswellCLCDOTPLOT20down <- dotplot(go_down_DM_Dome_V_TranswellCLC, showCategory = 20, font.size = 17) + ggtitle("Top 20 Down-regulated GO Biological Processes DM Dome Vs Transwell") + theme(plot.title = element_text(size = 25, face = "bold"), axis.text = element_text(size = 16),  # Increase axis text size
           legend.text = element_text(size = 16),  # Increase legend text size
           legend.title = element_text(size = 16))
ggsave("Top 20 Down-regulated GO Biological Processes DM Dome Vs Transwell.png",DM_Dome_V_TranswellCLCDOTPLOT20down , width = 16, height = 16) 




```
##  Uninfected vs wildtype
```
# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)


# Convert Log2 fold change and FDR p-value to numeric
Updated_Uninfected_VS_Wildtypesheet1$`Log2 fold change` <- as.numeric(as.character(Updated_Uninfected_VS_Wildtypesheet1$`Log2 fold change`))
Updated_Uninfected_VS_Wildtypesheet1$`FDR p-value` <- as.numeric(as.character(Updated_Uninfected_VS_Wildtypesheet1$`FDR p-value`))

# Remove rows with NA values
data_clean_Uninfected_VS_Wildtype <- Updated_Uninfected_VS_Wildtypesheet1[!is.na(Updated_Uninfected_VS_Wildtypesheet1$`Log2 fold change`) & !is.na(Updated_Uninfected_VS_Wildtypesheet1$`FDR p-value`), ]

# Create a named vector of log2 fold changes
gene_list_Uninfected_VS_Wildtype <- setNames(data_clean_Uninfected_VS_Wildtype$`Log2 fold change`, data_clean_Uninfected_VS_Wildtype$Name)

# Convert gene symbols to Entrez IDs
gene_ids_Uninfected_VS_Wildtype <- bitr(names(gene_list_Uninfected_VS_Wildtype), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Keep only the genes that were successfully mapped to Entrez IDs
gene_list_Uninfected_VS_Wildtype <- gene_list_Uninfected_VS_Wildtype[names(gene_list_Uninfected_VS_Wildtype) %in% gene_ids_SPI1$SYMBOL]
gene_ids_Uninfected_VS_Wildtype <- gene_ids_Uninfected_VS_Wildtype[match(names(gene_list_Uninfected_VS_Wildtype), gene_ids_Uninfected_VS_Wildtype$SYMBOL),]

# Create a named vector of log2 fold changes with Entrez IDs as names
gene_list_entrez_Uninfected_VS_Wildtype <- setNames(gene_list_Uninfected_VS_Wildtype, gene_ids_Uninfected_VS_Wildtype$ENTREZID)

# Perform GO enrichment analysis
go_result_Uninfected_VS_Wildtype <- enrichGO(gene = names(gene_list_entrez_Uninfected_VS_Wildtype),
                      universe = names(gene_list_entrez_Uninfected_VS_Wildtype),
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)

# Convert the result to a data frame
print(class(go_result_Uninfected_VS_Wildtype))

go_df_Uninfected_VS_Wildtype <- as.data.frame(go_result_Uninfected_VS_Wildtype)
if (!is.data.frame(go_result_Uninfected_VS_Wildtype)) {
  go_df_alt_Uninfected_VS_Wildtype <- go_result@result
  print(class(go_df_alt_Uninfected_VS_Wildtype))
  print(str(go_df_alt_Uninfected_VS_Wildtype))
}

# Display the top 20 enriched pathways
print(head(go_df_alt_Uninfected_VS_Wildtype, 20))

# Create a dot plot of enriched pathways
dotplot(go_df_alt_Uninfected_VS_Wildtype, showCategory = 20) +
  ggtitle("Top 20 Enriched GO Biological Processes SPI1 VS Uninfected")

ggsave("go_enrichment_dotplot.png", width = 12, height = 10)

# Create an enrichment map
emapplot(ggo_df_alt, showCategory = 20) +
  ggtitle("Enrichment Map of GO Biological Processes")

ggsave("go_enrichment_emapplot.png", width = 12, height = 10)

# Separate upregulated and downregulated genes
up_genes_Uninfected_VS_Wildtype <- names(gene_list_Uninfected_VS_Wildtype[gene_list_Uninfected_VS_Wildtype > 0])
down_genes_Uninfected_VS_Wildtype <- names(gene_list_Uninfected_VS_Wildtype[gene_list_Uninfected_VS_Wildtype < 0])


# Convert gene symbols to Entrez IDs using a more robust method
entrez_ids_GOUP_Uninfected_VS_Wildtype <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = up_genes_Uninfected_VS_Wildtype,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL")
# Convert gene symbols to Entrez IDs using a more robust method
entrez_ids_GODOWN_Uninfected_VS_Wildtype <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = down_genes_Uninfected_VS_Wildtype,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL")

# Perform GO enrichment analysis for upregulated genes
go_up_Uninfected_VS_Wildtype <- enrichPathway(gene =  entrez_ids_GOUP_Uninfected_VS_Wildtype$ENTREZID, 
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

# Perform GO enrichment analysis for downregulated genes
go_down_Uninfected_VS_Wildtype <-  enrichPathway(gene = entrez_ids_GODOWN_Uninfected_VS_Wildtype$ENTREZID,
                    pAdjustMethod = "BH",
                    readable = TRUE)
Uninfected_VS_Wildtype_DOTPLOT20Down <- dotplot(go_down_Uninfected_VS_Wildtype, showCategory = 20, font.size = 20) + ggtitle("Top 20 Down-regulated GO Biological Processes Uninfected VS Wildtype") + theme(plot.title = element_text(size = 25, face = "bold"), axis.text = element_text(size = 16),  # Increase axis text size
           legend.text = element_text(size = 16),  # Increase legend text size
           legend.title = element_text(size = 16))
ggsave("Top 20 Down-regulated GO Biological Processes Uninfected VS Wildtype.png", Uninfected_VS_Wildtype_DOTPLOT20Down , width = 16, height = 16)
# Check the result
print(class(go_down))
print(str(go_down))
# Create dot plots for upregulated and downregulated pathways

Uninfected_VS_Wildtype_DOTPLOT20 <- dotplot(go_up_Uninfected_VS_Wildtype, showCategory = 20, font.size = 20) + ggtitle("Top 20 Upregulated GO Biological Processes Uninfected VS Wildtype") + theme(plot.title = element_text(size = 25, face = "bold"), axis.text = element_text(size = 16),  # Increase axis text size
           legend.text = element_text(size = 16),  # Increase legend text size
           legend.title = element_text(size = 16))
ggsave("Top 20 Upregulated GO Biological Processes Uninfected VS Wildtype.png", Uninfected_VS_Wildtype_DOTPLOT20 , width = 16, height = 16)
 

``` 
