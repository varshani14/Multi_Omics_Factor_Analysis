library(MOFA2)
library(dplyr)
library(ggplot2)
MOFA_B <- load_model(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/MOFA_PDAC_v2_model.hdf5"
)

MOFA_B
factor5_variance <- get_variance_explained(
  MOFA_B,
  factors = 5
)

factor5_variance
factor5_contribution <- data.frame(
  Omics = names(factor5_variance$r2_per_factor$group1["Factor5", ]),
  Variance_Explained = as.numeric(
    factor5_variance$r2_per_factor$group1["Factor5", ]
  )
)

factor5_contribution
p_factor5_contribution <- ggplot(
  factor5_contribution,
  aes(
    x = reorder(Omics, Variance_Explained),
    y = Variance_Explained,
    fill = Omics
  )
) +
  geom_col(
    width = 0.7,
    color = "grey25"
  ) +
  geom_text(
    aes(
      label = sprintf("%.2f%%", Variance_Explained)
    ),
    hjust = -0.1,
    size = 4
  ) +
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "MOFA Factor 5: Variance Explained across Omics Views",
    subtitle = "Contribution of each molecular layer to Factor 5",
    x = "Omics view",
    y = "Variance explained (%)",
    fill = "Omics view"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 17
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.position = "none"
  )

ggsave(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/Downstream_analysis/Factor5_Omics_Contribution.png",
  p_factor5_contribution,
  width = 10,
  height = 7,
  dpi = 300
)

p_factor5_contribution
# ============================================================
# FACTOR 5 — CNV FEATURE CONTRIBUTIONS
# ============================================================

# Extract Factor 5 weights for the CNV view
factor5_cnv_weights <- get_weights(
  MOFA_B,
  views = "CNV",
  factors = "Factor5",
  as.data.frame = TRUE
)

# Inspect the result
head(factor5_cnv_weights)
# ============================================================
# RANK CNV FEATURES DRIVING FACTOR 5
# ============================================================

factor5_cnv_ranked <- factor5_cnv_weights %>%
  mutate(abs_weight = abs(value)) %>%
  arrange(desc(abs_weight))

# Show the 20 strongest CNV features
head(factor5_cnv_ranked, 20)
# ============================================================
# CHECK FACTOR 5 CNV WEIGHT DISTRIBUTION
# ============================================================

# How many unique weight values are there?
length(unique(factor5_cnv_ranked$value))

# Show the strongest positive CNV weights
factor5_cnv_ranked %>%
  filter(value > 0) %>%
  arrange(desc(value)) %>%
  select(feature, value) %>%
  head(20)

# Show the strongest negative CNV weights
factor5_cnv_ranked %>%
  filter(value < 0) %>%
  arrange(value) %>%
  select(feature, value) %>%
  head(20)
# Check how many unique Factor 5 CNV weight values exist
length(unique(factor5_cnv_ranked$value))
length(unique(factor5_cnv_ranked$value))

# positive weights
...

# negative weights
...
# ============================================================
# FACTOR 5 — TOP CNV FEATURES DRIVING THE FACTOR
# ============================================================

# Select strongest positive and negative CNV features
top_cnv_positive <- factor5_cnv_weights %>%
  filter(value > 0) %>%
  arrange(desc(value)) %>%
  slice_head(n = 15)

top_cnv_negative <- factor5_cnv_weights %>%
  filter(value < 0) %>%
  arrange(value) %>%
  slice_head(n = 15)

# Combine them
top_cnv_features <- bind_rows(
  top_cnv_positive,
  top_cnv_negative
) %>%
  mutate(
    direction = ifelse(value > 0, "Positive", "Negative"),
    feature = gsub("_CNV$", "", feature)
  ) %>%
  arrange(value) %>%
  mutate(
    feature = factor(feature, levels = feature)
  )

# Save the table
write.csv(
  top_cnv_features,
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/Downstream_analysis/Factor5_Top_CNV_Features.csv",
  row.names = FALSE
)

# Plot
p_cnv <- ggplot(
  top_cnv_features,
  aes(x = value, y = feature, color = direction)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = value,
      y = feature,
      yend = feature
    ),
    linewidth = 1.2
  ) +
  geom_point(
    size = 4
  ) +
  scale_color_manual(
    values = c(
      "Positive" = "#D95F02",
      "Negative" = "#1B9E77"
    )
  ) +
  labs(
    title = "MOFA Factor 5: Strongest CNV Contributors",
    subtitle = "Top positive and negative CNV feature weights",
    x = "Factor 5 CNV weight",
    y = "CNV feature",
    color = "Contribution"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18
    ),
    plot.subtitle = element_text(
      size = 12
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 9
    ),
    legend.title = element_text(
      face = "bold"
    ),
    panel.grid.major.y = element_blank()
  )

# Display plot
p_cnv

# Save high-resolution figure
ggsave(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/Downstream_analysis/Factor5_Top_CNV_Features.png",
  p_cnv,
  width = 10,
  height = 8,
  dpi = 300
)
# ============================================================
# FACTOR 5 CNV GENE ANNOTATION
# ============================================================

# Required annotation packages
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("AnnotationDbi")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Extract clean Ensembl IDs from CNV feature names
# ------------------------------------------------------------

factor5_cnv_annotated <- factor5_cnv_ranked %>%
  mutate(
    ensembl_id = sub("_CNV$", "", feature),
    ensembl_id = sub("\\..*$", "", ensembl_id)
  )

# ------------------------------------------------------------
# 2. Map Ensembl IDs to gene symbols
# ------------------------------------------------------------

gene_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(factor5_cnv_annotated$ensembl_id),
  columns = c("ENSEMBL", "SYMBOL", "GENENAME"),
  keytype = "ENSEMBL"
)

# Remove duplicate mappings
gene_annotation <- gene_annotation %>%
  distinct(ENSEMBL, .keep_all = TRUE)

# ------------------------------------------------------------
# 3. Add gene annotation to Factor 5 CNV weights
# ------------------------------------------------------------

factor5_cnv_annotated <- factor5_cnv_annotated %>%
  left_join(
    gene_annotation,
    by = c("ensembl_id" = "ENSEMBL")
  ) %>%
  mutate(
    gene_label = ifelse(
      is.na(SYMBOL) | SYMBOL == "",
      ensembl_id,
      SYMBOL
    )
  )

# ------------------------------------------------------------
# 4. Save complete annotated CNV table
# ------------------------------------------------------------

write.csv(
  factor5_cnv_annotated,
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/Downstream_analysis/Factor5_CNV_Annotated.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Show strongest annotated CNV contributors
# ------------------------------------------------------------

factor5_cnv_top_annotated <- factor5_cnv_annotated %>%
  arrange(desc(abs_weight)) %>%
  dplyr::select(
    feature,
    ensembl_id,
    SYMBOL,
    GENENAME,
    value,
    abs_weight
  ) %>%
  head(40)

factor5_cnv_top_annotated

# ------------------------------------------------------------
# 6. Save top annotated contributors
# ------------------------------------------------------------

write.csv(
  factor5_cnv_top_annotated,
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/Downstream_analysis/Factor5_Top_CNV_Annotated.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Prepare plot: strongest positive and negative genes
# ------------------------------------------------------------

plot_data <- factor5_cnv_annotated %>%
  filter(!is.na(gene_label)) %>%
  arrange(desc(abs_weight)) %>%
  group_by(value > 0) %>%
  slice_head(n = 15) %>%
  ungroup() %>%
  mutate(
    contribution = ifelse(value > 0, "Positive", "Negative")
  )

# Order genes by weight
plot_data <- plot_data %>%
  arrange(value) %>%
  mutate(
    gene_label = factor(gene_label, levels = gene_label)
  )

# ------------------------------------------------------------
# 8. Create publication-style plot
# ------------------------------------------------------------

p_factor5_cnv_genes <- ggplot(
  plot_data,
  aes(
    x = value,
    y = gene_label,
    color = contribution
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = value,
      y = gene_label,
      yend = gene_label
    ),
    linewidth = 1
  ) +
  geom_point(
    size = 3
  ) +
  scale_color_manual(
    values = c(
      "Negative" = "#1B9E77",
      "Positive" = "#D95F02"
    )
  ) +
  labs(
    title = "MOFA Factor 5: Strongest CNV-Associated Genes",
    subtitle = "Top positive and negative CNV feature weights",
    x = "Factor 5 CNV weight",
    y = "Gene",
    color = "Contribution"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 9),
    legend.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

p_factor5_cnv_genes

# ------------------------------------------------------------
# 9. Save plot
# ------------------------------------------------------------

ggsave(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/Downstream_analysis/Factor5_CNV_Annotated_Genes.png",
  p_factor5_cnv_genes,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 10. Quick summary
# ------------------------------------------------------------

cat(
  "\nFactor 5 CNV annotation completed.\n",
  "Total CNV features:", nrow(factor5_cnv_annotated), "\n",
  "Annotated with gene symbol:",
  sum(!is.na(factor5_cnv_annotated$SYMBOL)), "\n",
  "Not mapped:",
  sum(is.na(factor5_cnv_annotated$SYMBOL)), "\n"
)
# ============================================================
# BIOLOGICAL INTERPRETATION — FACTOR 5
# STEP 1: PREPARE GENE LIST
# ============================================================

factor5_gene_list <- factor5_cnv_annotated %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  distinct(SYMBOL)

# Check the number of genes
nrow(factor5_gene_list)

# Show first 20 genes
head(factor5_gene_list, 20)
# ============================================================
# STEP 2: LOAD ENRICHMENT PACKAGE
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)
BiocManager::install("GO.db", ask = FALSE, update = FALSE)
library(clusterProfiler)
library(org.Hs.eg.db)
# ============================================================
# STEP 3: CONVERT GENE SYMBOLS TO ENTREZ IDs
# ============================================================

factor5_entrez <- bitr(
  factor5_gene_list$SYMBOL,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Check conversion
head(factor5_entrez)

# Number successfully converted
nrow(factor5_entrez)
# ============================================================
# STEP 4: GO BIOLOGICAL PROCESS ENRICHMENT
# ============================================================

factor5_go <- enrichGO(
  gene = factor5_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# Show the results
head(as.data.frame(factor5_go), 20)
# ============================================================
# STEP 4: PLOT FACTOR 5 GO ENRICHMENT
# ============================================================

factor5_go_plot <- dotplot(
  factor5_go,
  showCategory = 15
) +
  ggtitle("MOFA Factor 5: GO Biological Process Enrichment")

factor5_go_plot

ggsave(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/Downstream_analysis/Factor5_GO_Enrichment.png",
  factor5_go_plot,
  width = 10,
  height = 8,
  dpi = 300
)
# ============================================================
# STEP 5: GO ENRICHMENT BAR PLOT
# ============================================================

factor5_go_barplot <- barplot(
  factor5_go,
  showCategory = 15,
  title = "MOFA Factor 5: GO Biological Process Enrichment"
)

factor5_go_barplot

ggsave(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/Downstream_analysis/Factor5_GO_Barplot.png",
  factor5_go_barplot,
  width = 10,
  height = 8,
  dpi = 300
)
# ============================================================
# STEP 6: KEGG PATHWAY ENRICHMENT
# ============================================================

factor5_kegg <- enrichKEGG(
  gene = factor5_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# Show KEGG results
head(as.data.frame(factor5_kegg), 15)
# Install the KEGG category data required by clusterProfiler
BiocManager::install("KEGGREST", ask = FALSE, update = FALSE)

# Reload clusterProfiler
library(clusterProfiler)

# Test KEGG again
factor5_kegg <- enrichKEGG(
  gene = factor5_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)
# Download the missing KEGG category data
dir.create(
  "C:/Users/varsh/AppData/Local/clusterProfiler",
  recursive = TRUE,
  showWarnings = FALSE
)

download.file(
  "https://yulab-smu.top/clusterProfiler/kegg_category.rda",
  "C:/Users/varsh/AppData/Local/clusterProfiler/kegg_category.rda",
  mode = "wb"
)
packageVersion("clusterProfiler")
factor5_kegg <- enrichKEGG(
  gene = factor5_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05,
  use_internal_data = FALSE
)
# ============================================================
# KEGG ENRICHMENT - BYPASS CATEGORY CACHE
# ============================================================

kegg_data <- download_KEGG(
  species = "hsa",
  keggType = "KEGG",
  keyType = "kegg"
)

# Perform enrichment using the downloaded KEGG pathways
factor5_kegg <- enricher(
  gene = factor5_entrez$ENTREZID,
  TERM2GENE = kegg_data$KEGGPATHID2EXTID,
  TERM2NAME = kegg_data$KEGGPATHID2NAME,
  pvalueCutoff = 0.05
)

# Show results
head(as.data.frame(factor5_kegg), 15)
# ============================================================
# STEP 7: KEGG ENRICHMENT DOT PLOT
# ============================================================

factor5_kegg_plot <- dotplot(
  factor5_kegg,
  showCategory = 15,
  title = "MOFA Factor 5: KEGG Pathway Enrichment"
)

factor5_kegg_plot
ggsave(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/Downstream_analysis/Factor5_KEGG_Dotplot.png",
  factor5_kegg_plot,
  width = 10,
  height = 8,
  dpi = 300
)
# ============================================================
# STEP 1: PREPARE RANKED GENE LIST FOR GSEA
# ============================================================

factor5_gsea_genes <- factor5_cnv_annotated %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(
    weight = mean(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(weight))

# Check the ranked gene list
head(factor5_gsea_genes, 20)

# Number of genes
nrow(factor5_gsea_genes)
# ============================================================
# STEP 2: CREATE RANKED GENE VECTOR FOR GSEA
# ============================================================

factor5_gene_rank <- factor5_gsea_genes$weight
names(factor5_gene_rank) <- factor5_gsea_genes$SYMBOL

# Remove any missing or non-finite values
factor5_gene_rank <- factor5_gene_rank[
  !is.na(factor5_gene_rank) &
    is.finite(factor5_gene_rank)
]

# Sort from strongest positive to strongest negative
factor5_gene_rank <- sort(factor5_gene_rank, decreasing = TRUE)

# Check
head(factor5_gene_rank, 10)
tail(factor5_gene_rank, 10)
length(factor5_gene_rank)
# ============================================================
# STEP 3: GSEA - GO BIOLOGICAL PROCESS
# ============================================================

factor5_gsea_go <- gseGO(
  geneList = factor5_gene_rank,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

# Check GSEA results
head(as.data.frame(factor5_gsea_go), 10)
# ============================================================
# STEP 4: GSEA GO DOT PLOT
# ============================================================

factor5_gsea_dotplot <- dotplot(
  factor5_gsea_go,
  showCategory = 15,
  title = "MOFA Factor 5: GSEA GO Biological Processes"
) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 10)
  )

factor5_gsea_dotplot
# ============================================================
# STEP 4: GSEA GO DOT PLOT
# ============================================================

factor5_gsea_dotplot <- dotplot(
  factor5_gsea_go,
  showCategory = 15,
  title = "MOFA Factor 5: GSEA GO Biological Processes"
) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 10)
  )

factor5_gsea_dotplot

# Save GSEA GO dot plot
ggsave(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/Downstream_analysis/Factor5_GSEA_GO_Dotplot.png",
  factor5_gsea_dotplot,
  width = 10,
  height = 8,
  dpi = 300
)
# ============================================================
# STEP 5: GSEA ENRICHMENT PLOT
# ============================================================

factor5_gsea_enrichment <- gseaplot2(
  factor5_gsea_go,
  geneSetID = c("GO:0016064", "GO:0031424"),
  title = "MOFA Factor 5: GSEA Enrichment",
  subplots = 1:3
)

factor5_gsea_enrichment
library(enrichplot)
factor5_gsea_enrichment <- enrichplot::gseaplot2(
  factor5_gsea_go,
  geneSetID = c("GO:0016064", "GO:0031424"),
  title = "MOFA Factor 5: GSEA Enrichment",
  subplots = 1:3
)

factor5_gsea_enrichment
ggsave(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/Downstream_analysis/Factor5_GSEA_Enrichment.png",
  factor5_gsea_enrichment,
  width = 10,
  height = 8,
  dpi = 300
)
# ============================================================
# FINAL: SAVE BIOLOGICAL ENRICHMENT RESULTS
# ============================================================

# Save GO enrichment results
write.csv(
  as.data.frame(factor5_go),
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/Downstream_analysis/Factor5_GO_Enrichment_Results.csv",
  row.names = FALSE
)

# Save KEGG enrichment results
write.csv(
  as.data.frame(factor5_kegg),
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/Downstream_analysis/Factor5_KEGG_Enrichment_Results.csv",
  row.names = FALSE
)

# Save GSEA GO results
write.csv(
  as.data.frame(factor5_gsea_go),
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/Downstream_analysis/Factor5_GSEA_GO_Results.csv",
  row.names = FALSE
)

cat("All biological enrichment result tables saved successfully.\n")
