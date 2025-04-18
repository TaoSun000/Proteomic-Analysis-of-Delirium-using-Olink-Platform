# -------------------------------------------------------------
# Project: Proteomic Analysis of Delirium Post-Cardiac Surgery
# Platform: Olink Proteomics
# Date: March 20, 2024
# Author: Tao Sun
# -------------------------------------------------------------

# Clear workspace and set up environment
rm(list = ls())

# -------------------------------------------------------------
# Step 1. Install and Load Required Packages
# -------------------------------------------------------------
required_packages <- c(
  "OlinkAnalyze", "dplyr", "tidyr", "stringr", "ggplot2",
  "clusterProfiler", "org.Hs.eg.db", "pheatmap", "ggrepel", "survival", "pROC"
)

install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

install_if_missing(required_packages)

# -------------------------------------------------------------
# Step 2. Load and Preprocess Proteomics Data
# -------------------------------------------------------------
batch1 <- read_NPX("PATH_TO_BATCH1_CSV")
batch2 <- read_NPX("PATH_TO_BATCH2_CSV")

status_updated <- read.csv("PATH_TO_STATUS_UPDATED_CSV", stringsAsFactors = FALSE)

# Combine batches
proteomics_data <- bind_rows(batch1, batch2) %>%
  filter(!str_detect(Sample_Type, "CONTROL")) %>%
  mutate(across(where(is.character), str_trim))

# Prepare metadata
status_long <- status_updated %>%
  pivot_longer(cols = c(Case.PID, Control.PID), names_to = "Case_Control", values_to = "PID") %>%
  pivot_longer(cols = c(Case.Specimen.IDs, Control.specimen.IDs), names_to = "Specimen_Type", values_to = "SampleID") %>%
  mutate(Case_Control = ifelse(str_detect(Case_Control, "Case"), "Case", "Control"))

proteomics_data <- proteomics_data %>%
  left_join(status_long, by = "SampleID") %>%
  filter(!is.na(Case_Control))

# -------------------------------------------------------------
# Step 3. Paired t-Test (Benjamini-Hochberg Adjusted)
# -------------------------------------------------------------
ttest_results <- olink_ttest(proteomics_data, variable = "Case_Control", pair_id = "pair") %>%
  arrange(Adjusted_pval)

significant_proteins <- ttest_results %>% filter(Adjusted_pval < 0.05)
write.csv(ttest_results, "paired_ttest_all.csv", row.names = FALSE)
write.csv(significant_proteins, "paired_ttest_significant.csv", row.names = FALSE)

# -------------------------------------------------------------
# Step 4. Visualization: Volcano Plot
# -------------------------------------------------------------
top20_proteins <- significant_proteins %>% slice_head(n = 20)

volcano_plot <- olink_volcano_plot(p.val_tbl = ttest_results, x_lab = "log2 Fold Change",
                                   olinkid_list = top20_proteins$OlinkID)
ggsave("top20_volcanoplot.pdf", plot = volcano_plot, width = 10, height = 8)

# -------------------------------------------------------------
# Step 5. Visualization: Boxplots for Top 5 Proteins
# -------------------------------------------------------------
top5_proteins <- top20_proteins %>% slice_head(n = 5) %>% pull(OlinkID)

boxplot_data <- proteomics_data %>% filter(OlinkID %in% top5_proteins)

ggplot(boxplot_data, aes(x = Assay, y = NPX, fill = Case_Control)) +
  geom_boxplot() + theme_bw() +
  labs(title = "Top 5 Significant Proteins", x = "Protein Assay", y = "NPX")
ggsave("top5_protein_boxplot.pdf", width = 10, height = 8)

# -------------------------------------------------------------
# Step 6. PCA Analysis
# -------------------------------------------------------------
pca_plot <- olink_pca_plot(proteomics_data, color_g = "Case_Control")
ggsave("pca_analysis.pdf", plot = pca_plot, width = 10, height = 8)

# -------------------------------------------------------------
# Step 7. Gene Set Enrichment Analysis (GSEA)
# -------------------------------------------------------------
gsea_results <- olink_pathway_enrichment(data = proteomics_data, test_results = ttest_results)
write.csv(gsea_results, "GSEA_results.csv", row.names = FALSE)

heatmap_gsea <- olink_pathway_heatmap(enrich_results = gsea_results,
                                      test_results = ttest_results, number_of_terms = 20)
ggsave("GSEA_heatmap.pdf", plot = heatmap_gsea, width = 10, height = 12)

# -------------------------------------------------------------
# Step 8. Over-Representation Analysis (ORA)
# -------------------------------------------------------------
ora_results <- olink_pathway_enrichment(data = proteomics_data, test_results = ttest_results, method = "ORA")
write.csv(ora_results, "ORA_results.csv", row.names = FALSE)

ora_heatmap <- olink_pathway_heatmap(enrich_results = ora_results,
                                     test_results = ttest_results, method = "ORA", number_of_terms = 20)
ggsave("ORA_heatmap.pdf", plot = ora_heatmap, width = 10, height = 12)

# -------------------------------------------------------------
# Step 9. Conditional Logistic Regression for Top Protein
# -------------------------------------------------------------
top_protein_data <- proteomics_data %>% filter(OlinkID == "FKBP1B") %>%
  pivot_wider(names_from = "Assay", values_from = "NPX") %>%
  na.omit()

clogit_model <- clogit(as.numeric(Case_Control == "Case") ~ FKBP1B + strata(pair), data = top_protein_data)
summary(clogit_model)

# ROC Analysis
roc_analysis <- roc(top_protein_data$Case_Control, predict(clogit_model, type = "risk"))
pdf("ROC_curve_FKBP1B.pdf")
plot(roc_analysis, main = paste("ROC Curve for FKBP1B", round(roc_analysis$auc, 2)))
dev.off()

# -------------------------------------------------------------
# End of Analysis
# -------------------------------------------------------------
