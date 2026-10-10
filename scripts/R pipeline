# Name: Kainoa Martinez
# Title: Linking MIDAS Data with Biological Databases: An R-Based Pipeline for Protein-Metabolite Discovery

# ================================================================================

# Setting up the environment - Install these packages & run libraries
install.packages("patchwork")

library(readxl)
library(stringr)

midas_stats = read.delim("midas_out_stats.tsv", header = TRUE, stringsAsFactors = FALSE)
midas_metab_pos = read_xlsx("20241121_protein_MIDAS_cHILIC_positive_metadata.xlsx")
midas_metab_neg = read_xlsx("20241121_protein_MIDAS_pHILIC_negative_metadata.xlsx")

# Are the names the same?
all(names(midas_metab_pos) == names(midas_metab_neg))

# Combine 
midas_combined = rbind(midas_metab_neg, midas_metab_pos)

hmdb = read.delim("parsed_hmdb_metabolites.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
hmdb_sub = subset(hmdb, hmdb$accession %in% midas_combined$HMDB_ID)


# midas_combined_hmdb = merge(midas_combined, hmdb_sub, by.x = "HMDB_ID", by.y = "accession", all.x = TRUE)
# treat midas_stats as the reference and then use what is in it to check the other tables

# unique_metabolites = unique(midas_stats$Metabolite_Name)
for (i in 1:nrow(midas_stats)) {
# for (i in 1:100) {
# start with the top of ms stats
current_metabolite = midas_stats$Metabolite_Name[i]
# find the corresponding HMDB ID
midas_combined_current_metabolite = subset(midas_metab_neg, Metabolite_Name == current_metabolite)
# midas_combined_current_metabolite_hmdb = midas_combined_current_metabolite$HMDB_ID
current_metabolite_hmdb_annotation = unique(midas_combined_current_metabolite$HMDB_ID)
if (length(current_metabolite_hmdb_annotation) > 1) {
    print("length(current_metabolite_hmdb_annotation) > 1")
    print(i)
}

hmdb_subset_from_metabolite = subset(hmdb_sub, accession == current_metabolite_hmdb_annotation)

if (length(hmdb_subset_from_metabolite$protein_associations..protein..uniprot_id) > 0 ) {
midas_stats$hmdb_interaction[i] = hmdb_subset_from_metabolite$protein_associations..protein..uniprot_id
} else {
    midas_stats$hmdb_interaction[i] = NA
    }
}

midas_stats_single_uniprot = subset(midas_stats, !str_detect(midas_stats$uniprot_names, " "))
midas_stats_single_uniprot = subset(midas_stats_single_uniprot, !str_detect(midas_stats_single_uniprot$uniprot_names, ";"))
midas_stats_single_uniprot = subset(midas_stats_single_uniprot, !is.na(midas_stats_single_uniprot$hmdb_interaction) & midas_stats_single_uniprot$hmdb_interaction != "")

midas_stats_single_uniprot$known_interaction = NA
for (i in 1:nrow(midas_stats_single_uniprot)) {
    current_uniprot = midas_stats_single_uniprot$uniprot_ids[i]
    current_known = unlist(str_split(midas_stats_single_uniprot$hmdb_interaction[i], "\\|\\|"))
    midas_stats_single_uniprot$known_interaction[i] = ifelse(current_uniprot %in% current_known, "known", "not known")
# } 

}
write.table(midas_stats_single_uniprot, "midas_stats_single_uniprot_2025-07-26.txt", sep = "\t", row.names = F)


midas_stats_single_uniprot$metabolite_uniprot = paste0(midas_stats_single_uniprot$Metabolite_Name, "_", midas_stats_single_uniprot$uniprot_ids)

midas_stats_single_uniprot_sig = subset(midas_stats_single_uniprot, midas_stats_single_uniprot$log2_fold_change > 1 & midas_stats_single_uniprot$q_value < 0.01)

table(midas_stats_single_uniprot_sig$known_interaction) # plot this

sort(table(midas_stats_single_uniprot_sig$Metabolite_Name)) # plot this
# ================================================================================

# ANALYSIS AND VISUALIZATION

# Barplot  ---- 
# Count known vs novel interactions (barplot)
pmi_status_counts <- table(midas_stats_single_uniprot_sig$known_interaction)
pmi_status_df <- as.data.frame(pmi_status_counts)
colnames(pmi_status_df) <- c("PMI_status", "Count")

# Plot bar chart
library(ggplot2)

ggplot(pmi_status_df, aes(x = PMI_status, y = Count, fill = PMI_status)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +  # Add count labels above bars
  scale_fill_manual(values = c("known" = "forestgreen", "not known" = "#e31a1c")) +
  scale_y_continuous(limits = c(0, 90)) +  # Set y-axis from 0 to 100
  labs(
    title = "Known vs Novel Protein–Metabolite Interactions",
    x = "Interaction Type",
    y = "Number of PMIs",
    fill = "PMI Status"
  ) +
  theme_minimal(base_size = 14)

# Barplot ----
#Combined Barplot showing top significant proteins and metabolites
library(dplyr)

# Count top 10 metabolites
top_metabolites <- midas_stats_single_uniprot_sig %>%
  count(Metabolite_Name, sort = TRUE) %>%
  slice_max(n, n = 10) %>%
  mutate(Category = "Metabolite")

# Count top 10 proteins
top_proteins <- midas_stats_single_uniprot_sig %>%
  count(uniprot_names, sort = TRUE) %>%
  slice_max(n, n = 10) %>%
  rename(Metabolite_Name = uniprot_names) %>%  # Rename for shared column name
  mutate(Category = "Protein")

# Combine both into one dataframe
library(ggplot2)
library(grid)

top_combined <- bind_rows(top_metabolites, top_proteins)

ggplot(top_combined, aes(x = reorder(Metabolite_Name, n), y = n, fill = Category)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = n), hjust = -0.1, size = 3.5) +
  facet_wrap(~ Category, scales = "free_y") +
  labs(
    title = "Top 10 Significant Proteins and Metabolites",
    x = "",
    y = "Number of Interactions"
  ) +
  scale_fill_manual(values = c("Metabolite" = "firebrick", "Protein" = "steelblue")) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10, margin = margin(r = 5)),
    legend.position = "none",
    panel.spacing = unit(2, "lines"),
    strip.placement = "outside",
    strip.text = element_text(margin = margin(b = 8))
  ) +
  coord_flip()


# PCA proteins  ----
library(tidyverse)
library(ggrepel)

# Step 1: Prepare data with composite score
adjusted_data = midas_stats %>%
    filter(
        !is.na(uniprot_names),
        !is.na(unique_metabolite_name),
        is.finite(log2_fold_change),
        is.finite(p_value)
    ) %>%
    mutate(
        p_value = if_else(p_value == 0, min(p_value[p_value > 0], na.rm = TRUE), p_value),
        composite_score = log2_fold_change * -log10(p_value)
    ) %>%
    select(uniprot_names, unique_metabolite_name, composite_score)

# Step 2: Deduplicate protein–metabolite pairs
dedup_data = adjusted_data %>%
    group_by(uniprot_names, unique_metabolite_name) %>%
    slice_max(order_by = abs(composite_score), n = 1, with_ties = FALSE) %>%
    ungroup()

# Step 3: Convert to protein × metabolite matrix
score_matrix = dedup_data %>%
    pivot_wider(names_from = unique_metabolite_name, values_from = composite_score) %>%
    column_to_rownames("uniprot_names")

# Step 4: Impute missing values (fill with 0)
score_matrix[is.na(score_matrix)] = 0

# Step 5: Perform PCA
pca_result = prcomp(score_matrix, scale. = TRUE)
explained_var = pca_result$sdev^2 / sum(pca_result$sdev^2)
pc1_var = round(explained_var[1] * 100, 1)
pc2_var = round(explained_var[2] * 100, 1)

# Step 6: Prepare PCA dataframe
pca_df = as_tibble(pca_result$x, rownames = "uniprot_names") %>%
    mutate(dist_from_origin = sqrt(PC1^2 + PC2^2))

# Step 7: Label top 10% farthest from origin
labeled_df = pca_df %>%
    filter(dist_from_origin >= quantile(dist_from_origin, 0.9))

# Step 8: Plot
ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(color = "steelblue") +
    geom_text_repel(
        data = labeled_df,
        aes(label = uniprot_names),
        size = 3,
        max.overlaps = Inf,
        segment.color = "gray50",
        segment.size = 0.5,
        box.padding = 0.5,
        point.padding = 0.3,
        force = 1,
        force_pull = 1
    ) +
    theme_minimal() +
    labs(
        title = "PCA of Protein–Metabolite Composite Scores",
        x = paste0("PC1 (", pc1_var, "%)"),
        y = paste0("PC2 (", pc2_var, "%)")
    )

# PCA metabolites ----
library(tidyverse)
library(ggrepel)

# Step 1: Prepare data with composite score
met_data = midas_stats %>%
    filter(
        !is.na(uniprot_names),
        !is.na(Metabolite_Name),
        is.finite(log2_fold_change),
        is.finite(p_value)
    ) %>%
    mutate(
        p_value = if_else(p_value == 0, min(p_value[p_value > 0], na.rm = TRUE), p_value),
        composite_score = log2_fold_change * -log10(p_value)
    ) %>%
    select(Metabolite_Name, uniprot_names, composite_score)

# Step 2: Deduplicate metabolite–protein pairs
met_dedup = met_data %>%
    group_by(Metabolite_Name, uniprot_names) %>%
    slice_max(order_by = abs(composite_score), n = 1, with_ties = FALSE) %>%
    ungroup()

# Step 3: Pivot to metabolite × protein matrix
met_matrix = met_dedup %>%
    pivot_wider(names_from = uniprot_names, values_from = composite_score) %>%
    column_to_rownames("Metabolite_Name")

# Step 4: Impute missing values
met_matrix[is.na(met_matrix)] = 0

# Step 5: PCA
met_pca = prcomp(met_matrix, scale. = TRUE)
met_var = met_pca$sdev^2 / sum(met_pca$sdev^2)
pc1_met = round(met_var[1] * 100, 1)
pc2_met = round(met_var[2] * 100, 1)

# Step 6: Prepare PCA dataframe
met_pca_df = as_tibble(met_pca$x, rownames = "Metabolite_Name")

# Step 7: Label metabolites with PC1 > 3 or < -3
met_labeled = met_pca_df %>%
    filter(PC1 > 4 | PC1 < -4)

# Step 8: Plot
ggplot(met_pca_df, aes(x = PC1, y = PC2)) +
    geom_point(color = "firebrick", alpha = 0.6) +
    geom_text_repel(
        data = met_labeled,
        aes(label = Metabolite_Name),
        size = 3,
        max.overlaps = 50,
        segment.color = "gray50",
        segment.size = 0.3,
        box.padding = 0.3,
        point.padding = 0.2,
        force = 2,
        force_pull = 1
    ) +
    theme_minimal() +
    labs(
        title = "PCA of Metabolite–Protein Composite Scores",
        x = paste0("PC1 (", pc1_met, "%)"),
        y = paste0("PC2 (", pc2_met, "%)")
    )


# Volcano Plot  ----

# Step 1: Get the protein (UniProt ID) with the most significant metabolite interactions
top_protein <- midas_stats_single_uniprot_sig %>%
  count(uniprot_ids, sort = TRUE) %>%
  slice(1) %>%
  pull(uniprot_ids)

# Step 1.5: Lookup gene symbol for top protein
top_gene_symbol <- midas_stats_single_uniprot %>%
  filter(uniprot_ids == top_protein) %>%
  pull(uniprot_names) %>%
  unique() %>%
  .[1]  # In case of duplicates, just take the first

# Step 2: Use the full dataset (not just significant) for volcano data
protein_data <- midas_stats_single_uniprot %>%
  filter(uniprot_ids == top_protein) %>%
  mutate(
    negLogQ = -log10(q_value),
    Significant = ifelse(q_value < 0.01 & abs(log2_fold_change) > 1, "Significant", "Not Significant")
  )

# Step 3: Plot
ggplot(protein_data, aes(x = log2_fold_change, y = negLogQ)) +
  geom_point(aes(color = Significant), alpha = 0.8, size = 2.5) +
  geom_text_repel(
    data = protein_data %>% filter(Significant == "Significant"),
    aes(label = Metabolite_Name),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Significant" = "firebrick", "Not Significant" = "grey60")) +
  labs(
    title = paste0("Volcano Plot: Metabolite Interactions with ", top_gene_symbol),
    x = "Log2 Fold Change",
    y = "-Log10(q-value)",
    color = "Significance"
  ) +
  theme_minimal(base_size = 14)
