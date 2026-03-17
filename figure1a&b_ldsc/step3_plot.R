# ==============================================================================
# Step 6: Compute p-values for each cell-type × GWAS trait combination
#
# Reads all Stratified LDSC output (.results) files from Step 5 and derives
# two p-value metrics from the Coefficient_z-score of the L2_0 annotation:
#
#   One-tailed p:  P = 1 - pnorm(z)
#     → Tests whether the cell type makes a *positive* contribution to the
#       GWAS trait (lenient; directional enrichment only).
#
#   Two-tailed p:  P = 2 * pnorm(-|z|)
#     → Tests enrichment in *either* direction (stricter; recommended for
#       unbiased reporting).
# ==============================================================================
library(stringr)
library(tidyr)
library(ggplot2)
library(tidyverse)   # Loads dplyr, readr, purrr, tibble, etc.

# FIX: Read .results files from the directory where Step 5 actually wrote them.
#      The original script incorrectly pointed to /ldsc/ldsc/ instead of
#      the z-result output directory.
results_dir <- "/data/h01005/hc/0.1ldsc-gwas/z-result/Blauwendraat_Parkinsons_sexdiff_2021_Sumstats"
files <- list.files(results_dir, pattern = "\\.results$", full.names = TRUE)

# Helper: load all .results files and extract cell-type & trait labels from
# the filename, then keep only the baseline annotation row (Category == "L2_0").
load_results <- function(files) {
  # FIX: replaced deprecated data_frame() with tibble()
  tibble(filename = files) %>%
    mutate(file_contents = map(filename, read_tsv, show_col_types = FALSE)) %>%
    mutate(
      # FIX: escape literal dots in gsub patterns so they are not treated as
      #      "match any character" wildcards in the regex engine.
      makenames = gsub("\\.bed_tissue_dir\\.results$",           "", basename(filename)),
      makenames = gsub("\\.bed_continuous_tissue_dir\\.results$","", makenames)
    ) %>%
    unnest(file_contents) %>%
    filter(Category == "L2_0") %>%   # Retain only the cell-type annotation row
    mutate(
      # Trait: everything before the first underscore in the filename stem
      Trait     = sub("_.*", "", makenames),
      # Cell_Type: everything after the first underscore (leading _ stripped)
      Cell_Type = gsub("^_", "", str_extract(makenames, "_.*"))
    ) %>%
    dplyr::select(Trait, Cell_Type, Enrichment, Enrichment_std_error,
                  Enrichment_p, `Coefficient_z-score`)
}

base_data <- load_results(files)

# --- One-tailed p-value (positive enrichment only) ---
d_one <- base_data %>%
  mutate(P = 1 - pnorm(`Coefficient_z-score`)) %>%    # P(Z > z): right-tail probability
  dplyr::select(-`Coefficient_z-score`) %>%
  arrange(Enrichment_p)

# --- Two-tailed p-value (enrichment in either direction) ---
d_two <- base_data %>%
  mutate(P = 2 * pnorm(-abs(`Coefficient_z-score`))) %>%  # P(|Z| > |z|): both tails
  dplyr::select(-`Coefficient_z-score`) %>%
  arrange(Enrichment_p)

write_tsv(d_one, file = "/data/h01005/hc/ldsc/ldsc/cell_types_GWAS_pvalues_one.txt")
write_tsv(d_two, file = "/data/h01005/hc/ldsc/ldsc/cell_types_GWAS_pvalues_two.txt")


# ==============================================================================
# Step 7: Visualise enrichment p-values — Figure 1A (scRNA) & 1B (scATAC)
#
# Raw p-values are –log10-transformed so that more significant results plot
# taller. Four Parkinson's disease GWAS datasets are compared side-by-side
# for each brain cell type.
#
# Dataset key:
#   ieu_b_7          → IEU Open GWAS PD summary statistics
#   LARGE_PD_...hg38 → LARGE-PD / TRACTOR EUR meta-analysis (hg38)
#   GP2_EUR_ONLY     → GP2 EUR-only GWAS
#   Global_PD_...    → Global Parkinson's Genetics Program 2023 trans-ancestry
# ==============================================================================

# Shared dataset column names and display labels (defined once to avoid repetition)
dataset_cols <- c(
  "ieu_b_7",
  "LARGE_PD_PAPER_TRACTOR_EUR_Meta_analysis_random_LARGE_PD_P1_P2_hg38",
  "GP2_EUR_ONLY",
  "Global_Parkinson_Genetics_Program_2023_GWAS_trans_ancestry"
)
dataset_labels <- c(
  "IEU b-7",
  "LARGE-PD TRACTOR EUR (hg38)",
  "GP2 EUR only",
  "Global PD Genetics Program 2023"
)
dataset_colors <- c("#1f78b4", "#33a02c", "pink", "#f57c6e")

# Helper: convert a wide data frame to plot-ready long format.
#   1. Apply –log10 transformation to all dataset columns at once.
#   2. Pivot to long format.
#   3. Deduplicate (keep max value per cell-type × dataset combination).
#   4. Set ordered factors so ggplot respects the intended axis order.
prepare_plot_data <- function(df) {
  df %>%
    # FIX: use mutate(across()) instead of transforming each column separately,
    #      which is both more concise and less error-prone when columns change.
    mutate(across(all_of(dataset_cols), ~ -log10(.x))) %>%
    pivot_longer(
      cols      = all_of(dataset_cols),
      names_to  = "Dataset",
      values_to = "Value"
    ) %>%
    # Deduplicate: if the same cell-type × dataset pair appears more than once,
    # retain the row with the highest –log10(p) value.
    group_by(celltype, Dataset) %>%
    summarise(Value = max(Value), .groups = "drop") %>%
    mutate(
      # Set factor levels so the x-axis and legend follow a consistent order
      celltype = factor(celltype, levels = sort(unique(celltype)), ordered = TRUE),
      Dataset  = factor(Dataset,  levels = dataset_cols,           ordered = TRUE)
    ) %>%
    # Sort rows to match factor level order (required for correct dodge alignment)
    arrange(celltype, Dataset)
}

# Helper: produce the grouped bar chart from a prepared long-format data frame.
make_barplot <- function(plot_data_clean, title) {
  ggplot(plot_data_clean, aes(x = celltype, y = Value, fill = Dataset)) +
    # Grouped bars; width = 0.7 leaves visible gaps between cell-type groups
    geom_col(position = "dodge", width = 0.7) +
    # Label each bar with its rounded –log10(p) value, positioned above the bar
    geom_text(
      aes(label = round(Value, 3)),
      position = position_dodge(width = 0.7),
      vjust = -0.5, size = 3
    ) +
    scale_fill_manual(values = dataset_colors, labels = dataset_labels) +
    labs(
      title = title,
      x     = "Cell Type",
      y     = expression(-log[10](p-value))
    ) +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
      legend.title = element_blank(),
      plot.title   = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    # Add 10 % headroom above the tallest bar so labels are not clipped
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
}


# --- Figure 1A: scRNA-seq cell types ---
scrna_raw <- data.frame(
  celltype = c("astrocyte", "fibroblast", "microglial_cell",
               "neuron", "oligodendrocyte", "oligodendrocyte_precursor_cell"),
  ieu_b_7 = c(0.46841380097643104, 0.4599130283932974, 0.5849195730947503,
              0.5041866058564755,  0.36373987360744664, 0.7110644408746453),
  LARGE_PD_PAPER_TRACTOR_EUR_Meta_analysis_random_LARGE_PD_P1_P2_hg38 =
    c(0.5454709150587502, 0.8365303856443889, 0.022398566647692908,
      0.45625740013938076, 0.004665446797138428, 0.2837346315062562),
  GP2_EUR_ONLY =
    c(0.6485323235785059, 0.9755694009841253, 0.03600845652966844,
      0.3378194085179831, 0.0350128676757504,  0.1752248710979334),
  Global_Parkinson_Genetics_Program_2023_GWAS_trans_ancestry =
    c(0.2019960321583696, 0.7289435573082299, 0.04621304519585989,
      0.7803115727872117, 0.13805993242476333, 0.1016093796315225)
)

make_barplot(prepare_plot_data(scrna_raw),
             "scRNA-seq: Cell-Type Enrichment Across PD GWAS Datasets")


# --- Figure 1B: scATAC-seq cell types ---
scatac_raw <- data.frame(
  celltype = c("astrocyte", "microglial_cell", "neuron",
               "oligodendrocyte", "oligodendrocyte_precursor_cell"),
  ieu_b_7 =
    c(0.2069535199204191, 0.35434086362761463, 0.2803916920577668,
      0.39677006108164237, 0.4822421736131449),
  LARGE_PD_PAPER_TRACTOR_EUR_Meta_analysis_random_LARGE_PD_P1_P2_hg38 =
    c(0.09134056287402958, 0.80436976630058,   0.6718499714043197,
      0.0190978322629749,  0.18286522849836762),
  GP2_EUR_ONLY =
    c(0.8329032069409021, 0.6563379317056216, 0.032279641914275126,
      0.6979193239835468, 0.6571429023425331),
  Global_Parkinson_Genetics_Program_2023_GWAS_trans_ancestry =
    c(0.15446340636008493, 0.02110552006675248, 0.01636324551245405,
      0.07445546398896186, 0.061631625319569894)
)

make_barplot(prepare_plot_data(scatac_raw),
             "scATAC-seq: Cell-Type Enrichment Across PD GWAS Datasets")