#===========================================================
# STEP 1: Select significant scATAC-seq peaks per cell type and generate
#         BED files for LDSC (LD Score Regression) annotation input
#
# Dataset : GSE147672 (human brain scATAC-seq, 24 Seurat clusters)
# Input   : per-cluster narrowPeak files in data_dir/
# Output  : one BED file per cell type containing the top 20% most
#           significant peaks, ranked by qValue (-log10 scale)
#===========================================================

# --- Load required libraries ---
library(data.table)  
library(purrr)       
library(dplyr)       
library(readr)       

#===========================================================
# 1. Define cluster-to-cell-type mapping
#    24 Seurat clusters manually annotated into 5 broad brain cell types
#    based on marker gene expression in the original GSE147672 study
#===========================================================

cell_type_map <- list(
  "Neurons"          = c(1, 2, 3, 4, 5, 6, 7, 11, 12, 18),
  "OPCs"             = 8:10,
  "Astrocytes"       = 13:17,
  "Oligodendrocytes" = 19:23,
  "Microglia"        = 24
)

# Build reverse lookup table: Cluster ID (integer) -> Cell Type label (string)
# rep() pairs each cell type name with all cluster IDs it contains
cluster_mapping <- data.table(
  Cluster_id = unlist(cell_type_map),
  Cell_Type  = rep(names(cell_type_map), times = sapply(cell_type_map, length))
)

#===========================================================
# 2. Define file paths
#    output_file: path for the optional full combined peak table (see Step 3)
#===========================================================

data_dir    <- "/data/h01005/hc/GSE147672/peak/"
output_file <- "/data/h01005/hc/GSE147672/combined_peaks_with_celltypes.tsv.gz"

#===========================================================
# 3. Read and combine all 24 cluster-level narrowPeak files
#
#    narrowPeak column definitions (ENCODE standard, 10 columns, no header):
#      V1  chrom       — chromosome name (e.g., chr1)
#      V2  start       — peak start coordinate (0-based)
#      V3  end         — peak end coordinate
#      V4  name        — peak identifier
#      V5  score       — integer score (0-1000)
#      V6  strand      — strand orientation (+/-/.)
#      V7  signalValue — fold-change or signal intensity
#      V8  pValue      — -log10(p-value); -1 if not computed
#      V9  qValue      — -log10(q-value, BH-FDR adjusted); -1 if not computed
#      V10 peak        — offset (bp) from start to peak summit
#===========================================================

combined_data <- rbindlist(
  lapply(1:24, function(cluster_num) {
    
    # Construct the file path for this cluster
    file_path <- file.path(
      data_dir,
      paste0("Cluster", cluster_num, ".overlap.optimal.narrowPeak.gz")
    )
    
    # Read narrowPeak file; assign temporary names V1-V10 (file has no header line)
    dt <- fread(
      file_path,
      header    = FALSE,
      sep       = "\t",
      col.names = paste0("V", 1:10)   # Temporary names; replaced in Step 4
    )
    
    # Look up the cell type label for this cluster via the reverse mapping table
    dt[, Cell_Type := cluster_mapping[Cluster_id == cluster_num, Cell_Type]]
    
    # Record the originating cluster number for downstream traceability
    dt[, Cluster := cluster_num]
    
    return(dt)
  }),
  fill = TRUE   # Pad with NA if any file unexpectedly has fewer columns
)


#===========================================================
# 4. Assign standard narrowPeak column names
#    Replaces opaque V1-V10 names with field names from the ENCODE spec,
#    plus the two appended columns (cell_type, cluster)
#    After this step, qValue (formerly V9) is confirmed as -log10(q-value)
#===========================================================

colnames(combined_data) <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "signalValue", "pValue", "qValue", "summit",
  "cell_type", "cluster"
)

#===========================================================
# 5. Select the top 20% most significant peaks per cell type
#
#    qValue = -log10(BH-adjusted q-value):
#      HIGHER value -> smaller q-value -> MORE significant
#      -> slice_max correctly retains the most significant peaks
#
#    prop = 0.2 retains floor(n x 0.2) rows per group.
#    with_ties = FALSE: drop boundary ties to return the exact top 20%
#===========================================================

significant_peaks <- combined_data %>%
  group_by(cell_type) %>%
  slice_max(
    order_by  = qValue,   # -log10(q): higher = more significant
    prop      = 0.2,      # Retain top 20% most significant peaks per cell type
    with_ties = FALSE     # Drop boundary ties for an exact 20% cutoff
  ) %>%
  ungroup() %>%
  as.data.frame()

#===========================================================
# 6. Prepare LDSC input table
#    Retain only genomic coordinates, peak name, and cell type label
#    dplyr:: prefix is required: data.table also exports select(),
#    causing a namespace conflict without the explicit prefix
#===========================================================

ldsc_input <- significant_peaks %>%
  dplyr::select(chrom, start, end, name, cell_type)

#===========================================================
# 7. Write one BED file per cell type
#
#    BED format (tab-delimited, no header): chrom | start | end | name
#    Directory "0.2" reflects the 20% significance threshold used above
#    make.names() converts cell-type labels to valid filenames:
#      spaces and special characters are replaced with "."
#      e.g., "Ex Neuron" -> "Ex.Neuron"
#===========================================================

dir.create(
  "/data/h01005/hc/GSE147672/ldsc/bed/0.2",
  recursive    = TRUE,
  showWarnings = FALSE
)

ldsc_input %>%
  split(.$cell_type) %>%
  walk(~ write_tsv(
    .x %>% dplyr::select(chrom, start, end, name),
    file      = paste0(
      "/data/h01005/hc/GSE147672/ldsc/bed/0.2/",
      make.names(unique(.x$cell_type)),   # Cell-type label -> valid filename
      ".bed"
    ),
    col_names = FALSE                     # BED format requires no header line
  ))