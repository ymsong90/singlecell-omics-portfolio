################################################################################
# Step 01: Debarcoding and Normalization
#
# Purpose: Process multiplexed CyTOF data with bead-based normalization
# Dataset: NMIBC bladder cancer CyTOF
#
# Input:
#   - ./rawdata/250904/*.FCS (raw CyTOF FCS files)
#   - Barcode key (defined in script)
#
# Output:
#   - ./data/sample_*.fcs (debarcoded individual samples)
#   - Quality control plots (beads, scatter, yields)
#
# Author: YMS
# Date: 2025
################################################################################

library(CATALYST)
library(cowplot)
library(flowCore)
library(ggplot2)
library(SingleCellExperiment)

# Create output directory
if (!dir.exists("./data")) {
    dir.create("./data", recursive = TRUE)
}
if (!dir.exists("./figure/01_debarcoding")) {
    dir.create("./figure/01_debarcoding", recursive = TRUE)
}

################################################################################
# 1. Load Raw CyTOF Data
################################################################################

# Specify raw FCS file paths
raw_data <- c("./rawdata/250904/25090A_bladder_01.FCS")

# Prepare data into SingleCellExperiment object
sce <- prepData(raw_data)

# Check sample information
table(sce$sample_id)
names(int_colData(sce))

################################################################################
# 2. Bead-Based Normalization
################################################################################

# NOTE: CyTOF data requires normalization using DVS beads
# This corrects for signal drift during long acquisition runs

res <- normCytof(
    sce,
    beads = "dvs",
    k = 50,
    assays = c("counts", "exprs"),
    overwrite = FALSE
)

# Calculate bead and removed event statistics
n <- ncol(sce)
ns <- c(ncol(res$beads), ncol(res$removed))

qc_stats <- data.frame(
    check.names = FALSE,
    "#" = c(ns[1], ns[2]),
    "%" = 100 * c(ns[1]/n, ns[2]/n),
    row.names = c("beads", "removed")
)

print(qc_stats)

# Extract normalized data (beads and doublets removed)
sce <- res$data
assayNames(sce)

################################################################################
# 3. Quality Control Plots
################################################################################

# Bead vs DNA scatter plot
bead_scatter <- res$scatter
ggsave(
    filename = "./figure/01_debarcoding/bead_dna_scatter.png",
    plot = bead_scatter,
    width = 8,
    height = 6,
    dpi = 300
)

# Smoothed bead intensities over time
bead_lines <- res$lines
ggsave(
    filename = "./figure/01_debarcoding/bead_intensity_timeline.png",
    plot = bead_lines,
    width = 10,
    height = 6,
    dpi = 300
)

################################################################################
# 4. Barcode Assignment
################################################################################

# NOTE: Sample key defines the barcode combinations
# Rows = samples, Columns = barcode channels
# 1 = barcode present, 0 = barcode absent

sample_key <- data.frame(
    '89'  = c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0),
    '194' = c(1,1,1,1,0,0,0,0,1,1,1,1,0,0,0),
    '195' = c(1,0,0,0,1,1,1,0,1,1,0,0,1,1,0),
    '196' = c(0,1,0,0,1,0,0,1,1,0,1,0,1,0,1),
    '198' = c(0,0,1,0,0,1,0,1,0,1,0,1,0,1,1),
    '113' = c(0,0,0,1,0,0,1,0,0,0,1,1,1,1,1),
    check.names = FALSE
)

colnames(sample_key) <- c(89, 194, 195, 196, 198, 113)

# Assign preliminary barcodes
sce <- assignPrelim(sce, sample_key)

# Check barcode channels
barcode_channels <- rownames(sce)[rowData(sce)$is_bc]
print(paste("Barcode channels:", paste(barcode_channels, collapse = ", ")))

# Check barcode assignment distribution
table(sce$bc_id)

################################################################################
# 5. Estimate and Apply Separation Cutoffs
################################################################################

# NOTE: Separation cutoffs distinguish positive from negative barcode signals
# Can use global cutoff or population-specific cutoffs

# Estimate separation cutoffs
sce <- estCutoffs(sce)
print(metadata(sce)$sep_cutoffs)

# Visualize yield estimates
yield_plot <- plotYields(sce, which = c(0, "1"))
ggsave(
    filename = "./figure/01_debarcoding/barcode_yields.png",
    plot = yield_plot,
    width = 8,
    height = 6,
    dpi = 300
)

# Apply population-specific cutoffs
sce_specific <- applyCutoffs(sce)

# Apply global cutoff for comparison
sce_global <- applyCutoffs(sce, sep_cutoffs = 0.35)

# Compare yields
yield_comparison <- c(
    specific = mean(sce_specific$bc_id != 0),
    global = mean(sce_global$bc_id != 0)
)
print(yield_comparison)

# NOTE: Proceed with population-specific filtering for better accuracy
sce <- sce_specific

################################################################################
# 6. Quality Assessment of Debarcoding
################################################################################

# Plot example events for unassigned and specific barcode populations
events_plot <- plotEvents(sce, which = c(0, "15"), n = 25)
ggsave(
    filename = "./figure/01_debarcoding/example_events.png",
    plot = events_plot,
    width = 10,
    height = 8,
    dpi = 300
)

# Mahalanobis distance plot for barcode population
mahal_plot <- plotMahal(sce, which = "1")
ggsave(
    filename = "./figure/01_debarcoding/mahalanobis_distance.png",
    plot = mahal_plot,
    width = 8,
    height = 6,
    dpi = 300
)

# Scatter plot of marker channels (QC check)
marker_scatter <- plotScatter(sce, chs = c("CD45", "CD45_e"))
ggsave(
    filename = "./figure/01_debarcoding/cd45_scatter.png",
    plot = marker_scatter,
    width = 8,
    height = 6,
    dpi = 300
)

################################################################################
# 7. Export Debarcoded Samples
################################################################################

# Convert SCE to flowSet, split by barcode
fs <- sce2fcs(sce, split_by = "bc_id")

# Get sample identifiers
ids <- fsApply(fs, identifier)
output_dir <- "./data"

# Export each sample as individual FCS file
for (id in ids) {
    ff <- fs[[id]]
    fn <- sprintf("sample_%s.fcs", id)
    fn <- file.path(output_dir, fn)
    write.FCS(ff, fn)
}

print(paste("Exported", length(ids), "debarcoded samples to", output_dir))

# Save debarcoded SCE object
save(sce, file = "./data/sce_debarcoded.RData")

