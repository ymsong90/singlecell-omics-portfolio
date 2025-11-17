################################################################################
# IMC Analysis Pipeline
# Step 04: Quality Control and Filtering
# 
# Description:
#   - DBSCAN-based outlier detection on UMAP/t-SNE coordinates
#   - k-NN distance plot for parameter selection
#   - Grid search for optimal DBSCAN parameters
#   - Remove noise clusters and spatial outliers
#   - Optional: downsample cells while preserving proportions
#   - Update metadata from sample_metadata.csv
#   - Remove "Junk" cell annotations
#
# Input:
#   - ./data/spe_annotated.rds (from Step 03)
#   - ./data/steinbock/sample_metadata.csv
#
# Output:
#   - ./data/spe_filtered.rds
#   - DBSCAN QC plots
#
# Author: YMS
# Date: 2024-11-17
################################################################################

# Load required packages
library(SpatialExperiment)
library(SingleCellExperiment)
library(dbscan)
library(ggplot2)
library(cowplot)
library(viridis)
library(dittoSeq)
library(dplyr)
library(scales)

# Set parameters
QC_PARAMS <- list(
    dbscan_eps            = 0.16,         # DBSCAN epsilon (adjust based on k-NN plot)
    dbscan_minpts         = 35,           # DBSCAN minimum points
    knn_k                 = 20,           # k for k-NN distance plot
    extra_drop_clusters   = c(7, 8, 9, 10), # Additional clusters to remove
    downsample_target     = 10000,        # Target cell number (0 = no downsampling)
    seed                  = 900223
)

# Create output directory
if (!dir.exists("./figure/04_qc")) {
    dir.create("./figure/04_qc", recursive = TRUE)
}

################################################################################
# 1. Load Annotated Data
################################################################################

spe <- readRDS("./data/spe_annotated.rds")

################################################################################
# 2. Select Reduction for DBSCAN
################################################################################

# NOTE: Use Harmony-corrected UMAP if available, otherwise use standard UMAP
umap_name <- if ("UMAP_harmonyCorrected" %in% names(reducedDims(spe))) {
    "UMAP_harmonyCorrected"
} else if ("UMAP" %in% names(reducedDims(spe))) {
    "UMAP"
} else {
    stop("No UMAP coordinates found in reducedDims")
}

umap <- as.matrix(reducedDims(spe)[[umap_name]])
colnames(umap) <- c("UMAP1", "UMAP2")

################################################################################
# 3. k-NN Distance Plot for eps Selection
################################################################################

# IMPORTANT: k-NN distance plot helps determine optimal eps parameter
# Look for "elbow" point where distance rapidly increases
kn <- QC_PARAMS$knn_k
d <- kNNdist(umap, k = kn)
dS <- sort(d)

# Calculate suggested eps values
eps_p95 <- as.numeric(quantile(dS, 0.95))  # 95th percentile
dd <- diff(dS)
idx <- which.max(diff(dd))
eps_elbow <- dS[idx]  # Elbow point

# Plot k-NN distances
png(
    filename = "./figure/04_qc/01_kNN_distance_plot.png",
    width = 10,
    height = 7,
    units = "in",
    res = 300
)
par(mar = c(4.2, 4.2, 2, 1), las = 1, mgp = c(2.2, 0.7, 0))
plot(dS, type = "l",
     xlab = "Points (sorted by k-NN distance)",
     ylab = sprintf("%d-NN distance", kn),
     main = "k-NN Distance Plot for eps Selection")
abline(h = eps_p95, col = "red", lty = 2)
text(x = length(dS) * 0.8, y = eps_p95,
     labels = sprintf("p95 ≈ %.3f", eps_p95),
     pos = 3, col = "red")
abline(h = eps_elbow, col = "blue", lty = 3)
text(x = length(dS) * 0.5, y = eps_elbow,
     labels = sprintf("elbow ≈ %.3f", eps_elbow),
     pos = 3, col = "blue")
dev.off()

################################################################################
# 4. Grid Search for Optimal DBSCAN Parameters (Optional)
################################################################################

# NOTE: This step helps find optimal eps/minPts combination
# Adjust grid ranges based on your data
eps_grid <- seq(0.12, 0.20, by = 0.02)
minpts_grid <- c(25, 35, 50)

res <- expand.grid(eps = eps_grid, minPts = minpts_grid)
res$noise_frac <- NA_real_
res$n_clusters <- NA_integer_

for (i in seq_len(nrow(res))) {
    fit <- dbscan(umap, eps = res$eps[i], minPts = res$minPts[i])
    tab <- table(fit$cluster)
    res$noise_frac[i] <- ifelse(
        "0" %in% names(tab),
        as.numeric(tab["0"]) / nrow(umap),
        0
    )
    res$n_clusters[i] <- length(setdiff(names(tab), "0"))
}

# Save grid search results
write.csv(
    res[order(res$noise_frac, decreasing = TRUE), ],
    file = "./figure/04_qc/dbscan_grid_search.csv",
    row.names = FALSE
)

################################################################################
# 5. Run DBSCAN with Selected Parameters
################################################################################

set.seed(QC_PARAMS$seed)

db <- dbscan(
    umap,
    eps = QC_PARAMS$dbscan_eps,
    minPts = QC_PARAMS$dbscan_minpts
)

# Store DBSCAN labels (0 = noise)
spe$umap_dbscan <- factor(db$cluster)

# Summary
cat(sprintf("\nDBSCAN Results (eps=%.2f, minPts=%d):\n",
            QC_PARAMS$dbscan_eps, QC_PARAMS$dbscan_minpts))
print(table(spe$umap_dbscan))

################################################################################
# 6. Visualize DBSCAN Results
################################################################################

df_all <- data.frame(
    UMAP1 = umap[, 1],
    UMAP2 = umap[, 2],
    dbscan_cluster = spe$umap_dbscan,
    cell_id = colnames(spe)
)

# DBSCAN clusters
p_dbscan <- ggplot(df_all, aes(UMAP1, UMAP2, color = dbscan_cluster)) +
    geom_point(size = 0.4, alpha = 0.7) +
    scale_color_viridis_d() +
    ggtitle(sprintf("DBSCAN (eps=%.2f, minPts=%d)",
                    QC_PARAMS$dbscan_eps, QC_PARAMS$dbscan_minpts)) +
    theme_classic()

# Faceted view of each cluster
p_facet <- ggplot(df_all, aes(UMAP1, UMAP2)) +
    geom_point(data = df_all, color = "grey90", size = 0.1) +
    geom_point(aes(color = dbscan_cluster), size = 0.5) +
    facet_wrap(~dbscan_cluster, ncol = 5) +
    scale_color_viridis_d() +
    ggtitle("Individual DBSCAN Clusters") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none")

ggsave(
    filename = "./figure/04_qc/02_DBSCAN_clusters.png",
    plot = plot_grid(p_dbscan, p_facet, ncol = 1, rel_heights = c(1, 2)),
    width = 14,
    height = 16,
    dpi = 300
)

################################################################################
# 7. Filter Outliers
################################################################################

# Remove noise (cluster 0) and specified outlier clusters
cl_drop <- c("0", as.character(QC_PARAMS$extra_drop_clusters))

# Convert factor to character for comparison
is_outlier <- as.character(spe$umap_dbscan) %in% cl_drop
keep <- !is_outlier

cat(sprintf("\nFiltering outliers:\n"))
cat(sprintf("  Kept: %d cells (%.1f%%)\n",
            sum(keep), 100 * sum(keep) / ncol(spe)))
cat(sprintf("  Removed: %d cells (%.1f%%) [clusters: %s]\n",
            sum(!keep), 100 * sum(!keep) / ncol(spe),
            paste(cl_drop, collapse = ", ")))

# Visualization: before/after filtering
df_all$keep <- keep
df_all$display_cluster <- ifelse(keep, 
                                 as.character(df_all$dbscan_cluster),
                                 "DROPPED")

p_before_after <- ggplot(df_all, aes(UMAP1, UMAP2)) +
    geom_point(data = df_all[!keep, ], color = "red", size = 0.5, alpha = 0.8) +
    geom_point(data = df_all[keep, ], 
               aes(color = dbscan_cluster), size = 0.3, alpha = 0.6) +
    scale_color_viridis_d(name = "Kept clusters") +
    ggtitle(sprintf("Red = Dropped clusters {%s}", 
                    paste(cl_drop, collapse = ", "))) +
    theme_classic(base_size = 12)

p_kept_only <- ggplot(df_all[keep, ], 
                      aes(UMAP1, UMAP2, color = dbscan_cluster)) +
    geom_point(size = 0.4, alpha = 0.7) +
    scale_color_viridis_d() +
    ggtitle(sprintf("Kept cells only (n=%d)", sum(keep))) +
    theme_classic(base_size = 12)

ggsave(
    filename = "./figure/04_qc/03_filtering_comparison.png",
    plot = plot_grid(p_before_after, p_kept_only, ncol = 2),
    width = 16,
    height = 7,
    dpi = 300
)

# Apply filter
spe <- spe[, keep, drop = FALSE]

# Clean up DBSCAN labels
if ("umap_dbscan" %in% colnames(colData(spe))) {
    spe$umap_dbscan <- droplevels(spe$umap_dbscan)
}

################################################################################
# 8. Remove Junk Cells
################################################################################

# NOTE: Remove cells annotated as "Junk" during manual annotation
if ("cluster_celltype" %in% colnames(colData(spe))) {
    keep_notjunk <- spe$cluster_celltype != "Junk"
    
    cat(sprintf("\nRemoving Junk cells:\n"))
    cat(sprintf("  Junk cells: %d (%.1f%%)\n",
                sum(!keep_notjunk), 
                100 * sum(!keep_notjunk) / ncol(spe)))
    
    spe <- spe[, keep_notjunk, drop = FALSE]
    spe$cluster_celltype <- droplevels(spe$cluster_celltype)
}

################################################################################
# 9. Update Metadata (group, block)
################################################################################

# IMPORTANT: Re-read metadata and update group/block assignments
# This is useful if metadata was changed after initial preprocessing
meta <- read.csv("./data/steinbock/sample_metadata.csv", 
                 stringsAsFactors = FALSE)

stopifnot(all(c("sample_id", "group", "block") %in% colnames(meta)))

# Map updated metadata
idx <- match(spe$sample_id, meta$sample_id)

if (anyNA(idx)) {
    missing_ids <- unique(as.character(spe$sample_id[is.na(idx)]))
    warning(sprintf("Metadata에 없는 sample_id: %s",
                    paste(missing_ids, collapse = ", ")))
}

# Update group and block
grp_chr <- meta$group[idx]
blk_chr <- meta$block[idx]

spe$group <- factor(grp_chr, levels = unique(meta$group))
spe$block <- factor(blk_chr, levels = unique(as.character(meta$block)))

# Update color palettes
cv <- metadata(spe)$color_vectors

# Group colors
base_group_cols <- c(Wt = "#4e79a7", Ko = "#e15759", Panc = "#59a14f")
levs_group <- levels(spe$group)
missing_group <- setdiff(levs_group, names(base_group_cols))

if (length(missing_group) > 0) {
    base_group_cols <- c(
        base_group_cols,
        setNames(hue_pal()(length(missing_group)), missing_group)
    )
}
cv$group <- base_group_cols[levs_group]
names(cv$group) <- levs_group

# Block colors
levs_block <- levels(spe$block)
cv$block <- setNames(hue_pal()(length(levs_block)), levs_block)

metadata(spe)$color_vectors <- cv

cat("\nUpdated metadata:\n")
print(table(spe$group, useNA = "ifany"))
print(table(spe$block, useNA = "ifany"))

################################################################################
# 10. Optional: Downsample Cells
################################################################################

if (QC_PARAMS$downsample_target > 0 && 
    ncol(spe) > QC_PARAMS$downsample_target) {
    
    cat(sprintf("\nDownsampling to %d cells...\n", 
                QC_PARAMS$downsample_target))
    
    # Extract celltype safely
    ct_col <- colData(spe)$cluster_celltype
    celltype_vec <- if (inherits(ct_col, "Rle")) {
        as.character(as.vector(ct_col))
    } else if (is.list(ct_col) || inherits(ct_col, "List")) {
        as.character(unlist(ct_col, use.names = FALSE))
    } else {
        as.character(ct_col)
    }
    
    celltype_vec[is.na(celltype_vec)] <- "_NA_"
    
    # Calculate proportions
    tab_all <- table(celltype_vec)
    celltype_props <- data.frame(
        celltype = names(tab_all),
        n = as.integer(tab_all),
        stringsAsFactors = FALSE
    )
    
    celltype_props$prop <- celltype_props$n / sum(celltype_props$n)
    celltype_props$target_n <- as.integer(
        round(celltype_props$prop * QC_PARAMS$downsample_target)
    )
    
    # Adjust for rounding errors
    diff <- QC_PARAMS$downsample_target - sum(celltype_props$target_n)
    if (diff != 0 && nrow(celltype_props) > 0) {
        ord <- order(-celltype_props$prop)
        i <- 1
        while (diff != 0) {
            idx_adj <- ord[((i - 1) %% nrow(celltype_props)) + 1]
            celltype_props$target_n[idx_adj] <- 
                celltype_props$target_n[idx_adj] + sign(diff)
            diff <- QC_PARAMS$downsample_target - sum(celltype_props$target_n)
            i <- i + 1
        }
    }
    
    # Sample cells per celltype
    meta_cells <- data.frame(
        cell_id = seq_len(ncol(spe)),
        celltype = celltype_vec,
        stringsAsFactors = FALSE
    )
    
    idx_by_ct <- split(meta_cells$cell_id, meta_cells$celltype)
    
    set.seed(QC_PARAMS$seed)
    sample_ids <- unlist(
        mapply(function(ids, k) {
            if (length(ids) == 0 || k <= 0) {
                integer(0)
            } else {
                sample(ids, min(k, length(ids)))
            }
        }, idx_by_ct[celltype_props$celltype], 
        celltype_props$target_n, 
        SIMPLIFY = FALSE),
        use.names = FALSE
    )
    
    cat(sprintf("  Downsampled to %d cells\n", length(sample_ids)))
    
    spe <- spe[, sample_ids]
}

################################################################################
# 11. Final UMAP Visualization
################################################################################

p_final <- dittoDimPlot(
    spe,
    var = "cluster_celltype",
    reduction.use = umap_name,
    size = 0.7,
    do.label = FALSE
) +
    scale_color_manual(values = metadata(spe)$color_vectors$cluster_celltype) +
    theme_classic(base_size = 14) +
    theme_bw() +
    theme(
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        plot.title = element_blank(),
        panel.grid = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))

ggsave(
    filename = "./figure/04_qc/04_final_UMAP_after_QC.tiff",
    plot = p_final,
    width = 10,
    height = 8,
    dpi = 400,
    compression = "lzw"
)

################################################################################
# 12. Save Filtered Object
################################################################################

saveRDS(spe, "./data/spe_filtered.rds")

# Clean up
rm(umap, db, df_all, res)
gc()

# NOTE: Next step is 05_basic_visualization.R
