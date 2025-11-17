################################################################################
# Step 03: Clustering and Cell Type Annotation
#
# Purpose: Identify marker genes and annotate cell types
# Dataset: Mouse PORCN KO vs WT
#
# Input:
#   - ./data/porcn.combined.harmony.RData
#
# Output:
#   - ./data/porcn.combined.harmony_annotated.RData
#   - ./results/03_clustering/*.csv (marker genes)
#   - ./results/03_clustering/*.png (visualizations)
#
# Author: YMS
# Date: 2025
################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)

# Create output directory
if (!dir.exists("./results/03_clustering")) {
    dir.create("./results/03_clustering", recursive = TRUE)
}

################################################################################
# 1. Load Integrated Data
################################################################################

load("./data/porcn.combined.harmony.RData")

################################################################################
# 2. Find Cluster Markers
################################################################################

# NOTE: Join layers before FindAllMarkers (required for Seurat v5)
porcn.combined.harmony <- JoinLayers(object = porcn.combined.harmony)

# Find differentially expressed genes for each cluster
porcn.markers <- FindAllMarkers(
    porcn.combined.harmony,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = FALSE
)

write.csv(
    porcn.markers,
    file = "./results/03_clustering/cluster_markers_all.csv",
    row.names = FALSE
)

# Get top markers per cluster
top_markers <- porcn.markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)

write.csv(
    top_markers,
    file = "./results/03_clustering/top_markers_per_cluster.csv",
    row.names = FALSE
)

################################################################################
# 3. Visualize Canonical Cell Type Markers
################################################################################

# Define canonical markers for major cell types
canonical_markers <- list(
    Tcell_NK = c('Ptprc', 'Cd3e', 'Cd3d', 'Cd4', 'Cd8a', 'Nkg7', 'Klrg1', 
                 'Sell', 'Prf1', 'Eomes', 'Trdc', 'Ccr7', 'Il7r', 'Cd28'),
    
    Epithelial = c('Krt19', 'Cdh1', 'Epcam', 'Gli3', 'Tff1', 'Tff2', 'Tff3', 
                   'Mki67', 'Top2a', 'Ido1', 'Ido2', 'Msln'),
    
    Myeloid = c('Ptprc', 'Cd68', 'Adgre1', 'Itgam', 'Cd14', 'Mrc1', 'H2-Eb1', 
                'Batf3', 'S100a8', 'Ly6c2', 'Ly6g', 'Ccr7', 'Clec9a'),
    
    Stromal = c('Col1a2', 'Pdpn', 'Pdgfa', 'Dcn', 'Cdh11', 'Apoe', 
                'Pecam1', 'Cdh5'),
    
    Bcell = c('Cd79a', 'Cd19', 'Ms4a1')
)

# Generate feature plots
for (celltype in names(canonical_markers)) {
    features <- canonical_markers[[celltype]]
    
    p <- FeaturePlot(
        object    = porcn.combined.harmony,
        features  = features,
        cols      = c("grey", "blue"),
        reduction = "umap",
        pt.size   = 0.1,
        ncol      = 4
    )
    
    ggsave(
        filename = paste0("./results/03_clustering/FeaturePlot_", celltype, ".png"),
        plot     = p,
        width    = 16,
        height   = ceiling(length(features) / 4) * 4,
        dpi      = 300
    )
}

################################################################################
# 4. Cell Type Annotation
################################################################################

# Store original cluster assignments
porcn.combined.harmony[["UMAP_Clusters"]] <- Idents(object = porcn.combined.harmony)

# NOTE: This annotation mapping is based on marker expression analysis
# Adjust cluster numbers according to your actual data
# Example mapping (replace with your actual cluster assignments):

new.cluster.ids <- c(
    "0"  = "Mki67hiTop2ahi Cancer cell",
    "1"  = "Tffhi Cancer cell",
    "2"  = "CTL",
    "3"  = "Mki67intTop2aint Cancer cell",
    "4"  = "Gli3hi Cancer cell",
    "5"  = "Pdgfahi CAF",
    "6"  = "Mono/Mac",
    "7"  = "CD4 T cell",
    "8"  = "Dcnhi CAF",
    "9"  = "Mki67hiTop2Ahi CTL",
    "10" = "Idohi cancer cell",
    "11" = "Batf3hi DC",
    "12" = "Neutrophil",
    "13" = "Epi/Cancer cell",
    "14" = "NK cell",
    "15" = "Mki67hiTop2AhiBatf3hi DC",
    "16" = "Mki67hiTop2Ahi B cell",
    "17" = "Cdh11hi CAF",
    "18" = "Apoehi CAF",
    "19" = "Epcamhi B cell",
    "20" = "B cell",
    "21" = "Endothelial cell"
)

# Apply annotations
porcn.combined.harmony$NH_labels <- new.cluster.ids[as.character(Idents(porcn.combined.harmony))]
Idents(porcn.combined.harmony) <- "NH_labels"

################################################################################
# 5. Create Myeloid Subset and Subclustering
################################################################################

# NOTE: Myeloid cells require finer resolution
# Extract and re-cluster for subtype identification

Myeloid_subset <- subset(
    porcn.combined.harmony,
    idents = c("Mono/Mac", "Batf3hi DC", "Mki67hiTop2AhiBatf3hi DC", "Neutrophil")
)

# Re-run standard workflow on myeloid subset
Myeloid_subset <- NormalizeData(Myeloid_subset)
Myeloid_subset <- FindVariableFeatures(Myeloid_subset)
Myeloid_subset <- ScaleData(Myeloid_subset)
Myeloid_subset <- RunPCA(Myeloid_subset)
Myeloid_subset <- RunUMAP(Myeloid_subset, dims = 1:30)

# Higher resolution clustering for myeloid subtypes
Myeloid_subset <- FindNeighbors(Myeloid_subset, dims = 1:30)
Myeloid_subset <- FindClusters(Myeloid_subset, resolution = 0.5)

# Annotate myeloid subclusters
# NOTE: Based on marker expression (Trem2, Selenop, Hexb, S100a8, etc.)
myeloid_annotations <- c(
    "0" = "Classical Monocyte",
    "1" = "Activated Monocyte",
    "2" = "Trem2+ Macrophage",
    "3" = "Selenop+ Macrophage",
    "4" = "Hexb+ Macrophage"
)

Myeloid_subset$Myeloid_labels <- myeloid_annotations[as.character(Idents(Myeloid_subset))]

# Filter low-quality cells
Myeloid_subset_fil <- subset(
    Myeloid_subset,
    subset = nFeature_RNA > 500 & nCount_RNA > 1000
)

################################################################################
# 6. Create Complete_Labels (Merged Cell Types)
################################################################################

# Create simplified cell type labels by merging subtypes
porcn.combined.harmony$Complete_Labels <- case_when(
    # B-cell subtypes → B cell
    porcn.combined.harmony$NH_labels %in% c(
        "B cell", "Epcamhi B cell", "Mki67hiTop2Ahi B cell"
    ) ~ "B cell",
    
    # CTL subtypes → CD8 T cell
    porcn.combined.harmony$NH_labels %in% c(
        "CTL", "Mki67hiTop2Ahi CTL"
    ) ~ "CD8 T cell",
    
    # DC subtypes → DC
    porcn.combined.harmony$NH_labels %in% c(
        "Batf3hi DC", "Mki67hiTop2AhiBatf3hi DC"
    ) ~ "DC",
    
    # Keep other labels as is
    TRUE ~ as.character(porcn.combined.harmony$NH_labels)
)

# Set factor levels for plotting order
porcn.combined.harmony$Complete_Labels <- factor(
    porcn.combined.harmony$Complete_Labels,
    levels = c(
        "Basal cancer cell", "Classical cancer cell",
        "CD8 T cell", "NK cell", "CD4 T cell",
        "CAF",
        "DC", "Mono/Mac", "Neutrophil",
        "B cell",
        "Endothelial cell"
    )
)

################################################################################
# 7. Generate Annotated UMAP
################################################################################

# UMAP by cell type
umap_annotated <- DimPlot(
    porcn.combined.harmony,
    reduction = "umap",
    label = TRUE,
    pt.size = 0.5,
    label.size = 3
) +
    ggtitle("Annotated Cell Types")

ggsave(
    filename = "./results/03_clustering/UMAP_annotated.png",
    plot     = umap_annotated,
    width    = 12,
    height   = 10,
    dpi      = 300
)

# UMAP split by condition
umap_split <- DimPlot(
    porcn.combined.harmony,
    reduction = "umap",
    split.by  = "ID",
    label     = TRUE,
    pt.size   = 0.3
) +
    ggtitle("Cell Types by Condition")

ggsave(
    filename = "./results/03_clustering/UMAP_split_condition.png",
    plot     = umap_split,
    width    = 16,
    height   = 7,
    dpi      = 300
)

# Myeloid subset UMAP
if (exists("Myeloid_subset_fil")) {
    myeloid_umap <- DimPlot(
        Myeloid_subset_fil,
        reduction = "umap",
        group.by  = "Myeloid_labels",
        label     = TRUE,
        pt.size   = 1
    ) +
        ggtitle("Myeloid Cell Subtypes")
    
    ggsave(
        filename = "./results/03_clustering/UMAP_myeloid_subtypes.png",
        plot     = myeloid_umap,
        width    = 10,
        height   = 8,
        dpi      = 300
    )
}

################################################################################
# 8. Save Annotated Objects
################################################################################

save(porcn.combined.harmony, 
     file = "./data/porcn.combined.harmony_annotated.RData")

if (exists("Myeloid_subset_fil")) {
    save(Myeloid_subset_fil, 
         file = "./data/Myeloid_subset_filtered.RData")
}

# Clean up
rm(porcn.markers, top_markers, canonical_markers)
if (exists("Myeloid_subset")) rm(Myeloid_subset)
gc()

################################################################################
# End of Step 03
################################################################################
