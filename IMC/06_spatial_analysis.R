################################################################################
# IMC Analysis Pipeline
# Step 06: Spatial Distance Analysis
# 
# Description:
#   - Define cell type categories (Epithelial, Immune, Myeloid)
#   - Region-based cell-cell distance analysis
#   - Calculate minimum distances between cell type pairs
#   - Four interaction categories:
#     * Immune → Epithelial (immune infiltration to tumor)
#     * Immune → Myeloid (T cell proximity to myeloid)
#     * Myeloid → Immune (myeloid proximity to T cells)
#     * Myeloid → Epithelial (myeloid infiltration to tumor)
#   - Statistical testing (Wilcoxon rank-sum test: Wt vs Ko)
#   - Multiple testing correction (Benjamini-Hochberg FDR)
#
# Input:
#   - ./data/spe_filtered.rds (from Step 04)
#
# Output:
#   - ./figure/06_spatial_analysis/region_distances_raw.csv
#   - ./figure/06_spatial_analysis/region_distances_statistics.csv
#   - Preliminary visualization plots
#
# Author: YMS
# Date: 2024-11-17
################################################################################

# Load required packages
library(imcRtools)
library(SpatialExperiment)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(patchwork)

# Set parameters
SPATIAL_PARAMS <- list(
    groups_to_compare = c("Wt", "Ko"),   # 비교할 그룹 (필수)
    min_samples       = 2,                # 통계 검정을 위한 최소 샘플 수
    min_cells         = 5,                # 거리 계산을 위한 최소 세포 수
    seed              = 900223
)

# Create output directory
if (!dir.exists("./figure/06_spatial_analysis")) {
    dir.create("./figure/06_spatial_analysis", recursive = TRUE)
}

################################################################################
# 1. Load Filtered Data
################################################################################

spe <- readRDS("./data/spe_filtered.rds")

# IMPORTANT: 비교할 그룹만 필터링
spe_region <- spe[, spe$group %in% SPATIAL_PARAMS$groups_to_compare]

# NOTE: region 컬럼이 없으면 에러 발생 (lisaClust 등으로 사전 분석 필요)
if (!"region" %in% colnames(colData(spe_region))) {
    stop("Column 'region' not found. Run spatial regionalization first (e.g., lisaClust).")
}

################################################################################
# 2. Define Cell Type Categories
################################################################################

# IMPORTANT: 본인의 cell type annotation에 맞게 수정 필요
all_celltypes <- unique(spe_region$cluster_celltype)

# Epithelial cells: "Epithelial" 문자열 포함
epithelial_types <- all_celltypes[grepl("Epithelial cell", all_celltypes)]

# Immune T cells: CTL, CD4 T cell 등
immune_types <- c("CTL", "CD4 T cell")
immune_types <- immune_types[immune_types %in% all_celltypes]

# Myeloid cells: SC1-SC35 또는 "Myeloid" 포함
# NOTE: Myeloid subcluster가 SC1, SC2 형식이면 첫 번째 조건 사용
if (any(grepl("^SC[0-9]+$", all_celltypes))) {
    myeloid_types <- all_celltypes[grepl("^SC[0-9]+$", all_celltypes)]
} else {
    myeloid_types <- all_celltypes[grepl("Myeloid|Monocyte|Macrophage", all_celltypes)]
}
myeloid_types <- sort(myeloid_types)

################################################################################
# 3. Define Analysis Combinations
################################################################################

# IMPORTANT: 4가지 카테고리의 cell-cell interaction 정의
analysis_combinations <- list()

# Category 1: Immune → Epithelial (면역세포가 종양으로 침투)
for (immune in immune_types) {
    for (epi in epithelial_types) {
        analysis_combinations[[length(analysis_combinations) + 1]] <- 
            list(from = immune, to = epi, category = "Immune->Epithelial")
    }
}

# Category 2: Immune → Myeloid (T세포와 골수세포 근접도)
for (immune in immune_types) {
    for (myeloid in myeloid_types) {
        analysis_combinations[[length(analysis_combinations) + 1]] <- 
            list(from = immune, to = myeloid, category = "Immune->Myeloid")
    }
}

# Category 3: Myeloid → Immune (골수세포와 T세포 근접도)
for (myeloid in myeloid_types) {
    for (immune in immune_types) {
        analysis_combinations[[length(analysis_combinations) + 1]] <- 
            list(from = myeloid, to = immune, category = "Myeloid->Immune")
    }
}

# Category 4: Myeloid → Epithelial (골수세포가 종양으로 침투)
for (myeloid in myeloid_types) {
    for (epi in epithelial_types) {
        analysis_combinations[[length(analysis_combinations) + 1]] <- 
            list(from = myeloid, to = epi, category = "Myeloid->Epithelial")
    }
}

################################################################################
# 4. Region-Based Distance Calculation
################################################################################

# NOTE: 대용량 데이터의 경우 수 시간 소요 가능
# 진행상황은 내부적으로 추적됨

regions <- sort(unique(na.omit(spe_region$region)))
within_results <- list()

start_time <- Sys.time()

for (reg in regions) {
    
    # 현재 region의 세포 인덱스
    region_cells <- which(spe_region$region == reg)
    available_types <- unique(spe_region$cluster_celltype[region_cells])
    
    for (combo in analysis_combinations) {
        from_type <- combo$from
        to_type <- combo$to
        category <- combo$category
        
        # IMPORTANT: 두 cell type이 모두 해당 region에 존재하는지 확인
        if (!from_type %in% available_types || !to_type %in% available_types) {
            next
        }
        
        # 각 cell type의 인덱스 추출
        to_cells_idx <- which(
            spe_region$region == reg & 
            spe_region$cluster_celltype == to_type
        )
        from_cells_idx <- which(
            spe_region$region == reg & 
            spe_region$cluster_celltype == from_type
        )
        
        # NOTE: 세포 수가 너무 적으면 건너뜀 (최소 5개)
        if (length(to_cells_idx) < SPATIAL_PARAMS$min_cells || 
            length(from_cells_idx) < SPATIAL_PARAMS$min_cells) {
            next
        }
        
        tryCatch({
            # Create logical vector for target cells
            to_cells_logical <- rep(FALSE, ncol(spe_region))
            to_cells_logical[to_cells_idx] <- TRUE
            
            # IMPORTANT: minDistToCells() 함수로 최소 거리 계산
            # 각 "from" 세포에서 가장 가까운 "to" 세포까지의 거리
            temp_spe <- minDistToCells(
                spe_region,
                x_cells = to_cells_logical,
                img_id = "sample_id",
                name = "temp_dist"
            )
            
            # 거리 데이터 추출 및 요약
            dist_data <- colData(temp_spe)[from_cells_idx, ] %>%
                as_tibble() %>%
                filter(!is.na(temp_dist)) %>%
                group_by(sample_id, group) %>%
                summarize(
                    n_cells = n(),
                    mean_dist = mean(temp_dist, na.rm = TRUE),      # 평균 거리
                    median_dist = median(temp_dist, na.rm = TRUE),  # 중앙값
                    sd_dist = sd(temp_dist, na.rm = TRUE),          # 표준편차
                    q25 = quantile(temp_dist, 0.25, na.rm = TRUE),  # 1사분위수
                    q75 = quantile(temp_dist, 0.75, na.rm = TRUE),  # 3사분위수
                    .groups = "drop"
                ) %>%
                mutate(
                    region = reg,
                    from_celltype = from_type,
                    to_celltype = to_type,
                    category = category
                )
            
            # 결과 저장 (고유 키 생성)
            if (nrow(dist_data) > 0) {
                key <- paste(reg, from_type, to_type, sep = "___")
                within_results[[key]] <- dist_data
            }
        }, error = function(e) {
            # 에러 발생 시 건너뜀
            NULL
        })
    }
}

end_time <- Sys.time()

################################################################################
# 5. Combine and Save Raw Results
################################################################################

if (length(within_results) > 0) {
    within_all <- bind_rows(within_results)
    
    # IMPORTANT: 원본 거리 측정 데이터 저장 (모든 후속 분석의 기초)
    write_csv(
        within_all,
        "./figure/06_spatial_analysis/region_distances_raw.csv"
    )
    
} else {
    stop("No distance results calculated. Check cell type definitions and region assignments.")
}

################################################################################
# 6. Statistical Testing
################################################################################

# 그룹별 샘플 수 확인
sample_counts <- within_all %>%
    group_by(region, from_celltype, to_celltype, group) %>%
    summarize(n_samples = n(), .groups = "drop") %>%
    pivot_wider(
        names_from = group,
        values_from = n_samples,
        values_fill = 0,
        names_prefix = "n_"
    )

# IMPORTANT: 통계 검정을 위한 최소 샘플 수 필터링
# 각 그룹에 최소 2개 이상의 샘플이 있어야 검정 가능
valid_pairs <- sample_counts %>%
    filter(if_all(starts_with("n_"), ~ . >= SPATIAL_PARAMS$min_samples))

if (nrow(valid_pairs) > 0) {
    
    # Wilcoxon rank-sum test (non-parametric test for two groups)
    stats_results <- within_all %>%
        semi_join(valid_pairs, by = c("region", "from_celltype", "to_celltype")) %>%
        group_by(region, from_celltype, to_celltype, category) %>%
        wilcox_test(mean_dist ~ group) %>%
        adjust_pvalue(method = "BH") %>%  # Benjamini-Hochberg FDR correction
        add_significance("p.adj") %>%
        arrange(p.adj)
    
    # IMPORTANT: 통계 결과 저장
    write_csv(
        stats_results,
        "./figure/06_spatial_analysis/region_distances_statistics.csv"
    )
    
}

################################################################################
# 7. Preliminary Visualization
################################################################################

# Category-wise comparison (카테고리별 전체 비교)
if (nrow(within_all) > 0) {
    
    category_summary <- within_all %>%
        group_by(category, group) %>%
        summarize(
            n = n(),
            mean_dist = mean(mean_dist, na.rm = TRUE),
            se_dist = sd(mean_dist, na.rm = TRUE) / sqrt(n()),  # Standard error
            .groups = "drop"
        )
    
    p_category <- ggplot(
        category_summary,
        aes(x = category, y = mean_dist, fill = group)
    ) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
        geom_errorbar(
            aes(ymin = mean_dist - se_dist, ymax = mean_dist + se_dist),
            position = position_dodge(0.9),
            width = 0.3
        ) +
        scale_fill_manual(
            values = c("Wt" = "#4DAF4A", "Ko" = "#E41A1C")
        ) +
        theme_classic(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
        labs(
            x = "",
            y = "Mean Distance (μm)",
            title = "Cell-Cell Distance by Category",
            fill = "Group"
        )
    
    ggsave(
        filename = "./figure/06_spatial_analysis/01_category_summary.tiff",
        plot = p_category,
        width = 11,
        height = 7,
        dpi = 400,
        compression = "lzw"
    )
}

# Top significant pairs boxplot (유의미한 결과가 있을 경우)
if (exists("stats_results") && nrow(stats_results) > 0) {
    
    # 유의미한 결과만 필터링
    sig_results <- stats_results %>% filter(p.adj < 0.05)
    
    if (nrow(sig_results) > 0) {
        n_plots <- min(12, nrow(sig_results))
        plot_list <- list()
        
        for (i in 1:n_plots) {
            data_sub <- within_all %>%
                filter(
                    region == sig_results$region[i],
                    from_celltype == sig_results$from_celltype[i],
                    to_celltype == sig_results$to_celltype[i]
                )
            
            p <- ggplot(data_sub, aes(x = group, y = mean_dist, fill = group)) +
                geom_boxplot(outlier.shape = NA, alpha = 0.7) +
                geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
                scale_fill_manual(values = c("Wt" = "#4DAF4A", "Ko" = "#E41A1C")) +
                theme_classic(base_size = 9) +
                labs(
                    title = sprintf("%s\n%s -> %s (p=%.1e)",
                                    sig_results$region[i],
                                    sig_results$from_celltype[i],
                                    sig_results$to_celltype[i],
                                    sig_results$p.adj[i]),
                    x = "",
                    y = "Distance (μm)"
                ) +
                stat_compare_means(
                    method = "wilcox.test",
                    label = "p.format",
                    size = 2.5
                ) +
                theme(
                    legend.position = "none",
                    plot.title = element_text(size = 8, face = "bold")
                )
            
            plot_list[[i]] <- p
        }
        
        p_combined <- wrap_plots(plot_list, ncol = 4)
        
        ggsave(
            filename = "./figure/06_spatial_analysis/02_significant_pairs_boxplot.tiff",
            plot = p_combined,
            width = 16,
            height = 12,
            dpi = 400,
            compression = "lzw"
        )
    }
}

################################################################################
# 8. Clean Up
################################################################################

rm(spe_region, within_results)
if (exists("temp_spe")) rm(temp_spe)
gc()

# NOTE: Next step is 07_advanced_visualization.R
