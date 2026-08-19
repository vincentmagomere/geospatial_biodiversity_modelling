# ============================================================
# ROBUST MULTI-SCENARIO SDM COMPARISON (6 SCENARIOS)
# GeoMAD (Raw / Indices) vs. TCI (Indices / Climate Integration)
# ============================================================

library(terra)
library(sf)
library(randomForest)
library(pROC)
library(dplyr)
library(ggplot2)
library(tidyr)

# Ensure reproducibility across runs
set.seed(42)

# ============================================================
# 1. INPUT CONFIGURATION & PATHS
# ============================================================

scenarios <- list(
  GeoMAD_Raw             = "/home/vincent/Development/geospatial-analysis/Biodiversity-Project/Shared_Functions/GeoMeadian/outputs/geomad_raw_kenya",
  GeoMad_Raw_Climate     = "/home/vincent/Development/geospatial-analysis/Biodiversity-Project/Shared_Functions/GeoMeadian/outputs/geomad_raw_climate_kenya",
  
  GeoMAD_Indices         = "/home/vincent/Development/geospatial-analysis/Biodiversity-Project/Shared_Functions/GeoMeadian/outputs/geomad_indices_kenya",
  GeoMAd_Indices_Climate = "/home/vincent/Development/geospatial-analysis/Biodiversity-Project/Shared_Functions/GeoMeadian/outputs/geomad_indices_climate_kenya",
  
  TCI_Indices            = "/home/vincent/Development/geospatial-analysis/Biodiversity-Project/Shared_Functions/GeoMeadian/outputs/TCI_indices_Kenya",
  TCI_Climate            = "/home/vincent/Development/geospatial-analysis/Biodiversity-Project/Shared_Functions/GeoMeadian/outputs/TCI_Climate_Kenya"
)

occurrence_file <- "/home/vincent/Development/geospatial-analysis/Biodiversity-Project/Shared_Functions/GeoMeadian/outputs/Shapefiles/Kenya_total_presence_points.csv"
aoi_file        <- "/home/vincent/Development/geospatial-analysis/Biodiversity-Project/Shared_Functions/GeoMeadian/outputs/Shapefiles/Kenya.shp"
output_dir      <- "/home/vincent/Development/geospatial-analysis/Biodiversity-Project/Shared_Functions/GeoMeadian/outputs/Scenario_Comparison_Results"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

# Helper to safely load and align rasters with mismatched extents/resolutions
load_and_align_predictors <- function(raster_files) {
  if (length(raster_files) == 1) {
    return(terra::rast(raster_files))
  }
  
  # Try loading all at once first (fastest if already matching)
  r_try <- try(terra::rast(raster_files), silent = TRUE)
  if (!inherits(r_try, "try-error")) {
    return(r_try)
  }
  
  cat("  Mismatch in extents/resolutions detected. Aligning and resampling layers...\n")
  
  # Load individually
  r_list <- lapply(raster_files, terra::rast)
  
  # Find reference raster: pick the one with the highest spatial resolution (smallest cell size)
  cell_sizes <- sapply(r_list, function(r) mean(terra::res(r)))
  ref_idx <- which.min(cell_sizes)
  ref <- r_list[[ref_idx]][[1]] # single-layer spatial template
  
  aligned_list <- list()
  for (i in seq_along(r_list)) {
    r <- r_list[[i]]
    
    # Check if geometry matches template grid
    if (!terra::compareGeom(r, ref, lyrs = FALSE, crs = TRUE, ext = TRUE, rowcol = TRUE, stopOnError = FALSE)) {
      # Reproject if CRS differs
      if (terra::crs(r) != terra::crs(ref)) {
        r <- terra::project(r, ref, method = "bilinear")
      }
      # Resample to match reference extent and grid resolution
      r <- terra::resample(r, ref, method = "bilinear")
    }
    aligned_list[[i]] <- r
  }
  
  return(terra::rast(aligned_list))
}

# Function to remove highly collinear variables (r > 0.95 cutoff)
remove_correlated <- function(cor_mat, cutoff = 0.95){
  cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
  remove <- c()
  repeat{
    max_cor <- suppressWarnings(max(abs(cor_mat), na.rm = TRUE))
    if(is.infinite(max_cor) || max_cor < cutoff) break
    idx <- which(abs(cor_mat) == max_cor, arr.ind = TRUE)[1, ]
    remove <- c(remove, colnames(cor_mat)[idx[2]])
    cor_mat[, idx[2]] <- NA
    cor_mat[idx[2], ] <- NA
  }
  return(unique(remove))
}

# Function to plot and save correlation matrix heatmap
save_correlation_heatmap <- function(cor_mat, scen_name, save_path) {
  cor_df <- as.data.frame(cor_mat) %>%
    mutate(Var1 = rownames(cor_mat)) %>%
    pivot_longer(-Var1, names_to = "Var2", values_to = "Correlation")
  
  p_cor <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Correlation, 2)), size = 2.2) +
    scale_fill_gradient2(low = "#2b5c8f", mid = "#ffffff", high = "#d95f02", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab", 
                         name = "Pearson\nCorrelation") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          panel.grid.major = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5)) +
    coord_fixed() +
    labs(title = paste0("Predictor Correlation Heatmap: ", scen_name),
         x = "", y = "")
  
  num_vars <- ncol(cor_mat)
  plot_dim <- max(6, min(16, num_vars * 0.55))
  
  ggsave(save_path, plot = p_cor, width = plot_dim + 1.5, height = plot_dim, dpi = 300)
}

# Function to calculate AUC, TSS, Sensitivity, Specificity
calculate_metrics <- function(labels, probs) {
  roc_obj <- pROC::roc(labels, probs, quiet = TRUE)
  auc_val <- pROC::auc(roc_obj)
  
  coords <- pROC::coords(roc_obj, "best", best.method = "youden")
  sens <- coords$sensitivity[1]
  spec <- coords$specificity[1]
  tss  <- sens + spec - 1
  
  return(list(
    auc = as.numeric(auc_val),
    tss = as.numeric(tss),
    sensitivity = as.numeric(sens),
    specificity = as.numeric(spec),
    roc_obj = roc_obj
  ))
}

# Spatial K-Fold assignment using k-means spatial partitioning
assign_spatial_folds <- function(pts, k = 5) {
  coords <- terra::crds(pts)
  km <- kmeans(coords, centers = k, nstart = 10)
  return(km$cluster)
}

# Predictor function wrapper for terra::predict
rf_fun <- function(model, data){
  predict(model, data, type = "prob")[, 2]
}

# ============================================================
# 3. LOAD AOI & OCCURRENCES (COMMON PREPARATION)
# ============================================================

cat("Loading spatial boundaries and occurrence points...\n")
kenya <- terra::vect(aoi_file)

occ <- read.csv(occurrence_file, stringsAsFactors = FALSE)

lon_col <- grep("lon|longitude|x", names(occ), ignore.case = TRUE, value = TRUE)[1]
lat_col <- grep("lat|latitude|y",  names(occ), ignore.case = TRUE, value = TRUE)[1]

if(is.na(lon_col) || is.na(lat_col)){
  stop("Longitude/Latitude columns not found in occurrence file.")
}

occ <- occ %>% dplyr::filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]]))
occ_vect <- terra::vect(occ, geom = c(lon_col, lat_col), crs = "EPSG:4326")

# ============================================================
# 4. SCENARIO EXECUTION LOOP
# ============================================================

scenario_results    <- list()
suitability_rasters <- list()
importance_list     <- list()
test_predictions    <- list()

for (scen_name in names(scenarios)) {
  cat("\n========================================================\n")
  cat(" RUNNING SCENARIO:", scen_name, "\n")
  cat("========================================================\n")
  
  pred_dir <- scenarios[[scen_name]]
  raster_files <- list.files(pred_dir, pattern = "\\.tif$", full.names = TRUE)
  
  if(length(raster_files) == 0){
    warning(paste("No raster files found for scenario:", scen_name, "- Skipping!"))
    next
  }
  
  # Load and auto-align rasters across different resolutions/extents
  predictors <- load_and_align_predictors(raster_files)
  names(predictors) <- make.names(names(predictors), unique = TRUE)
  
  # Align CRS, crop, and mask to Kenya AOI
  kenya_proj <- terra::project(kenya, predictors)
  predictors <- terra::crop(predictors, kenya_proj)
  predictors <- terra::mask(predictors, kenya_proj)
  
  # Correlation Analysis & Plot Generation
  cat("Sampling predictors and calculating correlation matrix...\n")
  sample_vals <- terra::spatSample(predictors, size = 10000, method = "random", na.rm = TRUE, as.df = TRUE)
  sample_vals <- sample_vals[complete.cases(sample_vals), ]
  cor_mat <- cor(sample_vals)
  
  # SAVE CORRELATION HEATMAP
  cor_plot_file <- file.path(output_dir, paste0("Correlation_Matrix_", scen_name, ".png"))
  cat("Saving correlation matrix heatmap to:", cor_plot_file, "\n")
  save_correlation_heatmap(cor_mat, scen_name, cor_plot_file)
  
  # Multicollinearity Removal (Strict 0.95 Threshold)
  cat("Removing collinear predictors (r > 0.95)...\n")
  drop_vars <- remove_correlated(cor_mat, cutoff = 0.95)
  
  if(length(drop_vars) > 0){
    cat("  Dropped variables:", paste(drop_vars, collapse = ", "), "\n")
    keep_vars <- setdiff(names(predictors), drop_vars)
    predictors <- terra::subset(predictors, keep_vars)
  }
  
  cat("  Remaining predictors count:", terra::nlyr(predictors), "\n")
  
  # Align Occurrences Spatially
  occ_scen <- terra::project(occ_vect, terra::crs(predictors))
  extracted_check <- terra::extract(predictors[[1]], occ_scen)
  occ_scen <- occ_scen[!is.na(extracted_check[, 2]), ]
  
  # Thin duplicate cells
  cells <- terra::cellFromXY(predictors[[1]], terra::crds(occ_scen))
  occ_scen <- occ_scen[!duplicated(cells), ]
  
  # Extract Presences
  pres_vals <- terra::extract(predictors, occ_scen)
  if("ID" %in% names(pres_vals)) pres_vals$ID <- NULL
  pres_vals$presence <- 1
  pres_coords <- terra::crds(occ_scen)
  
  # Generate Background Points (1:1 Ratio)
  n_pres <- nrow(pres_vals)
  bg_points <- terra::spatSample(predictors[[1]], size = n_pres * 3, method = "random", na.rm = TRUE, as.points = TRUE)
  
  bg_cells <- terra::cellFromXY(predictors[[1]], terra::crds(bg_points))
  pres_cells <- terra::cellFromXY(predictors[[1]], terra::crds(occ_scen))
  
  bg_points <- bg_points[!(bg_cells %in% pres_cells)]
  
  if(nrow(bg_points) > n_pres) {
    bg_points <- bg_points[sample(1:nrow(bg_points), n_pres), ]
  }
  
  # Extract Background Values
  bg_vals <- terra::extract(predictors, bg_points)
  if("ID" %in% names(bg_vals)) bg_vals$ID <- NULL
  bg_vals$presence <- 0
  bg_coords <- terra::crds(bg_points)
  
  # Combine Dataset
  sdm_data <- rbind(pres_vals, bg_vals)
  all_coords <- rbind(pres_coords, bg_coords)
  
  valid_rows <- complete.cases(sdm_data)
  sdm_data <- sdm_data[valid_rows, ]
  all_coords <- all_coords[valid_rows, ]
  
  # Assign Spatial Folds for Cross-Validation
  folds <- assign_spatial_folds(terra::vect(all_coords, crs = terra::crs(predictors)), k = 5)
  sdm_data$fold <- folds
  
  # 5-Fold Spatial Cross Validation Loop
  cat("Running 5-fold Spatial Cross-Validation...\n")
  cv_preds <- numeric(nrow(sdm_data))
  
  for(k in 1:5) {
    train_fold <- sdm_data %>% dplyr::filter(fold != k) %>% dplyr::select(-fold)
    test_fold  <- sdm_data %>% dplyr::filter(fold == k) %>% dplyr::select(-fold)
    
    rf_cv <- randomForest(as.factor(presence) ~ ., data = train_fold, ntree = 500, importance = FALSE)
    cv_preds[sdm_data$fold == k] <- predict(rf_cv, test_fold, type = "prob")[, 2]
  }
  
  # Evaluate Spatial CV Metrics
  metrics <- calculate_metrics(sdm_data$presence, cv_preds)
  cat(sprintf("Spatial CV AUC: %.4f | TSS: %.4f\n", metrics$auc, metrics$tss))
  
  # Train Final Full Model
  cat("Training Full Random Forest Model...\n")
  full_rf <- randomForest(as.factor(presence) ~ ., data = sdm_data %>% dplyr::select(-fold), ntree = 1000, importance = TRUE)
  
  # Predict Habitat Suitability Raster across Kenya
  cat("Predicting suitability map...\n")
  suitability <- terra::predict(predictors, full_rf, fun = rf_fun, na.rm = TRUE)
  names(suitability) <- paste0("Suitability_", scen_name)
  
  # Save Results Objects
  scenario_results[[scen_name]] <- list(
    AUC = metrics$auc,
    TSS = metrics$tss,
    Sensitivity = metrics$sensitivity,
    Specificity = metrics$specificity,
    ROC_Obj = metrics$roc_obj
  )
  
  test_predictions[[scen_name]] <- list(
    obs = sdm_data$presence,
    pred = cv_preds
  )
  
  suitability_rasters[[scen_name]] <- suitability
  
  # Extract Variable Importance
  imp_df <- data.frame(
    Scenario = scen_name,
    Variable = rownames(importance(full_rf)),
    MeanDecreaseGini = importance(full_rf)[, "MeanDecreaseGini"]
  )
  importance_list[[scen_name]] <- imp_df
  
  # Write Individual Raster Output
  terra::writeRaster(suitability, file.path(output_dir, paste0("Suitability_", scen_name, ".tif")), overwrite = TRUE)
}

# ============================================================
# 5. STATISTICAL COMPARISON & PAIRWISE DELONG TESTS
# ============================================================

cat("\n========================================================\n")
cat(" STATISTICAL COMPARISON & EVALUATION SUMMARY\n")
cat("========================================================\n")

summary_df <- data.frame(
  Scenario    = names(scenario_results),
  AUC         = sapply(scenario_results, function(x) x$AUC),
  TSS         = sapply(scenario_results, function(x) x$TSS),
  Sensitivity = sapply(scenario_results, function(x) x$Sensitivity),
  Specificity = sapply(scenario_results, function(x) x$Specificity)
)

print(summary_df)
write.csv(summary_df, file.path(output_dir, "Scenario_Performance_Metrics.csv"), row.names = FALSE)

# Pairwise DeLong Tests across all scenario pairs
cat("\nPerforming DeLong's Tests for AUC statistical significance...\n")
scen_pairs <- combn(names(scenario_results), 2, simplify = FALSE)

for(pair in scen_pairs) {
  m1_name <- pair[1]
  m2_name <- pair[2]
  
  roc1 <- scenario_results[[m1_name]]$ROC_Obj
  roc2 <- scenario_results[[m2_name]]$ROC_Obj
  
  test_res <- pROC::roc.test(roc1, roc2, method = "delong")
  
  p_val <- test_res$p.value
  sig <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "NS")))
  
  cat(sprintf("Comparison [%s vs %s]: p-value = %.5f (%s)\n", m1_name, m2_name, p_val, sig))
}

# ============================================================
# 6. VISUALIZATION & OUTPUT CHARTS
# ============================================================

# 1. Performance Metrics Bar Chart
metrics_long <- summary_df %>%
  tidyr::pivot_longer(cols = c("AUC", "TSS"), names_to = "Metric", values_to = "Value")

p_metrics <- ggplot(metrics_long, aes(x = Scenario, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  theme_minimal(base_size = 12) +
  labs(title = "Model Performance Comparison (Spatial CV)", y = "Score", x = "") +
  scale_fill_manual(values = c("AUC" = "#2b5c8f", "TSS" = "#d95f02")) +
  ylim(0, 1.0) +
  geom_text(aes(label = round(Value, 3)), position = position_dodge(0.8), vjust = -0.3, size = 3.2) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, face = "bold"))

ggsave(file.path(output_dir, "Scenario_Metrics_Comparison.png"), plot = p_metrics, width = 11, height = 6)

# 2. Variable Importance Plot across 6 Scenarios
all_imp <- do.call(rbind, importance_list)

p_imp <- ggplot(all_imp, aes(x = reorder(Variable, MeanDecreaseGini), y = MeanDecreaseGini, fill = Scenario)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~Scenario, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(title = "Predictor Importance (Mean Decrease Gini)", x = "Predictors", y = "Importance")

ggsave(file.path(output_dir, "Variable_Importance_Comparison.png"), plot = p_imp, width = 13, height = 11)

# 3. Spatial Alignment & Raster Stacking
cat("\nAligning raster extents across scenarios for comparison...\n")
template_raster <- suitability_rasters[[1]]

aligned_rasters <- list()
for(s_name in names(suitability_rasters)) {
  r <- suitability_rasters[[s_name]]
  if(!terra::compareGeom(r, template_raster, lyrs = FALSE, crs = TRUE, ext = TRUE, rowcol = TRUE, stopOnError = FALSE)) {
    cat(sprintf("  Resampling extent for scenario: %s...\n", s_name))
    r <- terra::resample(r, template_raster, method = "bilinear")
  }
  aligned_rasters[[s_name]] <- r
}

suit_stack <- terra::rast(aligned_rasters)

# Calculate difference map (GeoMAD Indices vs TCI Climate if available)
if("Suitability_GeoMAD_Indices" %in% names(suit_stack) && "Suitability_TCI_Climate" %in% names(suit_stack)) {
  diff_raster <- suit_stack[["Suitability_GeoMAD_Indices"]] - suit_stack[["Suitability_TCI_Climate"]]
  names(diff_raster) <- "Difference_GeoMAD_minus_TCI"
  terra::writeRaster(diff_raster, file.path(output_dir, "Difference_GeoMAD_vs_TCI.tif"), overwrite = TRUE)
}

# 4. Multi-Panel Suitability Map Grid
png(file.path(output_dir, "Habitat_Suitability_Maps.png"), width = 1500, height = 1200, res = 150)
par(mfrow = c(3, 3))

for(lyr_name in names(suit_stack)) {
  clean_title <- gsub("Suitability_", "", lyr_name)
  plot(suit_stack[[lyr_name]], main = clean_title, mar = c(2, 2, 2, 4))
}

if(exists("diff_raster")) {
  plot(diff_raster, main = "Diff (GeoMAD Ind - TCI)", col = hcl.colors(50, "Spectral"), mar = c(2, 2, 2, 4))
}

dev.off()

cat("\nAnalysis Complete! All 6 scenarios processed cleanly.\nOutputs saved to:\n", output_dir, "\n")