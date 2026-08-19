# ============================================================
# HABITAT SUITABILITY MODEL (RANDOM FOREST)
# ============================================================

library(terra)
library(sf)
library(randomForest)
library(pROC)
library(dplyr)

set.seed(42)

# ============================================================
# INPUTS
# ============================================================
# screnaio 1  - running raw geomad_raw predictors
predictor_dir_geomad_raw_kenya <- "/home/vincent/Development/geospatial-analysis/Biodiversity-Project/Shared_Functions/GeoMeadian/outputs/geomad_raw_kenya"

# screnaio 2  - running geomad indices predictors
predictor_dir_geomad_indices_kenya = "/home/vincent/Development/geospatial-analysis/Biodiversity-Project/Shared_Functions/GeoMeadian/outputs/geomad_indices_kenya"

# screnaio 3  - running tassled cap indices, humanfootprint and rainfall and temperature varliables predictors
predictor_dir_TCI_Climate_Kenya="/home/vincent/Development/geospatial-analysis/Biodiversity-Project/Shared_Functions/GeoMeadian/outputs/TCI_Climate_Kenya"


# ========================================================

predictor_dir = predictor_dir_geomad_raw_kenya
occurrence_file <- "/run/media/vincent/Extreme Pro/MyBiodiversity/Biodiversity/Resubmit/Optimized/data/Kenya_total_presence_points.csv"

output_dir <- "/home/vincent/Development/geospatial-analysis/Biodiversity-Project/Shared_Functions/GeoMeadian/outputs"

aoi_file <- "/run/media/vincent/Extreme Pro/Data/shapefiles/Kenya.shp"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD AOI
# ============================================================

cat("Loading Kenya boundary...\n")

kenya <- vect(aoi_file)

# ============================================================
# LOAD PREDICTORS
# ============================================================

cat("Loading predictor rasters...\n")

raster_files <- list.files(
  predictor_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

if(length(raster_files) == 0){
  stop("No raster files found.")
}

predictors <- rast(raster_files)

names(predictors) <- make.names(
  names(predictors),
  unique = TRUE
)

cat("Predictors loaded:", nlyr(predictors), "\n")

# ============================================================
# VERIFY CRS
# ============================================================

if(is.na(crs(predictors)) || crs(predictors) == ""){
  stop("Predictor rasters have no CRS.")
}

kenya <- project(kenya, predictors)

# ============================================================
# CROP TO KENYA
# ============================================================

predictors <- crop(predictors, kenya)
predictors <- mask(predictors, kenya)

# ============================================================
# REMOVE HIGH CORRELATION
# ============================================================

cat("Checking correlations...\n")

sample_vals <- spatSample(
  predictors,
  size = 10000,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE
)

sample_vals <- sample_vals[
  complete.cases(sample_vals),
]

cor_mat <- cor(sample_vals)

remove_correlated <- function(cor_mat, cutoff = 0.8){
  
  cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
  
  remove <- c()
  
  repeat{
    
    max_cor <- suppressWarnings(
      max(abs(cor_mat), na.rm = TRUE)
    )
    
    if(is.infinite(max_cor) || max_cor < cutoff)
      break
    
    idx <- which(
      abs(cor_mat) == max_cor,
      arr.ind = TRUE
    )[1, ]
    
    remove <- c(
      remove,
      colnames(cor_mat)[idx[2]]
    )
    
    cor_mat[, idx[2]] <- NA
    cor_mat[idx[2], ] <- NA
  }
  
  unique(remove)
}

drop_vars <- remove_correlated(cor_mat, 0.8)

if(length(drop_vars) > 0){
  
  cat("Removing correlated predictors:\n")
  print(drop_vars)
  
  predictors <- predictors[[!(
    names(predictors) %in% drop_vars
  )]]
}

cat("Remaining predictors:", nlyr(predictors), "\n")

# ============================================================
# LOAD OCCURRENCES
# ============================================================

cat("Loading occurrences...\n")

occ <- read.csv(
  occurrence_file,
  stringsAsFactors = FALSE
)

# Detect coordinate columns

lon_col <- grep(
  "lon|longitude|x",
  names(occ),
  ignore.case = TRUE,
  value = TRUE
)[1]

lat_col <- grep(
  "lat|latitude|y",
  names(occ),
  ignore.case = TRUE,
  value = TRUE
)[1]

if(is.na(lon_col) || is.na(lat_col)){
  stop("Longitude/Latitude columns not found.")
}

occ <- occ %>%
  filter(
    !is.na(.data[[lon_col]]),
    !is.na(.data[[lat_col]])
  )

# ============================================================
# CONVERT TO SPATVECTOR
# ============================================================

occ_vect <- vect(
  occ,
  geom = c(lon_col, lat_col),
  crs = "EPSG:4326"
)

occ_vect <- project(
  occ_vect,
  crs(predictors)
)

# Keep only points inside Kenya

occ_vect <- occ_vect[
  !is.na(
    extract(
      predictors[[1]],
      occ_vect
    )[,2]
  )
]

cat(
  "Occurrences in Kenya:",
  nrow(occ_vect),
  "\n"
)

# ============================================================
# REMOVE DUPLICATE CELLS
# ============================================================

cells <- cellFromXY(
  predictors[[1]],
  crds(occ_vect)
)

occ_vect <- occ_vect[
  !duplicated(cells)
]

cat(
  "Occurrences after thinning:",
  nrow(occ_vect),
  "\n"
)

# ============================================================
# PRESENCE DATA
# ============================================================

pres_vals <- extract(
  predictors,
  occ_vect
)

pres_vals <- pres_vals[, -1]

pres_vals$presence <- 1

# ============================================================
# BACKGROUND POINTS (1:1)
# ============================================================

cat("Generating background points...\n")

n_pres <- nrow(pres_vals)

bg_points <- spatSample(
  predictors[[1]],
  size = n_pres * 3,
  method = "random",
  na.rm = TRUE,
  as.points = TRUE
)

# Remove any point falling in same raster cell as presences

pres_cells <- cellFromXY(
  predictors[[1]],
  crds(occ_vect)
)

bg_cells <- cellFromXY(
  predictors[[1]],
  crds(bg_points)
)

bg_points <- bg_points[
  !(bg_cells %in% pres_cells)
]

bg_points <- bg_points[
  sample(
    1:nrow(bg_points),
    n_pres
  )
]

bg_vals <- extract(
  predictors,
  bg_points
)

bg_vals <- bg_vals[, -1]

bg_vals$presence <- 0

# ============================================================
# COMBINE DATA
# ============================================================

sdm_data <- rbind(
  pres_vals,
  bg_vals
)

sdm_data <- sdm_data[
  complete.cases(sdm_data),
]

cat(
  "Training samples:",
  nrow(sdm_data),
  "\n"
)

# ============================================================
# TRAIN / TEST SPLIT
# ============================================================

train_idx <- sample(
  seq_len(nrow(sdm_data)),
  size = round(
    0.8 * nrow(sdm_data)
  )
)

train <- sdm_data[train_idx, ]
test  <- sdm_data[-train_idx, ]

# ============================================================
# RANDOM FOREST
# ============================================================

cat("Training Random Forest...\n")

rf_model <- randomForest(
  
  as.factor(presence) ~ .,
  
  data = train,
  
  ntree = 1000,
  
  importance = TRUE
)

print(rf_model)

# ============================================================
# EVALUATION
# ============================================================

pred_test <- predict(
  rf_model,
  test,
  type = "prob"
)[,2]

auc <- pROC::auc(
  test$presence,
  pred_test
)

cat(
  "\nAUC =",
  round(auc, 4),
  "\n"
)

# ============================================================
# PREDICT HABITAT SUITABILITY
# ============================================================

cat("Predicting suitability raster...\n")

rf_fun <- function(model, data){
  
  predict(
    model,
    data,
    type = "prob"
  )[,2]
}

suitability <- terra::predict(
  predictors,
  rf_model,
  fun = rf_fun,
  na.rm = TRUE
)

names(suitability) <- "HabitatSuitability"

# ============================================================
# SAVE OUTPUT
# ============================================================

output_raster <- file.path(
  output_dir,
  "Kenya_Habitat_Suitability_RF.tif"
)

writeRaster(
  suitability,
  output_raster,
  overwrite = TRUE
)

cat(
  "\nSaved suitability raster:\n",
  output_raster,
  "\n"
)

# ============================================================
# VARIABLE IMPORTANCE
# ============================================================

importance_df <- data.frame(
  Variable = rownames(importance(rf_model)),
  Importance = importance(rf_model)[,1]
)

importance_df <- importance_df[
  order(
    importance_df$Importance,
    decreasing = TRUE
  ),
]

print(importance_df)

# ============================================================
# PLOT
# ============================================================

plot(
  suitability,
  main = paste0(
    "Habitat Suitability (AUC=",
    round(auc,3),
    ")"
  )
)

plot(
  kenya,
  add = TRUE
)