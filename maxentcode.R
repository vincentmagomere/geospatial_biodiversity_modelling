# ============================================================
# MAXENT SDM USING maxnet
# ============================================================

library(terra)
library(sf)
library(dplyr)
library(maxnet)
library(pROC)

set.seed(42)

# ============================================================
# INPUTS
# ============================================================

occurrence_file <- "/run/media/vincent/Extreme Pro/MyBiodiversity/Biodiversity/Resubmit/Optimized/data/Kenya_total_presence_points.csv"

predictor_dir <- "/home/vincent/Development/geospatial-analysis/Biodiversity-Project/Shared_Functions/GeoMeadian/outputs/Kenya_raw_1km"

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

cat("Loading predictors...\n")

raster_files <- list.files(
  predictor_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

predictors <- rast(raster_files)

names(predictors) <- make.names(
  names(predictors),
  unique = TRUE
)

kenya <- project(kenya, predictors)

predictors <- crop(predictors, kenya)
predictors <- mask(predictors, kenya)

cat("Predictors:", nlyr(predictors), "\n")

# ============================================================
# REMOVE HIGHLY CORRELATED VARIABLES
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

drop_vars <- remove_correlated(
  cor_mat,
  cutoff = 0.8
)

if(length(drop_vars) > 0){
  
  predictors <- predictors[[
    !(names(predictors) %in% drop_vars)
  ]]
  
  cat("Removed:\n")
  print(drop_vars)
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

occ_vect <- vect(
  occ,
  geom = c(lon_col, lat_col),
  crs = "EPSG:4326"
)

occ_vect <- project(
  occ_vect,
  crs(predictors)
)

# ============================================================
# REMOVE POINTS OUTSIDE PREDICTORS
# ============================================================

occ_vect <- occ_vect[
  !is.na(
    extract(
      predictors[[1]],
      occ_vect
    )[,2]
  )
]

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
  "Presence points:",
  nrow(occ_vect),
  "\n"
)

# ============================================================
# EXTRACT PRESENCE DATA
# ============================================================

pres_env <- extract(
  predictors,
  occ_vect
)

pres_env <- pres_env[, -1]

# ============================================================
# GENERATE 1:1 BACKGROUND
# ============================================================

cat("Generating background points...\n")

n_pres <- nrow(pres_env)

bg_points <- spatSample(
  predictors[[1]],
  size = n_pres * 3,
  method = "random",
  na.rm = TRUE,
  as.points = TRUE
)

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
    seq_len(nrow(bg_points)),
    n_pres
  )
]

bg_env <- extract(
  predictors,
  bg_points
)

bg_env <- bg_env[, -1]

# ============================================================
# BUILD MAXENT DATASET
# ============================================================

presence <- rep(1, nrow(pres_env))
background <- rep(0, nrow(bg_env))

response <- c(
  presence,
  background
)

predictor_data <- rbind(
  pres_env,
  bg_env
)

predictor_data <- predictor_data[
  complete.cases(predictor_data),
]

response <- response[
  complete.cases(
    rbind(pres_env, bg_env)
  )
]

cat(
  "Training records:",
  length(response),
  "\n"
)

# ============================================================
# TRAIN / TEST SPLIT
# ============================================================

train_idx <- sample(
  seq_len(length(response)),
  round(length(response) * 0.8)
)

train_x <- predictor_data[train_idx, ]
train_y <- response[train_idx]

test_x <- predictor_data[-train_idx, ]
test_y <- response[-train_idx]

# ============================================================
# FIT MAXENT
# ============================================================

cat("Training MaxEnt...\n")

mx_model <- maxnet(
  p = train_y,
  data = train_x,
  f = maxnet.formula(
    train_y,
    train_x,
    classes = "lqph"
  )
)

# ============================================================
# EVALUATION
# ============================================================

pred_test <- predict(
  mx_model,
  test_x,
  type = "cloglog"
)

auc <- pROC::auc(
  test_y,
  pred_test
)

cat(
  "\nAUC =",
  round(auc, 4),
  "\n"
)

# ============================================================
# PREDICT SUITABILITY MAP
# ============================================================

cat("Predicting habitat suitability...\n")

mx_fun <- function(model, data){
  
  predict(
    model,
    data,
    type = "cloglog"
  )
}

suitability <- terra::predict(
  predictors,
  mx_model,
  fun = mx_fun,
  na.rm = TRUE
)

names(suitability) <- "HabitatSuitability"

# ============================================================
# MASK TO KENYA
# ============================================================

suitability <- crop(
  suitability,
  kenya
)

suitability <- mask(
  suitability,
  kenya
)

# ============================================================
# SAVE OUTPUT
# ============================================================

output_raster <- file.path(
  output_dir,
  "Kenya_Habitat_Suitability_MaxEnt.tif"
)

writeRaster(
  suitability,
  output_raster,
  overwrite = TRUE
)

cat(
  "\nSaved:\n",
  output_raster,
  "\n"
)

# ============================================================
# VARIABLE CONTRIBUTIONS
# ============================================================

cat("\nModel coefficients:\n")

print(
  sort(
    abs(mx_model$betas),
    decreasing = TRUE
  )
)

# ============================================================
# PLOT
# ============================================================

# plot(
#   suitability,
#   main = paste(
#     "MaxEnt Habitat Suitability",
#     "\nAUC =",
#     round(auc, 3)
#   )
# )
# 
# plot(
#   kenya,
#   add = TRUE
# )