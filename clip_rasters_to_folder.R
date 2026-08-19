library(terra)
library(sf)

# ============================================================
# INPUTS
# ============================================================

input_folder <- "/run/media/vincent/Extreme Pro/Data/variables/Africa_MS7_250M/clipped_rasters"

output_folder <- "/home/vincent/Development/geospatial-analysis/Biodiversity-Project/Shared_Functions/GeoMeadian/outputs/TCI_Climate_Kenya"

kenya_shp <- "/run/media/vincent/Extreme Pro/Data/shapefiles/Kenya.shp"

# Create output folder if it doesn't exist
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD KENYA SHAPEFILE
# ============================================================

cat("Loading Kenya shapefile...\n")

kenya <- vect(kenya_shp)

# ============================================================
# FIND ALL RASTERS
# ============================================================

raster_files <- list.files(
  input_folder,
  pattern = "\\.(tif|img)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

cat("Found", length(raster_files), "raster(s)\n")

# ============================================================
# CLIP EACH RASTER
# ============================================================

for (i in seq_along(raster_files)) {
  
  raster_file <- raster_files[i]
  
  cat("\n[", i, "/", length(raster_files), "] Processing:",
      basename(raster_file), "\n")
  
  try({
    
    # Load raster
    r <- rast(raster_file)
    
    # Reproject shapefile if CRS differs
    kenya_proj <- if (!same.crs(r, kenya)) {
      project(kenya, crs(r))
    } else {
      kenya
    }
    
    # Crop first (faster)
    r_crop <- crop(r, kenya_proj)
    
    # Mask to exact Kenya boundary
    r_clip <- mask(r_crop, kenya_proj)
    
    # Output filename
    output_file <- file.path(
      output_folder,
      basename(raster_file)
    )
    
    # Save
    writeRaster(
      r_clip,
      output_file,
      overwrite = TRUE,
      gdal = c("COMPRESS=LZW")
    )
    
    cat("Saved:", output_file, "\n")
    
  }, silent = FALSE)
  
}

cat("\nFinished clipping all rasters.\n")