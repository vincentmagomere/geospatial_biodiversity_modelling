library(terra)

#-----------------------------------------
# Paths
#-----------------------------------------
africa_shp <- "data/shapefiles/Africa.shp"

input_dir <- "data/rasters"
output_dir <- "outputs/clipped_rasters"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#-----------------------------------------
# Read Africa boundary
#-----------------------------------------
africa <- vect(africa_shp)

#-----------------------------------------
# Find rasters
#-----------------------------------------
raster_files <- list.files(
  input_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

cat("Found", length(raster_files), "rasters\n")

#-----------------------------------------
# Loop through rasters
#-----------------------------------------
for (f in raster_files) {
  
  cat("Processing:", basename(f), "\n")
  
  r <- rast(f)
  
  africa_proj <- africa
  
  if (!same.crs(r, africa)) {
    africa_proj <- project(africa, crs(r))
  }
  
  r_crop <- crop(r, africa_proj)
  r_mask <- mask(r_crop, africa_proj)
  
  out_file <- file.path(
    output_dir,
    paste0(
      tools::file_path_sans_ext(basename(f)),
      "_Africa.tif"
    )
  )
  
  writeRaster(
    r_mask,
    out_file,
    overwrite = TRUE,
    wopt = list(
      gdal = c(
        "COMPRESS=LZW",
        "TILED=YES"
      )
    )
  )
  
  gc()
}

cat("Finished!\n")