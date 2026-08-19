library(terra)

# --------------------------------------------------
# SETTINGS
# --------------------------------------------------
raster_dir <- "/run/media/vincent/Extreme Pro/Data/variables/Africa_MS7_250M/intactness"

africa_shp <- "/run/media/vincent/Extreme Pro/FaithAhiono/Species occurrence Africa/Species_occurrence/Species_occurrence/Africa.shp"

output_file <- file.path(
  raster_dir,
  "Africa_Intactness_Mosaic.tif"
)

# Put Terra temp files on the external drive
terraOptions(
  tempdir = file.path(raster_dir, "terra_tmp"),
  memfrac = 0.8,
  progress = 1
)

dir.create(
  file.path(raster_dir, "terra_tmp"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(raster_dir, "clipped_tiles"),
  recursive = TRUE,
  showWarnings = FALSE
)

# --------------------------------------------------
# LOAD AFRICA
# --------------------------------------------------
africa <- vect(africa_shp)

# --------------------------------------------------
# PROCESS EACH TILE INDIVIDUALLY
# --------------------------------------------------
raster_files <- list.files(
  raster_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

clipped_files <- c()

for(i in seq_along(raster_files)) {
  
  cat("\nProcessing:", basename(raster_files[i]), "\n")
  
  r <- rast(raster_files[i])
  
  if (!same.crs(r, africa)) {
    africa_proj <- project(africa, crs(r))
  } else {
    africa_proj <- africa
  }
  
  # Skip tiles completely outside Africa
  if (is.null(intersect(ext(r), ext(africa_proj)))) {
    next
  }
  
  r_crop <- crop(r, africa_proj)
  r_mask <- mask(r_crop, africa_proj)
  
  outfile <- file.path(
    raster_dir,
    "clipped_tiles",
    paste0("tile_", i, ".tif")
  )
  
  writeRaster(
    r_mask,
    outfile,
    overwrite = TRUE,
    wopt = list(
      datatype = "FLT4S",
      gdal = c(
        "COMPRESS=DEFLATE",
        "BIGTIFF=YES"
      )
    )
  )
  
  clipped_files <- c(clipped_files, outfile)
  
  rm(r, r_crop, r_mask)
  gc()
}

# --------------------------------------------------
# MOSAIC CLIPPED TILES
# --------------------------------------------------
cat("\nBuilding final mosaic...\n")

r_list <- lapply(clipped_files, rast)

final <- do.call(
  mosaic,
  c(r_list, fun = "mean")
)

writeRaster(
  final,
  output_file,
  overwrite = TRUE,
  wopt = list(
    datatype = "FLT4S",
    gdal = c(
      "COMPRESS=DEFLATE",
      "BIGTIFF=YES"
    )
  )
)

cat("\nFinished\n")
cat("Saved to:\n", output_file, "\n")