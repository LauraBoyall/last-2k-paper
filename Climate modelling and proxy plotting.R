library(ggplot2) # used for plotting
library(ncdf4) # used for reading and extracting data from NetCDF model file
library(raster) # helpful package for plotting maps
library(dplyr) # helpful tool for filtering data
library(maps) # to create the map outlines
library(sf) # essential for spatial data handing
library(reshape2) # reshaping the spatial 3D data
library(scales)
library(reshape2)# Reshape to long-format data frame
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(akima) # interpolates the model


######### SORTING THE CLIMATE MODELS ##############
setwd("/Volumes/LaCie/PAGES/Hydro2k/")
iso <- nc_open("iCESM_d18O_850-1849.nc")
ncatt_get(iso, "time", "units")


process_model <- function(nc_path, 
                          var_name = "d18O", 
                          lat_name = "lat", 
                          lon_name = "lon", 
                          time_name = "time",
                          time_units = "months since 850",
                          model_name = NULL) {
  
  nc <- nc_open(nc_path)
  TS <- ncvar_get(nc, var_name)
  lat <- ncvar_get(nc, lat_name)
  lon <- ncvar_get(nc, lon_name)
  time <- ncvar_get(nc, time_name)
  
  # Adjust longitude
  lon <- ifelse(lon > 180, lon - 360, lon)
  
  # Extract years from time
  
  if (grepl("%Y%m", time_units)) {
    years <- floor(time / 100)
  } else if (grepl("months since", time_units)) {
    start_year <- as.numeric(gsub(".*since\\s+", "", time_units))
    years <- floor(start_year + time / 12)
  } else if (grepl("years since", time_units)) {
    start_year <- as.numeric(gsub(".*since\\s+", "", time_units))
    years <- floor(start_year + time)
  } else {
    warning("Unrecognized time format")
    years <- time
  }
  
  unique_years <- sort(unique(years))
  
  # Z-score
  TS_z <- TS
  for (i in 1:dim(TS)[1]) {
    for (j in 1:dim(TS)[2]) {
      x <- TS[i, j, ]
      if (all(is.na(x))) next
      TS_z[i, j, ] <- scale(x)[, 1]
    }
  }
  
  # Aggregate annually
  n_years <- length(unique_years)
  X_annual_z <- array(NA, dim = c(dim(TS)[1], dim(TS)[2], n_years))
  
  for (i in seq_along(unique_years)) {
    idx <- which(years == unique_years[i])
    X_annual_z[, , i] <- apply(TS_z[, , idx, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  }
  
  return(list(
    model = model_name,
    X_annual_z = X_annual_z,
    lon = lon,
    lat = lat,
    unique_years = unique_years
  ))
}
# iCESM (months since 850)
icesm <- process_model("iCESM_d18O_850-1849.nc", 
                       time_units = "months since 850",
                       model_name = "iCESM")

# GISS-E2-R (years since 0)
giss <- process_model("GISS-E2-R_d18O_850-1849.nc", 
                      time_units = "month as %Y%m.%f", 
                      model_name = "GISS")

# iHadCM3 (months since 850 but uses different variable names)
ihad <- process_model("iHadCM3_d18O_850-1849.nc", 
                      lat_name = "latitude", 
                      lon_name = "longitude", 
                      time_name = "t", 
                      time_units = "months since 850",
                      model_name = "iHadCM3")




mapping <- function(X_annual_z, lon, lat, unique_years, start, end, title = NULL, return_data = FALSE) {
  target_years <- which(unique_years >= start & unique_years <= end)
  
  avg_slice <- apply(X_annual_z[, , target_years, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  
  df <- reshape2::melt(avg_slice)
  names(df) <- c("lon_idx", "lat_idx", "value")
  df$lon <- lon[df$lon_idx]
  df$lat <- lat[df$lat_idx]
  df <- df[!is.na(df$value), ]
  
  if (return_data) return(df)
  if (is.null(title)) title <- paste0(start, "–", end, " CE")
  
  colors <- c(
    "#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0",
    "#ffffbf", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"
  )
  
  ggplot() +
    geom_raster(data = df, aes(x = lon, y = lat, fill = value)) +
    geom_sf(data = world, fill = NA, color = "black", size = 0.7) +
    scale_fill_gradientn(
      colours = c(
        "#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0",
        "#ffffbf", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"
      ),
      limits = c(-0.1, 0.1),
      oob = scales::squish,
      name = "Mean (Z-score)",
      guide = guide_colorbar(
        title.position = "bottom",
        title.hjust = 0.5,
        barwidth = unit(5, "cm")  # 👈 adjust as needed
      )) +
    coord_sf(xlim = c(-120, 70), ylim = c(20, 90), expand = FALSE) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      #plot.title = element_text(hjust = 0.5, size = 16)
    ) +
    
    guides(fill = guide_colorbar(title.position = "bottom", title.hjust = 0.5, barwidth = 15))
}




######### SORTING THE PROXY DATA ##################
######### MCA #################


setwd("~/Documents/Documents – Laura’s MacBook Pro/Hydroclimate/")
load("iso2k1_0_1.RData")
# Remove extraneous objects
rm(D, TS)

shapes = c("GlacierIce" = 20, "LakeSediment" = 17, "Speleothem" = 15, "MolluskShells" = 15,
           "MarineSediment" = 17, "Wood" = 18, "TerrestrialSediement" = 17)

mca_start = 1100
mca_end = 1250 

world <- ne_countries(scale = "medium", returnclass = "sf")

# Use geoChronR to extract primary TS
all_iso_ts = sTS[which(pullTsVariable(sTS, 
                                      variable = "paleoData_iso2kPrimaryTimeseries") == "TRUE")]

# Extract Atlantic series
all_iso_ts_green = all_iso_ts[which(pullTsVariable(all_iso_ts, 
                                                   variable = "geo_longitude") > -90)]
all_iso_ts_green = all_iso_ts_green[which(pullTsVariable(all_iso_ts_green, 
                                                         variable = "geo_longitude") < 50)]
all_iso_ts_green = all_iso_ts_green[which(pullTsVariable(all_iso_ts_green, 
                                                         variable = "geo_latitude") > 30)]
all_iso_ts_green = all_iso_ts_green[which(pullTsVariable(all_iso_ts_green, 
                                                         variable = "geo_latitude") < 85)]

mcaRecs = matrix(NA, length(all_iso_ts_green), 12) %>%
  set_colnames(c("record", "resolution", "duration", 
                 "archive", "lat", "lon", "infMat", "var", "interp", 
                 "season", "start_year", "end_year")) # creating column names in an empty matrix

for (i in 1:(length(all_iso_ts_green))) {
  
  year_vals <- all_iso_ts_green[[i]]$year
  data_vals <- all_iso_ts_green[[i]]$paleoData_values
  
  # Make sure year and value vectors are the same length
  if (length(year_vals) != length(data_vals)) {
    message(paste("Skipping record", i, "due to mismatched year/value lengths"))
    next
  }
  
  thisYearVec <- na.omit(year_vals)
  
  # Create test record
  testrec <- na.omit(data.frame(year = year_vals, val = data_vals))
  
  # Skip if time series doesn't overlap MCA
  if (min(testrec$year) > mca_end) { next }
  if (length(which(testrec$year < mca_end & testrec$year > mca_start)) < 3) { next }
  
  # Fill out record metadata
  mcaRecs[i, 1] = all_iso_ts_green[[i]]$paleoData_TSid
  mcaRecs[i, 2] = (max(thisYearVec) - min(thisYearVec)) / length(thisYearVec)
  mcaRecs[i, 3] = max(thisYearVec) - min(thisYearVec)
  mcaRecs[i, 4] = all_iso_ts_green[[i]]$archiveType
  mcaRecs[i, 5] = all_iso_ts_green[[i]]$geo_latitude
  mcaRecs[i, 6] = all_iso_ts_green[[i]]$geo_longitude
  mcaRecs[i, 7] = all_iso_ts_green[[i]]$paleoData_inferredMaterial
  mcaRecs[i, 8] = all_iso_ts_green[[i]]$paleoData_variableName
  mcaRecs[i, 9] = all_iso_ts_green[[i]]$isotopeInterpretation1_variableGroup
  mcaRecs[i,10] = if (!is.null(all_iso_ts_green[[i]]$isotopeInterpretation1_seasonality)) {
    all_iso_ts_green[[i]]$isotopeInterpretation1_seasonality
  } else { NA }
  mcaRecs[i,11] = min(thisYearVec)
  mcaRecs[i,12] = max(thisYearVec)
}

# Format columns
mcaRecs = as.data.frame((mcaRecs)) %>%
  mutate_at(c("resolution", "lat", "lon", "duration", "start_year"), as.numeric)

# Extract TS
mcaTS = all_iso_ts_green[!is.na(mcaRecs$record)]

# Remove NA
mcaRecs = mcaRecs[-which(is.na(mcaRecs$record)),]

mca_means = vector()
for(i in 1:length(mcaTS)){
  testrec = na.omit(data.frame(year = mcaTS[[i]]$year, 
                               val = mcaTS[[i]]$paleoData_values))
  # Trim to after 900CE
  include = which(testrec$year > 0)
  testrec = testrec[include,]
  # z-score
  testrec$val = scale(testrec$val)
  # mca subset
  testrec = testrec[which(testrec$year < mca_end & testrec$year > mca_start),]
  print(length(testrec$year))
  mca_means = append(mca_means, mean(testrec$val))
}
mcaRecs = cbind(mcaRecs, mca_means)

# Remove marine sediments
mcaRecs <- mcaRecs[mcaRecs$archive != "MarineSediment", ]

MCA_plot <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey40", size = 0.3) +
  geom_point(data = mcaRecs, aes(x = lon, y = lat, shape = archive), 
             color = "black", size = 4) +
  geom_point(data = mcaRecs, aes(x = lon, y = lat, shape = archive, color = mca_means), 
             size = 3) +
  scale_shape_manual(values = shapes) +
  scale_color_gradient2(
    name = "Mean (Z-score)",
    low = "blue", mid = "white", high = "red",
    limits = c(-1.5, 1.5),
    breaks = c(-1.5, 0, 1.5),
    labels = c("-1.5", "0", "1.5")
  )+
  coord_sf(xlim = c(-120, 70), ylim = c(20, 90), expand = FALSE) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  labs(title = "MCA (1100 - 1250 CE)", shape = NULL) +
  guides(
    color = guide_colorbar(title.position = "bottom", title.hjust = 0.5),
    shape = guide_legend(override.aes = list(size = 4))
  ) +
  # Your rectangles
  geom_rect(aes(xmin = -15, xmax = 40, ymin = 47, ymax = 85), 
            color = "black", fill = NA, lwd = 0.6, , linetype = "dashed") +
  geom_rect(aes(xmin = -15, xmax = 40, ymin = 30, ymax = 47), 
            color = "black", fill = NA, lwd = 0.6, linetype = "dashed") +
  geom_rect(aes(xmin = -90, xmax = -15, ymin = 50, ymax = 85), 
            color = "black", fill = NA, lwd = 0.6,, linetype = "dashed")




# Get model data for MCA period (e.g., from iCESM)

# Define time periods
mca_start <- 1100; mca_end <- 1250
lia_start <- 1650; lia_end <- 1850

# MCA plots
make_combined_plot(icesm, "iCESM", mcaRecs, "mca_means", mca_start, mca_end)
make_combined_plot(giss,  "GISS",  mcaRecs, "mca_means", mca_start, mca_end)
make_combined_plot(ihad,  "iHadCM3", mcaRecs, "mca_means", mca_start, mca_end)

# LIA plots (assuming you have created liaRecs and lia_means similarly)
make_combined_plot(icesm, "iCESM", liaRecs, "lia_means", lia_start, lia_end)
make_combined_plot(giss,  "GISS",  liaRecs, "lia_means", lia_start, lia_end)
make_combined_plot(ihad,  "iHadCM3", liaRecs, "lia_means", lia_start, lia_end)
 