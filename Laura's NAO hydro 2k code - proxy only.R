##### April 2025 #######
# Laura Boyall

##### February 23, 2025 #######
# Matt Jones

#Adding LALIA and PI time periods and changing some code based on Andrew's initial analyses

####Based on code from

# October 4, 2024 #
# Andrew Flaim
# aflaim@wustl.edu
# 
# Iso2k MCA vs. LIA comparisons continued from PAGES Potsdam NAO meeting.
# Adjusting previous analysis to include shorter records in the LIA subset,
# as well as changing regional subsets based on modeled NAO-d18O fingerprints.

#
#### Libraries ####
library(magrittr)
library(tidyverse)
library(geoChronR)
library(lubridate)
library(progress)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

#
#### Extracting Atlantic Series ####
# Import the country and continent borders from the maptools package

setwd("~/Documents/Documents – Laura’s MacBook Pro/Hydroclimate/")
load("iso2k1_0_1.RData")
# Remove extraneous objects
rm(D, TS)

shapes = c("GlacierIce" = 20, "LakeSediment" = 17, "Speleothem" = 15, "MolluskShells" = 15,
           "MarineSediment" = 17, "Wood" = 18, "TerrestrialSediement" = 17)


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

# Define intervals of interest
lalia_start = 535
lalia_end = 660

mca_start = 1100
mca_end = 1250

lia_start = 1650
lia_end = 1850

pi_start = 1850
pi_end = 2000




###############################################################################
###############################################################################
######################    LATE ANTIQUE LITTLE ICE AGE   #######################
###############################################################################
###############################################################################

laliaRecs = matrix(NA, length(all_iso_ts_green), 12) %>%
  set_colnames(c("record", "resolution", "duration", 
                 "archive", "lat", "lon", "infMat", "var", "interp", 
                 "season", "start_year", "end_year"))
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
  if (min(testrec$year) > lalia_end) { next }
  if (length(which(testrec$year < lalia_end & testrec$year > lalia_start)) < 3) { next }
  
  # Fill out record metadata
  laliaRecs[i, 1] = all_iso_ts_green[[i]]$paleoData_TSid
  laliaRecs[i, 2] = (max(thisYearVec) - min(thisYearVec)) / length(thisYearVec)
  laliaRecs[i, 3] = max(thisYearVec) - min(thisYearVec)
  laliaRecs[i, 4] = all_iso_ts_green[[i]]$archiveType
  laliaRecs[i, 5] = all_iso_ts_green[[i]]$geo_latitude
  laliaRecs[i, 6] = all_iso_ts_green[[i]]$geo_longitude
  laliaRecs[i, 7] = all_iso_ts_green[[i]]$paleoData_inferredMaterial
  laliaRecs[i, 8] = all_iso_ts_green[[i]]$paleoData_variableName
  laliaRecs[i, 9] = all_iso_ts_green[[i]]$isotopeInterpretation1_variableGroup
  laliaRecs[i,10] = if (!is.null(all_iso_ts_green[[i]]$isotopeInterpretation1_seasonality)) {
    all_iso_ts_green[[i]]$isotopeInterpretation1_seasonality
  } else { NA }
  laliaRecs[i,11] = min(thisYearVec)
  laliaRecs[i,12] = max(thisYearVec)
}
# Format columns
laliaRecs = as.data.frame((laliaRecs)) %>%
  mutate_at(c("resolution", "lat", "lon", "duration", "start_year"), as.numeric)

# Extract TS
laliaTS = all_iso_ts_green[!is.na(laliaRecs$record)]

# Remove NA
laliaRecs = laliaRecs[-which(is.na(laliaRecs$record)),]

# LALIA means
lalia_means = vector()
for(i in 1:length(laliaTS)){
  testrec = na.omit(data.frame(year = laliaTS[[i]]$year, 
                               val = laliaTS[[i]]$paleoData_values))
  # Trim to after 0CE
  include = which(testrec$year > 0)
  testrec = testrec[include,]
  # z-score
  testrec$val = scale(testrec$val)
  # lalia subset
  testrec = testrec[which(testrec$year < lalia_end & testrec$year > lalia_start),]
  lalia_means = append(lalia_means, mean(testrec$val,na.rm=TRUE))
}
laliaRecs = cbind(laliaRecs, lalia_means)

# Remove marine sediments
laliaRecs <- laliaRecs[laliaRecs$archive != "MarineSediment", ]

# plotting the LALIA
LALIA_plot <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey40", size = 0.3) +
  geom_point(data = laliaRecs, aes(x = lon, y = lat, shape = archive), 
             color = "black", size = 4) +
  geom_point(data = laliaRecs, aes(x = lon, y = lat, shape = archive, color = lalia_means), 
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
  labs(title = "LALIA (535 - 660 CE)", shape = NULL) +
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




######## removing temperature and effective moisture as only 3 of each here
laliaRecs_filt <- laliaRecs[!grepl("Temperature\\|", laliaRecs$interp) &
                              laliaRecs$interp != "EffectiveMoisture", ]



filtered_LALIA <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey40", size = 0.3) +
  geom_point(data = laliaRecs_filt, aes(x = lon, y = lat, shape = archive), 
             color = "black", size = 4) +
  geom_point(data = laliaRecs_filt, aes(x = lon, y = lat, shape = archive, color = lalia_means), 
             size = 3,) +
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
  labs(title = "LALIA (535 - 660 CE) ", shape = NULL) +
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





###############################################################################
###############################################################################
######################    MEDIEVAL CLIMATE ANOMALY   ##########################
###############################################################################
###############################################################################

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




######## removing temperature and effective moisture as only 3 of each here
mca_filt <- mcaRecs[!grepl("Temperature\\|", mcaRecs$interp) &
                              mcaRecs$interp != "EffectiveMoisture", ]



filtered_MCA <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey40", size = 0.3) +
  geom_point(data = mca_filt, aes(x = lon, y = lat, shape = archive), 
             color = "black", size = 4) +
  geom_point(data = mca_filt, aes(x = lon, y = lat, shape = archive, color = mca_means), 
             size = 3,) +
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



###############################################################################
###############################################################################
##########################    LITTLE ICE AGE   ################################
###############################################################################
###############################################################################

liaRecs = matrix(NA, length(all_iso_ts_green), 12) %>%
  set_colnames(c("record", "resolution", "duration", 
                 "archive", "lat", "lon", "infMat", "var", "interp", 
                 "season", "start_year", "end_year"))
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
  if (min(testrec$year) > lia_end) { next }
  if (length(which(testrec$year < lia_end & testrec$year > lia_start)) < 3) { next }
  
  # Fill out record metadata
  liaRecs[i, 1] = all_iso_ts_green[[i]]$paleoData_TSid
  liaRecs[i, 2] = (max(thisYearVec) - min(thisYearVec)) / length(thisYearVec)
  liaRecs[i, 3] = max(thisYearVec) - min(thisYearVec)
  liaRecs[i, 4] = all_iso_ts_green[[i]]$archiveType
  liaRecs[i, 5] = all_iso_ts_green[[i]]$geo_latitude
  liaRecs[i, 6] = all_iso_ts_green[[i]]$geo_longitude
  liaRecs[i, 7] = all_iso_ts_green[[i]]$paleoData_inferredMaterial
  liaRecs[i, 8] = all_iso_ts_green[[i]]$paleoData_variableName
  liaRecs[i, 9] = all_iso_ts_green[[i]]$isotopeInterpretation1_variableGroup
  liaRecs[i,10] = if (!is.null(all_iso_ts_green[[i]]$isotopeInterpretation1_seasonality)) {
    all_iso_ts_green[[i]]$isotopeInterpretation1_seasonality
  } else { NA }
  liaRecs[i,11] = min(thisYearVec)
  liaRecs[i,12] = max(thisYearVec)
}
# Format columns
liaRecs = as.data.frame((liaRecs)) %>%
  mutate_at(c("resolution", "lat", "lon", "duration", "start_year"), as.numeric)

# Extract TS
liaTS = all_iso_ts_green[!is.na(liaRecs$record)]

# Remove NA
liaRecs = liaRecs[-which(is.na(liaRecs$record)),]

lia_means = vector()
for(i in 1:length(liaTS)){
  testrec = na.omit(data.frame(year = liaTS[[i]]$year, 
                               val = liaTS[[i]]$paleoData_values))
  # Trim to after 900CE
  include = which(testrec$year > 0)
  testrec = testrec[include,]
  # z-score
  testrec$val = scale(testrec$val)
  # lia subset
  testrec = testrec[which(testrec$year < lia_end & testrec$year > lia_start),]
  lia_means = append(lia_means, mean(testrec$val))
}
liaRecs = cbind(liaRecs, lia_means)

# Remove marine sediments and corals
liaRecs <- liaRecs[liaRecs$archive != "MarineSediment", ]
liaRecs <- liaRecs[liaRecs$archive != "Coral", ]



LIA_plot <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey40", size = 0.3) +
  geom_point(data = liaRecs, aes(x = lon, y = lat, shape = archive), 
             color = "black", size = 4) +
  geom_point(data = liaRecs, aes(x = lon, y = lat, shape = archive, color = lia_means), 
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
  labs(title = "LIA (1650 - 1850 CE)", shape = NULL) +
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




######## removing temperature and effective moisture as only 3 of each here
lia_filt <- liaRecs[!grepl("Temperature\\|", liaRecs$interp) &
                      liaRecs$interp != "EffectiveMoisture", ]



filtered_LIA <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey40", size = 0.3) +
  geom_point(data = lia_filt, aes(x = lon, y = lat, shape = archive), 
             color = "black", size = 4) +
  geom_point(data = lia_filt, aes(x = lon, y = lat, shape = archive, color = lia_means), 
             size = 3,) +
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
  labs(title = "LIA (1650 - 1850 CE)", shape = NULL) +
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


###############################################################################
###############################################################################
###########################    PRE-INDUSTRIAL   ###############################
###############################################################################
###############################################################################

piRecs = matrix(NA, length(all_iso_ts_green), 12) %>%
  set_colnames(c("record", "resolution", "duration", 
                 "archive", "lat", "lon", "infMat", "var", "interp", 
                 "season", "start_year", "end_year"))
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
  if (min(testrec$year) > pi_end) { next }
  if (length(which(testrec$year < pi_end & testrec$year > pi_start)) < 3) { next }
  
  # Fill out record metadata
  piRecs[i, 1] = all_iso_ts_green[[i]]$paleoData_TSid
  piRecs[i, 2] = (max(thisYearVec) - min(thisYearVec)) / length(thisYearVec)
  piRecs[i, 3] = max(thisYearVec) - min(thisYearVec)
  piRecs[i, 4] = all_iso_ts_green[[i]]$archiveType
  piRecs[i, 5] = all_iso_ts_green[[i]]$geo_latitude
  piRecs[i, 6] = all_iso_ts_green[[i]]$geo_longitude
  piRecs[i, 7] = all_iso_ts_green[[i]]$paleoData_inferredMaterial
  piRecs[i, 8] = all_iso_ts_green[[i]]$paleoData_variableName
  piRecs[i, 9] = all_iso_ts_green[[i]]$isotopeInterpretation1_variableGroup
  piRecs[i,10] = if (!is.null(all_iso_ts_green[[i]]$isotopeInterpretation1_seasonality)) {
    all_iso_ts_green[[i]]$isotopeInterpretation1_seasonality
  } else { NA }
  piRecs[i,11] = min(thisYearVec)
  piRecs[i,12] = max(thisYearVec)
}
# Format columns
piRecs = as.data.frame((piRecs)) %>%
  mutate_at(c("resolution", "lat", "lon", "duration", "start_year"), as.numeric)

# Extract TS
piTS = all_iso_ts_green[!is.na(piRecs$record)]

# Remove NA
piRecs = piRecs[-which(is.na(piRecs$record)),]

pi_means = vector()
for(i in 1:length(piTS)){
  testrec = na.omit(data.frame(year = piTS[[i]]$year, 
                               val = piTS[[i]]$paleoData_values))
  # Trim to after 900CE
  include = which(testrec$year > 0)
  testrec = testrec[include,]
  # z-score
  testrec$val = scale(testrec$val)
  # lia subset
  testrec = testrec[which(testrec$year < pi_end & testrec$year > pi_start),]
  pi_means = append(pi_means, mean(testrec$val))
}
piRecs = cbind(piRecs, pi_means)

# Remove marine sediments and coral
piRecs <- piRecs[piRecs$archive != "MarineSediment", ]
piRecs <- piRecs[piRecs$archive != "Coral", ]


PI_plot <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey40", size = 0.3) +
  geom_point(data = piRecs, aes(x = lon, y = lat, shape = archive), 
             color = "black", size = 4) +
  geom_point(data = piRecs, aes(x = lon, y = lat, shape = archive, color = pi_means), 
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
  labs(title = "PI (1850 - 2000 CE)", shape = NULL) +
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




######## removing temperature and effective moisture as only 3 of each here
pi_filt <- piRecs[!grepl("Temperature\\|", piRecs$interp) &
                      piRecs$interp != "EffectiveMoisture", ]

filtered_pI <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey40", size = 0.3) +
  geom_point(data = pi_filt, aes(x = lon, y = lat, shape = archive), 
             color = "black", size = 4) +
  geom_point(data = pi_filt, aes(x = lon, y = lat, shape = archive, color = pi_means), 
             size = 3,) +
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
  labs(title = "PI ( 1850 - 2000 CE)", shape = NULL) +
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



###############################################################################
###############################################################################
########################    PLOTTING TOGETHER   ###############################
###############################################################################
###############################################################################
add_borders_and_spacing <- function(p) {
  p + 
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.margin = margin(10, 10, 10, 10)  # top, right, bottom, left
    )
}
LALIA_plot <- add_borders_and_spacing(LALIA_plot)
MCA_plot <- add_borders_and_spacing(MCA_plot)
LIA_plot <- add_borders_and_spacing(LIA_plot)
PI_plot <- add_borders_and_spacing(PI_plot)

plots <- (LALIA_plot| MCA_plot) / (LIA_plot | PI_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")




filtered_LALIA <- add_borders_and_spacing(filtered_LALIA)
filtered_MCA <- add_borders_and_spacing(filtered_MCA)
filtered_LIA <- add_borders_and_spacing(filtered_LIA)
filtered_pI <- add_borders_and_spacing(filtered_pI)

p_isotope_plots <- (filtered_LALIA | filtered_MCA) / (filtered_LIA | filtered_pI) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


###############################################################################
###############################################################################
############################    BOX PLOTS   ###################################
###############################################################################
###############################################################################

#### Regional boxplots of regional mean values ####
# Define individual regions:
# Greenland
lon_min = -90
lon_max = -15
lat_min = 50
lat_max = 85
# Filter Greenland records
mcaRecs_green = mcaRecs[which((mcaRecs$lat > lat_min) & 
                                (mcaRecs$lat < lat_max) & 
                                (mcaRecs$lon > lon_min) & 
                                (mcaRecs$lon < lon_max)),]
liaRecs_green = liaRecs[which((liaRecs$lat > lat_min) & 
                                (liaRecs$lat < lat_max) & 
                                (liaRecs$lon > lon_min) & 
                                (liaRecs$lon < lon_max)),]
laliaRecs_green = laliaRecs[which((laliaRecs$lat > lat_min) & 
                                    (laliaRecs$lat < lat_max) & 
                                    (laliaRecs$lon > lon_min) & 
                                    (laliaRecs$lon < lon_max)),]
piRecs_green = piRecs[which((piRecs$lat > lat_min) & 
                                    (piRecs$lat < lat_max) & 
                                    (piRecs$lon > lon_min) & 
                                    (piRecs$lon < lon_max)),]

# N.EU
lon_min = -15
lon_max = 40
lat_min = 50
lat_max = 85
# Filter N.EU records
mcaRecs_NEU = mcaRecs[which((mcaRecs$lat > lat_min) & 
                              (mcaRecs$lat < lat_max) & 
                              (mcaRecs$lon > lon_min) & 
                              (mcaRecs$lon < lon_max)),]
liaRecs_NEU = liaRecs[which((liaRecs$lat > lat_min) & 
                              (liaRecs$lat < lat_max) & 
                              (liaRecs$lon > lon_min) & 
                              (liaRecs$lon < lon_max)),]
laliaRecs_NEU = laliaRecs[which((laliaRecs$lat > lat_min) & 
                                  (laliaRecs$lat < lat_max) & 
                                  (laliaRecs$lon > lon_min) & 
                                  (laliaRecs$lon < lon_max)),]
piRecs_NEU = piRecs[which((piRecs$lat > lat_min) & 
                                  (piRecs$lat < lat_max) & 
                                  (piRecs$lon > lon_min) & 
                                  (piRecs$lon < lon_max)),]
# S. EU
lon_min = -15
lon_max = 40
lat_min = 30
lat_max = 50

# Filter N.EU records
mcaRecs_SEU = mcaRecs[which((mcaRecs$lat > lat_min) & 
                              (mcaRecs$lat < lat_max) & 
                              (mcaRecs$lon > lon_min) & 
                              (mcaRecs$lon < lon_max)),]
liaRecs_SEU = liaRecs[which((liaRecs$lat > lat_min) & 
                              (liaRecs$lat < lat_max) & 
                              (liaRecs$lon > lon_min) & 
                              (liaRecs$lon < lon_max)),]
laliaRecs_SEU = laliaRecs[which((laliaRecs$lat > lat_min) & 
                                  (laliaRecs$lat < lat_max) & 
                                  (laliaRecs$lon > lon_min) & 
                                  (laliaRecs$lon < lon_max)),]
piRecs_SEU = piRecs[which((piRecs$lat > lat_min) & 
                                  (piRecs$lat < lat_max) & 
                                  (piRecs$lon > lon_min) & 
                                  (piRecs$lon < lon_max)),]

#### Create boxplots for MCA and LIA and Hist ###
box.list = list("Grn LALIA" = laliaRecs_green$lalia_means,
                "Grn MCA" = mcaRecs_green$mca_means,
                "Grn LIA" = liaRecs_green$lia_means,
                "Grn PI" = piRecs_green$pi_means,
                "N.EU LALIA" = laliaRecs_NEU$lalia_means,
                "N.EU MCA" = mcaRecs_NEU$mca_means, 
                "N.EU LIA" = liaRecs_NEU$lia_means,
                "N.EU PI" = piRecs_NEU$pi_means,
                "S.EU LALIA" = laliaRecs_SEU$lalia_means,
                "S.EU MCA" = mcaRecs_SEU$mca_means,
                "S.EU LIA" = liaRecs_SEU$lia_means,
                "S.EU PI" = piRecs_SEU$pi_means)

boxplot(box.list, outline = FALSE, ylim = c(-2.0, 2.0), las = 2)  # las = 2 makes labels vertical

stripchart(box.list,
           method = "jitter", 
           col = c("black"), 
           pch = 16,
           vertical = T,
           add = T)

#### coloured

region_colors <- c(
  "Grn" = "#4C72B0",   
  "N.EU" = "#DD8452",   
  "S.EU" = "#55A868"    
)

box.names <- names(box.list)

# Extract region tag from each name (before first space)
box.regions <- sapply(strsplit(box.names, " "), function(x) x[1])

# Map region colors to box.list order
box.cols <- region_colors[box.regions]

# Plot boxplot with custom colors
boxplot(box.list,
        col = box.cols,
        outline = FALSE,
        ylim = c(-2.0, 2.0),
        las = 2,
        ylab="Z-Score")

# Add jittered points
stripchart(box.list,
           method = "jitter", 
           col = "black", 
           pch = 16,
           vertical = TRUE,
           add = TRUE)







# Extract model data for MCA period from iCESM
model_bg_df <- mapping(
  X_annual_z = icesm$X_annual_z,
  lon = icesm$lon,
  lat = icesm$lat,
  unique_years = icesm$unique_years,
  start = 1100,
  end = 1250,
  return_data = TRUE
)

# Plot combining model background and proxy data
MCA_combined_plot <- ggplot() +
  geom_raster(data = model_bg_df, aes(x = lon, y = lat, fill = value)) +
  geom_point(data = mcaRecs, aes(x = lon, y = lat, shape = archive), 
             color = "black", size = 4) +
  geom_point(data = mcaRecs, aes(x = lon, y = lat, shape = archive, color = mca_means), 
             size = 3) +
  geom_sf(data = world, fill = NA, color = "black", size = 0.3) +
  scale_shape_manual(values = shapes) +
  scale_fill_gradientn(
    name = "Model Mean (Z-score)",
    colours = c(
      "#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0",
      "#ffffbf", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"
    ),
    limits = c(-0.1, 0.1),
    oob = scales::squish
  ) +
  scale_color_gradient2(
    name = "Proxy Mean (Z-score)",
    low = "blue", mid = "white", high = "red",
    limits = c(-1.5, 1.5),
    breaks = c(-1.5, 0, 1.5),
    labels = c("-1.5", "0", "1.5")
  ) +
  coord_sf(xlim = c(-120, 70), ylim = c(20, 90), expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  labs(title = "MCA (1100–1250 CE)", shape = NULL) +
  guides(
    fill = guide_colorbar(title.position = "bottom", title.hjust = 0.5, barwidth = unit(5, "cm")),
    color = guide_colorbar(title.position = "bottom", title.hjust = 0.5, barwidth = unit(5, "cm")),
    shape = guide_legend(override.aes = list(size = 4))
  ) +
  geom_rect(aes(xmin = -15, xmax = 40, ymin = 47, ymax = 85), 
            color = "black", fill = NA, lwd = 0.6, linetype = "dashed") +
  geom_rect(aes(xmin = -15, xmax = 40, ymin = 30, ymax = 47), 
            color = "black", fill = NA, lwd = 0.6, linetype = "dashed") +
  geom_rect(aes(xmin = -90, xmax = -15, ymin = 50, ymax = 85), 
            color = "black", fill = NA, lwd = 0.6, linetype = "dashed")

# Show it
MCA_combined_plot



















