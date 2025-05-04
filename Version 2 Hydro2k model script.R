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

setwd("/Volumes/LaCie/PAGES/Hydro2k/")
iso <- nc_open("iCESM_d18O_850-1849.nc")
ncatt_get(iso, "time", "units")

TS <- ncvar_get(iso, "d18O") # 
lat <- ncvar_get(iso, "lat")
lon <- ncvar_get(iso, "lon")
# normalising / Z-score the data
TS_z <- TS  # copy to preserve original

for (i in 1:dim(TS)[1]) {     #looping over longitu
  for (j in 1:dim(TS)[2]) {   # looping over lat
    x <- TS[i, j, ]
    if (all(is.na(x))) next  # skip missing data
    TS_z[i, j, ] <- scale(x)[, 1]
  }
}


Time <- ncvar_get(iso, "time")
years <- floor(850 + Time / 12)
unique_years <- unique(years)

n_years <- length(unique_years)
X_annual_z <- array(NA, dim = c(dim(TS)[1], dim(TS)[2], n_years))  # [lon, lat, year]

for (i in seq_along(unique_years)) {
  idx <- which(years == unique_years[i])
  X_annual_z[, , i] <- apply(TS_z[, , idx, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
}


# A function to map the different time staps using NORMALISED data
world <- ne_countries(scale = "medium", returnclass = "sf")

mapping <- function(start, end, title = NULL) {
  target_years <- which(unique_years >= start & unique_years <= end)
  
  # Average across the year dimension
  avg_slice_z <- apply(X_annual_z[, , target_years, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  
  # Convert to dataframe
  df <- reshape2::melt(avg_slice_z)
  names(df) <- c("lon_idx", "lat_idx", "value")
  df$lon <- lon[df$lon_idx]
  df$lat <- lat[df$lat_idx]
  df$lon <- ifelse(df$lon > 180, df$lon - 360, df$lon)
  
  # Remove NAs for better rendering
  df <- df[!is.na(df$value), ]
  
  if (is.null(title)) {
    title <- paste0(start, "–", end, " CE")
  }
  
  colors <- c(
    "#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0",
    "#ffffbf", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"
  )
  
  ggplot() +
    geom_tile(data = df, aes(x = lon, y = lat, fill = value)) +
    geom_sf(data = world, fill = NA, color = "black", size = 0.7) +  # Make sure `world` is loaded
    scale_fill_gradientn(
      colours = colors,
      limits = c(-.1, .1),  # You can tighten this if values are close to 0
      oob = scales::squish,
      name = "Mean (Z-score)"
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
    labs(title = title) +
    guides(fill = guide_colorbar(title.position = "bottom", title.hjust = 0.5)) 
    
}









mapping <- function(start, end, title = NULL, return_data = FALSE) {
  target_years <- which(unique_years >= start & unique_years <= end)
  avg_slice_z <- apply(X_annual_z[, , target_years, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  
  df <- reshape2::melt(avg_slice_z)
  names(df) <- c("lon_idx", "lat_idx", "value")
  df$lon <- lon[df$lon_idx]
  df$lat <- lat[df$lat_idx]
  df$lon <- ifelse(df$lon > 180, df$lon - 360, df$lon)
  df <- df[!is.na(df$value), ]
  
  if (is.null(title)) title <- paste0(start, "–", end, " CE")
  
  if (return_data) return(df)
  
  colors <- c(
    "#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0",
    "#ffffbf", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"
  )
  
  ggplot() +
    geom_raster(data = df, aes(x = lon, y = lat, fill = value)) +
    geom_sf(data = world, fill = NA, color = "black", size = 0.7) +
    scale_fill_gradientn(
      colours = colors,
      limits = c(-0.1, 0.1),
      oob = scales::squish,
      name = "Mean (Z-score)"
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
    labs(title = title) +
    guides(
      fill = guide_colorbar(
        title.position = "bottom",
        title.hjust = 0.5,
        barwidth = 15  # Stretch width
      )
    )}






# Plot MCA and LIA
plot_mca <- mapping(1100, 1250, title = "MCA")
plot_lia <- mapping(1650, 1850, title = "LIA")

# Get data for difference plot
mca_df <- mapping(1100, 1250, return_data = TRUE)
lia_df <- mapping(1650, 1850, return_data = TRUE)

# Join and compute difference
diff_df <- merge(mca_df, lia_df, by = c("lon_idx", "lat_idx"), suffixes = c("_mca", "_lia"))
diff_df$value <- diff_df$value_mca - diff_df$value_lia
diff_df$lon <- diff_df$lon_mca
diff_df$lat <- diff_df$lat_mca
diff_df <- diff_df[, c("lon", "lat", "value")]

# Custom color palette for difference
library(scales)
library(patchwork)
diff_colors <- colorRampPalette(c("deeppink3", "white", "darkgreen"))(11)

# Plot the difference
plot_diff <- ggplot() +
  geom_tile(data = diff_df, aes(x = lon, y = lat, fill = value)) +
  geom_sf(data = world, fill = NA, color = "black", size = 0.7) +
  scale_fill_gradientn(
    colours = diff_colors,
    limits = c(-0.1, 0.1),
    oob = scales::squish,
    name = "MCA - LIA"
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
  labs(title = "MCA – LIA") +
  guides(
    fill = guide_colorbar(
      title.position = "bottom",
      title.hjust = 0.5,
      barwidth = 15  # Stretch width
    )
  )

(plot_mca | plot_lia | plot_diff) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")



####### Box plot
# Define regions again (short names for labels)
regions <- list(
  Grn = list(lon_min = -90, lon_max = -15, lat_min = 50, lat_max = 85),
  N.EU = list(lon_min = -15, lon_max = 40, lat_min = 50, lat_max = 85),
  S.EU = list(lon_min = -15, lon_max = 40, lat_min = 30, lat_max = 50)
)

years_mca <- which(unique_years >= 1100 & unique_years <= 1250)
years_lia <- which(unique_years >= 1650 & unique_years <= 1850)

lon_adj <- ifelse(lon > 180, lon - 360, lon)
box.list <- list()

for (reg in names(regions)) {
  r <- regions[[reg]]
  lon_idx <- which(lon_adj >= r$lon_min & lon_adj <= r$lon_max)
  lat_idx <- which(lat >= r$lat_min & lat <= r$lat_max)
  
  # MCA
  mca_vals <- as.vector(X_annual_z[lon_idx, lat_idx, years_mca])
  mca_vals <- mca_vals[!is.na(mca_vals)]
  box.list[[paste(reg, "MCA")]] <- mca_vals
  
  # LIA
  lia_vals <- as.vector(X_annual_z[lon_idx, lat_idx, years_lia])
  lia_vals <- lia_vals[!is.na(lia_vals)]
  box.list[[paste(reg, "LIA")]] <- lia_vals
}

region_colors <- c(
  "Grn" = "#4C72B0",   # blue
  "N.EU" = "#DD8452",  # orange
  "S.EU" = "#55A868"   # green
)

box.names <- names(box.list)

# Extract region tag from each name (e.g. "Grn" from "Grn MCA")
box.regions <- sapply(strsplit(box.names, " "), function(x) x[1])
box.cols <- region_colors[box.regions]
# Base R boxplot with custom colors
boxplot(box.list,
        col = box.cols,
        outline = FALSE,
        ylim = c(-0.75, .75),
        las = 2,  # vertical x-axis labels
        ylab = "Z-Score",
        cex.axis = 0.9,
        main = "MCA vs LIA across regions")















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



# Updated difference plot function without title
make_diff_plot <- function(mca_df, lia_df) {
  df_diff <- merge(mca_df, lia_df, by = c("lon_idx", "lat_idx"))
  df_diff$value <- df_diff$value.x - df_diff$value.y
  df_diff$lon <- df_diff$lon.x
  df_diff$lat <- df_diff$lat.x
  df_diff <- df_diff[, c("lon", "lat", "value")]
  
  ggplot() +
    geom_tile(data = df_diff, aes(x = lon, y = lat, fill = value)) +
    geom_sf(data = world, fill = NA, color = "black", size = 0.5) +
    scale_fill_gradientn(
      colours = colorRampPalette(c("purple", "white", "darkgreen"))(11),
      limits = c(-0.1, 0.1),
      oob = scales::squish,
      name = "MCA − LIA",
      guide = guide_colorbar(
        title.position = "bottom",
        title.hjust = 0.5,
        barwidth = unit(5, "cm")  # 👈 same as above
      )
    ) +
    coord_sf(xlim = c(-120, 70), ylim = c(20, 90), expand = FALSE) +
    theme_void()
}

# iCESM
p1_icesm <- mapping(icesm$X_annual_z, icesm$lon, icesm$lat, icesm$unique_years, 1100, 1250, title = NULL)
p2_icesm <- mapping(icesm$X_annual_z, icesm$lon, icesm$lat, icesm$unique_years, 1650, 1850, title = NULL)
df_mca_icesm <- mapping(icesm$X_annual_z, icesm$lon, icesm$lat, icesm$unique_years, 1100, 1250, return_data = TRUE)
df_lia_icesm <- mapping(icesm$X_annual_z, icesm$lon, icesm$lat, icesm$unique_years, 1650, 1850, return_data = TRUE)
p3_icesm <- make_diff_plot(df_mca_icesm, df_lia_icesm)

# GISS
p1_giss <- mapping(giss$X_annual_z, giss$lon, giss$lat, giss$unique_years, 1100, 1250, title = NULL)
p2_giss <- mapping(giss$X_annual_z, giss$lon, giss$lat, giss$unique_years, 1650, 1850, title = NULL)
df_mca_giss <- mapping(giss$X_annual_z, giss$lon, giss$lat, giss$unique_years, 1100, 1250, return_data = TRUE)
df_lia_giss <- mapping(giss$X_annual_z, giss$lon, giss$lat, giss$unique_years, 1650, 1850, return_data = TRUE)
p3_giss <- make_diff_plot(df_mca_giss, df_lia_giss)

# iHadCM3
p1_ihad <- mapping(ihad$X_annual_z, ihad$lon, ihad$lat, ihad$unique_years, 1100, 1250, title = NULL)
p2_ihad <- mapping(ihad$X_annual_z, ihad$lon, ihad$lat, ihad$unique_years, 1650, 1850, title = NULL)
df_mca_ihad <- mapping(ihad$X_annual_z, ihad$lon, ihad$lat, ihad$unique_years, 1100, 1250, return_data = TRUE)
df_lia_ihad <- mapping(ihad$X_annual_z, ihad$lon, ihad$lat, ihad$unique_years, 1650, 1850, return_data = TRUE)
p3_ihad <- make_diff_plot(df_mca_ihad, df_lia_ihad)

# Patch rows
row1 <- p1_icesm + p2_icesm + p3_icesm
row2 <- p1_giss  + p2_giss  + p3_giss
row3 <- p1_ihad  + p2_ihad  + p3_ihad

# Combine with shared legend
final_panel <- (row1 / row2 / row3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Add letter labels A–I to actual plots
final_panel <- final_panel + plot_annotation(tag_levels = "A")

# Show final panel
final_panel





######## box plots: 


boxplot_model <- function(X, lon, lat, unique_years, model_name) {
  # Define your regions
  regions <- list(
    Greenland = list(lon_min = -90, lon_max = -15, lat_min = 50, lat_max = 85),
    NorthernEurope = list(lon_min = -15, lon_max = 40, lat_min = 50, lat_max = 85),
    SouthernEurope = list(lon_min = -15, lon_max = 40, lat_min = 30, lat_max = 50)
  )
  
  # Time slices
  years_mca <- which(unique_years >= 1100 & unique_years <= 1250)
  years_lia <- which(unique_years >= 1650 & unique_years <= 1850)
  
  box.list <- list()
  
  for (reg in names(regions)) {
    r <- regions[[reg]]
    lon_idx <- which(lon >= r$lon_min & lon <= r$lon_max)
    lat_idx <- which(lat >= r$lat_min & lat <= r$lat_max)
    
    # MCA
    mca_vals <- as.vector(X[lon_idx, lat_idx, years_mca])
    mca_vals <- mca_vals[!is.na(mca_vals)]
    box.list[[paste(reg, "MCA", model_name)]] <- mca_vals
    
    # LIA
    lia_vals <- as.vector(X[lon_idx, lat_idx, years_lia])
    lia_vals <- lia_vals[!is.na(lia_vals)]
    box.list[[paste(reg, "LIA", model_name)]] <- lia_vals
  }
  
  return(box.list)
}




box_icesm <- boxplot_model(icesm$X_annual_z, icesm$lon, icesm$lat, icesm$unique_years, "iCESM")
box_giss <- boxplot_model(giss$X_annual_z, giss$lon, giss$lat, giss$unique_years, "GISS")
box_ihad <- boxplot_model(ihad$X_annual_z, ihad$lon, ihad$lat, ihad$unique_years, "iHadCM3")

box.all <- c(box_icesm, box_giss, box_ihad)


# Set custom region colors
region_colors <- c(
  "Greenland" = "#4C72B0",   # blue
  "NorthernEurope" = "#DD8452",  # orange
  "SouthernEurope" = "#55A868"   # green
)

# Get box names
box.names <- names(box.all)

# Extract region tag from each name (e.g. "Grn" from "Grn MCA iCESM")
box.regions <- sapply(strsplit(box.names, " "), function(x) x[1])
box.cols <- region_colors[box.regions]


boxplot(box.all,
        col = box.cols,
        outline = FALSE,
        las = 2,  # vertical axis labels
        ylab = "Z-score",
        cex.axis = 0.8,
        main = "MCA vs LIA by Region and Model")



###### grouping by regions; 

df_box <- do.call(rbind, lapply(names(box.all), function(nm) {
  vals <- box.all[[nm]]
  parts <- unlist(strsplit(nm, " "))  # e.g. "Grn MCA iCESM"
  data.frame(
    value = vals,
    Region = parts[1],
    Period = parts[2],
    Model = parts[3],
    stringsAsFactors = FALSE
  )
}))

# Factor levels for consistent ordering
df_box$Region <- factor(df_box$Region, levels = c("Greenland", "NorthernEurope", "SouthernEurope"))
df_box$Period <- factor(df_box$Period, levels = c("MCA", "LIA"))
df_box$Model <- factor(df_box$Model, levels = c("iCESM", "GISS", "iHadCM3"))

# Create full group labels
df_box$PlotGroup <- interaction(df_box$Region, df_box$Period, df_box$Model, sep = "_")

# Set factor order manually for grouped plotting
df_box$PlotGroup <- factor(df_box$PlotGroup,
                           levels = with(df_box, unique(PlotGroup[order(Region, Period, Model)])))



# Step 1: Get the correct group levels
group_levels <- levels(df_box$PlotGroup)

# Step 2: Extract region name from each level
box.regions <- sapply(strsplit(group_levels, "_"), function(x) x[1])

# Step 3: Check if any region names are misspelled
if (any(!box.regions %in% names(region_colors))) {
  warning("Some regions in plot group levels do not match your region_colors keys.")
}

# Step 4: Assign region-based colors per box (ordered!)
box.cols <- unname(region_colors[box.regions])
boxplot(value ~ PlotGroup,
        data = df_box,
        col = box.cols,  # <- this now works correctly
        outline = FALSE,
        las = 2,
        ylab = "Z-score",
        xlab = "",
        main = "MCA vs LIA by Region and Model",
        cex.axis = 0.75)



#### plotting just MCA


# Modify PlotGroup for custom x-axis labels (without Region)
df_box$PlotGroup_NoRegion <- interaction(df_box$Period, df_box$Model, sep = "_")

# Filter the data for MCA (MCA data only)
df_box_mca <- df_box[df_box$Period == "MCA", ]

# Filter the data for LIA (LIA data only)
df_box_lia <- df_box[df_box$Period == "LIA", ]

# Plot for MCA only
boxplot(value ~ PlotGroup,
        data = df_box_mca,
        col = box.cols,  # Use the same color scheme
        outline = FALSE,
        las = 2,  # vertical labels
        ylab = "Z-score",
        xlab = "",
        main = "MCA by Region and Model",
        cex.axis = 0.75)

# Custom x-axis labels for MCA (only Model and Period, no Region)
# Use the levels of PlotGroup_NoRegion for labels (MCA, LIA, etc.)
axis(1, at = 1:length(unique(df_box_mca$PlotGroup)), 
     labels = levels(df_box_mca$PlotGroup_NoRegion), las = 2)

# Plot for LIA only
boxplot(value ~ PlotGroup,
        data = df_box_lia,
        col = box.cols,  # Use the same color scheme
        outline = FALSE,
        las = 2,  # vertical labels
        ylab = "Z-score",
        xlab = "",
        main = "LIA by Region and Model",
        cex.axis = 0.75)

# Custom x-axis labels for LIA (only Model and Period, no Region)
# Use the levels of PlotGroup_NoRegion for labels (MCA, LIA, etc.)
axis(1, at = 1:length(unique(df_box_lia$PlotGroup)), 
     labels = levels(df_box_lia$PlotGroup_NoRegion), las = 2)




####### Now I want to plot the proxy data onto the maps 
# Get model data for MCA period (e.g., from iCESM)
model_bg_df <- mapping(icesm$X_annual_z, 
                       icesm$lon, 
                       icesm$lat, 
                       icesm$unique_years, 
                       1100, 1250, 
                       return_data = TRUE)






