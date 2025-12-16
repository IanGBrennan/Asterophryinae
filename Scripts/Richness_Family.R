# let's try: https://github.com/RS-eco/rasterSp

if(!"rasterSp" %in% installed.packages()[,"Package"]) remotes::install_github("RS-eco/rasterSp", build_vignettes = T)

library(rasterSp)

filedir <- "/Users/ianbrennan/Downloads/ANURA/"

# Convert shape files into rasters and save to file
rasterizeRange(dsn=paste0(filedir, "MICROHYLIDAE.shp"),
               id="sci_name",
               resolution=0.25, save=TRUE, touches=T,
               seasonal=c(1,2), origin=1, presence=c(1,2), 
               path=paste0(filedir, "/SpeciesData/"))

# Calculate the species richness from the rasters we generated in the previous step
sr_micros <- calcSR(species_names=microhylidae$sci_name, path=paste0(filedir, "/SpeciesData/"))
raster::plot(sr_micros)

# set the map color scale
map.cols <- rev(viridis::inferno(n=50))
ggmap2(sr_micros, name="richness", split=TRUE, ncol=1, country=T, colours=map.cols)

# Load countries and coastline sf data
countries <- ne_countries(scale = "small")
coast <- ne_coastline(scale = "small")

# Plot the continental outline if you'd like
ggplot() +
  geom_sf(data = countries, colour = "lightgrey", fill = "lightgrey") +
  geom_sf(data = coast, colour = "grey50") +
  theme_void()
