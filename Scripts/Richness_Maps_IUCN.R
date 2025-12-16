# packages
library(tidyverse)  # Data wrangling
#library(here)       # Safe paths
library(sf)         # Spatial features

################################################################################

# download IUCN species range maps for anuran amphibians: https://www.iucnredlist.org/resources/spatial-data-download
frog1 <- read_sf('/Users/ianbrennan/Downloads/ANURA/ANURA_PART1.shp')
frog2 <- read_sf('/Users/ianbrennan/Downloads/ANURA/ANURA_PART2.shp')

microhylidae <- rbind(dplyr::filter(frog1, family == "MICROHYLIDAE"),
                      dplyr::filter(frog2, family == "MICROHYLIDAE"))

# identify the genera of interest (all Asterophryinae)
astero.genera <- c("Aphantophryne","Asterophrys","Austrochaperina","Barygenys",
                   "Callulops","Choerophryne","Cophixalus","Copiula","Gastrophrynoides",
                   "Hylophorbus","Mantophryne","Oninia","Oreophryne","Paedophryne",
                   "Siamophryne","Sphenophryne","Vietnamophryne","Xenorhina",
                   "Glyphoglossus")

# downsample to just the genera we're interested in
asters <- dplyr::filter(microhylidae, genus %in% astero.genera)

# we can inspect this by plotting a single taxon
ggplot() +
  geom_sf(data=dplyr::filter(asters, sci_name == "Cophixalus saxatilis"), fill="blue")

# load the world maps
library(rnaturalearth)
worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                         returnclass = 'sf')
ggplot() + geom_sf(data = worldmap) + theme_classic() + coord_sf(crs = "+proj=robin")

# subset to the countries we're interested in
oz.ng <- worldmap[worldmap$name %in% c('Australia', 'Papua New Guinea', 'Indonesia', 'Solomon Is.'),]
# plot this to inspect
ggplot() + geom_sf(data = oz.ng) + theme_bw()

# trim down our region of interest
region = c(xmin = 128, ymin = -25, xmax = 155, ymax = 0)
sahul <- st_crop(oz.ng, st_bbox(region))
ggplot() + geom_sf(data = sahul)

# create grid
sahul_grid <- st_make_grid(sahul,
                        what = "polygons",
                        cellsize = 0.25,
                        square = FALSE,
                        flat_topped = TRUE)

# subset to grid cells that are within land
keep_hexes <- st_intersects(sahul_grid, sahul)
keep_hexes <- as.data.frame(keep_hexes)$row.id
sahul_grid <- sahul_grid[keep_hexes]

# check that the density of grid cells is appropriate
ggplot() +
  geom_sf(data = sahul_grid) +
  geom_sf(data = sahul, fill = NA)


# load the rangeBuilder package
library(rangeBuilder)

# move the polygons out of a spatialPolygonsDataframe and into a list
poly.list <- split(asters, seq(nrow(asters)))
names(poly.list) <- asters$sci_name

# test by plotting a single species atop our grid
ggplot() + 
  geom_sf(data = dplyr::filter(asters, sci_name == "Austrochaperina gracilipes"), fill="blue") +
  geom_sf(data = sahul_grid, fill=NA) +
  #geom_sf(data = sahul_grid[1025], fill="red") +
  geom_sf(data = sahul_grid[1131], fill="pink") +
  theme_classic()


# now the big step, run a loop that determines intersection
# between each taxon and each hex cell, this can take a while
richness <- NULL; rich.coph <- NULL; rich.aust <- NULL
#for(k in 1:length(sahul_grid)){
for(k in 1:length(sahul_grid)){
  print(paste("assessing hex:",k))
  curr.hex <- sahul_grid[k]
  
#  res <- unlist(lapply(poly.list, function(x) st_intersects(x, curr.hex, sparse=F)))
  res <- unlist(parallel::mclapply(poly.list, function(x) st_intersects(x, curr.hex, sparse=F), mc.cores=8))
  spp <- sum(res)
  richness <- append(richness, spp)
  
  res.coph <- res[grep("Cophixalus",names(res))]
  spp.coph <- sum(res.coph)
  rich.coph <- append(rich.coph, spp.coph)
  
  res.aust <- res[grep("Austrochaperina",names(res))]
  spp.aust <- sum(res.aust)
  rich.aust <- append(rich.aust, spp.aust)

}

# save the data as a dataframe of the richness information per grid cell
species.richness <- data.frame(Asterophryinae = richness,
                               Cophixalus = rich.coph,
                               Austrochaperina = rich.aust)
save(species.richness, file="Data/Richness_Maps_IUCN.RData")

# make an empty frame for the asteroprhyine richness
micros <- st_as_sf(sahul_grid)
micros$count <- richness

# plot it all together
micro_map <- ggplot() +
  geom_sf(data = micros, mapping = aes(fill = count), alpha = 1, color = NA) +
  #scale_fill_distiller(direction = 1, palette = "Greens", limits=c(0,32)) +
  scale_fill_viridis_c(direction=-1, option="rocket") +
  # add map
  geom_sf(data = sahul, color = "black", fill = NA) +
  # crop map
  coord_sf(xlim = c(129, 155), ylim = c(-20, 0)) +
  theme_classic() + theme(legend.position = "none")


# map the Cophixalus richness to the grid
cophs <- st_as_sf(sahul_grid)
cophs$count <- rich.coph

coph_map <- ggplot() +
  geom_sf(data = cophs, mapping = aes(fill = count), alpha = 1, color = NA) +
  scale_fill_distiller(direction = 1, palette = "Blues", limits=c(0,7)) +
  #scale_fill_viridis_c(direction=-1, option="rocket") +
  # add map
  geom_sf(data = sahul, color = "black", fill = NA) +
  # crop map
  coord_sf(xlim = c(129, 155), ylim = c(-20, 0)) +
  theme_classic() + theme(legend.position = "none")

# map the Austrochaperina richness to the grid
austs <- st_as_sf(sahul_grid)
austs$count <- rich.aust; max(rich.aust)

aust_map <- ggplot() +
  geom_sf(data = austs, mapping = aes(fill = count), alpha = 1, color = NA) +
  scale_fill_distiller(direction = 1, palette = "Greens", limits=c(0,6)) +
  #scale_fill_viridis_c(direction=-1, option="rocket") +
  # add map
  geom_sf(data = sahul, color = "black", fill = NA) +
  # crop map
  coord_sf(xlim = c(129, 155), ylim = c(-20, 0)) +
  theme_classic() + theme(legend.position = "none")


library(patchwork)

micro_map / coph_map / aust_map

(coph_map |  aust_map) /
  micro_map /
  (coph_map |  aust_map)

##############################################


