# Metadata
This directory contains all files with sample information produced from sampling, sequencing, and analysis, as well as sample lists for ease of coding. Files are organized into subdirectories based on the step of the project they came from.


## Field Sampling

| Filename                                | Data Source               | Description                                  |
| --------------------------------------- | ------------------------- | -------------------------------------------- |
| 40samples_KaseyUF_updateDec2020.txt.kmz | Rebecca Jones, U Tasmania | Google Earth file with provenance of samples |
| 40samples_KaseyUF_updateDec2020.xlsx    | Rebecca Jones, U Tasmania | Spreadsheet of provenance and genotype data for samples |

Plotted coordinates of samples in `R`. Used [Geocomputation with R](https://r.geocompx.org/adv-map) and [R for Geospatial Sciences](https://jsimkins2.github.io/geog473-673/spatial-plots-with-ggplot2.html).
```R
library(sf)
library(tmap)
library(ggplot2)
library(ozmaps)
library(grid)

# import location data for samples
loc_tab <- read.csv("C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/01.field_sampling/40samples_KaseyUF_updateDec2020.csv")
loc_tab <- loc_tab[,c("lat", "long", "sp_pop")]
# remove individuals for which there isn't recorded location data
loc_tab <- loc_tab[which(loc_tab$lat != "." & loc_tab$long != "."),]
loc_points <- st_as_sf(loc_tab, coords = c("long", "lat"), crs = st_crs(4326))

# FIGURE 1A
# extract shape of Tasmania and square around Hobart
tas_shape <- ozmaps::ozmap_states[6,]
tas_reg <- st_bbox(c(xmin = 147.00, xmax = 147.75, ymin = -43.17, ymax = -42.5), crs = st_crs(4326))
# make background for map with Tasmania landform
tas_map <- tm_shape(tas_shape, bbox = tas_reg) + tm_fill("white") + tm_borders() + tm_layout(bg.color = "#79B9B4")
tas_map <- tas_map + tm_graticules(labels.size = c(0.7)) + tm_scale_bar(breaks = c(0, 10, 20), text.size = 1, position = c("left", "top"))
# add points for samples with color and shape coding by sampling group
tas_map <- tas_map + tm_shape(loc_points) + tm_symbols(shape = "sp_pop", shapes = c("glob_pure" = 22, "glob_MR" = 21, "cord_MR" = 24), col = "sp_pop", palette = c("glob_pure" = "#1D91C0", "glob_MR" = "#454545", "cord_MR" = "#FED976"), border.col = "black", legend.col.show = FALSE, legend.shape.show = FALSE)
# add legend
tas_map <- tas_map + tm_add_legend(type = "symbol", labels = c("globulus (ref.)", "globulus (intr.)", "cordata"), title = "Sample Group", col = c("#1D91C0", "#454545", "#FED976"), shape = c(22,21,24)) + tm_layout(legend.outside = TRUE, legend.outside.position = "left", legend.outside.size = c(0.2, 0.3), legend.frame = TRUE, legend.text.size = 0.8, legend.title.size = 1.1, legend.bg.color = "white")
# add inset map with context location in Tasmania
inset_map <- tm_shape(tas_shape) + tm_fill("white") + tm_borders() + tm_shape(st_as_sfc(tas_reg)) + tm_borders(col = "black", lwd = 3) + tm_layout(bg.color = "#57939b")
# normalize aspect ratios (function copied from https://r.geocompx.org/adv-map)
norm_dim = function(obj){
    bbox = st_bbox(obj)
    width = bbox[["xmax"]] - bbox[["xmin"]]
    height = bbox[["ymax"]] - bbox[["ymin"]]
    w = width / max(width, height)
    h = height / max(width, height)
    return(unit(c(w, h), "snpc"))
}
main_dim = norm_dim(tas_reg)
inset_dim = norm_dim(tas_shape)
# create viewports
main_view <- viewport(width = main_dim[1] * 1.1, height = main_dim[2] * 1.1, x = unit(1, "npc"), y = unit(0.15, "npc"), just = c("right", "bottom"))
inset_view <- viewport(width = inset_dim[1] * 0.225, height = inset_dim[2] * 0.225, x = unit(0.8825,"npc"), y = unit(0.05, "npc"), just = c("right", "bottom"))
# plot maps
grid.newpage()
print(tas_map, vp = main_view)
pushViewport(main_view)
print(inset_map, vp = inset_view)

# FIGURE 1B (zoomed into Meehan Range)
# define Meehan Range box
mr_reg <- st_bbox(c(xmin = 147.38, xmax = 147.4175, ymin = -42.8575, ymax = -42.825), crs = st_crs(4326))
# add elements to Meehan Range map
mr_map <- tm_shape(tas_shape, bbox = mr_reg) + tm_fill("white") + tm_borders() + tm_layout(bg.color = "#79B9B4") # landmasss
mr_map <- mr_map + tm_graticules(labels.size = c(0.7)) + tm_scale_bar(breaks = c(0, 0.5, 1), text.size = 1, position = c("left", "top")) # axis labels and scale bar
# filter point sets by sample group
glob_tab <- loc_tab[which(loc_tab$sp_pop == "glob_MR"),]
glob_points <- st_as_sf(glob_tab, coords = c("long", "lat"), crs = st_crs(4326))
cord_tab <- loc_tab[which(loc_tab$sp_pop == "cord_MR"),]
cord_points <- st_as_sf(cord_tab, coords = c("long", "lat"), crs = st_crs(4326))
# plot sample groups separately so E. cordata is on top
mr_map <- mr_map + tm_shape(glob_points) + tm_symbols(shape = 21, col = "#454545", border.col = "black", legend.col.show = FALSE, legend.shape.show = FALSE) # E. globulus points
mr_map <- mr_map + tm_shape(cord_points) + tm_symbols(shape = 24, col = "#FED976", border.col = "black", legend.col.show = FALSE, legend.shape.show = FALSE)
# no legend since I'll be placing this next to the zoomed-out map
# add inset map of SW Tasmania with sample points plotted for context
tas_inset_map <- tm_shape(tas_shape, bbox = tas_reg) + tm_fill("white") + tm_borders() + tm_layout(bg.color = "#79B9B4")
tas_inset_map <- tas_inset_map + tm_shape(loc_points) + tm_symbols(size = 0.25, shape = "sp_pop", shapes = c("glob_pure" = 22, "glob_MR" = 21, "cord_MR" = 24), col = "sp_pop", palette = c("glob_pure" = "#1D91C0", "glob_MR" = "#999999", "cord_MR" = "#FED976"), legend.col.show = FALSE, legend.shape.show = FALSE)
tas_inset_map <- tas_inset_map + tm_shape(st_as_sfc(mr_reg)) + tm_borders(col = "black", lwd = 2)
# normalize aspect ratios
mr_dim <- norm_dim(mr_reg)
# create viewports
mr_view <- viewport(width = mr_dim[1], height = mr_dim[2])
tas_inset_view <- viewport(width = main_dim[1] * 0.35, height = main_dim[2] * 0.35, x = unit(0.9,"npc"), y = unit(0.06, "npc"), just = c("right", "bottom"))
# plot maps
grid.newpage()
print(mr_map, vp = mr_view)
pushViewport(mr_view)
print(tas_inset_map, vp = tas_inset_view)
```


## RAPiD Sequencing

| Filename                                | Data Source               | Description                                  |
| --------------------------------------- | ------------------------- | -------------------------------------------- |
| UFL_014801_Normal_Plate_Layout_22-02-21_GOOD.xlsx | Ariane Gélinas Marion, U Tasmania | Plate layout of samples sent to RAPiD Genomics for sequencing |
| UFL_014801_SampleSheet_june08.csv  | RAPiD Genomics | Table with sample codes and barcodes for run 2 of sequencing |
| UFL_014801_SampleSheet_may04.csv   | RAPiD Genomics | Table with sample codes and barcodes for run 1 of sequencing |


## Sequence Analysis

| Filename                                | Data Source               | Description                                  |
| --------------------------------------- | ------------------------- | -------------------------------------------- |
| chromosome_list.txt  | _Eucalyptus grandis_ genome | List of chromosome names in reference genome for code looping |
| sample_ids.txt          | sample_sequencing_metadata_all.xlsx | List of RAPiD ID numbers for each sample sequenced |
| sample_sequencing_metadata.csv | sample_sequencing_metadata_all.xlsx | Table of compiled metadata from all sources for each sequencing run file provided by RAPiD |
| sample_sequencing_metadata_all.xlsx | Filenames and headers of FASTQ files sent by RAPiD, sample sheets sent by RAPiD | Metadata for each FASTQ file, manually extracted and compiled. Files are split by sample, run, and direction. |
| sample_spp_table.csv                    | sample_sequencing_metadata_all.xlsx, 40samples_KaseyUF_updateDec2020.xlsx | Table translating RAPiD sample codes to U Tasmania accession numbers and species designation |
| seq_ids.txt | sample_sequencing_metadata_all.xlsx | List of all IDs unique to sample and sequencing run for code looping |

### Columns of note in sample_sequencing_metadata_all.xlsx
* **Filename**: The raw name of the FASTQ file sent by RAPiD Genomics
* **Sample**: RAPiD's ID for each _sample_, each of which was subject to two runs of sequencing.
* **RunSample**: RAPiD's ID for each _run_, unique for each sample and run, should identify each individual library prep.
* **Customer_Code**: Original accession numbers for each sample in U Tasmania data
