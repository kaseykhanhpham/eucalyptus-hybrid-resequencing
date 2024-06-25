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
loc_tab <- loc_tab[,c("RJgeno", "lat", "long", "sp_pop", "sumDArTorAFLP", "JLA.")]
# remove individuals for which there isn't recorded location data
loc_tab <- loc_tab[which(loc_tab$lat != "." & loc_tab$long != "."),]
loc_points <- st_as_sf(loc_tab, coords = c("long", "lat"), crs = st_crs(4326))

# configure column for JLA+ AND nuclear marker combos
intr_markers <- rep(21, nrow(loc_points))
names(intr_markers) <- loc_points$RJgeno
# just JLA+
intr_markers[which(loc_tab$JLA. == "S83")] <- rep(24, length(which(loc_tab$JLA. == "S83")))
# just nuclear markers
intr_markers[which(loc_tab$sumDArTorAFLP > 0)] <- rep(22, length(which(loc_tab$sumDArTorAFLP > 0)))
# both
both_indices <- which(loc_tab$JLA. == "S83" & loc_tab$sumDArTorAFLP > 0)
intr_markers[both_indices] <- rep(23, length(both_indices))
# labels for legend
intr_markers_leg <- unique(intr_markers)
intr_markers_lab <- c("S83 chloroplast haplotype", "DArT or AFLP nuclear", "S83 & nuclear", "None")

# FIGURE 1A
# extract shape of Tasmania and square around Hobart
tas_shape <- ozmaps::ozmap_states[6,]
tas_reg <- st_bbox(c(xmin = 147.00, xmax = 147.75, ymin = -43.17, ymax = -42.5),
                   crs = st_crs(4326)) 
# make background for map with Tasmania landform
tas_map <- tm_shape(tas_shape, bbox = tas_reg) + tm_fill("white") +
           tm_borders(lwd = 2, col = "white") +
           tm_layout(bg.color = "#7898D1")
tas_map <- tas_map + tm_graticules(labels.size = 3.2, lwd = 5, col = "#8A9062") +
           tm_scale_bar(breaks = c(0, 10, 20), text.size = 4, lwd = 5, 
                        position = c("left", "top"))

# add points for samples with color and shape coding by sampling group
tas_map <- tas_map + tm_shape(loc_points) +
           tm_symbols(size = 14,
                      border.lwd = 3.5,
                      shape = "RJgeno", 
                      shapes = intr_markers,
                      col = "sp_pop",
                      palette = c("glob_pure" = "#044075", "glob_MR" = "#13BDD7", "cord_MR" = "#FFCB3D"),
                      border.col = "black",
                      legend.col.show = FALSE, legend.shape.show = FALSE)

# add inset map with context location in Tasmania
inset_map <- tm_shape(tas_shape) + tm_fill("white") + tm_borders(col = "black", lwd = 3.5) +
             tm_shape(st_as_sfc(tas_reg)) + tm_borders(col = "black", lwd = 5) +
             tm_layout(bg.color = "#7898D1")

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
main_view <- viewport(width = main_dim[1], height = main_dim[2])
inset_view <- viewport(width = inset_dim[1] * 0.35, height = inset_dim[2] * 0.35,
                       x = unit(0.9, "npc"), y = unit(0.06, "npc"),
                       just = c("right", "bottom"))
# plot maps
png("sampling_map_zoomed_out.png", width = 2000, height = 2400)
grid.newpage()
print(tas_map, vp = main_view)
pushViewport(main_view)
print(inset_map, vp = inset_view)
dev.off()

# FIGURE 1B (zoomed into Meehan Range)
# filter point sets by sample group
glob_tab <- loc_tab[which(loc_tab$sp_pop == "glob_MR"),]
glob_points <- st_as_sf(glob_tab, coords = c("long", "lat"), crs = st_crs(4326))
cord_tab <- loc_tab[which(loc_tab$sp_pop == "cord_MR"),]
cord_points <- st_as_sf(cord_tab, coords = c("long", "lat"), crs = st_crs(4326))

# define Meehan Range box
mr_reg <- st_bbox(c(xmin = 147.38, xmax = 147.4175,
                    ymin = -42.8575, ymax = -42.825), crs = st_crs(4326))
# add elements to Meehan Range map
mr_map <- tm_shape(tas_shape, bbox = mr_reg) + tm_fill("white") +
          tm_borders() + tm_layout(bg.color = "#7898D1") # landmasss
mr_map <- mr_map + tm_graticules(labels.size = 2.8, lwd = 5, col = "#B2A691") +
                   tm_scale_bar(breaks = c(0, 0.5, 1), text.size = 4, lwd = 5, 
                                position = c("left", "top")) # axis labels and scale bar

# plot sample groups separately so E. cordata is on top
# E. globulus points
mr_map <- mr_map + tm_shape(glob_points) +
          tm_symbols(size = 14,
                     shape = "RJgeno", shapes = intr_markers,
                     col = "#13BDD7",
                     border.col = "black", border.lwd = 3.5,
                     legend.col.show = FALSE, legend.shape.show = FALSE)
 # E. cordata points
mr_map <- mr_map + tm_shape(cord_points) +
          tm_symbols(size = 14,
                     shape = "RJgeno", shapes = intr_markers,
                     col = "#FFCB3D",
                     border.col = "black", border.lwd = 3.5,
                     legend.col.show = FALSE, legend.shape.show = FALSE)
# no legend since I'll be placing this next to the zoomed-out map
# add inset map of SW Tasmania with sample points plotted for context
tas_inset_map <- tm_shape(tas_shape, bbox = tas_reg) + tm_fill("white") +
                 tm_borders() + tm_layout(bg.color = "#7898D1")
tas_inset_map <- tas_inset_map + tm_shape(loc_points) +
                 tm_symbols(size = 3,
                            shape = "RJgeno", shapes = intr_markers,
                            col = "sp_pop",
                            palette = c("glob_pure" = "#044075", "glob_MR" = "#13BDD7", "cord_MR" = "#FFCB3D"),
                            border.col = "#B2B793", border.lwd = 1,
                            legend.col.show = FALSE,
                            legend.shape.show = FALSE)
tas_inset_map <- tas_inset_map + tm_shape(st_as_sfc(mr_reg)) +
                 tm_borders(col = "black", lwd = 2)
# normalize aspect ratios
mr_dim <- norm_dim(mr_reg)
# create viewports
mr_view <- viewport(width = mr_dim[1], height = mr_dim[2])
tas_inset_view <- viewport(width = main_dim[1] * 0.35, height = main_dim[2] * 0.35,
                           x = unit(0.9,"npc"), y = unit(0.06, "npc"),
                           just = c("right", "bottom"))
# plot maps
png("sampling_map_zoomed_in.png", width = 2000, height = 2400)
grid.newpage()
print(mr_map, vp = mr_view)
pushViewport(mr_view)
print(tas_inset_map, vp = tas_inset_view)
dev.off()

## FIGURE 1 LEGEND
# add legend for sample group (color) and introgression evidence (shape)
tas_map_leg <- tas_map + tm_add_legend(type = "fill",
                                       border.col = "black", border.lwd = 3.5,
                                       labels = c("globulus ('pure')", "globulus (MR)", "cordata"),
                                       title = "Sample Group",
                                       col = c("#044075", "#13BDD7", "#FFCB3D")) +
                         tm_add_legend(size = 7.5,
                                       type = "symbol",
                                       labels = intr_markers_lab,
                                       shape = intr_markers_leg,
                                       title = "Evidence of Introgression",
                                       col = "black") +
                         tm_layout(legend.only = TRUE,
                                   legend.frame = "black", legend.frame.lwd = 4,
                                   legend.text.size = 7.5, legend.title.size = 9,
                                   legend.bg.color = "white")
png("sampling_map_legend.png", width = 2400, height = 2400)
print(tas_map_leg)
dev.off()
```

Atlas of Living Australia Euc occurrence map session: https://spatial.ala.org.au/?ss=1713500431205


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
