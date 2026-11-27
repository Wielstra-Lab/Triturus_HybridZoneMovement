library(sf)
library(rnaturalearth)
library(dplyr)

library(ggplot2)

setwd("D:/Analyses/Triturus/Isolines")

species <- st_read("TivaTana.shp")

files <- list.files("D:/Analyses/Triturus/Isolines/isolines_0_5_TanaMiss/",
                    pattern = "\\.shp$",
                    full.names = TRUE)

for (f in files) {
  
  iso <- st_read(f, quiet = TRUE)
  
  # Match CRS
  if (st_crs(iso) != st_crs(species)) {
    iso <- st_transform(iso, st_crs(species))
  }
  
  # Clip
  iso_clip <- st_intersection(iso, species)
  
  # Output filename
  out <- file.path(
    "D:/Analyses/Triturus/Isolines/isolines_0_5_TanaMiss/clipped",
    paste0(tools::file_path_sans_ext(basename(f)), "_clipped.shp")
  )
  
  st_write(iso_clip, out, delete_dsn = TRUE, quiet = TRUE)
}

############################################################################################
# Folder containing all shapefiles
foldermac <- "D:/Analyses/Triturus/Isolines/isolines_0_5_TmacMiss/clipped"
folderana <- "D:/Analyses/Triturus/Isolines/isolines_0_5_TanaMiss/clipped"
foldercri <- "D:/Analyses/Triturus/Isolines/isolines_0_5_TcriMiss/clipped"

# List all .shp files
shp_filesmac <- list.files(foldermac,
                        pattern = "\\.shp$",
                        full.names = TRUE)
shp_filesana <- list.files(folderana,
                           pattern = "\\.shp$",
                           full.names = TRUE)
shp_filescri <- list.files(foldercri,
                           pattern = "\\.shp$",
                           full.names = TRUE)

# Read and combine
isolinesmac <- lapply(shp_filesmac, st_read, quiet = TRUE) |>
  bind_rows()
isolinesana <- lapply(shp_filesana, st_read, quiet = TRUE) |>
  bind_rows()
isolinescri <- lapply(shp_filescri, st_read, quiet = TRUE) |>
  bind_rows()

world <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")

average_mac <- st_read("D:/Analyses/Triturus/Isolines/isolines_0_5_TmacMiss/clipped/0_5_isoline_Average_clipped.shp")
average_ana <- st_read("D:/Analyses/Triturus/Isolines/isolines_0_5_TanaMiss/clipped/0_5_isoline_Average_clipped.shp")
average_cri <- st_read("D:/Analyses/Triturus/Isolines/isolines_0_5_TcriMiss/clipped/0_5_isoline_Average_clipped.shp")

polygonmac <- st_read("D:/Analyses/Triturus/Isolines/Triturus_Shapefiles_Unclipped/disTmac.shp")
polygonana <- st_read("D:/Analyses/Triturus/Isolines/Triturus_Shapefiles_Unclipped/disTkar_central.shp")
polygoncri <- st_read("D:/Analyses/Triturus/Isolines/Triturus_Shapefiles_Unclipped/disTcri.shp")
polygoniva <- st_read("D:/Analyses/Triturus/Isolines/Triturus_Shapefiles_Unclipped/disTkar_west.shp")

###################
### Per species ###
###################
mac <- ggplot() +
  geom_sf(data = world,
          fill = "grey98",
          color = "grey60",
          linewidth = 0.2) +
  geom_sf(data = polygonmac,
          color = "grey50",
          linetype = 2,
          linewidth = 0.3) +
  geom_sf(data = polygoniva,
          color = "grey50",
          linetype = 2,
          linewidth = 0.3) +
  geom_sf(data = isolinesmac,
          color = "mediumorchid1",
          linewidth = 0.05,
          alpha = 0.3) +
  geom_sf(data = average_mac,
          color = "mediumorchid3",
          linewidth = 0.5) +
  coord_sf(xlim = c(17, 35), 
           ylim = c(36, 46), 
           expand = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

ggsave(
  filename = "TivaTmac_Clipped.png",
  plot = mac,
  width = 10,
  height = 6,
  dpi = 600
)

ana <- ggplot() +
  geom_sf(data = world,
          fill = "grey98",
          color = "grey60",
          linewidth = 0.2) +
  geom_sf(data = polygonana,
          color = "grey50",
          linetype = 2,
          linewidth = 0.3) +
  geom_sf(data = polygoniva,
          color = "grey50",
          linetype = 2,
          linewidth = 0.3) +
  geom_sf(data = isolinesana,
          color = "steelblue2",
          linewidth = 0.05,
          alpha = 0.3) +
  geom_sf(data = average_ana,
          color = "steelblue3",
          linewidth = 0.5) +
  coord_sf(xlim = c(17, 35), 
           ylim = c(36, 46), 
           expand = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

ggsave(
  filename = "TivaTana_Clipped.png",
  plot = ana,
  width = 10,
  height = 6,
  dpi = 600
)

cri <- ggplot() +
  geom_sf(data = world,
          fill = "grey98",
          color = "grey60",
          linewidth = 0.2) +
  geom_sf(data = polygoncri,
          color = "grey50",
          linetype = 2,
          linewidth = 0.3) +
  geom_sf(data = polygoniva,
          color = "grey50",
          linetype = 2,
          linewidth = 0.3) +
  geom_sf(data = isolinescri,
          color = "firebrick1",
          linewidth = 0.05,
          alpha = 0.3) +
  geom_sf(data = average_cri,
          color = "firebrick3",
          linewidth = 0.5) +
  coord_sf(xlim = c(17, 35), 
           ylim = c(36, 46), 
           expand = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

ggsave(
  filename = "TivaTcri_Clipped.png",
  plot = cri,
  width = 10,
  height = 6,
  dpi = 600
)

###################
### All species ###
###################
all <- ggplot() +
  geom_sf(data = world,
          fill = "grey98",
          color = "grey60",
          linewidth = 0.2) +
  geom_sf(data = polygonmac,
          color = "grey50",
          linetype = 2,
          linewidth = 0.3) +
  geom_sf(data = polygonana,
          color = "grey50",
          linetype = 2,
          linewidth = 0.3) +
  geom_sf(data = polygoncri,
          color = "grey50",
          linetype = 2,
          linewidth = 0.3) +
  geom_sf(data = polygoniva,
          color = "grey50",
          linetype = 2,
          linewidth = 0.3) +
  geom_sf(data = isolinesmac,
          color = "mediumorchid1",
          linewidth = 0.05,
          alpha = 0.3) +
  geom_sf(data = average_mac,
          color = "mediumorchid3",
          linewidth = 0.5) +
  geom_sf(data = isolinesana,
          color = "steelblue2",
          linewidth = 0.05,
          alpha = 0.3) +
  geom_sf(data = average_ana,
          color = "steelblue3",
          linewidth = 0.5) +
  geom_sf(data = isolinescri,
          color = "firebrick1",
          linewidth = 0.05,
          alpha = 0.3) +
  geom_sf(data = average_cri,
          color = "firebrick3",
          linewidth = 0.5) +
  coord_sf(xlim = c(17, 35), 
           ylim = c(36, 46), 
           expand = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

ggsave(
  filename = "All_Clipped.png",
  plot = all,
  width = 10,
  height = 6,
  dpi = 600
)
###########
# Spare code

ggplot() +
  geom_sf(data = world, fill = "grey90") +
  geom_sf(data = isolines, color = "blue", linewidth = 0.05) +
  geom_sf(
    data = average_line,
    color = "black",
    linewidth = 1
  ) +
  coord_sf(
    xlim = c(17, 35),
    ylim = c(36, 45),
    expand = FALSE
  ) +
  theme_minimal()

ggplot() +
  geom_sf(data = world,
          fill = "grey95",
          color = "grey60",
          linewidth = 0.2) +
  geom_sf(data = isolines,
          color = "steelblue",
          linewidth = 0.05,
          alpha = 0.3) +
  geom_sf(data = average_line,
          color = "black",
          linewidth = 0.8) +
  coord_sf(
    xlim = c(17, 35),
    ylim = c(36, 45),
    expand = FALSE
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p <- ggplot() +
  geom_sf(data = world,
          fill = "grey95",
          color = "grey60",
          linewidth = 0.2) +
  geom_sf(data = isolines,
          color = "steelblue",
          linewidth = 0.05,
          alpha = 0.3) +
  geom_sf(data = average_line,
          color = "black",
          linewidth = 0.8) +
  coord_sf(
    xlim = c(17, 35),
    ylim = c(36, 45),
    expand = FALSE
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  )
