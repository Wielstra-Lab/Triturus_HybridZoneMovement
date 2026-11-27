################################################################################
library(vcfR)
library(dplyr)
library(tidyr)
library(raster)
library(sf)
library(akima)
library(ggplot2)
library(viridis)
library(data.table)

setwd("D:/Analyses/Triturus/Isolines")

# Read in your vcf, subsetted to 1 diagnostic SNP per marker with max 2 alleles per site, optional missing data filter
vcf <- read.vcfR("TivaTmac.subset.vcf")

# Read in sample metadata
samples <- read.csv("coords_TivaTmac.csv", stringsAsFactors = FALSE)

# Extract genotype matrix (SNPs x samples)
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

# Save genotype matrix to file, then change encoding of REF SNP, and change to hybrid index(based on REF samples, change it so REF sample REF allele is always 0, alt allele is always 1 and heterozygotes are 0.5)
write.csv(gt, file = "TivaTmacgt.csv")

# Read updated genotype matrix and flip for plotting
snp_data <- read.csv("TivaTmacgtHIMiss.csv", stringsAsFactors = FALSE, check.names =  FALSE)
snp_dataT <- transpose(snp_data)
rownames(snp_dataT) <- colnames(snp_data)
colnames(snp_dataT) <- snp_dataT[1,]
snp_dataT <- snp_dataT[-1, ] 

snp_dataT$latitude <- samples$latitude
snp_dataT$longitude <- samples$longitude

snp_dataT <- snp_dataT %>% mutate_if(is.character, as.numeric)

### Average 0.5 isolines over all snps
interp_avg <- interp(x = snp_dataT$longitude, y = snp_dataT$latitude, z = snp_dataT$Average, 
                     duplicate = "mean", linear = TRUE, extrap = FALSE,
                     nx = 200,
                     ny = 200)

#Per SNP
save_0_5_isoline_per_snp <- function(snp_col, samples_df, output_dir = "isolines_0_5_TmacLTET") {
  # Filter out missing data for coordinates and SNP value
  snp_dataT <- samples_df %>%
    filter(!is.na(longitude), !is.na(latitude), !is.na(.data[[snp_col]]))
  
  if (nrow(snp_dataT) < 3) {
    message(paste("Not enough data points for", snp_col, "- skipping"))
    return(NULL)
  }
  
  # Interpolate SNP values over spatial grid
  interp_res <- interp(x = snp_dataT$longitude, y = snp_dataT$latitude, z = snp_dataT[[snp_col]],
                       duplicate = "mean", linear = TRUE, extrap = FALSE,
                       nx = 200,
                       ny = 200)
  
  # Extract contour lines at level 0.5
  contour_list <- contourLines(x = interp_res$x, y = interp_res$y, z = interp_res$z, levels = 0.5)
  
  if (length(contour_list) == 0) {
    message(paste("No 0.5 contour found for", snp_col))
    return(NULL)
  }
  
  # Convert contour lines to sf objects
  sf_lines <- lapply(seq_along(contour_list), function(i) {
    cl <- contour_list[[i]]
    coords <- cbind(cl$x, cl$y)
    sf::st_linestring(coords) %>%
      sf::st_sfc() %>%
      sf::st_sf(level = cl$level)
  })
  
  sf_contours <- do.call(rbind, sf_lines)
  sf_contours <- sf::st_set_crs(sf_contours, 4326)  # WGS84
  
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  shp_path <- file.path(output_dir, paste0("0_5_isoline_", snp_col, ".shp"))
  sf::st_write(sf_contours, shp_path, delete_layer = TRUE, quiet = TRUE)
  
  message(paste("Saved 0.5 isoline shapefile for", snp_col, "to", shp_path))
}

# Assuming your dataframe with coordinates and SNP data is called samples_df
# and SNP columns are all except longitude and latitude
snp_cols <- setdiff(colnames(snp_dataT), c("longitude", "latitude"))

for (snp_col in snp_cols) {
  save_0_5_isoline_per_snp(snp_col, snp_dataT, output_dir = "isolines_0_5_TmacLTET")
}

### Admixture isolines
adm_data <- read.csv("Tana_admixture_results_aligned.csv", stringsAsFactors = FALSE, check.names =  FALSE)

adm_dataT <- transpose(adm_data)
rownames(adm_dataT) <- colnames(adm_data)
colnames(adm_dataT) <- adm_dataT[1,]
adm_dataT <- adm_dataT[-1, ] 

adm_dataT <- adm_dataT %>% mutate_if(is.character, as.numeric)

#Per SNP
save_0_5_isoline_per_snp <- function(snp_col, samples_df, output_dir = "ADMisolines_0_5_Tana") {
  # Filter out missing data for coordinates and SNP value
  adm_dataT <- samples_df %>%
    filter(!is.na(longitude), !is.na(latitude), !is.na(.data[[snp_col]]))
  
  if (nrow(adm_dataT) < 3) {
    message(paste("Not enough data points for", snp_col, "- skipping"))
    return(NULL)
  }
  
  # Interpolate SNP values over spatial grid
  interp_res <- interp(x = adm_dataT$longitude, y = adm_dataT$latitude, z = adm_dataT[[snp_col]],
                       duplicate = "mean", linear = TRUE, extrap = FALSE)
  
  # Extract contour lines at level 0.5
  contour_list <- contourLines(x = interp_res$x, y = interp_res$y, z = interp_res$z, levels = 0.5)
  
  if (length(contour_list) == 0) {
    message(paste("No 0.5 contour found for", snp_col))
    return(NULL)
  }
  
  # Convert contour lines to sf objects
  sf_lines <- lapply(seq_along(contour_list), function(i) {
    cl <- contour_list[[i]]
    coords <- cbind(cl$x, cl$y)
    sf::st_linestring(coords) %>%
      sf::st_sfc() %>%
      sf::st_sf(level = cl$level)
  })
  
  sf_contours <- do.call(rbind, sf_lines)
  sf_contours <- sf::st_set_crs(sf_contours, 4326)  # WGS84
  
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  shp_path <- file.path(output_dir, paste0("0_5_isoline_", snp_col, ".shp"))
  sf::st_write(sf_contours, shp_path, delete_layer = TRUE, quiet = TRUE)
  
  message(paste("Saved 0.5 isoline shapefile for", snp_col, "to", shp_path))
}

# Assuming your dataframe with coordinates and SNP data is called samples_df
# and SNP columns are all except longitude and latitude
snp_cols <- setdiff(colnames(adm_dataT), c("longitude", "latitude"))

for (snp_col in snp_cols) {
  save_0_5_isoline_per_snp(snp_col, adm_dataT, output_dir = "ADMisolines_0_5_Tana")
}
