
library(rgdal)
library(ggplot2)
library(dplyr)
library(ggmap)

# getting the coordinates
location_coord <- read.table("Locations_coordinates.txt", header = T)
str(location_coord)
location_coord %>% head
coord <- location_coord[, c(1, 7:8)]
coord %>% head
rownames(coord) <- coord[,1]


colnames(coord)[1] <- "Population"

# getting the keys_azurea to do a left-join
keys_azurea <- read.table("keys_azurea.txt", header = T, fileEncoding="UTF-8-BOM")
str(keys_azurea)
levels(keys_azurea$Population)


# after a few checks, everything seems to be fine
# notice that I have finally fixed the teles pires and itauba in the population names

setdiff(coord$Population, keys_azurea$Population)


# Now will add the coordinates to the keys_azurea data.frame

keys_geo_azurea <- dplyr::left_join(keys_azurea, coord , by = "Population") 
rownames(keys_geo_azurea) <- keys_geo_azurea$DNASample
head(keys_geo_azurea)
dim(keys_geo_azurea)

keys_geo_azurea_to_remove <- subset(keys_geo_azurea, 	
                     DNASample == "1000520" | # ok
                     DNASample == "1000558"  | # ok
                     DNASample == "1000588" | # ok
                     DNASample == "1000590" | # ok
                     DNASample == "1000624"  | # ok
                     DNASample == "1000639"  | # ok
                     DNASample == "218"     | # ok
                     DNASample == "230"     | # ok new
                     DNASample == "2311"    | # ok
                     DNASample == "2307"    | # ok new
                     DNASample == "2303"     | # ok new
                     DNASample == "4359"     | # ok
                     DNASample == "35256"     | # ok
                     DNASample == "35291B"    |
                     DNASample == "35292B"    |
                     DNASample == "P3_RO_EAZU_EX_1") # ok

dim(keys_geo_azurea_to_remove)

keys_geo_azurea <- keys_geo_azurea %>% dplyr::anti_join(keys_geo_azurea_to_remove)
dim(keys_geo_azurea)
head(keys_geo_azurea)
rownames(keys_geo_azurea) <- keys_geo_azurea$DNASample

  library(geosphere)
  library(vegan)
  library(sp)
  library(geosphere)
  
  geo_dist_azurea <- as.matrix(distm(keys_geo_azurea[, c("long", "lat") ], keys_geo_azurea[, c("long", "lat") ], fun = distHaversine))/1000 # unit here is meters
  head(geo_dist_azurea)
  rownames(geo_dist_azurea) <- keys_geo_azurea$DNASample
  colnames(geo_dist_azurea) <- keys_geo_azurea$DNASample
  keys_geo_azurea %>% head
  
  head(geo_dist_azurea)
  
  PVE_azurea_PGAbfp_HC_50miss_indmiss_qualOdp_fqbd_noindels_minDp30Pl30.dt1 <- read.table("PVE_azurea_PGAbfp_HC_50miss_indmiss_qualOdp_fqbd_noindels_minDp30Pl30.distance1Table.DpPl30.tbl", 
             header = T, check.names=FALSE, row.names = 1)
  
  mantel(geo_dist_azurea, PVE_azurea_PGAbfp_HC_50miss_indmiss_qualOdp_fqbd_noindels_minDp30Pl30.dt1)
  
  
  
  
  
  keys_crassipes <- read.table("keys_crassipes.txt", header = T, fileEncoding="UTF-8-BOM")
  str(keys_crassipes)
  levels(keys_crassipes$Population)
  
  
  # after a few checks, everything seems to be fine
  # notice that I have finally fixed the teles pires and itauba in the population names
  
  setdiff(coord$Population, keys_crassipes$Population)
  
  
  # Now will add the coordinates to the keys_crassipes data.frame
  
  keys_geo_crassipes <- dplyr::left_join(keys_crassipes, coord , by = "Population") 
  rownames(keys_geo_crassipes) <- keys_geo_crassipes$DNASample
  head(keys_geo_crassipes)
  dim(keys_geo_crassipes)
  
  keys_geo_crassipes_to_remove <- subset(keys_geo_crassipes, 	
                                         DNASample == "1000520" | # ok
                                           DNASample == "1000588" | # ok
                                           DNASample == "218"     | # ok
                                           DNASample == "2311")

  dim(keys_geo_crassipes_to_remove)
  
  keys_geo_crassipes <- keys_geo_crassipes %>% dplyr::anti_join(keys_geo_crassipes_to_remove)
  dim(keys_geo_crassipes)
  head(keys_geo_crassipes)
  rownames(keys_geo_crassipes) <- keys_geo_crassipes$DNASample
  
  library(geosphere)
  library(vegan)
  library(sp)
  library(geosphere)
  
  geo_dist_crassipes <- as.matrix(distm(keys_geo_crassipes[, c("long", "lat") ], keys_geo_crassipes[, c("long", "lat") ], fun = distHaversine))/1000 # unit here is meters
  head(geo_dist_crassipes)
  rownames(geo_dist_crassipes) <- keys_geo_crassipes$DNASample
  colnames(geo_dist_crassipes) <- keys_geo_crassipes$DNASample
  keys_geo_crassipes %>% head
  
  head(geo_dist_crassipes)
  
  PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_fqbd_noindels_minDp30Pl30.dt1 <- read.table("PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_fqbd_noindels_minDp30Pl30.distance1Table.DpPl30.tbl", 
                                                                                             header = T, check.names=FALSE, row.names = 1)
  
  mantel(geo_dist_crassipes, PVE_crassipes_PGAbfp_HC_50miss_indmiss_qualOdp_fqbd_noindels_minDp30Pl30.dt1)
  
  