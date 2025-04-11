#' Assign habitat patches (polygons) to geographical coordinates (points)
#'
#' @param df A data frame containing at least an UTM Easting (named 'UTM33E')
#' and an UTM Northing (named 'UTM33N') column, with geographical coordinates in
#' UTM zone 33 (EPSG 25833).
#' @return The original data frame with an extra column with the assigned patch number
#' @export
#' @examples
#'      getPatch(df)

getPatch <- function(df) {
  if(!require(sf)) { print("Package 'sf' not installed"); break; }

  # inform user
  cat("  - Associating habitat patch information to coordinates...")
  
  # convert data frame to an sf object
  datapoints <- st_as_sf(df, coords = c("UTM33E", "UTM33N"))
  
  # set projection (UTM33, EPSG 25833) explicitly
  datapoints <- st_set_crs(datapoints, 25833)
  
  # read provided shapefile with mountain patch polygons
  omr <- st_read("maps/fjellpatch_inter.shp", quiet = T, options = "ENCODING=WINDOWS-1252")
  
  # st_intersect will use row order of polygons to match points within, so make
  # sure that patch ID's ('fid') are ordered according to row numbers
  omr <- omr[order(omr$fid),]
  rownames(omr) <- NULL

  # convert map attributes to data frame
  # - bookkeeping for associating mountain (management) areas and habitat patches
  patches <- as.data.frame(omr)
  patches <- patches %>% dplyr::select(fid, Mnt_area)
  patches$id <- as.integer(rownames(patches))

  # extract the patch number (fid) for each geographical point
  areaIntersect <- st_intersects(datapoints, omr, sparse = T)
  
  # convert NULL values to NA (for data consistency, otherwise NA's are removed)
  areaIntersect[sapply(areaIntersect, function(x) length(x)==0L)] <- NA
  
  # retain only one habitat patch per data point (use first match in border cases)
  areaIntersect <- sapply(areaIntersect, "[[", 1)
  
  # convert to data frame and merge with mountain (management) area names (by row ID)
  areaIntersect <- as.data.frame(areaIntersect)
  colnames(areaIntersect)[1] <- "id"
  areaIntersect <- dplyr::left_join(areaIntersect, patches, by = "id")
  
  # assign the intersected patches to the original data frame and return
  df$patch <- areaIntersect$fid
  df$patchArea <- areaIntersect$Mnt_area
  df$patch <- as.character(df$patch)

  cat("Done.\n")
    
  return(df)
}