#' Creates a list of unique dens based on the formatted den control data.
#' 
#' @param kontroller a data frame with wrangled den control data
#' @keywords arctic fox unique dens data frame utm
#' @export
#' @return a data frame with unique dens
#' @examples 
#'   uniqueDens(df = df)

uniqueDens <- function(df = control_data) {
  # keep only unique dens based on the last control year (with the most up to date coordinates)
  unike.hi <- as.data.frame(df)
  unike.hi = unike.hi[with(unike.hi, ave(year, loc_code, FUN=max)==year),]
  
  # keep only relevant columns
  unike.hi <- subset(unike.hi, select = c("loc_code",'loc_name','mnt_area','kom_no','kom','fylke','region','UTM33E','UTM33N'))
  
  # sort the list ascending by den location code
  unike.hi <- unike.hi[order(unike.hi$loc_code),]
  
  # remove any spaces from the den location name
  unike.hi$loc_name <- trimws(unike.hi$loc_name)

  # remove any spaces from the mountain (management) area name
  unike.hi$mnt_area <- trimws(unike.hi$mnt_area)
  
  # remove any dens with missing geospatial coordinates
  unike.hi <- subset(unike.hi, !is.na(unike.hi$UTM33E))
  unike.hi <- subset(unike.hi, !is.na(unike.hi$UTM33N))
  
  return(unike.hi)
}