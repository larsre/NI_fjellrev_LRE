#' Calculate the number of litters for all viable spatial units
#' Spatial Units:
#'    municipality (NI-level),
#'    patch (habitat fragment level),
#'    mountain area (sub-population/management level)
#' @param df A data frame with wrangled den control data.
#' @param export Save the output to CSV file (one per spatial unit)
#' @return A list of data frames with litter sizes, one data frame per spatial unit
#' @export
#' @examples
#'      litterSize(df = control_data, export = TRUE)

litterSize <- function(df, export = TRUE) {
  # set globals
  options(dplyr.summarise.inform = FALSE)
  returnList <- list()

  # inform user  
  cat("  Calculating number of litters:\n")
  cat("  - Per municipality...")
  
  
  ### Municipality
  
  ## Number of litters per unique den-year and municipality
  lt.dy <- df %>%
    group_by(year, kom, loc_code) %>%
    summarise(yngling = sum(yngling, na.rm = T), valper = sum(valper, na.rm = T))
  
  # retain only records of successful breedings
  lt.dy <- subset(lt.dy, lt.dy$yngling > 0)
  
  ## Number of litters per year and municipality
  # create a vector of unique municipalities
  koms <- names(table(lt.dy$kom))
  
  # filter the complete control data set on municipalities with breeding
  lt.cd.koms <- control_data[control_data$kom %in% koms,]
  colnames(lt.cd.koms)[2] <- "cyear"
  lt.cd.koms <- as.data.frame(lt.cd.koms)
  lt.dy3 <- lt.cd.koms %>%
    group_by(cyear, kom) %>%
    summarise(yngling = sum(yngling, na.rm = T), valper = sum(valper, na.rm = T)) %>%
    as.data.frame()
  
  # fill in missing years with zeros
  # NOTE: input must be a data.frame (NOT a tibble) for 'complete' function to work properly
  lt.dy3 <- lt.dy3 %>%
    tidyr::complete(kom, cyear = full_seq(cyear, 1),
                    fill = list(yngling = 0, valper = 0)) %>%
    as.data.frame()
  returnList$municipality <- lt.dy3[order(lt.dy3$kom, lt.dy3$cyear),]

  cat(" Done.\n")
  
  # export data if flag is TRUE
  if(export) {
    cat("  - [Export] Writing output to file\n")
    data.table::fwrite(returnList$municipality, file = "output/table_litters_per_year-kommune.csv", sep = ";", dec = ",", row.names = FALSE, bom = F)
  }
  

  ### Habitat fragment (patch)
  cat("  - Per habitat...")
  
  # remove samples missing patch association
  df2 <- subset(df, !is.na(df$patch))
  
  ## Number of litters per unique den-year and patch area
  lt.dy <- df2 %>%
    group_by(year, patch, loc_code) %>%
    summarise(yngling = sum(yngling, na.rm = T), valper = sum(valper, na.rm = T))
  
  # retain only records of successful breedings
  lt.dy <- subset(lt.dy, lt.dy$yngling > 0)
  
  ## Number of litters per year and municipality
  
  # create a vector of unique patch areas
  ptch <- names(table(lt.dy$patch))
  
  # filter the complete control data set on municipalities with breeding
  lt.cd.ptch <- control_data[control_data$patch %in% ptch,]
  colnames(lt.cd.ptch)[2] <- "cyear"
  lt.cd.ptch <- as.data.frame(lt.cd.ptch)
  lt.dy3 <- lt.cd.ptch %>%
    group_by(cyear, patch) %>%
    summarise(yngling = sum(yngling, na.rm = T), valper = sum(valper, na.rm = T)) %>%
    as.data.frame()
  
  # fill in missing years with zeros
  # NOTE: input must be a data.frame (NOT a tibble) for 'complete' function to work properly
  lt.dy3 <- lt.dy3 %>%
    tidyr::complete(patch, cyear = full_seq(cyear, 1),
                    fill = list(yngling = 0, valper = 0)) %>%
    as.data.frame()
  returnList$patch <- lt.dy3[order(lt.dy3$patch, lt.dy3$cyear),]
  
  cat(" Done.\n")
  
  # export data if flag is TRUE
  if(export) {
    cat("  - [Export] Writing output to file\n")
    data.table::fwrite(returnList$patch, file = "output/table_litters_per_year-patch.csv", sep = ";", dec = ",", row.names = FALSE, bom = F)
  }
  
  
  ### Mountain area (management/sub-population)
  ### (derived from the sub-population the habitat patch belongs to)
  cat("  - Per sub-population...")
  
  ## Number of litters per unique den-year and sub-population
  lt.dy <- df %>%
    group_by(year, mnt_area, loc_code) %>%
    summarise(yngling = sum(yngling, na.rm = T), valper = sum(valper, na.rm = T))
  
  # retain only records of successful breedings
  lt.dy <- subset(lt.dy, lt.dy$yngling > 0)
  
  ## Number of litters per year and municipality
  # create a vector of unique patch areas
  mnta <- names(table(lt.dy$mnt_area))
  
  # filter the complete control data set on municipalities with breeding
  lt.cd.mnta <- control_data[control_data$mnt_area %in% mnta,]
  colnames(lt.cd.mnta)[2] <- "cyear"
  lt.cd.mnta <- as.data.frame(lt.cd.mnta)
  lt.dy3 <- lt.cd.mnta %>%
    group_by(cyear, mnt_area) %>%
    summarise(yngling = sum(yngling, na.rm = T), valper = sum(valper, na.rm = T)) %>%
    as.data.frame()
  
  # fill in missing years with zeros
  # NOTE: input must be a data.frame (NOT a tibble) for 'complete' function to work properly
  lt.dy3 <- lt.dy3 %>%
    tidyr::complete(mnt_area, cyear = full_seq(cyear, 1),
                    fill = list(yngling = 0, valper = 0)) %>%
    as.data.frame()
  returnList$area <- lt.dy3[order(lt.dy3$mnt_area, lt.dy3$cyear),]
  
  cat(" Done.\n")
  
  # export data if flag is TRUE
  if(export) {
    cat("  - [Export] Writing output to file\n")
    data.table::fwrite(returnList$area, file = "output/table_litters_per_year-area.csv", sep = ";", dec = ",", row.names = FALSE, bom = F)
  }
  
  return(returnList)
}