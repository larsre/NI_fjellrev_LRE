#' Wrangle output from models and combine data
#' @param control_data a data frame with wrangled den controls
#' @param dna a data frame with wrangled dna data
#' @param cmrResults a list with output from the closed capture-recapture models
#' @param litterList a list with the number of Arctic fox litters per spatial unit per year
#' @param strata the spatial unit to project population estimates to. Can be
#'  either 'Kommune' (municipality; NI default), 'Patch' (habitat fragment), or
#'  'Area' (management/sub-population).
#' @param export if TRUE, saves the wrangled output to a CSV file
#' @export
#' @examples 
#'   wrangleOutput(control_data, dna, cmrResults, litterList, strata = "Kommune", export = FALSE)

wrangleOutput <- function(control_data, cmrResults, litterList, strata = "Kommune", export = FALSE) {
  
  # inform user of progress
  cat("  Wrangling CMR output data...")
  
  # create a data frame with connections between possible strata
  strata_df <- control_data %>%
    group_by(mnt_area, patch, kom) %>%
    distinct_at(vars(c(mnt_area, patch, kom)), .keep_all = F)
  strata_df <- subset(strata_df, !is.na(strata_df$patch))
  strata_df <- strata_df[order(strata_df$mnt_area, strata_df$kom),]
  
  # add regions
  strata_df$region <- "NA"
  strata_df <- within(strata_df, {
    region[mnt_area == "Artfjellet"] <- "Grp2"
    region[mnt_area == "Blåfjellet"] <- "Grp3"
    region[mnt_area == "Børgefjell"] <- "Grp3"
    region[mnt_area == "Finse"] <- "Grp5"
    region[mnt_area == "Forollhogna"] <- "Grp4"
    region[mnt_area == "Hardangervidda"] <- "Grp5"
    region[mnt_area == "Hestkjølen"] <- "Grp3"
    region[mnt_area == "Ifjordfjellet"] <- "Grp1"
    region[mnt_area == "Indre Troms"] <- "Grp1"
    region[mnt_area == "Junkeren"] <- "Grp2"
    region[mnt_area == "Kjølifjellet/Sylane"] <- "Grp3"
    region[mnt_area == "Knutshø"] <- "Grp4"
    region[mnt_area == "Reinheimen"] <- "Grp4"
    region[mnt_area == "Reisa nord"] <- "Grp1"
    region[mnt_area == "Reisa sør"] <- "Grp1"
    region[mnt_area == "Rondane"] <- "Grp4"
    region[mnt_area == "Saltfjellet"] <- "Grp2"
    region[mnt_area == "Skjækerfjellet"] <- "Grp3"
    region[mnt_area == "Snøhetta"] <- "Grp4"
    region[mnt_area == "Varangerhalvøya"] <- "Grp1"
    region[mnt_area == "Trollheimen"] <- "Grp4"
    region[mnt_area == "Anarjohka"] <- "Grp1"
    region[mnt_area == "Porsanger vest"] <- "Grp1"
    region[mnt_area == "Sitas"] <- "Grp1"
    region[mnt_area == "Andre områder Sør Norge (Nord)"] <- "Grp4"
    region[mnt_area == "Andre områder Sør Norge (Sør)"] <- "Grp5"
  })
  strata_df <- as.data.frame(strata_df)
  strata_df <- subset(strata_df, !strata_df$region=="NA")
  
  # use litter size table as a base for tidy data
  temp_df <- data.frame()
  if(identical("Kommune", strata)) { #default
    temp_df <- litterList$municipality
    temp_df <- subset(temp_df, select = c("cyear","kom","yngling"))
    colnames(temp_df) <- c("year","kom","litters")
    temp_strata <- strata_df %>% distinct_at(vars(c(kom,region)), .keep_all = F)
    
    # remove one municipality that appears in two regions
    temp_strata <- subset(temp_strata, !(temp_strata$kom=="Holtålen" & temp_strata$region=="Grp4"))
    temp_df <- left_join(temp_df, temp_strata, by = "kom")
    colnames(temp_df)[2] <- "strata"
  }
  
  if(identical("Patch", strata)) {
    temp_df <- litterList$patch
    temp_df <- subset(temp_df, select = c("cyear","patch","yngling"))
    colnames(temp_df) <- c("year","patch","litters")
    temp_strata <- strata_df %>% distinct_at(vars(c(patch,region)), .keep_all = F)
    temp_df <- left_join(temp_df, temp_strata, by = "patch")
    colnames(temp_df)[2] <- "strata"
  }
  
  if(identical("Area", strata)) {
    temp_df <- litterList$area
    temp_df <- subset(temp_df, select = c("cyear","mnt_area","yngling"))
    colnames(temp_df) <- c("year","mnt_area","litters")
    temp_strata <- strata_df %>% distinct_at(vars(c(mnt_area,region)), .keep_all = F)
    temp_df <- left_join(temp_df, temp_strata, by = "mnt_area")
    colnames(temp_df)[2] <- "strata"
  }
  
  # calculate the sum of litters per spatial unit
  temp_litt <- temp_df %>% group_by(region, year) %>% summarise(litt_reg_sum = sum(litters, na.rm = T))
  temp_df <- left_join(temp_df, temp_litt, by = c("region","year"))
  
  # calculate the proportion of litters per spatial unit out of the regional total per year
  temp_df$prop_litt <- round(temp_df$litters / temp_df$litt_reg_sum, 3)
  temp_df$prop_litt <- ifelse(is.nan(temp_df$prop_litt), 0, temp_df$prop_litt) # for new calculation

  # calculate 3-year moving average of spatial unit litter proportion
  temp_df <- get.movAvg.props(temp_df)
  
  # calculate 3-year moving average of regional litter sums
  temp_df <- get.movAvg.litters(temp_df)
  
  
  ## organize model estimates
  temp_modres <- data.frame()
  for(i in 1:length(cmrResults)) {
    cur.name <- names(cmrResults)[i] #regional group
    cur.years <- as.numeric(as.character(cmrResults[[i]]$years)) #years
    cur.nhat <- cmrResults[[i]]$estimated.popsize$estimate #estimated regional population size
    cur.mna <- cmrResults[[i]]$minimum.number.alive #minimum number alive (known individuals) per region
    
    temptempmod <- data.frame(year = cur.years,
                              region = cur.name,
                              mna = as.vector(cur.mna),
                              nhat = cur.nhat)
    temp_modres <- rbind(temp_modres, temptempmod)
  }
  
  # calculate 3-year moving average of regional population estimates (nhat)
  temp_modres <- get.movAvg.nhat(temp_modres)
  temp_df <- left_join(temp_df, temp_modres, by = c("year","region"))
  
  # inform user
  cat(" Done.\n")
  
  # export data if flag is TRUE
  if(export) {
    cat("  - Export flag is set. Writing wrangled table to 'output/table_pre-estimation.csv'.\n")
    data.table::fwrite(temp_df, file = "output/table_pre-estimation.csv", sep = ";", dec = ",", row.names = FALSE, bom = F)
  }
  
  return(temp_df)
}