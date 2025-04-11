#' Format den control data export from ROVBASE.
#' @param df a data frame with the complete and untouched records from ROVBASE (all years).
#' @export
#' @examples 
#'   wrangleControls(df = df)

wrangleControls <- function(df, min.year, max.year, refs) {
  if(!require(dplyr)) { print("Package 'dplyr' not installed."); break; }
  
  # check that data is provided
  if(is.null(df) | nrow(df) == 0 | missing(df)) {
    stop("\nData frame is empty, NULL or not provided.")
  }
  
  # inform user of progress
  cat("  Wrangling den control data:\n")
  cat("  - Organizing data by control year (October 1st to September 30th)...")
  
  # format date/time
  df$control_date <- parse_date_time(df$Kontrolldato, c("%Y-%m-%d"), exact = T)
  
  # define the control year from October 1st to September 30th
  df$control_year <- year(df$control_date)
  df$kontrollaar <- NA
  for(i in min(df$control_year, na.rm = T)+1:max(df$control_year, na.rm = T)) {
    df$kontrollaar[df$control_year == (i-1) & month(df$control_date) > 9] <- i
    df$kontrollaar[df$control_year == i & month(df$control_date) < 10] <- i
  }
  df$control_year <- df$kontrollaar
  df <- df %>% select(-kontrollaar)
  
  # control date is missing for a few den controls - set control_year = År
  df$control_year <- ifelse(is.na(df$control_year) & !is.na(df$År), df$År, df$control_year)
  cat(" Done.\n")
  
  cat("  - Formatting and filtering data...")
  # sort descending by control date (last control date first)
  df <- df[order(df$control_date, decreasing = T),]
  
  # extract den location code (for later reference/join)
  df$lok_kode <- sapply(df$Lokalitet, function(x) strsplit(x, split = " ", fixed = T)[[1]][1])
  
  # format municipality codes
  df$Komnr <- sapply(df$Kommunenummer, function(x) strsplit(x, split = "-", fixed = T)[[1]][2])
  
  # remove space before/after municipality codes
  df$Komnr <- as.integer(trimws(df$Komnr))
  
  # format municipality names
  df$Kom <- sapply(df$Kommune, function(x) strsplit(x, split = " (N)", fixed = T)[[1]][1])
  
  # format county names
  df$Fylk <- sapply(df$Fylke, function(x) strsplit(x, split = " (N)", fixed = T)[[1]][1])
  
  # add country
  temp.country <- strsplit(df$Fylke, " ", fixed = T)
  temp.country <- unlist(lapply(temp.country, tail, n = 1L))
  df$Land <- ifelse(temp.country == "(N)", "Norway",
                             ifelse(temp.country == "(S)", "Sweden", "Finland"))
  
  # filter only Norwegian dens
  df <- subset(df, df$Land == "Norway")
  
  # add region (for use in abundance estimation)
  df$region <- NA
  df <- within(df, {
    region[Fjellområde == "Finse" | Fjellområde == "Hardangervidda" | Fjellområde == "Andre områder Sør Norge (Sør)" | Fjellområde == "Setesdalsheiene"] <- "Sør-Norge Sør"
    region[Fjellområde == "Ottadalen nord" | Fjellområde == "Snøhetta" | Fjellområde == "Knutshø" | Fjellområde == "Forollhogna" | Fjellområde == "Trollheimen" | Fjellområde == "Rondane" | Fjellområde == "Andre områder Sør Norge (Nord)"] <- "Sør-Norge Nord"
    region[Fjellområde == "Blåfjellet" | Fjellområde == "Hestkjølen" | Fjellområde == "Skjækerfjellet" | Fjellområde == "Børgefjell" | Fjellområde == "Kjølifjellet/Sylane"] <- "Midt-Norge"
    region[Fjellområde == "Artfjellet" | Fjellområde == "Junkeren" | Fjellområde == "Saltfjellet"] <- "Nordland"
    region[Fjellområde == "Indre Troms" | Fjellområde == "Reisa sør" | Fjellområde == "Reisa nord" | Fjellområde == "Ifjordfjellet" | Fjellområde == "Varangerhalvøya" | Fjellområde == "Porsanger vest" | Fjellområde == "Sitas" | Fjellområde == "Anarjohka" | Fjellområde == "Andre områder"] <- "Nord-Norge"
  })
  
  # nomenclature: change "Ottadalen nord" to "Reinheimen"
  df$Fjellområde[df$Fjellområde == "Ottadalen nord"] <- "Reinheimen"
  
  
  ## Filter columns
  df2 <- subset(df, select = c('RovbaseID', 'control_year', 'lok_kode', 'Lokalitet', 'Fjellområde', 'Komnr', 'Kom', 'Fylk', 'region', 'Art', 'Dokumentert', 'Antatt sikker', 'Valper sett', 'Øst (UTM33/SWEREF99 TM)', 'Nord (UTM33/SWEREF99 TM)'))
  
  # change column names
  colnames(df2) <- c('rbID', 'cyear', 'loc_code', 'loc_name', 'mntarea', 'mun_no', 'mun', 'county', 'region', 'species', 'dok', 'ant', 'pups', 'UTM33E', 'UTM33N')
  
  # drop if missing information for a den locality
  df3 <- df2[is.na(df2$loc_code)!="TRUE",]
  
  # sort by den location and descending control year
  df3$cyear <- as.integer(df3$cyear)
  df3 <- df3[with(df3, order(loc_code, -cyear)), ]
  cat(" Done.\n")
  
  # filter unique dens (latest control) and control year
  # NOTE: use max(table) for species, matching 'Fjellrev'
  cat("  - Grouping data by den location and control year...")
  
  # convert character to integer
  df3$pups <- as.integer(df3$pups)
  df3$UTM33E <- as.integer(df3$UTM33E)
  df3$UTM33N <- as.integer(df3$UTM33N)
  
  # remove controls without UTM coordinates
  df3 <- subset(df3, !is.na(df3$UTM33E))

  # group data by den location and control year
  df4 <- df3 %>%
    group_by(loc_code, cyear) %>%
    mutate(rbID = names(which.max(table(rbID))),
           loc_name = names(which.max(table(loc_name))),
           mntarea = names(which.max(table(mntarea))),
           mun_no = max(mun_no),
           mun = names(which.max(table(mun))),
           county = names(which.max(table(county))),
           region = names(which.max(table(region))),
           species = ifelse(all(is.na(species)), names(which.max(table(species, useNA = "ifany"))), names(which.max(table(species)))),
           yngling = ifelse(any(dok == "Ja"), 1, 0),
           pups = ifelse(length(pups) == 1 & is.na(pups), NA, max(pups, na.rm = T)),
           UTM33E = max(UTM33E, na.rm = T),
           UTM33N = max(UTM33N, na.rm = T)) %>%
    ungroup() %>%
    distinct_at(vars(c(loc_code, cyear)), .keep_all = T)
  cat(" Done.\n")
  
  # filter data by min.year and max.year
  usable_minyear <- NA
  usable_maxyear <- NA
  if(!missing(min.year) & !is.na(min.year) & is.numeric(min.year) & (min.year > 1900) & (min.year < 2100)) {
    usable_minyear <- min.year
  } else {
    cat(paste0("  min.year is either not set or not in a valid format (", min.year, "). Including all years.\n"))
  }
  
  if(!missing(max.year) & !is.na(max.year) & is.numeric(max.year) & (max.year > 1900) & (max.year < 2100)) {
    usable_maxyear <- max.year
  } else {
    cat(paste0("  max.year is either not set or not in a valid format (", max.year, "). Including all years.\n"))
  }
  
  if(!is.na(usable_minyear) & !is.na(usable_maxyear)) {
    df4 <- subset(df4, df4$cyear >= usable_minyear & df4$cyear <= usable_maxyear)
  }
  
  cat("  - Conforming municipalities to NIdb...")
  
  # sort data
  df <- subset(df4, select = c('rbID','cyear','loc_code','loc_name','mntarea','mun_no','mun','county','region','species','UTM33E','UTM33N','yngling','pups'))
  colnames(df) <- c('rovbaseID','year','loc_code','loc_name','mnt_area','kom_no','kom','fylke','region','species','UTM33E','UTM33N','yngling','valper')
  
  # conform municipality names to the Nature Index database
  df <- within(df, {
    kom[kom == "Deatnu - Tana"] <- "Deatnu Tana"
    kom[kom == "Gáivuotna - Kåfjord - Kaivuono"] <- "Gáivuotna Kåfjord"
    kom[kom == "Guovdageaidnu - Kautokeino"] <- "Guovdageaidnu Kautokeino"
    kom[kom == "Hammerfest - Hámmerfeasta"] <- "Hammerfest"
    kom[kom == "Kárášjohka - Karasjok"] <- "Karasjohka Karasjok"
    kom[kom == "Nordreisa - Ráisa - Raisi"] <- "Nordreisa"
    kom[kom == "Os"] <- "Os (Hedmark)"
    kom[kom == "Porsanger - Porsá?gu - Porsanki"] <- "Porsanger Porsángu Porsanki"
    kom[kom == "Porsanger - Porsángu - Porsanki"] <- "Porsanger Porsángu Porsanki"
    kom[kom == "Rosse - Røros"] <- "Røros"
    kom[kom == "Raarvihke - Røyrvik"] <- "Røyrvik"
    kom[kom == "Snåase - Snåsa"] <- "Snåsa"
    kom[kom == "Storfjord - Omasvuotna - Omasvuono"] <- "Storfjord"
    kom[kom == "Unjárga - Nesseby"] <- "Unjargga Nesseby"
    kom[kom == "Aarborte - Hattfjelldal"]  <- "Hattfjelldal"
  })
  
  # remove controls not in any NIdb municipalities
  df <- df[df$kom %in% refs$Area.Name,]
  
  # assign NIdb area ID's
  aID <- refs[,c(1:2)]
  df <- left_join(df, aID, by = c("kom" = "Area.Name"))
  
  cat(" Done.\n")
  
  return(df)
}