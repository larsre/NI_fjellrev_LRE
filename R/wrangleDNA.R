#' Format DNA data export from ROVBASE.
#' @param df a data frame with the complete and untouched records from ROVBASE (all years).
#' @export
#' @examples 
#'   wrangleDNA(df = df)

wrangleDNA <- function(df, min.year, max.year, refs) {
  if(!require(dplyr)) { print("Package 'dplyr' not installed."); break; }
  
  # check that data is provided
  if(is.null(df) | nrow(df) == 0 | missing(df)) {
    stop("\nData frame is empty, NULL or not provided.")
  }
  
  # inform the user
  cat("  Wrangling DNA sample data:\n")
  cat("  - Organizing data by control year (October 1st to September 30th)...")
  
  # set proper encoding
  Encoding(df$Fjellområde) <- "UTF-8"
  
  # remove rows where sampling date is missing (only 3 older samples)
  df <- subset(df, !is.na(df$Funnetdato))
  
  # define the control year from October 1st to September 30th
  df$control_year <- year(df$Funnetdato)
  df$kontrollaar <- NA
  for(i in min(df$control_year, na.rm = T)+1:max(df$control_year, na.rm = T)) {
    df$kontrollaar[df$control_year == (i-1) & month(df$Funnetdato) > 9] <- i
    df$kontrollaar[df$control_year == i & month(df$Funnetdato) < 10] <- i
  }
  df$control_year <- df$kontrollaar
  df <- df %>% select(-kontrollaar)

  # keep only winter samples (between October 1st and May 31st)
  df <- subset(df, month(df$Funnetdato) > 9 | month(df$Funnetdato) < 6)
  
  # keep only data from years specified
  df <- subset(df, df$control_year >= min.year & df$control_year <= max.year)
  
  cat(" Done.\n")
  
  cat("  - Removing species not identified as Arctic fox...")
  
  # keep only rows where "Art (Analyse)" = Fjellrev
  df <- subset(df, df$`Art (Analyse)` == "Fjellrev")
  
  cat(" Done.\n")
  
  cat("  - Removing samples not identified to Arctic fox individual...")
  
  # keep only rows with identified individuals
  df <- subset(df,!is.na(df$Individ))
  
  # fix nomenclature of individuals - use only fox name code
  df$Individ <- sapply(df$Individ, function(x) strsplit(x, split = " ", fixed = T)[[1]][1])
  
  # remove trailing and leading white spaces from fox name code
  df$Individ <- trimws(df$Individ)
  
  cat(" Done.\n")
  
  cat("  - Removing non-Norwegian samples...")
  
  # remove non-Norwegian data samples
  nonNor <- c(
    "Dalarnas län (S)",
    "Gotlands län (S)",
    "Jämtlands län (S)",
    "Lapin lääni (F)",
    "Norrbottens län (S)",
    "Västerbottens län (S)"
  )
  df <- df[!df$Fylke %in% nonNor, ]
  
  cat(" Done.\n")
  
  cat("  - Conforming municipalities to NIdb...")
  
  # format municipality codes
  df$Komnr <- sapply(df$Kommunenummer, function(x) strsplit(x, split = "-", fixed = T)[[1]][2])
  
  # remove space before/after municipality codes
  df$Komnr <- as.integer(trimws(df$Komnr))
  
  # format municipality names
  df$Kom <- sapply(df$Kommune, function(x) strsplit(x, split = " (N)", fixed = T)[[1]][1])
  
  # format county names
  df$Fylk <- sapply(df$Fylke, function(x) strsplit(x, split = " (N)", fixed = T)[[1]][1])
  
  # conform municipality names to the Nature Index database
  df <- within(df, {
    Kom[Kom == "Deatnu - Tana"] <- "Deatnu Tana"
    Kom[Kom == "Guovdageaidnu - Kautokeino"] <- "Guovdageaidnu Kautokeino"
    Kom[Kom == "Os"] <- "Os (Hedmark)"
    Kom[Kom == "Rosse - Røros"] <- "Røros"
    Kom[Kom == "Raarvihke - Røyrvik"] <- "Røyrvik"
    Kom[Kom == "Storfjord - Omasvuotna - Omasvuono"] <- "Storfjord"
    Kom[Kom == "Unjárga - Nesseby"] <- "Unjargga Nesseby"
    Kom[Kom == "Aarborte - Hattfjelldal"]  <- "Hattfjelldal"
  })
  
  # remove samples not in any NIdb municipalities
  df <- df[df$Kom %in% refs$Area.Name,]
  
  # assign NIdb area ID's
  aID <- refs[,c(1:2)]
  df <- left_join(df, aID, by = c("Kom" = "Area.Name"))
  
  cat(" Done.\n")
  
  # return only relevant variables
  df <- subset(df, select = c("Strekkode (Prøve)", "Funnetdato", "control_year", "Individ",
                              "Øst (UTM33/SWEREF99 TM)","Nord (UTM33/SWEREF99 TM)", "Fjellområde",
                              "Lokalitet","Komnr","Kom","Fylk","Area_ID"))
  colnames(df) <- c("sampleID","date","year","individ","UTM33E","UTM33N","mnt_area","loc_name","kom_no","kom","fylke","Area_ID")
  
  return(df)
}
