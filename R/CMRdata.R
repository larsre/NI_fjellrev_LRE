#' Prepare DNA data for CMR models
#' @param df a data frame with wrangled DNA data.
#' @export
#' @examples 
#'   CMRdata(df = dna)

CMRdata <- function(df = dna, unique_dens = unique_dens) {
  # set globals
  if(!require(lubridate)) { print("Package 'lubridate' not installed"); break; }
  data_filtered <- FALSE
  
  # inform user
  cat("  - Preparing DNA data for capture histories...")
  
  
  ## Prepare DNA data for capture histories
  
  # extract collection month of sample
  df$month <- month(df$date)
  
  # filter only winter months (January to May), as this is the intensive sampling period
  df <- subset(df, df$month >= 1 & df$month <= 5)
  
  # define 2-month blocks (3 blocks total)
  df$month[df$month==2] = 1
  df$month[df$month==4] = 3
  
  # match DNA samples to nearest den
  df$dist <- "NA"
  df$mnt_area2 <- "NA"
  df$loc_name2 <- "NA"
  df$UTM33E <- as.numeric(df$UTM33E)
  df$UTM33N <- as.numeric(df$UTM33N)
  
  for(i in 1:nrow(df)) {
    check = unique_dens
    check$nordDNA = df$UTM33N[i]
    check$eastDNA = df$UTM33E[i]
    check$diffnord = check$UTM33N - check$nordDNA
    check$diffeast = check$UTM33E - check$eastDNA
    check$distance = ((check$diffnord^2)+(check$diffeast^2))^0.5
    check$distance = round(check$distance,0)
    check = check[order(check$distance),]
    df$dist[i] = check$distance[1]
    df$mnt_area2[i] = check$mnt_area[1]
    df$loc_name2[i] = check$loc_code[1]
  }
  df$dist = as.numeric(df$dist)
  
  # censor all records that are far from den sites (> 5 km; approx. 1.5 % of samples)
  df = df[df$dist<=5000,]
  
  # define regional blocks for abundance estimation, to ensure enough samples
  # NOTE: blocks are defined per sub-population (all municipalities and patches
  # are already assigned to a sub-population)
  df$region <- "NA"
  df <- within(df, {
    region[df$mnt_area2 == "Artfjellet"] <- "Grp2"
    region[df$mnt_area2 == "Blåfjellet"] <- "Grp3"
    region[df$mnt_area2 == "Børgefjell"] <- "Grp3"
    region[df$mnt_area2 == "Finse"] <- "Grp5"
    region[df$mnt_area2 == "Forollhogna"] <- "Grp4"
    region[df$mnt_area2 == "Hardangervidda"] <- "Grp5"
    region[df$mnt_area2 == "Hestkjølen"] <- "Grp3"
    region[df$mnt_area2 == "Ifjordfjellet"] <- "Grp1"
    region[df$mnt_area2 == "Indre Troms"] <- "Grp1"
    region[df$mnt_area2 == "Junkeren"] <- "Grp2"
    region[df$mnt_area2 == "Kjølifjellet/Sylane"] <- "Grp3"
    region[df$mnt_area2 == "Knutshø"] <- "Grp4"
    region[df$mnt_area2 == "Reinheimen"] <- "Grp4"
    region[df$mnt_area2 == "Reisa nord"] <- "Grp1"
    region[df$mnt_area2 == "Reisa sør"] <- "Grp1"
    region[df$mnt_area2 == "Rondane"] <- "Grp4"
    region[df$mnt_area2 == "Saltfjellet"] <- "Grp2"
    region[df$mnt_area2 == "Skjækerfjellet"] <- "Grp3"
    region[df$mnt_area2 == "Snøhetta"] <- "Grp4"
    region[df$mnt_area2 == "Varangerhalvøya"] <- "Grp1"
    region[df$mnt_area2 == "Trollheimen"] <- "Grp4"
    region[df$mnt_area2 == "Anarjohka"] <- "Grp1"
    region[df$mnt_area2 == "Porsanger vest"] <- "Grp1"
    region[df$mnt_area2 == "Sitas"] <- "Grp1"
  })
  
  ## Build encounter histories based on strata
  ch <- subset(df, select = c('region', 'mnt_area2', 'individ', 'year', 'month'))

  # detections
  ch$det = 1
  
  # take the maximum value for each combination of den, year and month
  ch.max <- aggregate(ch$det, by = list(ch$region, ch$mnt_area2, ch$individ, ch$year, ch$month), max)
  colnames(ch.max) <- c('region', 'mtnarea2', 'individ', 'year', 'month', 'maxcode')
  
  # create a cross table with columns for month 
  ch.table = spread(data=ch.max, key=month, value=maxcode, fill=NA)
  
  # replace NA values with dots for not sampled
  ch.table[is.na(ch.table)] <- 0
  colnames(ch.table) = c('region', 'strata', 'individ', 'year', 'janfeb', 'marapr', 'maymay')
  ch = ch.table[, 5:7]
  
  # squish detections into a string
  ch.table$ch = apply(ch,1,paste,collapse="")
  ch.table$freq = 1

  cat(" Done.")
  
  return(ch.table)
}
