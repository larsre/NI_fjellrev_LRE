#' Calculate the 3-year average of Arctic fox litter size per region
#' @param df a data frame with formatted output data
#' @export
#' @examples 
#'   get.movAvg.litters(df)

get.movAvg.litters <- function(df = temp_df) {
  
  # get the unique sum of litters per year per region
  unique.region <- df %>% group_by(region, year) %>% summarise(litt_reg_sum = mean(litt_reg_sum))
  
  # iterate unique sum of litters to calculate 3-year moving averages per region
  returndata <- NULL
  for(i in 1:length(unique(unique.region$region))) {
    movavg.litt <- NA
    cur.region <- sort(unique(unique.region$region))[i]
    cur.litt <- subset(unique.region, unique.region$region == cur.region) #subset current region
    cur.vect <- as.numeric(as.vector(cur.litt$litt_reg_sum)) #vectorize
    for(j in 2:length(cur.vect)-1) {
      cur.avg <- mean(c(cur.vect[j-1],cur.vect[j],cur.vect[j+1]), na.rm = T)
      movavg.litt[j-1] <- cur.avg
    }
    movavg.litt <- round(movavg.litt, 1)
    movavg.litt <- c(NA,movavg.litt,NA)
    cur.litt$litt3yr <- movavg.litt
    returndata <- rbind(returndata, cur.litt)
  }
  
  # format data
  returndata <- subset(returndata, select = c("region","year","litt3yr"))
  returndata <- left_join(df, returndata, by = c("region","year"))
  
  return(returndata)
}
