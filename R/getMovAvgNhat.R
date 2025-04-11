#' Calculate 3-year average of estimated Arctic fox abundance per strata
#' @param df a data frame with formatted output data
#' @export
#' @examples 
#'   get.movAvg.nhat(df)

get.movAvg.nhat <- function(df = temp_modres) {
  
  # iterate over annual population estimates per region to calculate 3-year moving averages
  returndata <- NULL
  for(i in 1:length(unique(df$region))) {
    movavg.nhat <- NA
    cur.region <- sort(unique(df$region))[i]
    cur.nhat <- subset(df, df$region == cur.region) #subset current region
    cur.vect <- as.numeric(as.vector(cur.nhat$nhat)) #vectorize
    for(j in 2:length(cur.vect)-1) {
      cur.avg <- mean(c(cur.vect[j-1],cur.vect[j],cur.vect[j+1]), na.rm = T)
      movavg.nhat[j-1] <- cur.avg
    }
    movavg.nhat <- round(movavg.nhat, 1)
    movavg.nhat <- c(NA,movavg.nhat,NA)
    cur.nhat$nhat3yr <- movavg.nhat
    returndata <- rbind(returndata, cur.nhat)
  }
  
  # format data
  returndata <- subset(returndata, select = c("region","year","nhat3yr"))
  returndata <- left_join(df, returndata, by = c("region","year"))
  
  return(returndata)
}