#' Calculate 3-year average proportion of Arctic fox litter size per spatial unit
#' @param df a data frame with formatted output data
#' @export
#' @examples 
#'   get.movAvg.litters(df)

get.movAvg.props <- function(df = temp_df) {
  
  # set globals
  returndata <- NULL
  
  # calculate the moving average of the proportion of litter size ('prop_litt')
  for(i in 1:length(unique(df$strata))) {
    movavg.props <- NA
    cur.strata <- sort(unique(df$strata))[i]
    cur.props <- subset(df, df$strata == cur.strata) #subset current strata
    cur.vect <- as.numeric(as.vector(cur.props$prop_litt)) #vectorize
    for(j in 2:length(cur.vect)-1) {
      cur.avg <- mean(c(cur.vect[j-1],cur.vect[j],cur.vect[j+1]), na.rm = T)
      movavg.props[j-1] <- cur.avg
    }
    movavg.props <- round(movavg.props, 3)
    movavg.props <- c(NA,movavg.props,NA)
    cur.props$litt_prop_3yr <- movavg.props
    returndata <- rbind(returndata, cur.props)
  }
  
  # format data
  returndata <- subset(returndata, select = c("strata","year","litt_prop_3yr"))
  returndata <- left_join(df, returndata, by = c("strata","year"))
  
  return(returndata)
}
