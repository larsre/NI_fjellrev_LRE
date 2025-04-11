#' Estimates the population abundance of Arctic fox per spatial unit and year
#' @param df a data frame with formatted output data from CMR models
#' @param frTable a table with the downloaded values for the current indicator
#' from the Nature Index database
#' @param export if TRUE, saves the estimates to a CSV file
#' @param show_existing_values displays the existing NI-values in the NIdb
#' alongside the new estimates
#' @export
#' @examples 
#'   wrangleOutput(df = tidy.output, export = FALSE, show_existing_values = FALSE)

estimateAbundance <- function(df = tidy.output, frTable = fjellrevTable, export = FALSE, show_existing_values = FALSE) {
  
  cat("  Estimating Arctic fox abundance per spatial unit:\n")
  cat("  - Calculating correction factor...")
  
  # aggregate regional means for correction factor estimation
  reg_grouped <- df %>% group_by(region, year) %>% summarise(litt3yr = mean(litt3yr), nhat3yr = mean(nhat3yr))
  
  # set up 5-year blocks for correction factor estimation (horrible code, but so little time...)
  cf_blocks <- 1
  counter <- 1
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
  for(i in 1:length(unique(reg_grouped$year))) {
    if((i != 1) & (!is.wholenumber(i/5))) {
      cf_blocks <- c(cf_blocks, counter)
    } else if((i != 1) & (is.wholenumber(i/5))) {
      cf_blocks <- c(cf_blocks, counter)
      counter <- counter + 1
    }
  }
  cf_blocks <- c(rep(cf_blocks, 5))
  reg_grouped$cf_blocks <- cf_blocks

  # estimate the correction factor (linear model) per unique region and 5-year block
  cf.values <- data.frame()
  for(i in 1:length(unique(reg_grouped$region))) {
    cur.reg.name <- sort(unique(reg_grouped$region))[i]
    cur.reg.df <- subset(reg_grouped, reg_grouped$region == cur.reg.name)
    for(j in 1:length(unique(cur.reg.df$cf_blocks))) {
      # cumulative inclusion of blocks
      cur.block.df <- subset(cur.reg.df, cur.reg.df$cf_blocks <= j)
      # linear model
      cur.mod <- lm(nhat3yr ~ litt3yr, data = cur.block.df)
      # get coefficients from model
      cur.stats <- as.data.frame(summary(cur.mod)$coefficients)
      # get intercept
      cur.ic <- cur.stats$Estimate[1]
      # get intercept SE
      cur.ic.se <- cur.stats$`Std. Error`[1]
      # get correction factor (slope)
      cur.cf <- cur.stats$Estimate[2]
      # get correction factor SE
      cur.cf.se <- cur.stats$`Std. Error`[2]
      temp.val <- data.frame(region = cur.reg.name,
                             block = j,
                             intercept = cur.ic,
                             intercept_se = cur.ic.se,
                             cf = cur.cf,
                             cf_se = cur.cf.se)
      cf.values <- rbind(cf.values, temp.val)
    }
  }
 
  # rearrange input table
  df <- subset(df, select = c("year","strata","litters",
                              "region","litt_prop_3yr","litt3yr",
                              "nhat","nhat3yr","mna"))
  
  # add year-blocks to input table
  block_table <- subset(reg_grouped,
                        select = c("region","year","cf_blocks"))
  df <- left_join(df,
                  block_table,
                  by = c("region","year"))
  
  # add values from correction factor model to input table
  df <- left_join(df,
                  cf.values,
                  by = c("region" = "region",
                         "cf_blocks" = "block"))

  # inform user  
  cat(" Done.\n")
  
  cat("  - Calculating abundance...")
  
  ## calculate abundance per strata based on the number of confirmed litters
  ##   that each spatial unit contributed per year
  ## all parameters except 'litters' are regional estimates
  
  # abundance
  df$abundance <- (df$litt_prop_3yr * df$intercept) +
                          ((df$litters) * df$cf)
  # limit lower bounds to zero
  df$abundance <- ifelse(df$abundance < 0, 0, df$abundance)
  
  # lower confidence limit
  df$abun_lower <- ((df$litt_prop_3yr * df$intercept) -
                            (df$litt_prop_3yr * df$intercept_se * 1.96)) + 
                           ((df$litters * df$cf) - 
                              (df$litters * df$cf_se * 1.96))
  # limit lower bounds to zero
  df$abun_lower <- ifelse(df$abun_lower < 0, 0, df$abun_lower)
  
  # upper confidence limit
  df$abun_upper <- ((df$litt_prop_3yr * df$intercept) + 
                            (df$litt_prop_3yr * df$intercept_se * 1.96)) +
                           ((df$litters * df$cf) + 
                              (df$litters * df$cf_se * 1.96))

  # set the last year (t) to the same value as t-1, as we operate with 3-year moving averages
  last_year <- max(df$year)
  temp_df <- subset(df, df$year == last_year - 1)
  temp_df$year <- last_year
  df <- subset(df, df$year != last_year)
  df <- rbind(df, temp_df)
  df <- df[order(df$strata, df$year),]
  
  # remove first year (NA)
  first_year <- min(df$year)
  df <- subset(df, df$year != first_year)
  
  # inform user  
  cat(" Done.\n")
  
  cat("  - Calculating custom distribution...")
  
  # estimate the custom distribution per strata and year
  # 1. estimate the minimum number alive per spatial unit (strata) and year
  #    NOTE: this is calculated as the regional known number of individuals (mna)
  #    multiplied by the proportional contribution to the population by each
  #    strata in each year (litt_prop_3yr)
  df$mna_strata <- df$litt_prop_3yr * df$mna
  # 2. as there are no real beta estimates on the spatial unit level, we log-
  #    transform the final abundance estimate and use this as logmean for the
  #    custom distribution. Non-transformable estimates (e.g. zeros) are set to NA.
  #    We first extract the necessary data for the transformation.
  df_templog <- subset(df, select = c("year","strata","mna_strata","abundance","abun_lower","abun_upper"))
  # 3. the standard deviation is calculated as the log-transformed distance
  #    between the upper and lower 95 % confidence limits divided by the width
  #    of the 95 % standard errors (2 * 1.96 = 3.92), then multiplied by the
  #    square root of the estimated sample size per spatial unit and year.
  #    Several values for the lower confidence limit has previously been rounded
  #    to 0.0. In these cases, only the upper confidence limit is used
  df_templog$sd <- ifelse(df_templog$abun_lower == 0,
                          (sqrt(df_templog$mna_strata) * (df_templog$abun_upper / 3.92)),
                          (sqrt(df_templog$mna_strata) * ((df_templog$abun_upper - df_templog$abun_lower) / 3.92))
                          )
  # 4. we use the 'normal2Lognormal' function from the NIcalc package to
  #    transform the abundance estimates (mean) and calculated standard deviations.
  df_temp <- NIcalc::normal2Lognormal(mean = df_templog$abundance, 
                                      sd = df_templog$sd)
  df$meanlog <- df_temp$mean
  df$sdlog <- df_temp$sd
  df$meanlog[is.infinite(df$meanlog)] <- NA
  df$meanlog[is.nan(df$meanlog)] <- NA
  df$sdlog[is.infinite(df$sdlog)] <- NA
  df$sdlog[is.nan(df$sdlog)] <- NA
  # 5. non-transformable estimates reflect no estimate on Arctic fox abundance.
  #    Thus, make sure that estimated abundance, lower and upper CI are set to NA.
  df$abundance <- ifelse(is.na(df$meanlog), NA, round(df$abundance, 2))
  df$abun_lower <- ifelse(is.na(df$meanlog), NA, round(df$abun_lower, 2))
  df$abun_upper <- ifelse(is.na(df$meanlog), NA, round(df$abun_upper, 2))
  
  # inform user
  cat(" Done.\n")
  
  ## Finalize
  # if flag is TRUE, add the currently stored values in NIdb alongside the
  # new ones in the export table
  if(show_existing_values) {
    cat("  - show_existing_values flag is set. Adding existing values from NIdb to output table.\n")
    cur_data <- subset(frTable, select = c("yearName","areaName","verdi","nedre_Kvartil","ovre_Kvartil"))
    cur_data <- subset(cur_data, !cur_data$yearName %in% c("1950","1990","2000","Referanseverdi"))
    cur_data$yearName <- as.numeric(cur_data$yearName)
    colnames(cur_data) <- c("year","strata","NI_abun","NI_abun_lower","NI_abun_upper")
    df <- left_join(df, cur_data, by = c("year","strata"))
  }
  
  # if flag is TRUE, export the abundance estimates to a text file
  if(export) {
    cat("  - export flag is set. Writing estimates to 'output/table_estimated_abundance.csv'.\n")
    data.table::fwrite(df,
                       file = "output/table_estimated_abundance.csv",
                       sep = ";",
                       dec = ",",
                       row.names = FALSE,
                       bom = F,
                       encoding = "native")
  }
  
  return(df)
}
