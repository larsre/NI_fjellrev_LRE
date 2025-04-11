#' Run CMR models
#' @param df a data frame with encounter histories for CMR analysis
#' @export
#' @examples 
#'   capture.models.dna(df)

#function to run models per block
capture.models.dna <- function(df = data.frame()) {
  # set globals
  if(!require(RMark)) { print("Package 'RMark' not installed"); break; }
  returnlist <- list()
  
  # inform user
  cat("  |- Running DNA-based capture-recapture models for strata...\n")

  # run per region
  for(i in 1:length(unique(df$region))) {
    
    # subset current region
    cur.df <- subset(df, df$region == sort(unique(df$region))[i])
    cur.df$individ <- factor(cur.df$individ)
    cur.df$year <- factor(cur.df$year)

    temp.individs <- length(levels(cur.df$individ))
    temp.foxyears <- dim(cur.df)[1]

    #process data
    df.proc = process.data(cur.df, model = 'Closed', groups = c("year"))

    #closed model function
    run.closed = function() {
      
      #define parameter models
      pdotshared = list(formula = ~1, share = TRUE)
      ptimeshared = list(formula = ~time, share = TRUE)
      f0dot = list(formula = ~1)
      f0year = list(formula = ~year)

      #capture closed models
      closed.M0 = mark(df.proc, model = "Closed", model.parameters = list(p = pdotshared, f0 = f0year), brief = T, )
      closed.Mt = mark(df.proc, model = "Closed", model.parameters = list(p = ptimeshared, f0 = f0year), brief = T)

      #return model table and list
      return(collect.models())
    }

    # run models
    model.results = run.closed()
    
    # derived N (population abundance) per year on regional level
    temp.N <- model.results$closed.Mt$results$derived$`N Population Size`
    
    # unique foxes per year (i.e. minimum number alive)
    temp.mna <- table(cur.df$year)
    
    # sample size
    temp.sN <- model.results$closed.Mt$results$n
    
    # detection
    temp.p <- model.results$closed.Mt$results$real
    
    # formatted estimates of N
    temp.N.est <- round(as.numeric(temp.N$estimate), 1)
    
    # median conversion factor for estimated vs. minimum population size
    temp.Nmna.median <- round(median(as.numeric(temp.N$estimate/temp.mna)),2)
    
    # 95% confidence limits of conversion factor
    temp.Nmna.quant <- round(quantile(as.numeric(temp.N$estimate/temp.mna), c(0.025,0.975)),2)
    
    # beta values
    temp.beta <- model.results$closed.Mt$results$beta
    
    temp.returnlist <- list("years" = sort(unique(cur.df$year)),
                            "total.individs" = temp.individs,
                            "total.foxyears" = temp.foxyears,
                            "sample.size" = temp.sN,
                            "estimated.popsize" = temp.N,
                            "minimum.number.alive" = temp.mna,
                            "detection" = temp.p,
                            "est.formatted" = temp.N.est,
                            "median.N.mna" = temp.Nmna.median,
                            "quantile.N.mna" = temp.Nmna.quant,
                            "beta" = temp.beta)

    returnlist[[i]] <- temp.returnlist
    names(returnlist)[[i]] <- cur.df$region[1]
  }

  return(returnlist)
}