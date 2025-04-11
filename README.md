# NI_fjellrev_LRE
Workflow for the estimation of Arctic fox abundance per municipality for the Nature Index for Norway

## Preface
This repository contains the scripts and (some) data for estimating Arctic fox abundance using a closed capture-recapture model based on DNA sampling, and for further calculating municipality-level abundances for entry into the Nature Index for Norway database (NIdb). DNA samples are collected through the Arctic fox monitoring programme (Ulvund et al. 2023; NINA Report 2344). The base model was developed by Brett Sandercock in 2018-19, and estimates from the model are updated annually by the Arctic Fox Research Group to inform management in both Norway, Sweden and Finland.

The model and workflow has been prepared with functionality to estimate annual population abundance of Arctic fox on several spatial scales. In addition to municipality level (NIdb default), abundance can be estimated on the management level (sub-population) and habitat level (natural mountain fragments).

All functions (in 'R' folder) have been conformed to the roxygen2 commenting standard, and function steps explained. See the Appendix for a list of functions and their purpose. This 'Readme' detailes the workflow of the R-script.

## Switches
Limit data by time interval (min year to max year). Estimates are provided as 3-year moving averages, so the first and last
year will be censored. The last year will be set equal to last year-1.
```{r switches, echo=TRUE, message=FALSE}
min.year <- 2007
max.year <- 2024
```

## Import data
Data is downloaded from ROVBASE (Norwegian Environment Agency). There is (currently) no API available for ROVBASE, thus data must be downloaded via the web interface by one of the very few people with access. The raw data is then imported here. Necessary tables from ROVBASE are the Arctic fox den controls (complete record) and Arctic fox DNA samples (complete record).
```{r import_data, echo=TRUE, message=FALSE}
# DEN CONTROLS
inputControls <- as.data.frame(
  read_excel("data/eksport_NSF_kontroller_rovbase_2024.xlsx",
                                          sheet = "Rapport",
                                          trim_ws = T,
                                          col_types = c(rep("text",51),
                                                        "date",
                                                        rep("text",34)))
  )

# DNA DATA
inputDNA <- as.data.frame(
  read_excel("data/eksport_NSF_DNA_rovbase_2024.xlsx",
                                     sheet = "Rapport",
                                     trim_ws = T,
                                     col_types = c(rep("text",7),
                                                   "date",
                                                   rep("text",51)))
  )

# Download existing data from the Nature Index database
# NOTE: Requires login with username and password obtained from NIdb administrators
indicators <- c("Fjellrev")
currentData <- downloadData_NIdb(indicators = indicators, save = F, save_path = getwd())
fjellrevTable <- currentData$Fjellrev$indicatorValues

# create a NIdb reference table for area ID's and names (for data filtering)
ni_refs <- subset(fjellrevTable, fjellrevTable$yearName == "Referanseverdi", select = c("areaId","areaName"))
colnames(ni_refs) <- c("Area_ID","Area.Name")
```

## Wrangle data
Data must be properly formatted and structured. Data sets are limited to the min and max years set through the switches above. Proper area ID codes from NIdb are added to each municipality. Only area ID codes with at least one controlled Arctic fox den during the time interval (min year to max year) are retained.

Data on documented litters per den per year are derived from the den control data. This is used to estimate correction factors between estimated population abundance and the number of litters on a regional level. This factor is then used to calculate the abundance on smaller scales (municipality, sub-population or habitat).
```{r wrangle_data, echo=TRUE, message=FALSE}
## DEN CONTROLS
# control data
control_data <- wrangleControls(inputControls, min.year, max.year, ni_refs)

# unique dens (based on last control year)
unique_dens <- uniqueDens(control_data)

# assign patch number (habitat level) to control data
control_data <- getPatch(control_data)

## DNA
# dna data
dna <- wrangleDNA(inputDNA, min.year, max.year, ni_refs)

# assign patch number (habitat level) to DNA data
dna <- getPatch(dna)

## LITTER SIZE
# calculate the number of litters per strata (municipality, habitat patch and mountain area) per year
litterList <- litterSize(df = control_data, export = FALSE)
```

## Prepare and run the closed capture-recapture model
DNA data is formatted into presence/absence (encounter histories) per individual fox and year. Only foxes identified to individual by DNA are included, and limited to the more intensive winter sampling period (January to May). The data is then stratified by region (5 regional levels covering Norway). This is to ensure that the model has enough data to provide annual estimates with a tolerable level of precision.

Estimation is performed through a closed capture-recapture population model in RMark, utilizing only 'year' (time) as a covariate in the abundance and observation processes. The model is based on the Huggins (1989, 1991) closed capture model, providing derived estimates on population abundance based on the relationship between known and unknown individuals in the population. See https://sites.warnercnr.colostate.edu/gwhite/huggins-closed-captures-models/ for an in-depth explanation of the model. The model is run separately for each region, and the derived population estimates are retained for further calculations. Note that only 95% confidence intervals are available for the derived estimates, computed with a lognormal distribution.
```{r model_run, echo=TRUE, message=FALSE}
# Prepare data for estimating abundance
cmrData <- CMRdata(dna)

# Set start time
t.start <- Sys.time()

# Run CMR models to estimate population abundance
cmrResults <- capture.models.dna(cmrData)

# Total run time
Sys.time() - t.start

# Remove temporary files
cleanup(ask = F)
```

## Estimating Arctic fox abundance
The model results are combined into a tidy data frame. The spatial scale of the estimates must be defined here, using the 'strata' variable, either on municipality ("Kommune"), sub-population ("Area"), or habitat ("Patch") levels. This step also calculates the 3-year moving averages of contribution (proportion) of litters from each strata (e.g. municipality) as well as the regional litter averages.

The final population abundance per spatial scale is then calculated. In brief, the correction factor - i.e. linear relationship between the 3-year moving average of documented litters and 3-year moving average of estimated population abundance per region - is estimated, and then used to calculate the population abundance per strata. This is calculated as:
(the 3-year moving average proportion of litter contribution by the strata * intercept of the linear litter/abundance relationship on a regional level) + (number of litters contributed in a specific year by the strata * slope of the linear litter/abundance relationship on a regional level). See the 'estimateAbundance.R' function for further details.

Note that the last year of estimates is set to 'NA' due to the 3-year moving averages used in calculations. The estimateAbundance.R function changes this value to have the same value as the previous year.
```{r estimate_abundance, echo=TRUE, message=FALSE}
# Combine data and model output into a tidy list
tidy.output <- wrangleOutput(control_data,
                             cmrResults,
                             litterList,
                             strata = "Kommune",
                             export = FALSE)

# Estimate Arctic fox abundance using regional estimates,
# correction factors and litter size
est.abundance <- estimateAbundance(tidy.output,
								   fjellrevTable = fjellrevTable,
                                   export = FALSE,
                                   show_existing_values = FALSE)
```

## Data output
The main output file includes all data from the steps involved in calculating the final abundance estimates, including the regional levels, the 3-year moving averages of the number of litters (per strata and region), estimated regional abundances, the minimum number of foxes alive per region (mna), intercept and slope from the linear correction factor estimation, and the final calculated abundances per strata and year. A compacted example of the final data that is to be updated is shown below. In addition, the log-transformed estimates and standard deviance are calculated and added as a custom distribution function to each municipality and year (see estimateAbundance.R for further details).
```{r data_output, echo=TRUE, message=FALSE}
# Compact output for markdown presentation
present <- subset(est.abundance, select = c("year",
                                            "strata",
                                            "abundance",
                                            "abun_lower",
                                            "abun_upper"))

# Number of strata with estimated abundances
length(table(est.abundance$strata))

# Example of the final data set (compacted)
present[present$strata=="Dovre",]
```

## Data upload
The curated data is then conformed to the NIdb data structure and uploaded to the NI database.
```{r data_upload, echo=TRUE, message=FALSE}
# check for spatial unit/strata nomenclature consistency
all(unique(est.abundance$strata) %in% unique(fjellrevTable$areaName)) #should be TRUE

# make a copy of the current indicator value table for editing
editable <- currentData$Fjellrev

# save the copy of the current indicator value table
saveRDS(editable, "data/oldIndicatorValues_2025.rds")

# update the copy of the indicator table
for(i in 1:nrow(editable$indicatorValues)) {
  cur.row <- editable$indicatorValues[i,]
  cur.areaId <- cur.row$areaId
  cur.year <- cur.row$yearName
  
  # find the matching row in the estimated abundance table
  new.data <- subset(est.abundance, as.character(est.abundance$year) == cur.row$yearName & 
                       est.abundance$strata == cur.row$areaName)
  
  # if match = update
  if(nrow(new.data) > 0) {
    # make custom distribution if possible
    if(!is.na(new.data$meanlog)) {
      distObj <- NIcalc::makeDistribution(input = "logNormal",
                                          distParams = list(mean = new.data$meanlog,
                                                            sd = new.data$sdlog))

      # update table
      editable <- NIcalc::setIndicatorValues(indicatorData = editable,
                                             areaId = cur.areaId,
                                             years = as.integer(cur.year),
                                             est = new.data$abundance,
                                             lower = new.data$abun_lower,
                                             upper = new.data$abun_upper,
                                             distribution = distObj,
                                             datatype = 3,
                                             unitOfMeasurement = "Antall voksne/reproduserende individer")
    }
  }
}

# upload data
updatedIndicatorData <- currentData
updatedIndicatorData$Fjellrev <- editable
uploadData_NIdb(indicatorData = updatedIndicatorData)
```

## Appendix - Main functions
**wrangleControls.R** - Formats Arctic fox den control data imported from ROVBASE.  
**uniqueDens.R** - Creates a list of unique Arctic fox dens. Used by CMRdata.R.  
**getPatch.R** - Intersects den control and DNA geographic coordinates with habitat fragment polygons, for potential future estimation of abundance on habitat level.  
**wrangleDNA.R** - Formats Arctic fox DNA sample data imported from ROVBASE.  
**litterSize.R** - Calculate the number of documented Arctic fox litters for all viable spatial units.  
**CMRdata.R** - Format DNA data into yearly encounter histories for Arctic foxes identified from DNA. Calls uniqueDens.R.  
**CaptureModels.R** - Runs the CMR (capture-mark-recapture) model (capture.models.dna) for each regional level, to estimate Arctic fox abundance.  
**wrangleOutput.R** - Formats the output of the CMR models, and conforms the data into a tidy frame. Calls getMovAvgProps.R and getMovAvgLitters.R.  
**getMovAvgProps.R** - Calculates the 3-year moving average proportion of Arctic fox litter size per spatial unit.  
**getMovAvgLitters.R** - Calculates the 3-year moving average Arctic fox litter size per region.  
**getMovAvgNhat.R** - Calculates the 3-year moving average of estimated Arctic fox abundance per region.  
**estimateAbundance.R** - Calculates 3-year moving averages, and calculates abundance estimates for the specified spatial scale. Sorts and conforms data.  
**downloadData_NIdb.R** - Function to download specific indicator data from the Nature Index database.  
**uploadData_NIdb.R** - Function to upload updated indicator data to the Nature Index database.  
