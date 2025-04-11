# ------------------------------------------------------------ #
# Naturindeks Fjellrev - Lars Rød-Eriksen og Nina E. Eide 2025 #
#                                                              #
# This document details the workflow in estimating the nature  #
# index values for Arctic foxes (abundance) using a closed     #
# capture-recapture model from DNA samples. DNA samples are    #
# collected through the Arctic fox monitoring programme        #
# (Ulvund et al. 2023; NINA Report 2344).                      #
#                                                              #
# ------------------------------------------------------------ #

library(tidyverse)
library(readxl)
if(!("NIcalc" %in% installed.packages())){
  devtools::install_github("NINAnor/NIcalc", build_vignettes = T)
}
library(NIcalc)
library(dplyr)
library(magrittr)

# SETUP #
#-------#

## Source all functions in "R" folder
sourceDir <- function(path, trace = TRUE, ...) {
  cat("Loading functions:\n")
  cat("------------------\n")
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,"")
    source(file.path(path, nm), encoding = "UTF-8")
    if(trace) cat("\n")
  }
}
sourceDir('R')


## Set switches 

# limit data by time interval (min year to max year)
# NOTE: estimates are provided as 3-year moving averages. The first and last
# year will therefore be censored (last year will be set equal to last year-1).
min.year <- 2007
max.year <- 2024


# IMPORT DATA #
#-------------#
# NOTE: Data is downloaded from ROVBASE (Norwegian Environment Agency).
#       There is no API available for ROVBASE, thus data is downloaded
#       via the web interface by one of the very few people with access.
#       The raw data is then imported here.
# Den controls = a complete record of all controlled dens
# DNA = a complete record of all Arctic fox DNA samples

# DEN CONTROLS
inputControls <- as.data.frame(read_excel("data/eksport_NSF_kontroller_rovbase_2024.xlsx",
                                          sheet = "Rapport",
                                          trim_ws = T,
                                          col_types = c(rep("text",51),"date",rep("text",34))))

# DNA DATA
inputDNA <- as.data.frame(read_excel("data/eksport_NSF_DNA_rovbase_2024.xlsx",
                                     sheet = "Rapport",
                                     trim_ws = T,
                                     col_types = c(rep("text",7),"date",rep("text",51))))

# Download existing data from the Nature Index database
# NOTE: Requires login with username and password obtained from NIdb administrators
indicators <- c("Fjellrev")
currentData <- downloadData_NIdb(indicators = indicators, save = F, save_path = getwd())
fjellrevTable <- currentData$Fjellrev$indicatorValues

# create a NIdb reference table for area ID's and names (for data filtering)
ni_refs <- subset(fjellrevTable, fjellrevTable$yearName == "Referanseverdi", select = c("areaId","areaName"))
colnames(ni_refs) <- c("Area_ID","Area.Name")


# WRANGLE DEN CONTROL DATA #
#--------------------------#

# Format den control data - include only data from min.year to max.year
# Conform municipalities to NIdb names and area ID's
# NOTE: the wrangled data will only contain municipalities with at least 1 controlled Arctic fox den
control_data <- wrangleControls(inputControls, min.year, max.year, ni_refs)

# Unique dens (based on last control year)
# NOTE: this is used as reference for censoring dens while building
# encounter histories (CMRdata.R)
unique_dens <- uniqueDens(control_data)

# Assign patch number to control data
# NOTE: this is only used for (future) abundance estimation on habitat level
control_data <- getPatch(control_data)


# WRANGLE DNA DATA #
#------------------#

# Format DNA data - include only data from min.year to max.year
# Conform municipalities to NIdb names and area ID's
dna <- wrangleDNA(inputDNA, min.year, max.year, ni_refs)

# Assign patch number to DNA data
dna <- getPatch(dna)


# LITTER SIZE #
#-------------#

# Calculate the number of litters per municipality, patch and mountain area per year
# NOTE: returned data only includes spatial units with confirmed breedings
litterList <- litterSize(df = control_data, export = FALSE)


# PREPARE CAPTURE-RECAPTURE DATA #
#--------------------------------#

# Prepare data for estimating abundance for municipalities
cmrData <- CMRdata(dna, unique_dens)


# MODEL RUN #
#-----------#

# Run CMR models to estimate population abundance per municipality
t.start <- Sys.time()

cmrResults <- capture.models.dna(cmrData)

Sys.time() - t.start

cleanup(ask = F)


# ESTIMATE ARCTIC FOX ABUNDANCE #
#-------------------------------#

# Combine data and model output into a tidy list
# NOTE: Strata (spatial unit) must be defined here. Defaults to 'Kommune' (municipality; NI default).
#       Other strata are 'Patch' (habitat fragment) or 'Area' (sub-population).
tidy.output <- wrangleOutput(control_data, cmrResults, litterList, strata = "Kommune", export = FALSE)

# Estimate Arctic fox abundance using regional estimates, correction factors and litter size
# 'export = write output CSV
# 'show_existing_values = adds columns with current Naturindeks values for comparison in quality assessment
#  (currently only imports data from CSV - add functionality to read data directly from NI database)
est.abundance <- estimateAbundance(tidy.output, fjellrevTable, export = FALSE, show_existing_values = FALSE)


# PREPARE DATA FOR NI DATABASE AND UPDATE NI TABLE #
#--------------------------------------------------#

## Prepare new data
# last check for municipality nomenclature consistency
all(unique(est.abundance$strata) %in% unique(fjellrevTable$areaName)) #should be TRUE

# make a copy of the current indicator value table for editing
editable <- currentData$Fjellrev

# save the copy of the current indicator value table
#saveRDS(editable, "data/oldIndicatorValues_2025.rds")

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

# UPLOAD DATA TO NI DATABASE #
#----------------------------#
updatedIndicatorData <- currentData
updatedIndicatorData$Fjellrev <- editable

uploadData_NIdb(indicatorData = updatedIndicatorData)
