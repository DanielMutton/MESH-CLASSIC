#-------------------------------------------------------------------------------
# Adjust the variables here to configure the amber run. Currently there are 3
# types of runs based on the CLASSIC run types:
#       1) Global         
#       2) Canada-only      !! Uses a rotated grid !!
#       3) Site-level       !! Expects CLASSIC outputs in CSV, not NetCDF !!
#
# Also, the ability to compare 2 amber runs is configured by changing the 
# compareTwoRuns variable found below
#-------------------------------------------------------------------------------
rm(list = ls()) # For safety, this removes any and all R objects currently loaded
#-------------------------------------------------------------------------------
# Choose your CLASSIC model output data: "global", "canada", or "siteLevel"
modelOutputType <- "siteLevel"

# Set the CLASSIC model id string (what will be displayed on amber plots)
mod.id <- "FLUXNET"

# Set the git repository path (needed when referencing scripts and LaTeX assets)
amber.gitrepo.path <- ""

# Where the amber outputs will be stored (directory). Also acts as a working directory for the amber run
outputDir <- ""

# Where the CLASSIC NetCDF model output is stored (directory)
mod.path <- ""

# Where the CLASSIC CSV model output is stored (directory)  ??? Only necessary for siteLevel runs?
mod.csv.path <- ""

# Where the NetCDF reference data (ilamb) is stored (directory)
ref.path <- ""

# Where the CSV reference data (FluxNet) is stored (directory)
ref.csv.path <- ""

# Set your number of cores used to run amber in parallel:
numCores <- 2

# Define score weights used to calculate the overall score:
score.weights <- c(1, 2, 1, 1, 1) # S_bias, S_rmse, S_phase, S_iav, S_dist

# Define a time period over which model results should be averaged when comparing
# model outputs against reference data that do not have a time dimension (e.g. soil carbon)
timePeriod <- c("1980-01", "2017-12")

# Define a time period for plotting model outputs when not compared to observations
start.date.plot <- "1990-01"; 
end.date.plot <- "2017-12"

# Do you wish to convert PDF plots to PNG?
pdf2png <- TRUE

# Path to conda-managed CDO utility (only used in Canada runs)
amber.cdo <- ""

#-------------------------------------------------------------------------------
# Comparing two amber runs
#-------------------------------------------------------------------------------
# This toggles the comparison of two different amber runs (produces plot "scores_compare.pdf")
compareTwoRuns <- FALSE

# The next 4 lines will be ignored if compareTwoRuns is FALSE
if (compareTwoRuns == TRUE) {
  # baseline amber run directory you wish to compare against
  mod01.path <- outputDir
  # your current amber run directory
  mod02.path <- outputDir
  # how to label the baseline in amber plots
  mod01.id <- "baseline"
  # the current amber run labels will be the same as mod.id
  mod02.id <- mod.id
}
#-------------------------------------------------------------------------------

# Create (if necessary) and set the working directory
dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
setwd(outputDir)
print(paste("All outputs created by AMBER will be stored in", outputDir, sep = " "))

print("Configuration completed")

#-------------------------------------------------------------------------------
# No need to edit the following... unless more run types are developed
#
# Based on the chosen run type, copy the required scripts from the 
# gitlab repository folder to outputDir, then run the appropriate runScript
#-------------------------------------------------------------------------------
if (modelOutputType == "global") {
    runScript <- file.path(amber.gitrepo.path, "scripts", "runGlobal.R")
    preprocessScript <- file.path(amber.gitrepo.path, "scripts", "preprocessGlobal.R")
    evaluateScript <- file.path(amber.gitrepo.path, "scripts", "evaluateGlobal.R")
} else if (modelOutputType == "canada") {
    runScript <- file.path(amber.gitrepo.path, "scripts", "runCanada.R")
    preprocessScript <- file.path(amber.gitrepo.path, "scripts", "preprocessCanada.R")
    evaluateScript <- file.path(amber.gitrepo.path, "scripts", "evaluateCanada.R")
} else if (modelOutputType == "siteLevel") {
    runScript <- file.path(amber.gitrepo.path, "runSiteLevel.R")
    preprocessScript <- file.path(amber.gitrepo.path, "preprocessSiteLevel.R")
    evaluateScript <- file.path(amber.gitrepo.path, "evaluateSiteLevel.R")
} else {
    stop("modelOutputType must be one of 'global', 'canada', or 'siteLevel'")
    quit() # this line is probably overkill
}
# since we are already in outputDir, we can use relative paths safely
for (src in c(runScript, preprocessScript, evaluateScript))
{
    dst <- file.path(".", basename(src))
    file.copy(src, dst, overwrite=TRUE)
}
source(file.path(".", basename(runScript)))
