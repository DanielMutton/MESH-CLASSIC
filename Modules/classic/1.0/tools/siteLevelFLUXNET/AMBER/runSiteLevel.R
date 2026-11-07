#-------------------------------------------------------------------------------
# What does this script do?
#-------------------------------------------------------------------------------
# The script is divided in 6 main steps:
# (1) Delete any old files you may have
# (2) Load all required R packages
# (3) Preprocess selected variables
# (4) Run AMBER
# (5) Convert PDF output to PNG for website (optional)
# (6) Compile LaTeX file
#-------------------------------------------------------------------------------
# (1) Delete any old files you may have.
#-------------------------------------------------------------------------------
# Create (if doesn't exist) and Set the output directory as a working directory:
dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
# Create a new preprocess folder if it doesn't exist already
preprocessDir <- paste(outputDir, "preprocess", sep = "/")
dir.create(preprocessDir, recursive = TRUE, showWarnings = FALSE)
# Delete any old amber-related files under <outputDir>
# EXCEPT amber_methods.tex and our <preprocessDir>
to_delete <- c(
  # NetCDF
  list.files(path = ".", pattern = paste("(AGLBIO-|ALBS-|BURNT-|CSOIL-|GPP-|HFLS-|",
                                     "HFSS-|LAI-|MRRO-|NEE-|RECO-|RLS-|RNS-|",
                                     "RSS-|SNW-|score_|scoreinputs).*\\.nc$", sep = ""), full.names = TRUE),
  # PDF
  list.files(path = ".", pattern = paste("(AGLBIO-|ALBS-|BURNT-|CSOIL-|GPP-|HFLS-|",
                                     "HFSS-|LAI-|MRRO-|NEE-|RECO-|RLS-|RNS-|",
                                     "RSS-|SNW-|score_|scoreinputs|",
                                     "amberResults|-plotNc).*\\.pdf$", sep = ""), full.names = TRUE),
  # PNG
  list.files(path = ".", pattern = paste("(AGLBIO-|ALBS-|BURNT-|CSOIL-|GPP-|HFLS-|",
                                     "HFSS-|LAI-|MRRO-|NEE-|RECO-|RLS-|RNS-|",
                                     "RSS-|SNW-|score_|scoreinputs|-plotNc).*\\.png$", sep = ""), full.names = TRUE),
  # LATEX
  list.files(path = ".", pattern = "(score_|scoreinputs).*\\.tex$", full.names = TRUE),
  # Preprocessed inputs (if they exist, but leave the folder alone)
  list.files(path = "./preprocess", pattern = "(aglbio|nee|reco|rns)_(annually|monthly)\\.nc$", full.name = TRUE),
  # Miscellaneous tables and other stuff
  list.files(path = ".", pattern = "(RNS|RECO|NEE|HFSS|HFLS|GPP)_FLUXNET$", full.names = TRUE),
  list.files(path = ".", pattern = "(MRRO)_runoff", full.names = TRUE),
  list.files(path = ".", pattern = "(scoreinputs_|scorevalues_|allscorevalues-).*$", full.names = TRUE),
  list.files(path = ".", pattern = "amber_methods\\.*(aux|bbl|blg|log|out)", full.names = TRUE),
  list.files(path = ".", pattern = "reference.*\\.(csv|tex|bib)", full.names = TRUE)
)
for (f in to_delete) { unlink(f) }

#-------------------------------------------------------------------------------
# (2) Load all required packages
#-------------------------------------------------------------------------------
library(amber)
library(classInt)
library(doParallel)
library(foreach)
library(latex2exp)
library(ncdf4)
library(raster)
library(scico)
library(sp)
library(spatial.tools)
library(stats)
library(utils)
library(viridis)
library(xtable)

print("All R packages required to run AMBER have been loaded.")
#-------------------------------------------------------------------------------
# (3) Preprocess selected variables
#-------------------------------------------------------------------------------

# Preprocess selected CLASSIC outputs to make the data comparabale to reference data

source("preprocessSiteLevel.R")
print("preprocessSiteLevel.R: DONE")

#-------------------------------------------------------------------------------
# (4) Run AMBER
#-------------------------------------------------------------------------------

# record time when AMBER was started:
amber.start <- Sys.time()
# run the evaluation script:
source("evaluateSiteLevel.R")
print("evaluateSiteLevel.R: DONE")

#-------------------------------------------------------------------------------
# (5) Convert PDF output to PNG for website
#-------------------------------------------------------------------------------
# convert all PDF files to PNG files, which can be used in a web interface
if (pdf2png == TRUE) {system("for x in $(ls ./*.pdf); do pdftoppm $x $x -r 300 -png; done")}
# need to rename the new pngs to remove the ugly file extension suffix
to_rename <- list.files(path = ".", pattern = "\\.pdf-1.png$", full.names = TRUE)
for (src in to_rename) {
  dst <- file.path(".", gsub('\\.pdf-1.png$', '.png', basename(src)))
  file.rename(src, dst)
}

#-------------------------------------------------------------------------------
# (6) Compile LaTeX file
#-------------------------------------------------------------------------------

# Copy files related to the LaTeX document from amber.gitrepo.path to the working dir (outputPath)
to_copy <- grep(
  list.files(path = file.path(amber.gitrepo.path, "scripts", "assets"), full.names = TRUE),
  pattern = '*.R$', inv = TRUE, value = TRUE
)
for (src in to_copy) {
  dst <- file.path(".", basename(src))
  file.copy(src, dst, overwrite = TRUE)
}

# (a) Create reference_data_table.tex.
# This will update any potential changes in reference_data_table.csv

data <- read.csv("reference_data_tableGlobal.csv", header = TRUE, sep = ";")
data <- xtable(data)
caption(data) <- "Reference data"
print.xtable(data, type = "latex", label = "tab:ref_data", include.rownames = FALSE, file = "reference_data_table.tex",
             caption.placement = "top", sanitize.text.function = function(x) {x})

# (b) Compile LaTeX document amber_methods.tex
system("pdflatex amber_methods; bibtex amber_methods; pdflatex amber_methods; pdflatex amber_methods")

# (c) Add all plots to amber_methods.pdf.
# Plots are added in the order they were produced, except for the score plots, which are put first.
if (compareTwoRuns == FALSE) {system("pdfunite amber_methods.pdf score_without_massweighting.pdf $(ls -rt *-*.pdf) tmp.pdf")}
if (compareTwoRuns == TRUE) {system("pdfunite amber_methods.pdf scores_compare.pdf $(ls -rt *-*.pdf) tmp.pdf")}

# (d) Make a PDF outline using plot names (first page: n=13)
system('n=13 ; for i in $(ls -tr *-*.pdf) ; do echo "[/Page $n /Title ($i) /OUT pdfmark" ; n=$((n+1)) ; done > index.info')

# (e) Add outline to PDF document
system("gs -sDEVICE=pdfwrite -q -dBATCH -dNOPAUSE -sOutputFile=amberResults.pdf -dPDFSETTINGS=/prepress index.info -f tmp.pdf")

# (f) remove files that are no longer needed
system("rm tmp.pdf; rm index.info; rm amber_methods.pdf")
print("Your AMBER run has completed. Please find your results in the document amberResults.pdf.")

#-------------------------------------------------------------------------------
amber.end <- Sys.time()  # record time when AMBER has finished
duration <- amber.start - amber.end  # time duration
duration <- abs(duration)
#-------------------------------------------------------------------------------
print(duration)


