#-------------------------------------------------------------------------------
# Variable-specific settings
#-------------------------------------------------------------------------------
# Notes: a value is classified as an outlier if it is outside a given range this range is defined
# as the interquartile range multiplied by a factor this factor is called 'outlier.factor'
#-------------------------------------------------------------------------------
sites <- c('AU-Tum', 'CA-TPD', 'CZ-BK1', 'DK-Sor', 'FR-Fon', 'GH-Ank', 'IT-Tor',
           'PA-SPs', 'RU-Ha1', 'US-WCr', 'ZM-Mon', 'BR-Sa1', 'CG-Tch', 'DE-Kli', 'ES-LgS',
           'FR-Pue', 'IT-Lav', 'MY-PSO', 'RU-Che', 'RU-SkP', 'US-Wkg', 'CA-Qfo', 'CN-Dan',
           'DE-Tha', 'FI-Hyy', 'GF-Guy', 'IT-SRo', 'NL-Loo', 'RU-Fyo', 'SD-Dem', 'ZA-Kru')
#-------------------------------------------------------------------------------
# CARBON
#-------------------------------------------------------------------------------
# gpp (gross primary productivity):
long.name <- "Gross primary productivity"
mod.csv <- "gpp_monthly.csv"
ref.csv <- paste(ref.csv.path, "gpp_monthly_fluxnet.csv", sep = "/")
mod.id <- 'CLASSIC-Sitelevel' # define a model experiment ID
ref.id <- "Fluxnet"  # give reference dataset a name
unit.conv.mod <- 86400 * 1000  # optional unit conversion for model data
unit.conv.ref <- 86400 * 1000  # optional unit conversion for reference data
variable.unit <- "gC m$^{-2}$ day$^{-1}$"  # unit after conversion (LaTeX notation)
scores.fluxnet.site(long.name, mod.csv, mod.csv.path, ref.csv, mod.id, ref.id,
                    unit.conv.mod, unit.conv.ref, variable.unit, sites, outputDir = outputDir)
print(paste("Evaluation of", long.name, "for", mod.id, "and", ref.id, "completed.", sep = " "))
#-------------------------------------------------------------------------------
# reco (ecosystem respiration):
long.name <- "Ecosystem respiration"
mod.csv <- "reco_monthly.csv"
ref.csv <- paste(ref.csv.path, "reco_monthly_fluxnet.csv", sep = "/")
mod.id <- 'CLASSIC-Sitelevel' # define a model experiment ID
ref.id <- "Fluxnet"  # give reference dataset a name
unit.conv.mod <- 86400 * 1000  # optional unit conversion for model data
unit.conv.ref <- 86400 * 1000  # optional unit conversion for reference data
variable.unit <- "gC m$^{-2}$ day$^{-1}$"  # unit after conversion (LaTeX notation)
scores.fluxnet.site(long.name, mod.csv, mod.csv.path = preprocessDir, ref.csv, mod.id, ref.id,
                    unit.conv.mod, unit.conv.ref, variable.unit, sites, outputDir = outputDir)
print(paste("Evaluation of", long.name, "for", mod.id, "and", ref.id, "completed.", sep = " "))
#-------------------------------------------------------------------------------
# nee (net ecosystem exchange):
long.name <- "Net ecosystem exchange"
mod.csv <- "nee_monthly.csv"
ref.csv <- paste(ref.csv.path, "nee_monthly_fluxnet.csv", sep = "/")
mod.id <- 'CLASSIC-Sitelevel' # define a model experiment ID
ref.id <- "Fluxnet"  # give reference dataset a name
unit.conv.mod <- 86400 * 1000  # optional unit conversion for model data
unit.conv.ref <- (-1)*86400 * 1000  # optional unit conversion for reference data
variable.unit <- "gC m$^{-2}$ day$^{-1}$"  # unit after conversion (LaTeX notation)
scores.fluxnet.site(long.name, mod.csv, mod.csv.path = preprocessDir, ref.csv, mod.id, ref.id,
                    unit.conv.mod, unit.conv.ref, variable.unit, sites, outputDir = outputDir)
print(paste("Evaluation of", long.name, "for", mod.id, "and", ref.id, "completed.", sep = " "))
#-------------------------------------------------------------------------------
# RADIATION
#-------------------------------------------------------------------------------
# rns (net radidation):
long.name <- "Net surface radidation"
mod.csv <- "rns_monthly.csv"
ref.csv <- paste(ref.csv.path, "rns_monthly_fluxnet.csv", sep = "/")
mod.id <- 'CLASSIC-Sitelevel' # define a model experiment ID
ref.id <- "Fluxnet"  # give reference dataset a name
unit.conv.mod <- 1  # optional unit conversion
unit.conv.ref <- 1  # optional unit conversion
variable.unit <- "W m$^{-2}$"  # unit after conversion (LaTeX notation)
scores.fluxnet.site(long.name, mod.csv, mod.csv.path = preprocessDir, ref.csv, mod.id, ref.id,
                    unit.conv.mod, unit.conv.ref, variable.unit, sites, outputDir = outputDir)
print(paste("Evaluation of", long.name, "for", mod.id, "and", ref.id, "completed.", sep = " "))
#-------------------------------------------------------------------------------
# HEAT FLUXES
#-------------------------------------------------------------------------------
# hfls (latent heat flux):
long.name <- "Latent heat flux"
mod.csv <- "hfls_monthly.csv"
ref.csv <- paste(ref.csv.path, "hfls_monthly_fluxnet.csv", sep = "/")
mod.id <- 'CLASSIC-Sitelevel' # define a model experiment ID
ref.id <- "Fluxnet"  # give reference dataset a name
unit.conv.mod <- 1  # optional unit conversion for model data
unit.conv.ref <- 1  # optional unit conversion for reference data
variable.unit <- "W m$^{-2}$"  # unit after conversion (LaTeX notation)
scores.fluxnet.site(long.name, mod.csv, mod.csv.path, ref.csv, mod.id, ref.id,
                    unit.conv.mod, unit.conv.ref, variable.unit, sites, outputDir = outputDir)
print(paste("Evaluation of", long.name, "for", mod.id, "and", ref.id, "completed.", sep = " "))
#-------------------------------------------------------------------------------
# hfss (sensible heat flux):
long.name <- "Sensible heat flux"
mod.csv <- "hfss_monthly.csv"
ref.csv <- paste(ref.csv.path, "hfss_monthly_fluxnet.csv", sep = "/")
mod.id <- 'CLASSIC-Sitelevel' # define a model experiment ID
ref.id <- "Fluxnet"  # give reference dataset a name
unit.conv.mod <- 1  # optional unit conversion for model data
unit.conv.ref <- 1  # optional unit conversion for reference data
variable.unit <- "W m$^{-2}$"  # unit after conversion (LaTeX notation)
scores.fluxnet.site(long.name, mod.csv, mod.csv.path, ref.csv, mod.id, ref.id,
                    unit.conv.mod, unit.conv.ref, variable.unit, sites, outputDir = outputDir)
print(paste("Evaluation of", long.name, "for", mod.id, "and", ref.id, "completed.", sep = " "))
#-------------------------------------------------------------------------------
# summarize results in a table and a plot
scores.tables(plot.width = 6, plot.height = 3, myMargin = c(4, 7, 3, 4), outputDir = outputDir)
#-------------------------------------------------------------------------------'
