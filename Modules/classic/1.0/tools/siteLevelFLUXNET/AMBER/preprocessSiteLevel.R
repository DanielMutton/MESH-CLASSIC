# Preprocess selected CLASSIC outputs to make the data comparabale to reference data

#-------------------------------------------------------------------------------
# General settings
#-------------------------------------------------------------------------------
sites <- c('AU-Tum', 'CA-TPD', 'CZ-BK1', 'DK-Sor', 'FR-Fon', 'GH-Ank', 'IT-Tor',
           'PA-SPs', 'RU-Ha1', 'US-WCr', 'ZM-Mon', 'BR-Sa1', 'CG-Tch', 'DE-Kli', 'ES-LgS',
           'FR-Pue', 'IT-Lav', 'MY-PSO', 'RU-Che', 'RU-SkP', 'US-Wkg', 'CA-Qfo', 'CN-Dan',
           'DE-Tha', 'FI-Hyy', 'GF-Guy', 'IT-SRo', 'NL-Loo', 'RU-Fyo', 'SD-Dem', 'ZA-Kru')
#-------------------------------------------------------------------------------
for(i in 1:length(sites))
{
  sitename <- sites[i]
  csv.mod.gpp <- file.path(mod.csv.path, sitename, "csv", "gpp_monthly.csv")
  csv.mod.ra <- file.path(mod.csv.path, sitename, "csv", "ra_monthly.csv")
  csv.mod.rh <- file.path(mod.csv.path, sitename, "csv", "rh_monthly.csv")
  csv.mod.rss <- file.path(mod.csv.path, sitename, "csv", "rss_monthly.csv")
  csv.mod.rls <- file.path(mod.csv.path, sitename, "csv", "rls_monthly.csv")
  #-----------------------------------------------------------------------------
  # NEE
  #-----------------------------------------------------------------------------
  gpp <- read.csv(csv.mod.gpp)
  ra <- read.csv(csv.mod.ra)
  rh <- read.csv(csv.mod.rh)
  nee <- gpp[4] - ra[4] - rh[4]
  colnames(nee) <- "nee"
  latLonTime <- gpp[1:3]
  nee <- data.frame(latLonTime, nee)
  write.csv(nee, file = "nee_monthly.csv", row.names = FALSE)

  #-----------------------------------------------------------------------------
  # RECO
  #-----------------------------------------------------------------------------
  # ecosystem respirtation (reco) is not an outputfile in CLASSIC
  # this script preprocesses reco (reco = ra + rh)
  ra <- read.csv(csv.mod.ra)
  rh <- read.csv(csv.mod.rh)
  reco <- ra[4] + rh[4]
  colnames(reco) <- "reco"
  latLonTime <- ra[1:3]
  reco <- data.frame(latLonTime, reco)
  write.csv(reco, file = "reco_monthly.csv", row.names = FALSE)

  #-----------------------------------------------------------------------------
  # RNS
  #-----------------------------------------------------------------------------
  # RNS is not an outputfile in CLASSIC
  # this script preprocesses RNS (rns = rss + rls)
  rss <- read.csv(csv.mod.rss)
  rls <- read.csv(csv.mod.rls)
  rns <- rss[4] + rls[4]
  colnames(rns) <- "rns"
  latLonTime <- rss[1:3]
  rns <- data.frame(latLonTime, rns)
  write.csv(rns, file = "rns_monthly.csv", row.names = FALSE)

  #-----------------------------------------------------------------------------
  # move all files to the preprocess folder
  #-----------------------------------------------------------------------------
  preprocessDirSite <- file.path(preprocessDir, sitename, "csv")
  system(paste("mkdir -p", preprocessDirSite, sep = " "))
  system(paste("mv *.csv", preprocessDirSite, sep = " "))
}
