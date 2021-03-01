#CDS
#https://cran.r-project.org/web/packages/ecmwfr/vignettes/cds_vignette.html
#library(cds)
#ecmwfr
#mcera5

#UID78176
#API Key062c3e77-bcc8-4c56-8e72-4872e7a92be6

#-------
#https://www.erikkusch.com/post/krigr-mats/krigrworkshop/

#devtools::install_github("https://github.com/ErikKusch/KrigR")
library(KrigR)
library(raster)

Extent <- extent(c(14.8, 15.1, 50.1, 50.7)) # roughly the extent of Saxony

#Dir.StateExt <- file.path(Dir.Data, "State_Extent")
#dir.create(Dir.StateExt)

State_Raw <- download_ERA(
  Variable = "skin_temperature",
  DataSet = "era5-land",
  DateStart = "1995-07-20",
  DateStop = "1995-07-22",
  TResolution = "hour",
  TStep = 1,
  Extent = Extent,
  Dir = "/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/ERA5/",
  API_User = "78176",
  API_Key = "062c3e77-bcc8-4c56-8e72-4872e7a92be6"
)

plot(State_Raw)

Plot_Raw(State_Raw, Dates = c("20-07-1995", "21-07-1995"))

#https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Parameterlistings
#2m_temperature
#sea_surface_temperature, K
#skin_temperature
lake_total_layer_temperature

