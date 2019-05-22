#Loading marine data
#From http://www.marineheatwaves.org/tracker.html

# The three packages we will need
library(dplyr)
library(rerddap)
library(ncdf4)
library(tidyverse)
library(heatwaveR)

# The information for the NOAA OISST data
info(datasetid = "ncdc_oisst_v2_avhrr_by_time_zlev_lat_lon", url = "https://www.ncei.noaa.gov/erddap/")

# This function expects the user to provide it with two values 
# that match the time format of the target OISST dataset
OISST_sub <- function(times){
  oisst_res <- griddap(x = "ncdc_oisst_v2_avhrr_by_time_zlev_lat_lon", 
                       url = "https://www.ncei.noaa.gov/erddap/", 
                       time = times, 
                       depth = c(0, 0),
                       latitude = c(-40, -35),
                       longitude = c(15, 21),
                       fields = "sst")
}

#Download data
OISST1 <- OISST_sub(c("1981-09-01T00:00:00Z", "1990-12-31T00:00:00Z"))
OISST2 <- OISST_sub(c("1991-01-01T00:00:00Z", "1999-12-31T00:00:00Z"))
OISST3 <- OISST_sub(c("2000-01-01T00:00:00Z", "2008-12-31T00:00:00Z"))
OISST4 <- OISST_sub(c("2009-01-01T00:00:00Z", "2013-12-03T00:00:00Z"))
OISST5 <- OISST_sub(c("2014-01-01T00:00:00Z", "2018-12-03T00:00:00Z"))

#Preparing data
OISST_prep <- function(nc_file){
  
  # Open the NetCDF connection
  nc <- nc_open(nc_file$summary$filename)
  
  # Extract the SST values and add the lon/lat/time dimension names
  res <- ncvar_get(nc, varid = "sst")
  dimnames(res) <- list(lon = nc$dim$longitude$vals,
                        lat = nc$dim$latitude$vals,
                        t = nc$dim$time$vals)
  
  # Convert the data into a 'long' dataframe for use in the 'tidyverse' ecosystem
  res <- as.data.frame(reshape2::melt(res, value.name = "temp"), row.names = NULL) %>% 
    mutate(t = as.Date(as.POSIXct(t, origin = "1970-01-01 00:00:00")),
           temp = round(temp, 2))
  
  # Close the NetCDF connection and finish
  nc_close(nc)
  return(res)
}

# Prep the data
OISST1_prep <- OISST_prep(OISST1)
OISST2_prep <- OISST_prep(OISST2)
OISST3_prep <- OISST_prep(OISST3)
OISST4_prep <- OISST_prep(OISST4)
OISST5_prep <- OISST_prep(OISST5)

# Bind them together
OISST <- rbind(OISST1_prep, OISST2_prep, OISST3_prep, OISST4_prep, OISST5_prep)

# Save the data as an .Rda file
#saveRDS(OISST, file = "~/Desktop/OISST_vignette.Rda")
#OISST <- read_rds("~/Desktop/OISST_vignette.Rda")

#-------------------------
#Event detection

event_only <- function(df){
  # First calculate the climatologies
  clim <- ts2clm(data = df, climatologyPeriod = c("1982-01-01", "2011-01-01"))
  # Then the events
  event <- detect_event(data = clim)
  # Lastly we return only the event dataframe of results
  return(event$event)
}

#purrr method

system.time(
  # First we start by chosing the 'OISST' dataframe
  MHW_purrr <- OISST %>% 
    # Then we group the data by the 'lon' and 'lat' columns
    group_by(lon, lat) %>% 
    # After that we nest each lon/lat pixel into its own little dataframe
    # The column containing all of these little dataframes is named 'data'
    nest() %>% 
    # Next we 'map' the event_only() wrapper function we made to the little dataframes
    mutate(event = map(data, event_only)) %>% 
    # Lastly we get rid of the column that is not part of the final result 
    # before unnesting everything
    select(-data) %>% 
    unnest()
)


MHW_result<- MHW_purrr
#-------------------------
#Trend detection

# summarise the number of unique longitude, latitude and year combination:
event_freq <- MHW_result %>% 
  mutate(year = year(date_start)) %>% 
  group_by(lon, lat, year) %>% 
  summarise(n = n())
head(event_freq)

# create complete grid for merging with:
sst_grid <- OISST %>% 
  select(lon, lat, t) %>% 
  mutate(t = lubridate::year(t)) %>% 
  dplyr::rename(year = t) %>% 
  distinct()

# and merge:
OISST_n <- left_join(sst_grid, event_freq, by = c("lon", "lat", "year")) %>% 
  mutate(n = ifelse(is.na(n), 0, n))

lin_fun <- function(ev) {
  mod1 <- glm(n ~ year, family = poisson(link = "log"), data = ev)
  # extract slope coefficient and its p-value
  tr <- data.frame(slope = summary(mod1)$coefficients[2,1],
                   p = summary(mod1)$coefficients[2,4])
  return(tr)
}

OISST_nTrend <- OISST_n %>% 
  group_by(lon, lat) %>% 
  nest() %>% 
  mutate(res = map(data, lin_fun)) %>% 
  select(-data) %>%
  unnest() %>% 
  mutate(pval = cut(p, breaks = c(0, 0.001, 0.01, 0.05, 1)))
head(OISST_nTrend)

OISST_nTrend <- plyr::ddply(OISST_n, c("lon", "lat"), lin_fun)
OISST_nTrend$pval <- cut(OISST_nTrend$p, breaks = c(0, 0.001, 0.01, 0.05, 1))
head(OISST_nTrend)

#-------------------------
#Visualizing the results

# The base map
map_base <- ggplot2::fortify(maps::map(fill = TRUE, plot = FALSE)) %>% 
  dplyr::rename(lon = long)

map_slope <- ggplot(OISST_nTrend, aes(x = lon, y = lat)) +
  geom_rect(size = 0.2, fill = NA,
            aes(xmin = lon - 0.1, xmax = lon + 0.1, ymin = lat - 0.1, ymax = lat + 0.1,
                colour = pval)) +
  geom_raster(aes(fill = slope), interpolate = FALSE, alpha = 0.9) +
  scale_fill_gradient2(name = "count/year (slope)", high = "red", mid = "white",
                       low = "darkblue", midpoint = 0,
                       guide = guide_colourbar(direction = "horizontal",
                                               title.position = "top")) +
  scale_colour_manual(breaks = c("(0,0.001]", "(0.001,0.01]", "(0.01,0.05]", "(0.05,1]"),
                      values = c("firebrick1", "firebrick2", "firebrick3", "white"),
                      name = "p-value", guide = FALSE) +
  geom_polygon(data = map_base, aes(group = group), 
               colour = NA, fill = "grey80") +
  coord_fixed(ratio = 1, xlim = c(13.0, 23.0), ylim = c(-33, -42), expand = TRUE) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(legend.position = "bottom")

map_p <- ggplot(OISST_nTrend, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = pval), interpolate = FALSE) +
  scale_fill_manual(breaks = c("(0,0.001]", "(0.001,0.01]", "(0.01,0.05]",
                               "(0.05,0.1]", "(0.1,0.5]", "(0.5,1]"),
                    values = c("black", "grey20", "grey40",
                               "grey80", "grey90", "white"),
                    name = "p-value",
                    guide = guide_legend(direction = "horizontal",
                                         title.position = "top")) +
  geom_polygon(data = map_base, aes(group = group), 
               colour = NA, fill = "grey80") +
  coord_fixed(ratio = 1, xlim = c(13.0, 23.0), ylim = c(-33, -42), expand = TRUE) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(legend.position = "bottom")

library(ggpubr)
map_both <- ggpubr::ggarrange(map_slope, map_p, align = "hv")
map_both





