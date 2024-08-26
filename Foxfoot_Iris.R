#iris foxfoot code sample
#ifoxfoot@gmail.com
#august 25, 2024

#PURPOSE
################################################################################
#In this code script I will download PRISM climate data.
#Then I will extract climate variables for each USACE reservoir.
#Finally I will create some simple maps to check the data

#R PACKAGES
################################################################################
#first load R packages
library(here) #package for managing working directories
library(janitor) #package for cleaning column names
library(sf) #package for managing spatial data
library(dplyr) #package for wrangling data
library(prism) #package for downloading prism climate data
library(terra) #package for handling rasters
library(ggplot2) #for maps, plots

#DATA CLEANING
################################################################################
#load USACE reservoir data
#NOTE: CHANGE TO YOUR FILE PATH
#NOTE: Data downloaded from https://geospatial-usace.opendata.arcgis.com/datasets/03e322d7e89b48a9b48e9c3f4bcaf29e/explore
#and modified. modified data available on my github :)
res_raw <- read_sf("raw_data/USACE_Reservoirs/USACE_Reservoirs_Final.shp")

#clean names and geometries
res_clean <- res_raw %>% 
  janitor::clean_names() %>% #convert column names to lower case snake
  sf::st_make_valid()

#get centers of lakes
res_center <- st_centroid(res_clean) %>% 
  select(objectid)

#PRISM DOWNLOAD
################################################################################
#set download directory
prism_set_dl_dir(here("raw_data/prism_data"))

#create function to download, process, and extract PRISM data
process_prism <- function(cli_var, lake_center){
  
  #if prism data doesn't exist for all 12 months, download it
  if(length(prism_archive_subset(cli_var, "monthly normals", resolution = "800m")) != 12) {
    get_prism_normals(type = cli_var, resolution = "800m",
                      mon = 1:12, annual = FALSE, keepZip = FALSE)
  }
  
  #get subset of prism data by type
  subset <- prism_archive_subset(cli_var, "monthly normals", resolution = "800m")
  
  #create raster stack
  stack <- pd_stack(subset)
  
  #make sure lake center is same projection as prism data
  lake_center_reproj <- st_transform(lake_center, st_crs(stack))
  
  #extract values for each lake center
  extracted <- raster::extract(stack, lake_center_reproj, df= T)
  
  #return extracted df
  return(extracted)
  
}

#create list of variables we want from prism normals dataset
climate_vars <- c(
  "ppt", #precipitation
  "tmean", #mean temperature
  "tmin", #min temperature
  "tmax", #max temperature
  "tdmean", #mean dew point temp
  "vpdmin", #minimum vapor deficit
  "vpdmax" #max vapor deficit
)

#store new copy of res_clean for joining
res_with_prism <- res_clean

#loop through and join processed prism to lake data
#make sure you only run this once!!
for(i in climate_vars) {
  prism_dat <- process_prism(i, res_center) %>% select(-ID)
  res_with_prism <- cbind(as.data.frame(res_with_prism), prism_dat)
}

#SEASONAL SUMMARY
################################################################################
#condensing data by season
#winter = dec, jan, feb
#spring = mar, apr, may,
#summer = jun, jul, aug
#fall = sep, oct, nov

#function to condense data into seasonal summary
seasonal_summary <- function(x, season, data_type, fun) {
  season_months <- if(season == "winter") {c("12", "01", "02")} 
  else if(season == "spring") {c("03", "04", "05")} 
  else if(season == "summer") {c("06", "07", "08")}
  else if(season == "fall") {c("09", "10", "11")}
  
  ls <- c()
  
  if(data_type == "ppt") {
    for (i in 1:3) {
      ls[i] <- paste0("PRISM_", data_type, "_30yr_normal_800mM4_", season_months[i], "_bil")
    }
  } else{
    for (i in 1:3) {
      ls[i] <- paste0("PRISM_", data_type, "_30yr_normal_800mM5_", season_months[i], "_bil")
    }
  }
  
  
  new_col <- 
    if (fun == "add") {
      x[,ls[[1]]] + x[,ls[[2]]] + x[,ls[[3]]]
    }
  else if (fun == "mean") {
    (x[,ls[[1]]] + x[,ls[[2]]] + x[,ls[[3]]])/3
  }
  else if (fun == "max") {
    pmax(x[,ls[[1]]], x[,ls[[2]]],x[,ls[[3]]])
  }
  else if (fun == "min")  {
    pmin(x[,ls[[1]]], x[,ls[[2]]],x[,ls[[3]]])
  }
  
  as.data.frame(new_col)
  
  return(new_col)
}


#get total precip per season
seasonal_climate_vals <- res_with_prism %>% 
  #seasonal total precip
  mutate(winter_total_precip = seasonal_summary(res_with_prism, season = "winter", data_type = "ppt", fun = "add")) %>% 
  mutate(spring_total_precip = seasonal_summary(res_with_prism, season = "spring", data_type = "ppt", fun = "add")) %>%
  mutate(summer_total_precip = seasonal_summary(res_with_prism, season = "summer", data_type = "ppt", fun = "add")) %>% 
  mutate(fall_total_precip = seasonal_summary(res_with_prism, season = "fall", data_type = "ppt", fun = "add")) %>% 
  #seasonal mean temp
  mutate(winter_mean_temp = seasonal_summary(res_with_prism, season = "winter", data_type = "tmean", fun = "mean")) %>% 
  mutate(spring_mean_temp = seasonal_summary(res_with_prism, season = "spring", data_type = "tmean", fun = "mean")) %>% 
  mutate(summer_mean_temp = seasonal_summary(res_with_prism, season = "summer", data_type = "tmean", fun = "mean")) %>% 
  mutate(fall_mean_temp = seasonal_summary(res_with_prism, season = "fall", data_type = "tmean", fun = "mean")) %>% 
  #seasonal min temp
  mutate(winter_min_temp = seasonal_summary(res_with_prism, season = "winter", data_type = "tmin", fun = "min")) %>% 
  mutate(spring_min_temp = seasonal_summary(res_with_prism, season = "spring", data_type = "tmin", fun = "min")) %>% 
  mutate(summer_min_temp = seasonal_summary(res_with_prism, season = "summer", data_type = "tmin", fun = "min")) %>% 
  mutate(fall_min_temp = seasonal_summary(res_with_prism, season = "fall", data_type = "tmin", fun = "min")) %>% 
  #seasonal max temp
  mutate(winter_max_temp = seasonal_summary(res_with_prism, season = "winter", data_type = "tmax", fun = "max")) %>% 
  mutate(spring_max_temp = seasonal_summary(res_with_prism, season = "spring", data_type = "tmax", fun = "max")) %>% 
  mutate(summer_max_temp = seasonal_summary(res_with_prism, season = "summer", data_type = "tmax", fun = "max")) %>% 
  mutate(fall_max_temp = seasonal_summary(res_with_prism, season = "fall", data_type = "tmax", fun = "max")) %>% 
  #seasonal mean dew point
  mutate(winter_mean_dewp = seasonal_summary(res_with_prism, season = "winter", data_type = "tdmean", fun = "mean")) %>% 
  mutate(spring_mean_dewp = seasonal_summary(res_with_prism, season = "spring", data_type = "tdmean", fun = "mean")) %>% 
  mutate(summer_mean_dewp = seasonal_summary(res_with_prism, season = "summer", data_type = "tdmean", fun = "mean")) %>% 
  mutate(fall_mean_dewp = seasonal_summary(res_with_prism, season = "fall", data_type = "tdmean", fun = "mean")) %>% 
  #seasonal min vapor deficit
  mutate(winter_min_vpd= seasonal_summary(res_with_prism, season = "winter", data_type = "vpdmin", fun = "min")) %>% 
  mutate(spring_min_vpd = seasonal_summary(res_with_prism, season = "spring", data_type = "vpdmin", fun = "min")) %>% 
  mutate(summer_min_vpd = seasonal_summary(res_with_prism, season = "summer", data_type = "vpdmin", fun = "min")) %>% 
  mutate(fall_min_vpd = seasonal_summary(res_with_prism, season = "fall", data_type = "vpdmin", fun = "min")) %>% 
  #seasonal max vapor deficit
  mutate(winter_max_vpd= seasonal_summary(res_with_prism, season = "winter", data_type = "vpdmax", fun = "max")) %>% 
  mutate(spring_max_vpd = seasonal_summary(res_with_prism, season = "spring", data_type = "vpdmax", fun = "max")) %>% 
  mutate(summer_max_vpd = seasonal_summary(res_with_prism, season = "summer", data_type = "vpdmax", fun = "max")) %>% 
  mutate(fall_max_vpd = seasonal_summary(res_with_prism, season = "fall", data_type = "vpdmax", fun = "max"))

seasonal_climate <- seasonal_climate_vals %>% 
  select(!PRISM_ppt_30yr_normal_800mM4_01_bil:PRISM_vpdmax_30yr_normal_800mM5_12_bil)

#MAKE MAPS
################################################################################
#get state data
state_map <- map_data("state")

#get center points again
seasonal_temp_pts <- st_as_sf(seasonal_climate) %>% 
  st_centroid()

#define temp limits (for fixed color scales between plots)
min_temp <- min(seasonal_temp_pts$winter_mean_temp)
max_temp <- max(seasonal_temp_pts$summer_mean_temp)

#plot summer temps
ggplot() +
  geom_polygon(data = state_map, aes(x=long,y=lat,group=group), inherit.aes=F, 
               colour='black', fill=NA) +
  geom_sf(data = seasonal_temp_pts, aes(color = summer_mean_temp)) +
  scale_color_gradientn(colors = c("blue", "green", "orange", "red"),
                        limits = c(min_temp, max_temp)) + 
  labs(title = "Mean Summer Temp at USACE Reservoirs",
       color = "Temperature (°C)") +
  theme_minimal()

#plot winter temps
ggplot() +
  geom_polygon(data = state_map, aes(x=long,y=lat,group=group), inherit.aes=F, 
               colour='black', fill=NA) +
  geom_sf(data = seasonal_temp_pts, aes(color = winter_mean_temp)) +
  scale_color_gradientn(colors = c("blue", "green", "orange", "red"),
                        limits = c(min_temp, max_temp)) + 
  labs(title = "Mean Winter Temp at USACE Reservoirs",
       color = "Temperature (°C)") +
  theme_minimal()