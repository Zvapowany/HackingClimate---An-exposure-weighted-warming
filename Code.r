###############################################################################
###############################################################################
#
#
#     Project: An exposure weighted warming
#     Team: Marcin Gajewski, Salma Al-Zadjali, Yasna Palmeiro Silva
#     
#     Analysis of UKCP18 dataset for 
#       - Mean air temperature at 1.5m (°C)
#       - Maximum air temperature at 1.5m (°C)
#       - Minimum air temperature at 1.5m (°C)
#
#     In order to obtain meaningful results, we used monthly data.
#
#     Based on the data, we calculated temperature anomalies, considering
#     a historical baseline of 20 years (1980-2000). 
#     
#     In order to calculate population-weighted average anomalies, 
#     we used population estimations from 2000 onward. 
#
###############################################################################
###############################################################################

# Libraries we need
library(raster)   
library(sf)
library(sp)
library(ncdf4)
library(ggplot2)
library(rgdal)
library(rlang)
library(dplyr)
library(lubridate)
library(zoo)
library(stringr)
library(tidyverse) #include readr
library(ggthemes)
library(ggpubr)

###############################################################################
#   SHAPEFILES
###############################################################################

# We need: ukcp18-uk-land-country-hires 
uk_borders <- st_read("/nfs/cfs/home4/rejb/rejbypa/UKCP18_R/ukcp-spatial-files/spatial-files/ukcp18-uk-land-country-hires/ukcp18-uk-land-country-hires.shp")
as_Spatial(uk_borders)
plot(uk_borders)

uk_borders$geo_region # "Channel Islands", "England", "Isle of Man", "Northern Ireland", "Scotland", "Wales"

# In order to obtain some data for each region, we create subsets of polygons
Cha_borders <- uk_borders[uk_borders$geo_region == "Channel Islands",]
Eng_borders <- uk_borders[uk_borders$geo_region == "England",]
Isl_borders <- uk_borders[uk_borders$geo_region == "Isle of Man",]
NIr_borders <- uk_borders[uk_borders$geo_region == "Northern Ireland",]
Sco_borders <- uk_borders[uk_borders$geo_region == "Scotland",]
Wal_borders <- uk_borders[uk_borders$geo_region == "Wales",]

###############################################################################
#   POPULATION
###############################################################################

# Dataset: SSP5 population projections on the OSGB grid 
pop.nc <- nc_open("/nfs/cfs/home4/rejb/rejbypa/UKCP18_R/population/population-projections-ssp5-total.nc")
print(pop.nc)
names(pop.nc$var) #"crs" "ssp"

pop <- brick("/nfs/cfs/home4/rejb/rejbypa/UKCP18_R/population/population-projections-ssp5-total.nc")
pop   # This dataset has 11 years: 2000, 2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100.
      # We checked that resolution, extent, and CRS were the same for temperature projections. 
      # This dataset contains absolute number of people living per grid.

#################################################
# Manipulation & interpolation
#################################################

# Masking (M) and cropping according to uk_borders
popM <- mask(pop, uk_borders) %>% crop(uk_borders)

# Checking a time series of population
pop_TS <- c(cellStats(popM, sum, na.rm=T))
pop_TS <- ts(pop_TS)
plot(pop_TS)    # Population estimation grows in a approximately linear fashion.

# We can linearly interpolate years between given years
# formula for linear interpolation: P(T) = P1 + [(P2-P1)*(T-T1)/(T2-T1)]
# Interpolation is done decade by decade to preserve some different growth rates. 

UKpop2000_2100 <- brick(nrows=109, ncols=55, xmn=0, xmx=660000, ymn=-84000, ymx=1224000 , nl=81)

UKpop2000_2100[[1]] <- popM[[1]] #2000
UKpop2000_2100[[2]] <- popM[[1]] + (((popM[[2]] - popM[[1]]) * (2001-2000)) / (2010-2000))
UKpop2000_2100[[3]] <- popM[[1]] + (((popM[[2]] - popM[[1]]) * (2002-2000)) / (2010-2000))
UKpop2000_2100[[4]] <- popM[[1]] + (((popM[[2]] - popM[[1]]) * (2003-2000)) / (2010-2000))
UKpop2000_2100[[5]] <- popM[[1]] + (((popM[[2]] - popM[[1]]) * (2004-2000)) / (2010-2000))
UKpop2000_2100[[6]] <- popM[[1]] + (((popM[[2]] - popM[[1]]) * (2005-2000)) / (2010-2000))
UKpop2000_2100[[7]] <- popM[[1]] + (((popM[[2]] - popM[[1]]) * (2006-2000)) / (2010-2000))
UKpop2000_2100[[8]] <- popM[[1]] + (((popM[[2]] - popM[[1]]) * (2007-2000)) / (2010-2000))
UKpop2000_2100[[9]] <- popM[[1]] + (((popM[[2]] - popM[[1]]) * (2008-2000)) / (2010-2000))
UKpop2000_2100[[10]] <- popM[[1]] + (((popM[[2]] - popM[[1]]) * (2009-2000)) / (2010-2000))
UKpop2000_2100[[11]] <- popM[[2]] #2010
UKpop2000_2100[[12]] <- popM[[2]] + (((popM[[3]] - popM[[2]]) * (2011-2010)) / (2020-2010))
UKpop2000_2100[[13]] <- popM[[2]] + (((popM[[3]] - popM[[2]]) * (2012-2010)) / (2020-2010))
UKpop2000_2100[[14]] <- popM[[2]] + (((popM[[3]] - popM[[2]]) * (2013-2010)) / (2020-2010))
UKpop2000_2100[[15]] <- popM[[2]] + (((popM[[3]] - popM[[2]]) * (2014-2010)) / (2020-2010))
UKpop2000_2100[[16]] <- popM[[2]] + (((popM[[3]] - popM[[2]]) * (2015-2010)) / (2020-2010))
UKpop2000_2100[[17]] <- popM[[2]] + (((popM[[3]] - popM[[2]]) * (2016-2010)) / (2020-2010))
UKpop2000_2100[[18]] <- popM[[2]] + (((popM[[3]] - popM[[2]]) * (2017-2010)) / (2020-2010))
UKpop2000_2100[[19]] <- popM[[2]] + (((popM[[3]] - popM[[2]]) * (2018-2010)) / (2020-2010))
UKpop2000_2100[[20]] <- popM[[2]] + (((popM[[3]] - popM[[2]]) * (2019-2010)) / (2020-2010))
UKpop2000_2100[[21]] <- popM[[3]] #2020
UKpop2000_2100[[22]] <- popM[[3]] + (((popM[[4]] - popM[[3]]) * (2021-2020)) / (2030-2020))
UKpop2000_2100[[23]] <- popM[[3]] + (((popM[[4]] - popM[[3]]) * (2022-2020)) / (2030-2020))
UKpop2000_2100[[24]] <- popM[[3]] + (((popM[[4]] - popM[[3]]) * (2023-2020)) / (2030-2020))
UKpop2000_2100[[25]] <- popM[[3]] + (((popM[[4]] - popM[[3]]) * (2024-2020)) / (2030-2020))
UKpop2000_2100[[26]] <- popM[[3]] + (((popM[[4]] - popM[[3]]) * (2025-2020)) / (2030-2020))
UKpop2000_2100[[27]] <- popM[[3]] + (((popM[[4]] - popM[[3]]) * (2026-2020)) / (2030-2020))
UKpop2000_2100[[28]] <- popM[[3]] + (((popM[[4]] - popM[[3]]) * (2027-2020)) / (2030-2020))
UKpop2000_2100[[29]] <- popM[[3]] + (((popM[[4]] - popM[[3]]) * (2028-2020)) / (2030-2020))
UKpop2000_2100[[30]] <- popM[[3]] + (((popM[[4]] - popM[[3]]) * (2029-2020)) / (2030-2020))
UKpop2000_2100[[31]] <- popM[[4]] #2030
UKpop2000_2100[[32]] <- popM[[4]] + (((popM[[5]] - popM[[4]]) * (2031-2030)) / (2040-2030))
UKpop2000_2100[[33]] <- popM[[4]] + (((popM[[5]] - popM[[4]]) * (2032-2030)) / (2040-2030))
UKpop2000_2100[[34]] <- popM[[4]] + (((popM[[5]] - popM[[4]]) * (2033-2030)) / (2040-2030))
UKpop2000_2100[[35]] <- popM[[4]] + (((popM[[5]] - popM[[4]]) * (2034-2030)) / (2040-2030))
UKpop2000_2100[[36]] <- popM[[4]] + (((popM[[5]] - popM[[4]]) * (2035-2030)) / (2040-2030))
UKpop2000_2100[[37]] <- popM[[4]] + (((popM[[5]] - popM[[4]]) * (2036-2030)) / (2040-2030))
UKpop2000_2100[[38]] <- popM[[4]] + (((popM[[5]] - popM[[4]]) * (2037-2030)) / (2040-2030))
UKpop2000_2100[[39]] <- popM[[4]] + (((popM[[5]] - popM[[4]]) * (2038-2030)) / (2040-2030))
UKpop2000_2100[[40]] <- popM[[4]] + (((popM[[5]] - popM[[4]]) * (2039-2030)) / (2040-2030))
UKpop2000_2100[[41]] <- popM[[5]] #2040
UKpop2000_2100[[42]] <- popM[[5]] + (((popM[[6]] - popM[[5]]) * (2041-2040)) / (2050-2040))
UKpop2000_2100[[43]] <- popM[[5]] + (((popM[[6]] - popM[[5]]) * (2042-2040)) / (2050-2040))
UKpop2000_2100[[44]] <- popM[[5]] + (((popM[[6]] - popM[[5]]) * (2043-2040)) / (2050-2040))
UKpop2000_2100[[45]] <- popM[[5]] + (((popM[[6]] - popM[[5]]) * (2044-2040)) / (2050-2040))
UKpop2000_2100[[46]] <- popM[[5]] + (((popM[[6]] - popM[[5]]) * (2045-2040)) / (2050-2040))
UKpop2000_2100[[47]] <- popM[[5]] + (((popM[[6]] - popM[[5]]) * (2046-2040)) / (2050-2040))
UKpop2000_2100[[48]] <- popM[[5]] + (((popM[[6]] - popM[[5]]) * (2047-2040)) / (2050-2040))
UKpop2000_2100[[49]] <- popM[[5]] + (((popM[[6]] - popM[[5]]) * (2048-2040)) / (2050-2040))
UKpop2000_2100[[50]] <- popM[[5]] + (((popM[[6]] - popM[[5]]) * (2049-2040)) / (2050-2040))
UKpop2000_2100[[51]] <- popM[[6]] #2050
UKpop2000_2100[[52]] <- popM[[6]] + (((popM[[7]] - popM[[6]]) * (2051-2050)) / (2060-2050))
UKpop2000_2100[[53]] <- popM[[6]] + (((popM[[7]] - popM[[6]]) * (2052-2050)) / (2060-2050))
UKpop2000_2100[[54]] <- popM[[6]] + (((popM[[7]] - popM[[6]]) * (2053-2050)) / (2060-2050))
UKpop2000_2100[[55]] <- popM[[6]] + (((popM[[7]] - popM[[6]]) * (2054-2050)) / (2060-2050))
UKpop2000_2100[[56]] <- popM[[6]] + (((popM[[7]] - popM[[6]]) * (2055-2050)) / (2060-2050))
UKpop2000_2100[[57]] <- popM[[6]] + (((popM[[7]] - popM[[6]]) * (2056-2050)) / (2060-2050))
UKpop2000_2100[[58]] <- popM[[6]] + (((popM[[7]] - popM[[6]]) * (2057-2050)) / (2060-2050))
UKpop2000_2100[[59]] <- popM[[6]] + (((popM[[7]] - popM[[6]]) * (2058-2050)) / (2060-2050))
UKpop2000_2100[[60]] <- popM[[6]] + (((popM[[7]] - popM[[6]]) * (2059-2050)) / (2060-2050))
UKpop2000_2100[[61]] <- popM[[7]] #2060
UKpop2000_2100[[62]] <- popM[[7]] + (((popM[[8]] - popM[[7]]) * (2061-2060)) / (2070-2060))
UKpop2000_2100[[63]] <- popM[[7]] + (((popM[[8]] - popM[[7]]) * (2062-2060)) / (2070-2060))
UKpop2000_2100[[64]] <- popM[[7]] + (((popM[[8]] - popM[[7]]) * (2063-2060)) / (2070-2060))
UKpop2000_2100[[65]] <- popM[[7]] + (((popM[[8]] - popM[[7]]) * (2064-2060)) / (2070-2060))
UKpop2000_2100[[66]] <- popM[[7]] + (((popM[[8]] - popM[[7]]) * (2065-2060)) / (2070-2060))
UKpop2000_2100[[67]] <- popM[[7]] + (((popM[[8]] - popM[[7]]) * (2066-2060)) / (2070-2060))
UKpop2000_2100[[68]] <- popM[[7]] + (((popM[[8]] - popM[[7]]) * (2067-2060)) / (2070-2060))
UKpop2000_2100[[69]] <- popM[[7]] + (((popM[[8]] - popM[[7]]) * (2068-2060)) / (2070-2060))
UKpop2000_2100[[70]] <- popM[[7]] + (((popM[[8]] - popM[[7]]) * (2069-2060)) / (2070-2060))
UKpop2000_2100[[71]] <- popM[[8]] #2070
UKpop2000_2100[[72]] <- popM[[8]] + (((popM[[9]] - popM[[8]]) * (2071-2070)) / (2080-2070))
UKpop2000_2100[[73]] <- popM[[8]] + (((popM[[9]] - popM[[8]]) * (2072-2070)) / (2080-2070))
UKpop2000_2100[[74]] <- popM[[8]] + (((popM[[9]] - popM[[8]]) * (2073-2070)) / (2080-2070))
UKpop2000_2100[[75]] <- popM[[8]] + (((popM[[9]] - popM[[8]]) * (2074-2070)) / (2080-2070))
UKpop2000_2100[[76]] <- popM[[8]] + (((popM[[9]] - popM[[8]]) * (2075-2070)) / (2080-2070))
UKpop2000_2100[[77]] <- popM[[8]] + (((popM[[9]] - popM[[8]]) * (2076-2070)) / (2080-2070))
UKpop2000_2100[[78]] <- popM[[8]] + (((popM[[9]] - popM[[8]]) * (2077-2070)) / (2080-2070))
UKpop2000_2100[[79]] <- popM[[8]] + (((popM[[9]] - popM[[8]]) * (2078-2070)) / (2080-2070))
UKpop2000_2100[[80]] <- popM[[8]] + (((popM[[9]] - popM[[8]]) * (2079-2070)) / (2080-2070))
UKpop2000_2100[[81]] <- popM[[9]] #2080

years <- seq(as.Date("2000/01/01"), by = "year", length.out = 81)
UKpop2000_2100 <- setZ(UKpop2000_2100, years, "years")
names(UKpop2000_2100) <- years

#################################################
# Getting the proportion of people per grid
#################################################
# To get population weighted average, we need the fraction of the population
# that lives in each grid. So, we need to calculate D = pop in each grid / pop_TS

# UK
UKpop_TS2000_2100 <- c(cellStats(UKpop2000_2100, sum, na.rm=T)) # Total per layer

UKpopD2000_2100 <- brick(nrows=109, ncols=55, xmn=0, xmx=660000, ymn=-84000, ymx=1224000 , nl=81)

for (i in 1:nlayers(UKpop2000_2100)){
  UKpopD2000_2100[[i]] <- (UKpop2000_2100[[i]]  / UKpop_TS2000_2100[i] )
}   # UKpopD2000_2100 is a RasterBrick with the proportion of people living in each grid.


# Population by regions - To calculate regional exposures
# Channel Islands 
Cha_pop2000_2100 <- mask(UKpop2000_2100, Cha_borders) %>% crop(Cha_borders) #NAs (due to size?)

# England
Eng_pop2000_2100 <- mask(UKpop2000_2100, Eng_borders) %>% crop(Eng_borders)
Eng_pop_TS2000_2100 <- c(cellStats(Eng_pop2000_2100, sum, na.rm=T))
Eng_popD2000_2100 <- brick(nrows=55, ncols=48, xmn=84000, xmx=660000, ymn=0, ymx=660000, nl=81)
for (i in 1:nlayers(Eng_pop2000_2100)){
  Eng_popD2000_2100[[i]] <- (Eng_pop2000_2100[[i]]  / Eng_pop_TS2000_2100[i] )
} 

# Isle of Man
Isl_pop2000_2100 <- mask(UKpop2000_2100, Isl_borders) %>% crop(Isl_borders)
Isl_pop_TS2000_2100 <- c(cellStats(Isl_pop2000_2100, sum, na.rm=T))
Isl_popD2000_2100 <- brick(nrows=3, ncols=3, xmn=216000, xmx=252000, ymn=468000, ymx=504000, nl=81)
for (i in 1:nlayers(Isl_pop2000_2100)){
  Isl_popD2000_2100[[i]] <- (Isl_pop2000_2100[[i]]  / Isl_pop_TS2000_2100[i] )
} 

# Northen Ireland
NIr_pop2000_2100 <- mask(UKpop2000_2100, NIr_borders) %>% crop(NIr_borders)
NIr_pop_TS2000_2100 <- c(cellStats(NIr_pop2000_2100, sum, na.rm=T))
NIr_popD2000_2100 <- brick(nrows=12, ncols=15, xmn=0, xmx=180000, ymn=468000, ymx=612000, nl=81)
for (i in 1:nlayers(NIr_pop2000_2100)){
  NIr_popD2000_2100[[i]] <- (NIr_pop2000_2100[[i]]  / NIr_pop_TS2000_2100[i] )
} 

# Scotland
Sco_pop2000_2100 <- mask(UKpop2000_2100, Sco_borders) %>% crop(Sco_borders)
Sco_pop_TS2000_2100 <- c(cellStats(Sco_pop2000_2100, sum, na.rm=T))
Sco_popD2000_2100 <- brick(nrows=58, ncols=39, xmn=0, xmx=468000, ymn=528000, ymx=1224000, nl=81)
for (i in 1:nlayers(Sco_pop2000_2100)){
  Sco_popD2000_2100[[i]] <- (Sco_pop2000_2100[[i]]  / Sco_pop_TS2000_2100[i] )
} 

# Wales
Wal_pop2000_2100 <- mask(UKpop2000_2100, Wal_borders) %>% crop(Wal_borders)
Wal_pop_TS2000_2100 <- c(cellStats(Wal_pop2000_2100, sum, na.rm=T))
Wal_popD2000_2100 <- brick(nrows=19, ncols=17, xmn=156000, xmx=360000, ymn=168000, ymx=396000, nl=81)
for (i in 1:nlayers(Wal_pop2000_2100)){
  Wal_popD2000_2100[[i]] <- (Wal_pop2000_2100[[i]]  / Wal_pop_TS2000_2100[i] )
} 


###############################################################################
#   TEMPERATURE PROJECTIONS
###############################################################################

tasmin.nc <- nc_open("/nfs/cfs/home4/rejb/rejbypa/UKCP18_R/mon/tasmin_rcp85_land-rcm_uk_12km_01_mon_198012-208011.nc")
print(tasmin.nc)
names(tasmin.nc$var) #"tasmin"
tasmin <- brick("/nfs/cfs/home4/rejb/rejbypa/UKCP18_R/mon/tasmin_rcp85_land-rcm_uk_12km_01_mon_198012-208011.nc")
tasmin    # 1200 layers = 100 years

tasmax.nc <- nc_open("/nfs/cfs/home4/rejb/rejbypa/UKCP18_R/mon/tasmax_rcp85_land-rcm_uk_12km_01_mon_198012-208011.nc")
print(tasmax.nc)
names(tasmax.nc$var) #"tasmax"
tasmax <- brick("/nfs/cfs/home4/rejb/rejbypa/UKCP18_R/mon/tasmax_rcp85_land-rcm_uk_12km_01_mon_198012-208011.nc")

tas.nc <- nc_open("/nfs/cfs/home4/rejb/rejbypa/UKCP18_R/mon/tas_rcp85_land-rcm_uk_12km_01_mon_198012-208011.nc")
print(tas.nc)
names(tas.nc$var) #"tasmin"
tas <- brick("/nfs/cfs/home4/rejb/rejbypa/UKCP18_R/mon/tas_rcp85_land-rcm_uk_12km_01_mon_198012-208011.nc")

#################################################
# Masking and dates manipulations
#################################################

# Masking (M) and cropping some grids
tasminM <- mask(tasmin, uk_borders) %>% crop(uk_borders)
tasmaxM <- mask(tasmax, uk_borders) %>% crop(uk_borders)
tasM <- mask(tas, uk_borders) %>% crop(uk_borders)

# Standardized layers' dates - From .nc to rasterbrick there are some problems with layers' names. 
head(names(tasminM)) # "X1980.10.20" "X1980.11.19" "X1980.12.19" "X1981.01.18" "X1981.02.17" "X1981.03.19"
tail(names(tasminM)) # "X2078.11.15" "X2078.12.15" "X2079.01.14" "X2079.02.13" "X2079.03.15" "X2079.04.14"

time_layers <- seq(as.Date("1980/12/01"), as.Date("2080/11/30"), by = "month")

tasminM <- setZ(tasminM, time_layers, "months")
names(tasminM) <- time_layers

tasmaxM <- setZ(tasmaxM, time_layers, "months")
names(tasmaxM) <- time_layers

tasM <- setZ(tasM, time_layers, "months")
names(tasM) <- time_layers

#################################################
# Anomaly calculations
#################################################

# First, we need to calculate averages by month to get "historical" baseline - 1980 Dec (nl= 1) to 2000 Nov (nl= 240)
# We drop the layers we don't need (241:1200)
tasmin1980_2000 <- dropLayer(tasminM, c(241:1200))
tasmax1980_2000 <- dropLayer(tasmaxM, c(241:1200))
tas1980_2000 <- dropLayer(tasM, c(241:1200))

month_layers <- seq(as.Date("2000/12/1"), by = "month", length.out = 12)

# Using zApply, we have monthly means (of 20 Jans, 20 Febs...20 Decs). It stars from Dec-Nov
tasminSA1980_2000  <-  zApply(tasmin1980_2000, by=month, fun= mean, name = "months") 
tasminSA1980_2000  <- setZ(tasminSA1980_2000, month_layers, "months")
names(tasminSA1980_2000) <- month_layers 

tasmaxSA1980_2000  <-  zApply(tasmax1980_2000, by=month, fun= mean, name = "months") 
tasmaxSA1980_2000 <- setZ(tasmaxSA1980_2000, month_layers, "months")
names(tasmaxSA1980_2000) <- month_layers 

tasSA1980_2000  <-  zApply(tas1980_2000, by=month, fun= mean, name = "months")
tasSA1980_2000  <- setZ(tasSA1980_2000, month_layers, "months")
names(tasSA1980_2000) <- month_layers

# For current temperature: We're using from 2000 onward (Dec 1999 (nl= 229))
# We drop layers we do not need (Dec 1980 to Nov 1999)
tasmin2000_2080 <- dropLayer(tasminM, c(1:228))
tasmax2000_2080 <- dropLayer(tasmaxM, c(1:228))
tas2000_2080 <- dropLayer(tasM, c(1:228))

# Now, we need to subtract current temperature minus historical temperature by month
Dec_layers <- seq(as.Date("1999/12/01"), as.Date("2080/11/30"), by = "year")

TasMinAnom2000_2080 <-lapply(1:12, function(i) {
  m <- formatC(c(12,1:11), width = 2, flag= "0")[i]
  TasMinAnom <- (subset(tasmin2000_2080, which(substr(getZ(tasmin2000_2080),6,7) %in% c(m)))) - tasminSA1980_2000[[i]]
  if(m=="12") {Dec_layers <- as.Date(sapply(1999:2079, paste,m,"01",sep= "/"))} else {Dec_layers <- as.Date(sapply(2000:2080, paste,m,"01",sep= "/"))}
  names(TasMinAnom) <- Dec_layers
  return(TasMinAnom)
})

TasMinAnom2000_2080 <- stack(TasMinAnom2000_2080)

TasMaxAnom2000_2080 <-lapply(1:12, function(i) {
  m <- formatC(c(12,1:11), width = 2, flag= "0")[i]
  TasMaxAnom <- (subset(tasmax2000_2080, which(substr(getZ(tasmax2000_2080),6,7) %in% c(m)))) - tasmaxSA1980_2000[[i]]
  if(m=="12") {Dec_layers <- as.Date(sapply(1999:2079, paste,m,"01",sep= "/"))} else {Dec_layers <- as.Date(sapply(2000:2080, paste,m,"01",sep= "/"))}
  names(TasMaxAnom) <- Dec_layers
  return(TasMaxAnom)
})

TasMaxAnom2000_2080 <- stack(TasMaxAnom2000_2080)

TasAnom2000_2080 <-lapply(1:12, function(i) {
  m <- formatC(c(12,1:11), width = 2, flag= "0")[i]
  TasAnom <- (subset(tas2000_2080, which(substr(getZ(tas2000_2080),6,7) %in% c(m)))) - tasSA1980_2000[[i]]
  if(m=="12") {Dec_layers <- as.Date(sapply(1999:2079, paste,m,"01",sep= "/"))} else {Dec_layers <- as.Date(sapply(2000:2080, paste,m,"01",sep= "/"))}
  names(TasAnom) <- Dec_layers
  return(TasAnom)
})

TasAnom2000_2080 <- stack(TasAnom2000_2080)


Dec_layers <- seq(as.Date("1999/12/01"), as.Date("2080/11/30"), by = "year")
Jan_layers <- seq(as.Date("2000/01/01"), as.Date("2080/11/30"), by = "year")
Feb_layers <- seq(as.Date("2000/02/01"), as.Date("2080/11/30"), by = "year")
Mar_layers <- seq(as.Date("2000/03/01"), as.Date("2080/11/30"), by = "year")
Apr_layers <- seq(as.Date("2000/04/01"), as.Date("2080/11/30"), by = "year")
May_layers <- seq(as.Date("2000/05/01"), as.Date("2080/11/30"), by = "year")
Jun_layers <- seq(as.Date("2000/06/01"), as.Date("2080/11/30"), by = "year")
Jul_layers <- seq(as.Date("2000/07/01"), as.Date("2080/11/30"), by = "year")
Aug_layers <- seq(as.Date("2000/08/01"), as.Date("2080/11/30"), by = "year")
Sep_layers <- seq(as.Date("2000/09/01"), as.Date("2080/11/30"), by = "year")
Oct_layers <- seq(as.Date("2000/10/01"), as.Date("2080/11/30"), by = "year")
Nov_layers <- seq(as.Date("2000/11/01"), as.Date("2080/11/30"), by = "year")

Anom2000_2080_layers <- c(Dec_layers, Jan_layers, Feb_layers, Mar_layers, Apr_layers, May_layers,
                          Jun_layers, Jul_layers, Aug_layers, Sep_layers, Oct_layers, Nov_layers)

TasMinAnom2000_2080 <- setZ(TasMinAnom2000_2080 , Anom2000_2080_layers, name= "Date")
TasMaxAnom2000_2080 <- setZ(TasMaxAnom2000_2080, Anom2000_2080_layers, name= "Date")
TasAnom2000_2080 <- setZ(TasAnom2000_2080, Anom2000_2080_layers, name= "Date")

# Sorting raster layers by date
TasMinAnom2000_2080 <- TasMinAnom2000_2080[[order(getZ(TasMinAnom2000_2080))]]
TasMaxAnom2000_2080 <- TasMaxAnom2000_2080[[order(getZ(TasMaxAnom2000_2080))]]
TasAnom2000_2080 <- TasAnom2000_2080[[order(getZ(TasAnom2000_2080))]]

#################################################
# Area-weighted Average Monthly T Anomalies from 2000-2080
#################################################
UK_AWTasMinAnom2000_2080 <- c(cellStats(TasMinAnom2000_2080, mean, na.rm=T))
UK_AWTasMaxAnom2000_2080 <- c(cellStats(TasMaxAnom2000_2080, mean, na.rm=T))
UK_AWTasAnom2000_2080 <- c(cellStats(TasAnom2000_2080, mean, na.rm=T))

# England
Eng_TasMinAnom2000_2080  <- mask(TasMinAnom2000_2080, Eng_borders) %>% crop(Eng_borders)
Eng_TasMaxAnom2000_2080  <- mask(TasMaxAnom2000_2080, Eng_borders) %>% crop(Eng_borders)
Eng_TasAnom2000_2080  <- mask(TasAnom2000_2080, Eng_borders) %>% crop(Eng_borders)

Eng_AWTasMinAnom2000_2080 <- c(cellStats(Eng_TasMinAnom2000_2080, mean, na.rm=T))
Eng_AWTasMaxAnom2000_2080 <- c(cellStats(Eng_TasMaxAnom2000_2080, mean, na.rm=T))
Eng_AWTasAnom2000_2080 <- c(cellStats(Eng_TasAnom2000_2080, mean, na.rm=T))

# Isle of Man
Isl_TasMinAnom2000_2080  <- mask(TasMinAnom2000_2080, Isl_borders) %>% crop(Isl_borders)
Isl_TasMaxAnom2000_2080  <- mask(TasMaxAnom2000_2080, Isl_borders) %>% crop(Isl_borders)
Isl_TasAnom2000_2080  <- mask(TasAnom2000_2080, Isl_borders) %>% crop(Isl_borders)

Isl_AWTasMinAnom2000_2080 <- c(cellStats(Isl_TasMinAnom2000_2080, mean, na.rm=T))
Isl_AWTasMaxAnom2000_2080 <- c(cellStats(Isl_TasMaxAnom2000_2080, mean, na.rm=T))
Isl_AWTasAnom2000_2080 <- c(cellStats(Isl_TasAnom2000_2080, mean, na.rm=T))

# Northen Ireland
NIr_TasMinAnom2000_2080  <- mask(TasMinAnom2000_2080, NIr_borders) %>% crop(NIr_borders)
NIr_TasMaxAnom2000_2080  <- mask(TasMaxAnom2000_2080, NIr_borders) %>% crop(NIr_borders)
NIr_TasAnom2000_2080  <- mask(TasAnom2000_2080, NIr_borders) %>% crop(NIr_borders)

NIr_AWTasMinAnom2000_2080 <- c(cellStats(NIr_TasMinAnom2000_2080, mean, na.rm=T))
NIr_AWTasMaxAnom2000_2080 <- c(cellStats(NIr_TasMaxAnom2000_2080, mean, na.rm=T))
NIr_AWTasAnom2000_2080 <- c(cellStats(NIr_TasAnom2000_2080, mean, na.rm=T))

# Scotland
Sco_TasMinAnom2000_2080  <- mask(TasMinAnom2000_2080, Sco_borders) %>% crop(Sco_borders)
Sco_TasMaxAnom2000_2080  <- mask(TasMaxAnom2000_2080, Sco_borders) %>% crop(Sco_borders)
Sco_TasAnom2000_2080  <- mask(TasAnom2000_2080, Sco_borders) %>% crop(Sco_borders)

Sco_AWTasMinAnom2000_2080 <- c(cellStats(Sco_TasMinAnom2000_2080, mean, na.rm=T))
Sco_AWTasMaxAnom2000_2080 <- c(cellStats(Sco_TasMaxAnom2000_2080, mean, na.rm=T))
Sco_AWTasAnom2000_2080 <- c(cellStats(Sco_TasAnom2000_2080, mean, na.rm=T))

# Wales
Wal_TasMinAnom2000_2080  <- mask(TasMinAnom2000_2080, Wal_borders) %>% crop(Wal_borders)
Wal_TasMaxAnom2000_2080  <- mask(TasMaxAnom2000_2080, Wal_borders) %>% crop(Wal_borders)
Wal_TasAnom2000_2080  <- mask(TasAnom2000_2080, Wal_borders) %>% crop(Wal_borders)

Wal_AWTasMinAnom2000_2080 <- c(cellStats(Wal_TasMinAnom2000_2080, mean, na.rm=T))
Wal_AWTasMaxAnom2000_2080 <- c(cellStats(Wal_TasMaxAnom2000_2080, mean, na.rm=T))
Wal_AWTasAnom2000_2080 <- c(cellStats(Wal_TasAnom2000_2080, mean, na.rm=T))


#################################################
# Population-weighted Average Monthly T Anomalies
#################################################

# UK
UK_PWTasMinAnom2000_2080 <- brick(extend(UKpopD2000_2100, UKpopD2000_2100))
UK_PWTasMaxAnom2000_2080 <- brick(extend(UKpopD2000_2100, UKpopD2000_2100))
UK_PWTasAnom2000_2080 <- brick(extend(UKpopD2000_2100, UKpopD2000_2100))
for (i in 1:nlayers(UKpopD2000_2100)){
  for (j in 1:12) {
    UK_PWTasMinAnom2000_2080[[(12*(i-1))+j]] <- TasMinAnom2000_2080[[(12*(i-1))+j]] * UKpopD2000_2100[[i]]
    UK_PWTasMaxAnom2000_2080[[(12*(i-1))+j]] <- TasMaxAnom2000_2080[[(12*(i-1))+j]] * UKpopD2000_2100[[i]]
    UK_PWTasAnom2000_2080[[(12*(i-1))+j]] <- TasAnom2000_2080[[(12*(i-1))+j]] * UKpopD2000_2100[[i]]
  }
}

UKv_PWTasMinAnom2000_2080 <- c(cellStats(UK_PWTasMinAnom2000_2080, sum, na.rm=T))
UKv_PWTasMaxAnom2000_2080 <- c(cellStats(UK_PWTasMaxAnom2000_2080, sum, na.rm=T))
UKv_PWTasAnom2000_2080 <- c(cellStats(UK_PWTasAnom2000_2080, sum, na.rm=T))

# Eng
Eng_PWTasMinAnom2000_2080 <- brick(extend(Eng_popD2000_2100, Eng_popD2000_2100))
Eng_PWTasMaxAnom2000_2080 <- brick(extend(Eng_popD2000_2100, Eng_popD2000_2100))
Eng_PWTasAnom2000_2080 <- brick(extend(Eng_popD2000_2100, Eng_popD2000_2100))
for (i in 1:nlayers(Eng_popD2000_2100)){
  for (j in 1:12) {
    Eng_PWTasMinAnom2000_2080[[(12*(i-1))+j]] <- Eng_TasMinAnom2000_2080[[(12*(i-1))+j]] * Eng_popD2000_2100[[i]]
    Eng_PWTasMaxAnom2000_2080[[(12*(i-1))+j]] <- Eng_TasMaxAnom2000_2080[[(12*(i-1))+j]] * Eng_popD2000_2100[[i]]
    Eng_PWTasAnom2000_2080[[(12*(i-1))+j]] <- Eng_TasAnom2000_2080[[(12*(i-1))+j]] * Eng_popD2000_2100[[i]]
  }
}

Engv_PWTasMinAnom2000_2080 <- c(cellStats(Eng_PWTasMinAnom2000_2080, sum, na.rm=T))
Engv_PWTasMaxAnom2000_2080 <- c(cellStats(Eng_PWTasMaxAnom2000_2080, sum, na.rm=T))
Engv_PWTasAnom2000_2080 <- c(cellStats(Eng_PWTasAnom2000_2080, sum, na.rm=T))

# NIr
NIr_PWTasMinAnom2000_2080 <- brick(extend(NIr_popD2000_2100, NIr_popD2000_2100))
NIr_PWTasMaxAnom2000_2080 <- brick(extend(NIr_popD2000_2100, NIr_popD2000_2100))
NIr_PWTasAnom2000_2080 <- brick(extend(NIr_popD2000_2100, NIr_popD2000_2100))
for (i in 1:nlayers(NIr_popD2000_2100)){
  for (j in 1:12) {
    NIr_PWTasMinAnom2000_2080[[(12*(i-1))+j]] <- NIr_TasMinAnom2000_2080[[(12*(i-1))+j]] * NIr_popD2000_2100[[i]]
    NIr_PWTasMaxAnom2000_2080[[(12*(i-1))+j]] <- NIr_TasMaxAnom2000_2080[[(12*(i-1))+j]] * NIr_popD2000_2100[[i]]
    NIr_PWTasAnom2000_2080[[(12*(i-1))+j]] <- NIr_TasAnom2000_2080[[(12*(i-1))+j]] * NIr_popD2000_2100[[i]]
  }
}

NIrv_PWTasMinAnom2000_2080 <- c(cellStats(NIr_PWTasMinAnom2000_2080, sum, na.rm=T))
NIrv_PWTasMaxAnom2000_2080 <- c(cellStats(NIr_PWTasMaxAnom2000_2080, sum, na.rm=T))
NIrv_PWTasAnom2000_2080 <- c(cellStats(NIr_PWTasAnom2000_2080, sum, na.rm=T))

# Sco
Sco_PWTasMinAnom2000_2080 <- brick(extend(Sco_popD2000_2100, Sco_popD2000_2100))
Sco_PWTasMaxAnom2000_2080 <- brick(extend(Sco_popD2000_2100, Sco_popD2000_2100))
Sco_PWTasAnom2000_2080 <- brick(extend(Sco_popD2000_2100, Sco_popD2000_2100))
for (i in 1:nlayers(Sco_popD2000_2100)){
  for (j in 1:12) {
    Sco_PWTasMinAnom2000_2080[[(12*(i-1))+j]] <- Sco_TasMinAnom2000_2080[[(12*(i-1))+j]] * Sco_popD2000_2100[[i]]
    Sco_PWTasMaxAnom2000_2080[[(12*(i-1))+j]] <- Sco_TasMaxAnom2000_2080[[(12*(i-1))+j]] * Sco_popD2000_2100[[i]]
    Sco_PWTasAnom2000_2080[[(12*(i-1))+j]] <- Sco_TasAnom2000_2080[[(12*(i-1))+j]] * Sco_popD2000_2100[[i]]
  }
}

Scov_PWTasMinAnom2000_2080 <- c(cellStats(Sco_PWTasMinAnom2000_2080, sum, na.rm=T))
Scov_PWTasMaxAnom2000_2080 <- c(cellStats(Sco_PWTasMaxAnom2000_2080, sum, na.rm=T))
Scov_PWTasAnom2000_2080 <- c(cellStats(Sco_PWTasAnom2000_2080, sum, na.rm=T))


# Wal
Wal_PWTasMinAnom2000_2080 <- brick(extend(Wal_popD2000_2100, Wal_popD2000_2100))
Wal_PWTasMaxAnom2000_2080 <- brick(extend(Wal_popD2000_2100, Wal_popD2000_2100))
Wal_PWTasAnom2000_2080 <- brick(extend(Wal_popD2000_2100, Wal_popD2000_2100))
for (i in 1:nlayers(Wal_popD2000_2100)){
  for (j in 1:12) {
    Wal_PWTasMinAnom2000_2080[[(12*(i-1))+j]] <- Wal_TasMinAnom2000_2080[[(12*(i-1))+j]] * Wal_popD2000_2100[[i]]
    Wal_PWTasMaxAnom2000_2080[[(12*(i-1))+j]] <- Wal_TasMaxAnom2000_2080[[(12*(i-1))+j]] * Wal_popD2000_2100[[i]]
    Wal_PWTasAnom2000_2080[[(12*(i-1))+j]] <- Wal_TasAnom2000_2080[[(12*(i-1))+j]] * Wal_popD2000_2100[[i]]
  }
}

Walv_PWTasMinAnom2000_2080 <- c(cellStats(Wal_PWTasMinAnom2000_2080, sum, na.rm=T))
Walv_PWTasMaxAnom2000_2080 <- c(cellStats(Wal_PWTasMaxAnom2000_2080, sum, na.rm=T))
Walv_PWTasAnom2000_2080 <- c(cellStats(Wal_PWTasAnom2000_2080, sum, na.rm=T))


write.csv(as.data.frame(cbind(UK_AWTasMinAnom2000_2080, UK_AWTasMaxAnom2000_2080, UK_AWTasAnom2000_2080,
                              Eng_AWTasMinAnom2000_2080, Eng_AWTasMaxAnom2000_2080, Eng_AWTasAnom2000_2080, 
                              NIr_AWTasMinAnom2000_2080, NIr_AWTasMaxAnom2000_2080, NIr_AWTasAnom2000_2080,
                              Sco_AWTasMinAnom2000_2080, Sco_AWTasMaxAnom2000_2080, Sco_AWTasAnom2000_2080,
                              Wal_AWTasMinAnom2000_2080, Wal_AWTasMaxAnom2000_2080, Wal_AWTasAnom2000_2080,
                              UKv_PWTasMinAnom2000_2080, UKv_PWTasMaxAnom2000_2080, UKv_PWTasAnom2000_2080,
                              Engv_PWTasMinAnom2000_2080, Engv_PWTasMaxAnom2000_2080, Engv_PWTasAnom2000_2080, 
                              NIrv_PWTasMinAnom2000_2080, NIrv_PWTasMaxAnom2000_2080, NIrv_PWTasAnom2000_2080,
                              Scov_PWTasMinAnom2000_2080, Scov_PWTasMaxAnom2000_2080, Scov_PWTasAnom2000_2080,
                              Walv_PWTasMinAnom2000_2080, Walv_PWTasMaxAnom2000_2080, Walv_PWTasAnom2000_2080)), "AW_PW_MonthlyTAnomalies2000_2080.csv")



###############################################################################
#   VISUALISATIONS
###############################################################################

#################################################
# Area-weighted Average Monthly T Anomalies from 2000-2080
#################################################

AW_PW_MonthlyTAnomalies2000_2080 <- read.csv("AW_PW_MonthlyTAnomalies2000_2080.csv")
AW_PW_MonthlyTAnomalies2000_2080$date <- as.Date(seq(as.Date("1999/12/01"), as.Date("2080/11/01"), by = "month"))

#TasMin
AW_PW_MonthlyTAnomalies2000_2080$S_UK_AWTasMinAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$UK_AWTasMinAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
UK_TasMin <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, UK_AWTasMinAnom2000_2080, fill= S_UK_AWTasMinAnom2000_2080)) + 
  geom_bar(stat = "identity", show.legend = FALSE) + 
  #scale_x_date(date_breaks = "years", date_labels = "%Y") +
  scale_y_continuous(breaks = seq(-4, 10, 2)) +
  scale_fill_manual(values = c("#99000d", "#034e7b")) +
  labs(y = "Temperature anomaly (°C)", x = "", 
       title= "Monthly Minimum Temperature anomaly (°C) in the UK 2000-2080", subtitle= "compared to 1980-2000 average",
       caption= "Data: UKCP18") +
  theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Eng_AWTasMinAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Eng_AWTasMinAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Eng_TasMin <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Eng_AWTasMinAnom2000_2080, fill= S_Eng_AWTasMinAnom2000_2080)) + 
                geom_bar(stat = "identity", show.legend = FALSE) + 
                #scale_x_date(date_breaks = "years", date_labels = "%Y") +
                scale_y_continuous(breaks = seq(-4, 10, 2)) +
                scale_fill_manual(values = c("#99000d", "#034e7b")) +
                labs(y = "Temperature anomaly (°C)", x = "", 
                     title= "Monthly Minimum Temperature anomaly (°C) in England 2000-2080", subtitle= "compared to 1980-2000 average",
                     caption= "Data: UKCP18") +
                theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_NIr_AWTasMinAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$NIr_AWTasMinAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
NIr_TasMin <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, NIr_AWTasMinAnom2000_2080, fill= S_NIr_AWTasMinAnom2000_2080)) + 
                geom_bar(stat = "identity", show.legend = FALSE) + 
                scale_y_continuous(breaks = seq(-4, 10, 2)) +
                scale_fill_manual(values = c("#99000d", "#034e7b")) +
                labs(y = "Temperature anomaly (°C)", x = "", 
                     title= "Monthly Minimum Temperature anomaly (°C) in the Northern Ireland 2000-2080", subtitle= "compared to 1980-2000 average",
                     caption= "Data: UKCP18") +
                theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Sco_AWTasMinAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Sco_AWTasMinAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Sco_TasMin <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Sco_AWTasMinAnom2000_2080, fill= S_Sco_AWTasMinAnom2000_2080)) + 
                geom_bar(stat = "identity", show.legend = FALSE) + 
                scale_y_continuous(breaks = seq(-4, 10, 2)) +
                scale_fill_manual(values = c("#99000d", "#034e7b")) +
                labs(y = "Temperature anomaly (°C)", x = "", 
                     title= "Monthly Minimum Temperature anomaly (°C) in Scotland 2000-2080", subtitle= "compared to 1980-2000 average",
                     caption= "Data: UKCP18") +
                theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Wal_AWTasMinAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Wal_AWTasMinAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Wal_TasMin <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Wal_AWTasMinAnom2000_2080, fill= S_Wal_AWTasMinAnom2000_2080)) + 
                geom_bar(stat = "identity", show.legend = FALSE) + 
                scale_y_continuous(breaks = seq(-4, 10, 2)) +
                scale_fill_manual(values = c("#99000d", "#034e7b")) +
                labs(y = "Temperature anomaly (°C)", x = "", 
                     title= "Monthly Minimum Temperature anomaly (°C) in Wales 2000-2080", subtitle= "compared to 1980-2000 average",
                     caption= "Data: UKCP18") +
                theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                theme_hc())


p1 <- ggarrange(UK_TasMin, ncol=1, nrow=1)
p2 <- ggarrange(Eng_TasMin, NIr_TasMin, Sco_TasMin, Wal_TasMin, ncol=2, nrow=2)
p3 <- ggarrange(p1, p2, ncol=1, nrow= 2)


#TasMax
AW_PW_MonthlyTAnomalies2000_2080$S_UK_AWTasMaxAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$UK_AWTasMaxAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
UK_TasMax <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, UK_AWTasMaxAnom2000_2080, fill= S_UK_AWTasMaxAnom2000_2080)) + 
                geom_bar(stat = "identity", show.legend = FALSE) + 
                #scale_x_date(date_breaks = "years", date_labels = "%Y") +
                scale_y_continuous(breaks = seq(-4, 10, 2)) +
                scale_fill_manual(values = c("#99000d", "#034e7b")) +
                labs(y = "Temperature anomaly (°C)", x = "", 
                     title= "Monthly Maximum Temperature anomaly (°C) in the UK 2000-2080", subtitle= "compared to 1980-2000 average",
                     caption= "Data: UKCP18") +
                theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Eng_AWTasMaxAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Eng_AWTasMaxAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Eng_TasMax <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Eng_AWTasMaxAnom2000_2080, fill= S_Eng_AWTasMaxAnom2000_2080)) + 
                 geom_bar(stat = "identity", show.legend = FALSE) + 
                 #scale_x_date(date_breaks = "years", date_labels = "%Y") +
                 scale_y_continuous(breaks = seq(-4, 10, 2)) +
                 scale_fill_manual(values = c("#99000d", "#034e7b")) +
                 labs(y = "Temperature anomaly (°C)", x = "", 
                      title= "Monthly Maximum Temperature anomaly (°C) in England 2000-2080", subtitle= "compared to 1980-2000 average",
                      caption= "Data: UKCP18") +
                 theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                 theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_NIr_AWTasMaxAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$NIr_AWTasMaxAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
NIr_TasMax <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, NIr_AWTasMaxAnom2000_2080, fill= S_NIr_AWTasMaxAnom2000_2080)) + 
                 geom_bar(stat = "identity", show.legend = FALSE) + 
                 scale_y_continuous(breaks = seq(-4, 10, 2)) +
                 scale_fill_manual(values = c("#99000d", "#034e7b")) +
                 labs(y = "Temperature anomaly (°C)", x = "", 
                      title= "Monthly Maximum Temperature anomaly (°C) in the Northern Ireland 2000-2080", subtitle= "compared to 1980-2000 average",
                      caption= "Data: UKCP18") +
                 theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                 theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Sco_AWTasMaxAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Sco_AWTasMaxAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Sco_TasMax <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Sco_AWTasMaxAnom2000_2080, fill= S_Sco_AWTasMaxAnom2000_2080)) + 
                 geom_bar(stat = "identity", show.legend = FALSE) + 
                 scale_y_continuous(breaks = seq(-4, 10, 2)) +
                 scale_fill_manual(values = c("#99000d", "#034e7b")) +
                 labs(y = "Temperature anomaly (°C)", x = "", 
                      title= "Monthly Maximum Temperature anomaly (°C) in Scotland 2000-2080", subtitle= "compared to 1980-2000 average",
                      caption= "Data: UKCP18") +
                 theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                 theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Wal_AWTasMaxAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Wal_AWTasMaxAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Wal_TasMax <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Wal_AWTasMaxAnom2000_2080, fill= S_Wal_AWTasMaxAnom2000_2080)) + 
                 geom_bar(stat = "identity", show.legend = FALSE) + 
                 scale_y_continuous(breaks = seq(-4, 10, 2)) +
                 scale_fill_manual(values = c("#99000d", "#034e7b")) +
                 labs(y = "Temperature anomaly (°C)", x = "", 
                      title= "Monthly Maximum Temperature anomaly (°C) in Wales 2000-2080", subtitle= "compared to 1980-2000 average",
                      caption= "Data: UKCP18") +
                 theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                 theme_hc())


p4 <- ggarrange(UK_TasMax, ncol=1, nrow=1)
p5 <- ggarrange(Eng_TasMax, NIr_TasMax, Sco_TasMax, Wal_TasMax, ncol=2, nrow=2)
p6 <- ggarrange(p4, p5, ncol=1, nrow= 2)


#Tas
AW_PW_MonthlyTAnomalies2000_2080$S_UK_AWTasAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$UK_AWTasAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
UK_Tas <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, UK_AWTasAnom2000_2080, fill= S_UK_AWTasAnom2000_2080)) + 
                geom_bar(stat = "identity", show.legend = FALSE) + 
                #scale_x_date(date_breaks = "years", date_labels = "%Y") +
                scale_y_continuous(breaks = seq(-4, 10, 2)) +
                scale_fill_manual(values = c("#99000d", "#034e7b")) +
                labs(y = "Temperature anomaly (°C)", x = "", 
                     title= "Monthly Mean Temperature anomaly (°C) in the UK 2000-2080", subtitle= "compared to 1980-2000 average",
                     caption= "Data: UKCP18") +
                theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Eng_AWTasAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Eng_AWTasAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Eng_Tas <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Eng_AWTasAnom2000_2080, fill= S_Eng_AWTasAnom2000_2080)) + 
                 geom_bar(stat = "identity", show.legend = FALSE) + 
                 #scale_x_date(date_breaks = "years", date_labels = "%Y") +
                 scale_y_continuous(breaks = seq(-4, 10, 2)) +
                 scale_fill_manual(values = c("#99000d", "#034e7b")) +
                 labs(y = "Temperature anomaly (°C)", x = "", 
                      title= "Monthly Mean Temperature anomaly (°C) in England 2000-2080", subtitle= "compared to 1980-2000 average",
                      caption= "Data: UKCP18") +
                 theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                 theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_NIr_AWTasAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$NIr_AWTasAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
NIr_Tas <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, NIr_AWTasAnom2000_2080, fill= S_NIr_AWTasAnom2000_2080)) + 
                 geom_bar(stat = "identity", show.legend = FALSE) + 
                 scale_y_continuous(breaks = seq(-4, 10, 2)) +
                 scale_fill_manual(values = c("#99000d", "#034e7b")) +
                 labs(y = "Temperature anomaly (°C)", x = "", 
                      title= "Monthly Mean Temperature anomaly (°C) in the Northern Ireland 2000-2080", subtitle= "compared to 1980-2000 average",
                      caption= "Data: UKCP18") +
                 theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                 theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Sco_AWTasAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Sco_AWTasAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Sco_Tas <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Sco_AWTasAnom2000_2080, fill= S_Sco_AWTasAnom2000_2080)) + 
                 geom_bar(stat = "identity", show.legend = FALSE) + 
                 scale_y_continuous(breaks = seq(-4, 10, 2)) +
                 scale_fill_manual(values = c("#99000d", "#034e7b")) +
                 labs(y = "Temperature anomaly (°C)", x = "", 
                      title= "Monthly Mean Temperature anomaly (°C) in Scotland 2000-2080", subtitle= "compared to 1980-2000 average",
                      caption= "Data: UKCP18") +
                 theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                 theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Wal_AWTasAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Wal_AWTasAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Wal_Tas <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Wal_AWTasAnom2000_2080, fill= S_Wal_AWTasAnom2000_2080)) + 
                 geom_bar(stat = "identity", show.legend = FALSE) + 
                 scale_y_continuous(breaks = seq(-4, 10, 2)) +
                 scale_fill_manual(values = c("#99000d", "#034e7b")) +
                 labs(y = "Temperature anomaly (°C)", x = "", 
                      title= "Monthly Mean Temperature anomaly (°C) in Wales 2000-2080", subtitle= "compared to 1980-2000 average",
                      caption= "Data: UKCP18") +
                 theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                 theme_hc())


p7 <- ggarrange(UK_Tas, ncol=1, nrow=1)
p8 <- ggarrange(Eng_Tas, NIr_Tas, Sco_Tas, Wal_Tas, ncol=2, nrow=2)
p9 <- ggarrange(p7, p8, ncol=1, nrow= 2)


#################################################
# Population-weighted Average Monthly Anomalies
#################################################

#TasMin
AW_PW_MonthlyTAnomalies2000_2080$S_UK_PWTasMinAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$UKv_PWTasMinAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
UK_PWTasMin <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, UKv_PWTasMinAnom2000_2080, fill= S_UK_PWTasMinAnom2000_2080)) + 
                geom_bar(stat = "identity", show.legend = FALSE) + 
                #scale_x_date(date_breaks = "years", date_labels = "%Y") +
                scale_y_continuous(breaks = seq(-4, 10, 2)) +
                scale_fill_manual(values = c("#99000d", "#034e7b")) +
                labs(y = "Temperature anomaly (°C)", x = "", 
                     title= "Monthly Population-weighted Minimum Temperature anomaly (°C) in the UK 2000-2080", subtitle= "compared to 1980-2000 average",
                     caption= "Data: UKCP18") +
                theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Eng_PWTasMinAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Engv_PWTasMinAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Eng_PWTasMin <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Eng_AWTasMinAnom2000_2080, fill= S_Eng_PWTasMinAnom2000_2080)) + 
                 geom_bar(stat = "identity", show.legend = FALSE) + 
                 #scale_x_date(date_breaks = "years", date_labels = "%Y") +
                 scale_y_continuous(breaks = seq(-4, 10, 2)) +
                 scale_fill_manual(values = c("#99000d", "#034e7b")) +
                 labs(y = "Temperature anomaly (°C)", x = "", 
                      title= "Monthly Population-weighted Minimum Temperature anomaly (°C) in England 2000-2080", subtitle= "compared to 1980-2000 average",
                      caption= "Data: UKCP18") +
                 theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                 theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_NIr_PWTasMinAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$NIrv_PWTasMinAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
NIr_PWTasMin <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, NIrv_PWTasMinAnom2000_2080, fill= S_NIr_PWTasMinAnom2000_2080)) + 
                 geom_bar(stat = "identity", show.legend = FALSE) + 
                 scale_y_continuous(breaks = seq(-4, 10, 2)) +
                 scale_fill_manual(values = c("#99000d", "#034e7b")) +
                 labs(y = "Temperature anomaly (°C)", x = "", 
                      title= "Monthly Population-weighted Minimum Temperature anomaly (°C) in the Northern Ireland 2000-2080", subtitle= "compared to 1980-2000 average",
                      caption= "Data: UKCP18") +
                 theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                 theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Sco_PWTasMinAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Scov_PWTasMinAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Sco_PWTasMin <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Scov_PWTasMinAnom2000_2080, fill= S_Sco_PWTasMinAnom2000_2080)) + 
                 geom_bar(stat = "identity", show.legend = FALSE) + 
                 scale_y_continuous(breaks = seq(-4, 10, 2)) +
                 scale_fill_manual(values = c("#99000d", "#034e7b")) +
                 labs(y = "Temperature anomaly (°C)", x = "", 
                      title= "Monthly Population-weighted Minimum Temperature anomaly (°C) in Scotland 2000-2080", subtitle= "compared to 1980-2000 average",
                      caption= "Data: UKCP18") +
                 theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                 theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Wal_PWTasMinAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Walv_PWTasMinAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Wal_PWTasMin <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Walv_PWTasMinAnom2000_2080, fill= S_Wal_PWTasMinAnom2000_2080)) + 
                 geom_bar(stat = "identity", show.legend = FALSE) + 
                 scale_y_continuous(breaks = seq(-4, 10, 2)) +
                 scale_fill_manual(values = c("#99000d", "#034e7b")) +
                 labs(y = "Temperature anomaly (°C)", x = "", 
                      title= "Monthly Population-weighted Minimum Temperature anomaly (°C) in Wales 2000-2080", subtitle= "compared to 1980-2000 average",
                      caption= "Data: UKCP18") +
                 theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                 theme_hc())


p10 <- ggarrange(UK_PWTasMin, ncol=1, nrow=1)
p11 <- ggarrange(Eng_PWTasMin, NIr_PWTasMin, Sco_PWTasMin, Wal_PWTasMin, ncol=2, nrow=2)
p12 <- ggarrange(p10, p11, ncol=1, nrow= 2)



#TasMax
AW_PW_MonthlyTAnomalies2000_2080$S_UK_PWTasMaxAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$UKv_PWTasMaxAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
UK_PWTasMax <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, UKv_PWTasMaxAnom2000_2080, fill= S_UK_PWTasMaxAnom2000_2080)) + 
                geom_bar(stat = "identity", show.legend = FALSE) + 
                #scale_x_date(date_breaks = "years", date_labels = "%Y") +
                scale_y_continuous(breaks = seq(-4, 10, 2)) +
                scale_fill_manual(values = c("#99000d", "#034e7b")) +
                labs(y = "Temperature anomaly (°C)", x = "", 
                     title= "Monthly Population-weighted Maximum Temperature anomaly (°C) in the UK 2000-2080", subtitle= "compared to 1980-2000 average",
                     caption= "Data: UKCP18") +
                theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Eng_PWTasMaxAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Engv_PWTasMaxAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Eng_PWTasMax <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Engv_PWTasMaxAnom2000_2080, fill= S_Eng_PWTasMaxAnom2000_2080)) + 
                 geom_bar(stat = "identity", show.legend = FALSE) + 
                 #scale_x_date(date_breaks = "years", date_labels = "%Y") +
                 scale_y_continuous(breaks = seq(-4, 10, 2)) +
                 scale_fill_manual(values = c("#99000d", "#034e7b")) +
                 labs(y = "Temperature anomaly (°C)", x = "", 
                      title= "Monthly Population-weighted Maximum Temperature anomaly (°C) in England 2000-2080", subtitle= "compared to 1980-2000 average",
                      caption= "Data: UKCP18") +
                 theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                 theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_NIr_PWTasMaxAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$NIrv_PWTasMaxAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
NIr_PWTasMax <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, NIrv_PWTasMaxAnom2000_2080, fill= S_NIr_PWTasMaxAnom2000_2080)) + 
                 geom_bar(stat = "identity", show.legend = FALSE) + 
                 scale_y_continuous(breaks = seq(-4, 10, 2)) +
                 scale_fill_manual(values = c("#99000d", "#034e7b")) +
                 labs(y = "Temperature anomaly (°C)", x = "", 
                      title= "Monthly Population-weighted Maximum Temperature anomaly (°C) in the Northern Ireland 2000-2080", subtitle= "compared to 1980-2000 average",
                      caption= "Data: UKCP18") +
                 theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                 theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Sco_PWTasMaxAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Scov_PWTasMaxAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Sco_PWTasMax <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Scov_PWTasMaxAnom2000_2080, fill= S_Sco_PWTasMaxAnom2000_2080)) + 
                 geom_bar(stat = "identity", show.legend = FALSE) + 
                 scale_y_continuous(breaks = seq(-4, 10, 2)) +
                 scale_fill_manual(values = c("#99000d", "#034e7b")) +
                 labs(y = "Temperature anomaly (°C)", x = "", 
                      title= "Monthly Population-weighted Maximum Temperature anomaly (°C) in Scotland 2000-2080", subtitle= "compared to 1980-2000 average",
                      caption= "Data: UKCP18") +
                 theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                 theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Wal_PWTasMaxAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Walv_PWTasMaxAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Wal_PWTasMax <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Walv_PWTasMaxAnom2000_2080, fill= S_Wal_PWTasMaxAnom2000_2080)) + 
                 geom_bar(stat = "identity", show.legend = FALSE) + 
                 scale_y_continuous(breaks = seq(-4, 10, 2)) +
                 scale_fill_manual(values = c("#99000d", "#034e7b")) +
                 labs(y = "Temperature anomaly (°C)", x = "", 
                      title= "Monthly Population-weighted Maximum Temperature anomaly (°C) in Wales 2000-2080", subtitle= "compared to 1980-2000 average",
                      caption= "Data: UKCP18") +
                 theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
                 theme_hc())


p13 <- ggarrange(UK_PWTasMax, ncol=1, nrow=1)
p14 <- ggarrange(Eng_PWTasMax, NIr_PWTasMax, Sco_PWTasMax, Wal_PWTasMax, ncol=2, nrow=2)
p15 <- ggarrange(p13, p14, ncol=1, nrow= 2)


#Tas
AW_PW_MonthlyTAnomalies2000_2080$S_UK_PWTasAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$UKv_PWTasAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
UK_PWTas <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, UKv_PWTasAnom2000_2080, fill= S_UK_PWTasAnom2000_2080)) + 
             geom_bar(stat = "identity", show.legend = FALSE) + 
             #scale_x_date(date_breaks = "years", date_labels = "%Y") +
             scale_y_continuous(breaks = seq(-4, 10, 2)) +
             scale_fill_manual(values = c("#99000d", "#034e7b")) +
             labs(y = "Temperature anomaly (°C)", x = "", 
                  title= "Monthly Population-weighted Mean Temperature anomaly (°C) in the UK 2000-2080", subtitle= "compared to 1980-2000 average",
                  caption= "Data: UKCP18") +
             theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Eng_PWTasAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Engv_PWTasAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Eng_PWTas <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Engv_PWTasAnom2000_2080, fill= S_Eng_PWTasAnom2000_2080)) + 
              geom_bar(stat = "identity", show.legend = FALSE) + 
              #scale_x_date(date_breaks = "years", date_labels = "%Y") +
              scale_y_continuous(breaks = seq(-4, 10, 2)) +
              scale_fill_manual(values = c("#99000d", "#034e7b")) +
              labs(y = "Temperature anomaly (°C)", x = "", 
                   title= "Monthly Population-weighted Mean Temperature anomaly (°C) in England 2000-2080", subtitle= "compared to 1980-2000 average",
                   caption= "Data: UKCP18") +
              theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
              theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_NIr_PWTasAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$NIrv_PWTasAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
NIr_PWTas <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, NIrv_PWTasAnom2000_2080, fill= S_NIr_PWTasAnom2000_2080)) + 
              geom_bar(stat = "identity", show.legend = FALSE) + 
              scale_y_continuous(breaks = seq(-4, 10, 2)) +
              scale_fill_manual(values = c("#99000d", "#034e7b")) +
              labs(y = "Temperature anomaly (°C)", x = "", 
                   title= "Monthly Population-weighted Mean Temperature anomaly (°C) in the Northern Ireland 2000-2080", subtitle= "compared to 1980-2000 average",
                   caption= "Data: UKCP18") +
              theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
              theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Sco_PWTasAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Scov_PWTasAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Sco_PWTas <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Scov_PWTasAnom2000_2080, fill= S_Sco_PWTasAnom2000_2080)) + 
              geom_bar(stat = "identity", show.legend = FALSE) + 
              scale_y_continuous(breaks = seq(-4, 10, 2)) +
              scale_fill_manual(values = c("#99000d", "#034e7b")) +
              labs(y = "Temperature anomaly (°C)", x = "", 
                   title= "Monthly Population-weighted Mean Temperature anomaly (°C) in Scotland 2000-2080", subtitle= "compared to 1980-2000 average",
                   caption= "Data: UKCP18") +
              theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
              theme_hc())

AW_PW_MonthlyTAnomalies2000_2080$S_Wal_PWTasAnom2000_2080 <- ifelse(AW_PW_MonthlyTAnomalies2000_2080$Walv_PWTasAnom2000_2080 >0, "pos","neg") %>% factor(c("pos", "neg"))
Wal_PWTas <- (ggplot(AW_PW_MonthlyTAnomalies2000_2080, aes(date, Walv_PWTasAnom2000_2080, fill= S_Wal_PWTasAnom2000_2080)) + 
              geom_bar(stat = "identity", show.legend = FALSE) + 
              scale_y_continuous(breaks = seq(-4, 10, 2)) +
              scale_fill_manual(values = c("#99000d", "#034e7b")) +
              labs(y = "Temperature anomaly (°C)", x = "", 
                   title= "Monthly Population-weighted Mean Temperature anomaly (°C) in Wales 2000-2080", subtitle= "compared to 1980-2000 average",
                   caption= "Data: UKCP18") +
              theme(plot.title= element_text(size= 10), plot.subtitle = element_text(size=9)) +
              theme_hc())


p16 <- ggarrange(UK_PWTas, ncol=1, nrow=1)
p17 <- ggarrange(Eng_PWTas, NIr_PWTas, Sco_PWTas, Wal_PWTas, ncol=2, nrow=2)
p18 <- ggarrange(p16, p17, ncol=1, nrow= 2)


#################################################
# Time series plots for AW an PW
#################################################
data.frame(colnames(AW_PW_MonthlyTAnomalies2000_2080))

layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE))

AW_PW_MonthlyTAnomalies2000_2080$Year <- as.numeric(substr(AW_PW_MonthlyTAnomalies2000_2080$date,1,4))
AW_PW_MonthlyTAnomalies2000_2080$Month <- as.numeric(substr(AW_PW_MonthlyTAnomalies2000_2080$date,6,7))
AW_PW_MonthlyTAnomalies2000_2080$DecimalDate <- AW_PW_MonthlyTAnomalies2000_2080$Year + (AW_PW_MonthlyTAnomalies2000_2080$Month / 12) - 1/24

#TasMin
# UK
UK_TasMin.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$UK_AWTasMinAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
UK_PWTasMin.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$UKv_PWTasMinAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
matplot(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, AW_PW_MonthlyTAnomalies2000_2080[,c(2,17)], type="l", 
        xlab="Year", main= "Monthly Minimum Temperature anomaly (°C) in the UK 2000-2080", ylab=expression(degree * "C"), 
        col=c("lightblue", "salmon"), lwd=2, lty=1)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, UK_TasMin.decomp$time.series[,2], col="blue", lwd=3)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, UK_PWTasMin.decomp$time.series[,2], col="red", lwd=3)
legend("bottomright", legend= c("Area Weighted", "Population Weighted"), lwd= 2, col=c("lightblue", "salmon"), bg= ("white"), horiz= F)

# England
Eng_TasMin.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Eng_AWTasMinAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
Eng_PWTasMin.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Engv_PWTasMinAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
matplot(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, AW_PW_MonthlyTAnomalies2000_2080[,c(5,20)], type="l", 
        xlab="Year", main= "Monthly Minimum Temperature anomaly (°C) in the England 2000-2080", ylab=expression(degree * "C"), 
        col=c("lightblue", "salmon"), lwd=2, lty=1)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Eng_TasMin.decomp$time.series[,2], col="blue", lwd=3)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Eng_PWTasMin.decomp$time.series[,2], col="red", lwd=3)

# NIr
NIr_TasMin.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$NIr_AWTasMinAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
NIr_PWTasMin.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$NIrv_PWTasMinAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
matplot(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, AW_PW_MonthlyTAnomalies2000_2080[,c(8,23)], type="l", 
        xlab="Year", main= "Monthly Minimum Temperature anomaly (°C) in the N.I. 2000-2080", ylab=expression(degree * "C"), 
        col=c("lightblue", "salmon"), lwd=2, lty=1)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, NIr_TasMin.decomp$time.series[,2], col="blue", lwd=3)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, NIr_PWTasMin.decomp$time.series[,2], col="red", lwd=3)

# Scotland
Sco_TasMin.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Sco_AWTasMinAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
Sco_PWTasMin.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Scov_PWTasMinAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
matplot(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, AW_PW_MonthlyTAnomalies2000_2080[,c(11,26)], type="l", 
        xlab="Year", main= "Monthly Minimum Temperature anomaly (°C) in the Scotland 2000-2080", ylab=expression(degree * "C"), 
        col=c("lightblue", "salmon"), lwd=2, lty=1)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Sco_TasMin.decomp$time.series[,2], col="blue", lwd=3)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Sco_PWTasMin.decomp$time.series[,2], col="red", lwd=3)

# Wales
Wal_TasMin.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Wal_AWTasMinAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
Wal_PWTasMin.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Walv_PWTasMinAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
matplot(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, AW_PW_MonthlyTAnomalies2000_2080[,c(14,29)], type="l", 
        xlab="Year", main= "Monthly Minimum Temperature anomaly (°C) in the Wales 2000-2080", ylab=expression(degree * "C"), 
        col=c("lightblue", "salmon"), lwd=2, lty=1)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Wal_TasMin.decomp$time.series[,2], col="blue", lwd=3)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Wal_PWTasMin.decomp$time.series[,2], col="red", lwd=3)


#TasMax
# UK
UK_TasMax.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$UK_AWTasMaxAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
UK_PWTasMax.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$UKv_PWTasMaxAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
matplot(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, AW_PW_MonthlyTAnomalies2000_2080[,c(3,18)], type="l", 
        xlab="Year", main= "Monthly Maximum Temperature anomaly (°C) in the UK 2000-2080", ylab=expression(degree * "C"), 
        col=c("lightblue", "salmon"), lwd=2, lty=1)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, UK_TasMax.decomp$time.series[,2], col="blue", lwd=3)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, UK_PWTasMax.decomp$time.series[,2], col="red", lwd=3)
legend("bottomright", legend= c("Area Weighted", "Population Weighted"), lwd= 2, col=c("lightblue", "salmon"), bg= ("white"), horiz= F)

# England
Eng_TasMax.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Eng_AWTasMaxAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
Eng_PWTasMax.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Engv_PWTasMaxAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
matplot(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, AW_PW_MonthlyTAnomalies2000_2080[,c(6,21)], type="l", 
        xlab="Year", main= "Monthly Maximum Temperature anomaly (°C) in the England 2000-2080", ylab=expression(degree * "C"), 
        col=c("lightblue", "salmon"), lwd=2, lty=1)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Eng_TasMax.decomp$time.series[,2], col="blue", lwd=3)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Eng_PWTasMax.decomp$time.series[,2], col="red", lwd=3)

# NIr
NIr_TasMax.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$NIr_AWTasMaxAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
NIr_PWTasMax.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$NIrv_PWTasMaxAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
matplot(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, AW_PW_MonthlyTAnomalies2000_2080[,c(9,24)], type="l", 
        xlab="Year", main= "Monthly Maximum Temperature anomaly (°C) in the N.I. 2000-2080", ylab=expression(degree * "C"), 
        col=c("lightblue", "salmon"), lwd=2, lty=1)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, NIr_TasMax.decomp$time.series[,2], col="blue", lwd=3)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, NIr_PWTasMax.decomp$time.series[,2], col="red", lwd=3)

# Scotland
Sco_TasMax.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Sco_AWTasMaxAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
Sco_PWTasMax.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Scov_PWTasMaxAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
matplot(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, AW_PW_MonthlyTAnomalies2000_2080[,c(12,27)], type="l", 
        xlab="Year", main= "Monthly Maximum Temperature anomaly (°C) in the Scotland 2000-2080", ylab=expression(degree * "C"), 
        col=c("lightblue", "salmon"), lwd=2, lty=1)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Sco_TasMax.decomp$time.series[,2], col="blue", lwd=3)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Sco_PWTasMax.decomp$time.series[,2], col="red", lwd=3)

# Wales
Wal_TasMax.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Wal_AWTasMaxAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
Wal_PWTasMax.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Walv_PWTasMaxAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
matplot(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, AW_PW_MonthlyTAnomalies2000_2080[,c(15,30)], type="l", 
        xlab="Year", main= "Monthly Maximum Temperature anomaly (°C) in the Wales 2000-2080", ylab=expression(degree * "C"), 
        col=c("lightblue", "salmon"), lwd=2, lty=1)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Wal_TasMax.decomp$time.series[,2], col="blue", lwd=3)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Wal_PWTasMax.decomp$time.series[,2], col="red", lwd=3)


#Tas
# UK
UK_Tas.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$UK_AWTasAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
UK_PWTas.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$UKv_PWTasAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
matplot(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, AW_PW_MonthlyTAnomalies2000_2080[,c(4,19)], type="l", 
        xlab="Year", main= "Monthly Mean Temperature anomaly (°C) in the UK 2000-2080", ylab=expression(degree * "C"), 
        col=c("lightblue", "salmon"), lwd=2, lty=1)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, UK_Tas.decomp$time.series[,2], col="blue", lwd=3)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, UK_PWTas.decomp$time.series[,2], col="red", lwd=3)
legend("bottomright", legend= c("Area Weighted", "Population Weighted"), lwd= 2, col=c("lightblue", "salmon"), bg= ("white"), horiz= F)

# England
Eng_Tas.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Eng_AWTasAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
Eng_PWTas.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Engv_PWTasAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
matplot(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, AW_PW_MonthlyTAnomalies2000_2080[,c(7,22)], type="l", 
        xlab="Year", main= "Monthly Mean Temperature anomaly (°C) in the England 2000-2080", ylab=expression(degree * "C"), 
        col=c("lightblue", "salmon"), lwd=2, lty=1)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Eng_Tas.decomp$time.series[,2], col="blue", lwd=3)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Eng_PWTas.decomp$time.series[,2], col="red", lwd=3)

# NIr
NIr_Tas.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$NIr_AWTasAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
NIr_PWTas.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$NIrv_PWTasAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
matplot(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, AW_PW_MonthlyTAnomalies2000_2080[,c(10,25)], type="l", 
        xlab="Year", main= "Monthly MeanTemperature anomaly (°C) in the N.I. 2000-2080", ylab=expression(degree * "C"), 
        col=c("lightblue", "salmon"), lwd=2, lty=1)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, NIr_Tas.decomp$time.series[,2], col="blue", lwd=3)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, NIr_PWTas.decomp$time.series[,2], col="red", lwd=3)

# Scotland
Sco_Tas.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Sco_AWTasAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
Sco_PWTas.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Scov_PWTasAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
matplot(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, AW_PW_MonthlyTAnomalies2000_2080[,c(13,28)], type="l", 
        xlab="Year", main= "Monthly Mean Temperature anomaly (°C) in the Scotland 2000-2080", ylab=expression(degree * "C"), 
        col=c("lightblue", "salmon"), lwd=2, lty=1)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Sco_Tas.decomp$time.series[,2], col="blue", lwd=3)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Sco_PWTas.decomp$time.series[,2], col="red", lwd=3)

# Wales
Wal_Tas.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Wal_AWTasAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
Wal_PWTas.decomp <- 
  stl(ts(AW_PW_MonthlyTAnomalies2000_2080$Walv_PWTasAnom2000_2080, frequency=12, start = c(1999, 12)), 
      s.window=11, t.window=121)
matplot(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, AW_PW_MonthlyTAnomalies2000_2080[,c(16,31)], type="l", 
        xlab="Year", main= "Monthly Mean Temperature anomaly (°C) in the Wales 2000-2080", ylab=expression(degree * "C"), 
        col=c("lightblue", "salmon"), lwd=2, lty=1)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Wal_Tas.decomp$time.series[,2], col="blue", lwd=3)
lines(AW_PW_MonthlyTAnomalies2000_2080$DecimalDate, Wal_PWTas.decomp$time.series[,2], col="red", lwd=3)



#################################################
# GIF 
################################################


# download data
d <- read.table("http://data.giss.nasa.gov/gistemp/tabledata_v3/GLB.Ts.txt", 
                skip=7, header=TRUE, nrows=143, colClasses = "character")

# clean data
d <- d[-grep("Year", d$Year),]
rownames(d) <- d[,1]
d <- d[,-1]
d <- d[,-13:-19]
d[d=="****"] <- NA

# convert to numeric and celsius
for (i in 1:12) {
  d[,i] <- as.numeric(d[,i]) / 100
}	

# download seasonal adjustments
s <- read.table("http://data.giss.nasa.gov/gistemp/faq/merra2_seas_anom.txt", 
                skip=3, header=TRUE, colClasses = "character")
sa <- as.numeric(s$seas_anom)

# create colours from blue through to red
colours <- colorRampPalette(c("grey","blue","red"))(nrow(d))

# create initial plot
mmin <- -3
mmax <- 3
par(mar=c(2,2,2,0))
plot(1:12, d[1,]+sa, type="l", ylim=c(mmin,mmax), 
     yaxt="n", xaxt="n", bty="n")

# add axes
axis(side=1, at=1:12, 
     labels=c("Jan","Feb","Mar",
              "Apr","May","Jun",
              "Jul","Aug","Sep",
              "Oct","Nov","Dec"), 
     lwd=0, lwd.ticks=1, cex.axis=0.8)

axis(side=2, at=-3:3, labels=-3:3, 
     lwd=0, lwd.ticks=1, las=2, cex.axis=0.8)

# add title
title(main=expression(paste("Temperature ","Anomaly ","(",degree,"C)")), 
      adj=0, cex.main=0.8)
title(main="(Difference from 1980-2015 annual mean)", adj=0, cex.main=0.6, 
      line=0, font.main=1)

# add horizontal lines
for (j in -3:3) {
  lines(1:12, rep(j, 12), lty=3)
}

# add yearly temperature lines
for (j in 1:nrow(d)) {
  lines(1:12, as.numeric(d[j,])+sa, col=colours[j])
}





Assumes you have run the code above to create d and sa
#

doplot <- function(i) {
  
  x   <- d[i,] + sa
  col <- colours[i]
  
  mmin <- -3
  mmax <- 3
  par(mar=c(2.5,2.5,2,0))
  plot(1:12, x, type="l", ylim=c(mmin,mmax), 
       yaxt="n", xaxt="n", bty="n", col=col)
  axis(side=1, at=1:12, labels=c("Jan","Feb","Mar",
                                 "Apr","May","Jun",
                                 "Jul","Aug","Sep",
                                 "Oct","Nov","Dec"), 
       lwd=0, lwd.ticks=1, cex.axis=1.5)
  
  axis(side=2, at=-3:3, labels=-3:3, lwd=0, 
       lwd.ticks=1, las=2, cex.axis=1.5)
  
  title(main=expression(paste("Temperature ","Anomaly ","(",degree,"C)")), 
        adj=0, cex.main=1.5)
  title(main="(Difference from 1980-2015 annual mean)", 
        adj=0, cex.main=1, line=-0.5, font.main=1)
  
  # add horizontal lines
  for (j in -3:3) {
    lines(1:12, rep(j, 12), lty=3)
  }
  
  # add other years
  for (k in 1:i) {
    lines(1:12, d[k,]+sa, col=colours[k])
  }
  
  # add year label
  text(7, 0.8, rownames(d)[i], col=colours[i], font=2, cex=2)
  
}

# create plots/PNGs
for (j in 1:nrow(d)) {
  
  name <- paste("plot",sprintf("%03.f", j),".png", sep="")
  png(name, width=800, height=600)
  doplot(j)
  dev.off()
}

convert *.png -delay 3 -loop 1 globaltemp.gif
convert globaltemp.gif \( +clone -set delay 500 \) +swap +delete globaltempwpause.gif


###############################################################################
# IGNORE THIS SECTION
###############################################################################
#################################################
# Creating averages by seasons
#################################################
# Winter (December-January-February, DJF); 
# Spring (March-April-May, MAM); 
# Summer (June-July-August, JJA);
# Autumn (September-October-November, SON)

#Anom_layers2000_2080 <- seq(as.Date("1999/12/01"), as.Date("2080/11/30"), by = "month")
#
#AnomYears <- as.integer(format(Anom_layers2000_2080, "%Y"))
#AnomMonths <- as.integer(format(Anom_layers2000_2080, "%m"))
#
#n <- length(AnomMonths)
#mnt <- c(AnomMonths[-1], ifelse(AnomMonths[n] < 12, AnomMonths[n]+1, 1)) # move all months back one month  
#yrs <- c(AnomYears[-1],  ifelse(AnomYears[n] < 12, AnomYears[n], AnomYears[n]+1)) # move the years along
#yrs[972] <- 2080
#trims <- (mnt-1) %/% 3  # group by trimesters using integer division (or do: floor((mnt-1) / 3))
#trimnms <- c("DJF", "MAM", "JJA", "SON")[trims + 1] # get names instead of 0, 1, 2, 3
#yt <- paste(yrs, trimnms, sep="_")
#
#SeasTasMin2000_2080 <- stackApply(TasMinAnom2000_2080, indices=yt, fun=mean, na.rm=TRUE) 
#SeasTasMax2000_2080 <- stackApply(TasMaxAnom2000_2080, indices=yt, fun=mean, na.rm=TRUE) 
#SeasTas2000_2080 <- stackApply(TasAnom2000_2080, indices=yt, fun=mean, na.rm=TRUE) 
#
#plot(SeasTasMin2000_2080[[1]], main= "T Anomaly (C) - DJF 2000")
#
#
#write.csv(namesSeasons, "namesSeasons.csv", row.names = F)
#
#writeRaster(SeasTasMin2000_2080 , 
#            filename='SeasTasMin2000_2080.tif', 
#            format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
#
#writeRaster(SeasTasMax2000_2080 , 
#            filename='SeasTasMax2000_2080.tif', 
#            format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
#
#writeRaster(SeasTas2000_2080 , 
#            filename='SeasTas2000_2080.tif', 
#            format="GTiff", overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
#
## Per country: England, Scotland, Wales, North Ireland 
## Simple average per country
## Pop weighted average per country 
#
#write.csv(rasterToPoints(SeasTasMin2000_2080 ), "SeasTasMin2000_2080.csv")
#write.csv(rasterToPoints(SeasTasMax2000_2080 ), "SeasTasMax2000_2080.csv")
#write.csv(rasterToPoints(SeasTas2000_2080 ), "SeasTas2000_2080.csv")
#
## Calculate means by decade and season
## Definition of decades:
## 2000-2009; 2010-2019; 2020-2029; 2030-2039; 2040-2049; 2050-2059; 2060-2069; 2070-2079
#
## I need to delete last year (2080)
#SeasTasMin2000_2079 <- dropLayer(SeasTasMin2000_2080, c(321:324))
#SeasTasMax2000_2079 <- dropLayer(SeasTasMax2000_2080, c(321:324))
#SeasTas2000_2079 <- dropLayer(SeasTas2000_2080, c(321:324))
#
## TasMin
#DecSeasTasMin2000_2009_DJF<- calc(stack(SeasTasMin2000_2079[[1]], SeasTasMin2000_2079[[5]], SeasTasMin2000_2079[[9]], SeasTasMin2000_2079[[13]], SeasTasMin2000_2079[[17]], 
#                                        SeasTasMin2000_2079[[21]], SeasTasMin2000_2079[[25]], SeasTasMin2000_2079[[29]], SeasTasMin2000_2079[[33]], SeasTasMin2000_2079[[37]]), mean)
#DecSeasTasMin2010_2019_DJF<- calc(stack(SeasTasMin2000_2079[[41]], SeasTasMin2000_2079[[45]], SeasTasMin2000_2079[[49]], SeasTasMin2000_2079[[53]], SeasTasMin2000_2079[[57]], 
#                                        SeasTasMin2000_2079[[61]], SeasTasMin2000_2079[[65]], SeasTasMin2000_2079[[69]], SeasTasMin2000_2079[[73]], SeasTasMin2000_2079[[77]]), mean)
#DecSeasTasMin2020_2029_DJF<- calc(stack(SeasTasMin2000_2079[[81]], SeasTasMin2000_2079[[85]], SeasTasMin2000_2079[[89]], SeasTasMin2000_2079[[93]], SeasTasMin2000_2079[[97]], 
#                                        SeasTasMin2000_2079[[101]], SeasTasMin2000_2079[[105]], SeasTasMin2000_2079[[109]], SeasTasMin2000_2079[[113]], SeasTasMin2000_2079[[117]]), mean)
#DecSeasTasMin2030_2039_DJF<- calc(stack(SeasTasMin2000_2079[[121]], SeasTasMin2000_2079[[125]], SeasTasMin2000_2079[[129]], SeasTasMin2000_2079[[133]], SeasTasMin2000_2079[[137]], 
#                                        SeasTasMin2000_2079[[141]], SeasTasMin2000_2079[[145]], SeasTasMin2000_2079[[149]], SeasTasMin2000_2079[[153]], SeasTasMin2000_2079[[157]]), mean)
#DecSeasTasMin2040_2049_DJF<- calc(stack(SeasTasMin2000_2079[[161]], SeasTasMin2000_2079[[165]], SeasTasMin2000_2079[[169]], SeasTasMin2000_2079[[173]], SeasTasMin2000_2079[[177]], 
#                                        SeasTasMin2000_2079[[181]], SeasTasMin2000_2079[[185]], SeasTasMin2000_2079[[189]], SeasTasMin2000_2079[[193]], SeasTasMin2000_2079[[197]]), mean)
#DecSeasTasMin2050_2059_DJF<- calc(stack(SeasTasMin2000_2079[[201]], SeasTasMin2000_2079[[205]], SeasTasMin2000_2079[[209]], SeasTasMin2000_2079[[213]], SeasTasMin2000_2079[[217]], 
#                                        SeasTasMin2000_2079[[221]], SeasTasMin2000_2079[[225]], SeasTasMin2000_2079[[229]], SeasTasMin2000_2079[[233]], SeasTasMin2000_2079[[237]]), mean)
#DecSeasTasMin2060_2069_DJF<- calc(stack(SeasTasMin2000_2079[[241]], SeasTasMin2000_2079[[245]], SeasTasMin2000_2079[[249]], SeasTasMin2000_2079[[253]], SeasTasMin2000_2079[[257]], 
#                                        SeasTasMin2000_2079[[261]], SeasTasMin2000_2079[[265]], SeasTasMin2000_2079[[269]], SeasTasMin2000_2079[[273]], SeasTasMin2000_2079[[277]]), mean)
#DecSeasTasMin2070_2079_DJF<- calc(stack(SeasTasMin2000_2079[[281]], SeasTasMin2000_2079[[285]], SeasTasMin2000_2079[[289]], SeasTasMin2000_2079[[293]], SeasTasMin2000_2079[[297]], 
#                                        SeasTasMin2000_2079[[301]], SeasTasMin2000_2079[[305]], SeasTasMin2000_2079[[309]], SeasTasMin2000_2079[[313]], SeasTasMin2000_2079[[317]]), mean)
#
#DecSeasTasMin2000_2009_MAM<- calc(stack(SeasTasMin2000_2079[[2]], SeasTasMin2000_2079[[6]], SeasTasMin2000_2079[[10]], SeasTasMin2000_2079[[14]], SeasTasMin2000_2079[[18]], 
#                                        SeasTasMin2000_2079[[22]], SeasTasMin2000_2079[[26]], SeasTasMin2000_2079[[30]], SeasTasMin2000_2079[[34]], SeasTasMin2000_2079[[38]]), mean)
#DecSeasTasMin2010_2019_MAM<- calc(stack(SeasTasMin2000_2079[[42]], SeasTasMin2000_2079[[46]], SeasTasMin2000_2079[[50]], SeasTasMin2000_2079[[54]], SeasTasMin2000_2079[[58]], 
#                                        SeasTasMin2000_2079[[62]], SeasTasMin2000_2079[[66]], SeasTasMin2000_2079[[70]], SeasTasMin2000_2079[[74]], SeasTasMin2000_2079[[78]]), mean)
#DecSeasTasMin2020_2029_MAM<- calc(stack(SeasTasMin2000_2079[[82]], SeasTasMin2000_2079[[86]], SeasTasMin2000_2079[[90]], SeasTasMin2000_2079[[94]], SeasTasMin2000_2079[[98]], 
#                                        SeasTasMin2000_2079[[102]], SeasTasMin2000_2079[[106]], SeasTasMin2000_2079[[110]], SeasTasMin2000_2079[[114]], SeasTasMin2000_2079[[118]]), mean)
#DecSeasTasMin2030_2039_MAM<- calc(stack(SeasTasMin2000_2079[[122]], SeasTasMin2000_2079[[126]], SeasTasMin2000_2079[[130]], SeasTasMin2000_2079[[134]], SeasTasMin2000_2079[[138]], 
#                                        SeasTasMin2000_2079[[142]], SeasTasMin2000_2079[[146]], SeasTasMin2000_2079[[150]], SeasTasMin2000_2079[[154]], SeasTasMin2000_2079[[158]]), mean)
#DecSeasTasMin2040_2049_MAM<- calc(stack(SeasTasMin2000_2079[[162]], SeasTasMin2000_2079[[166]], SeasTasMin2000_2079[[170]], SeasTasMin2000_2079[[174]], SeasTasMin2000_2079[[178]], 
#                                        SeasTasMin2000_2079[[182]], SeasTasMin2000_2079[[186]], SeasTasMin2000_2079[[190]], SeasTasMin2000_2079[[194]], SeasTasMin2000_2079[[198]]), mean)
#DecSeasTasMin2050_2059_MAM<- calc(stack(SeasTasMin2000_2079[[202]], SeasTasMin2000_2079[[206]], SeasTasMin2000_2079[[210]], SeasTasMin2000_2079[[214]], SeasTasMin2000_2079[[218]], 
#                                        SeasTasMin2000_2079[[222]], SeasTasMin2000_2079[[226]], SeasTasMin2000_2079[[230]], SeasTasMin2000_2079[[234]], SeasTasMin2000_2079[[238]]), mean)
#DecSeasTasMin2060_2069_MAM<- calc(stack(SeasTasMin2000_2079[[242]], SeasTasMin2000_2079[[246]], SeasTasMin2000_2079[[250]], SeasTasMin2000_2079[[254]], SeasTasMin2000_2079[[258]], 
#                                        SeasTasMin2000_2079[[262]], SeasTasMin2000_2079[[266]], SeasTasMin2000_2079[[270]], SeasTasMin2000_2079[[274]], SeasTasMin2000_2079[[278]]), mean)
#DecSeasTasMin2070_2079_MAM<- calc(stack(SeasTasMin2000_2079[[282]], SeasTasMin2000_2079[[286]], SeasTasMin2000_2079[[290]], SeasTasMin2000_2079[[294]], SeasTasMin2000_2079[[298]], 
#                                        SeasTasMin2000_2079[[302]], SeasTasMin2000_2079[[306]], SeasTasMin2000_2079[[310]], SeasTasMin2000_2079[[314]], SeasTasMin2000_2079[[318]]), mean)
#
#DecSeasTasMin2000_2009_JJA<- calc(stack(SeasTasMin2000_2079[[3]], SeasTasMin2000_2079[[7]], SeasTasMin2000_2079[[11]], SeasTasMin2000_2079[[15]], SeasTasMin2000_2079[[19]], 
#                                        SeasTasMin2000_2079[[23]], SeasTasMin2000_2079[[27]], SeasTasMin2000_2079[[31]], SeasTasMin2000_2079[[35]], SeasTasMin2000_2079[[39]]), mean)
#DecSeasTasMin2010_2019_JJA<- calc(stack(SeasTasMin2000_2079[[43]], SeasTasMin2000_2079[[47]], SeasTasMin2000_2079[[51]], SeasTasMin2000_2079[[55]], SeasTasMin2000_2079[[59]], 
#                                        SeasTasMin2000_2079[[63]], SeasTasMin2000_2079[[67]], SeasTasMin2000_2079[[71]], SeasTasMin2000_2079[[75]], SeasTasMin2000_2079[[79]]), mean)
#DecSeasTasMin2020_2029_JJA<- calc(stack(SeasTasMin2000_2079[[83]], SeasTasMin2000_2079[[87]], SeasTasMin2000_2079[[91]], SeasTasMin2000_2079[[95]], SeasTasMin2000_2079[[99]], 
#                                        SeasTasMin2000_2079[[103]], SeasTasMin2000_2079[[107]], SeasTasMin2000_2079[[111]], SeasTasMin2000_2079[[115]], SeasTasMin2000_2079[[119]]), mean)
#DecSeasTasMin2030_2039_JJA<- calc(stack(SeasTasMin2000_2079[[123]], SeasTasMin2000_2079[[127]], SeasTasMin2000_2079[[131]], SeasTasMin2000_2079[[135]], SeasTasMin2000_2079[[139]], 
#                                        SeasTasMin2000_2079[[143]], SeasTasMin2000_2079[[147]], SeasTasMin2000_2079[[151]], SeasTasMin2000_2079[[155]], SeasTasMin2000_2079[[159]]), mean)
#DecSeasTasMin2040_2049_JJA<- calc(stack(SeasTasMin2000_2079[[163]], SeasTasMin2000_2079[[167]], SeasTasMin2000_2079[[171]], SeasTasMin2000_2079[[175]], SeasTasMin2000_2079[[179]], 
#                                        SeasTasMin2000_2079[[183]], SeasTasMin2000_2079[[187]], SeasTasMin2000_2079[[191]], SeasTasMin2000_2079[[195]], SeasTasMin2000_2079[[199]]), mean)
#DecSeasTasMin2050_2059_JJA<- calc(stack(SeasTasMin2000_2079[[203]], SeasTasMin2000_2079[[207]], SeasTasMin2000_2079[[211]], SeasTasMin2000_2079[[215]], SeasTasMin2000_2079[[219]], 
#                                        SeasTasMin2000_2079[[223]], SeasTasMin2000_2079[[227]], SeasTasMin2000_2079[[231]], SeasTasMin2000_2079[[235]], SeasTasMin2000_2079[[239]]), mean)
#DecSeasTasMin2060_2069_JJA<- calc(stack(SeasTasMin2000_2079[[243]], SeasTasMin2000_2079[[247]], SeasTasMin2000_2079[[251]], SeasTasMin2000_2079[[255]], SeasTasMin2000_2079[[259]], 
#                                        SeasTasMin2000_2079[[263]], SeasTasMin2000_2079[[267]], SeasTasMin2000_2079[[271]], SeasTasMin2000_2079[[275]], SeasTasMin2000_2079[[279]]), mean)
#DecSeasTasMin2070_2079_JJA<- calc(stack(SeasTasMin2000_2079[[283]], SeasTasMin2000_2079[[287]], SeasTasMin2000_2079[[291]], SeasTasMin2000_2079[[295]], SeasTasMin2000_2079[[299]], 
#                                        SeasTasMin2000_2079[[303]], SeasTasMin2000_2079[[307]], SeasTasMin2000_2079[[311]], SeasTasMin2000_2079[[315]], SeasTasMin2000_2079[[319]]), mean)
#
#DecSeasTasMin2000_2009_SON<- calc(stack(SeasTasMin2000_2079[[4]], SeasTasMin2000_2079[[8]], SeasTasMin2000_2079[[12]], SeasTasMin2000_2079[[16]], SeasTasMin2000_2079[[20]], 
#                                        SeasTasMin2000_2079[[24]], SeasTasMin2000_2079[[28]], SeasTasMin2000_2079[[32]], SeasTasMin2000_2079[[36]], SeasTasMin2000_2079[[40]]), mean)
#DecSeasTasMin2010_2019_SON<- calc(stack(SeasTasMin2000_2079[[44]], SeasTasMin2000_2079[[48]], SeasTasMin2000_2079[[52]], SeasTasMin2000_2079[[56]], SeasTasMin2000_2079[[60]], 
#                                        SeasTasMin2000_2079[[64]], SeasTasMin2000_2079[[68]], SeasTasMin2000_2079[[72]], SeasTasMin2000_2079[[76]], SeasTasMin2000_2079[[80]]), mean)
#DecSeasTasMin2020_2029_SON<- calc(stack(SeasTasMin2000_2079[[84]], SeasTasMin2000_2079[[88]], SeasTasMin2000_2079[[92]], SeasTasMin2000_2079[[96]], SeasTasMin2000_2079[[100]], 
#                                        SeasTasMin2000_2079[[104]], SeasTasMin2000_2079[[108]], SeasTasMin2000_2079[[112]], SeasTasMin2000_2079[[116]], SeasTasMin2000_2079[[120]]), mean)
#DecSeasTasMin2030_2039_SON<- calc(stack(SeasTasMin2000_2079[[124]], SeasTasMin2000_2079[[128]], SeasTasMin2000_2079[[132]], SeasTasMin2000_2079[[136]], SeasTasMin2000_2079[[140]], 
#                                        SeasTasMin2000_2079[[144]], SeasTasMin2000_2079[[148]], SeasTasMin2000_2079[[152]], SeasTasMin2000_2079[[156]], SeasTasMin2000_2079[[160]]), mean)
#DecSeasTasMin2040_2049_SON<- calc(stack(SeasTasMin2000_2079[[164]], SeasTasMin2000_2079[[168]], SeasTasMin2000_2079[[172]], SeasTasMin2000_2079[[176]], SeasTasMin2000_2079[[180]], 
#                                        SeasTasMin2000_2079[[184]], SeasTasMin2000_2079[[188]], SeasTasMin2000_2079[[192]], SeasTasMin2000_2079[[196]], SeasTasMin2000_2079[[200]]), mean)
#DecSeasTasMin2050_2059_SON<- calc(stack(SeasTasMin2000_2079[[204]], SeasTasMin2000_2079[[208]], SeasTasMin2000_2079[[212]], SeasTasMin2000_2079[[216]], SeasTasMin2000_2079[[220]], 
#                                        SeasTasMin2000_2079[[224]], SeasTasMin2000_2079[[228]], SeasTasMin2000_2079[[232]], SeasTasMin2000_2079[[236]], SeasTasMin2000_2079[[240]]), mean)
#DecSeasTasMin2060_2069_SON<- calc(stack(SeasTasMin2000_2079[[244]], SeasTasMin2000_2079[[248]], SeasTasMin2000_2079[[252]], SeasTasMin2000_2079[[256]], SeasTasMin2000_2079[[260]], 
#                                        SeasTasMin2000_2079[[264]], SeasTasMin2000_2079[[268]], SeasTasMin2000_2079[[272]], SeasTasMin2000_2079[[276]], SeasTasMin2000_2079[[280]]), mean)
#DecSeasTasMin2070_2079_SON<- calc(stack(SeasTasMin2000_2079[[284]], SeasTasMin2000_2079[[288]], SeasTasMin2000_2079[[292]], SeasTasMin2000_2079[[296]], SeasTasMin2000_2079[[300]], 
#                                        SeasTasMin2000_2079[[304]], SeasTasMin2000_2079[[308]], SeasTasMin2000_2079[[312]], SeasTasMin2000_2079[[316]], SeasTasMin2000_2079[[320]]), mean)
#
#DecSeasTasMin2000_2079 <- stack(DecSeasTasMin2000_2009_DJF, DecSeasTasMin2000_2009_MAM, DecSeasTasMin2000_2009_JJA, DecSeasTasMin2000_2009_SON,
#                                DecSeasTasMin2010_2019_DJF, DecSeasTasMin2010_2019_MAM, DecSeasTasMin2010_2019_JJA, DecSeasTasMin2010_2019_SON,
#                                DecSeasTasMin2020_2029_DJF, DecSeasTasMin2020_2029_MAM, DecSeasTasMin2020_2029_JJA, DecSeasTasMin2020_2029_SON,
#                                DecSeasTasMin2030_2039_DJF, DecSeasTasMin2030_2039_MAM, DecSeasTasMin2030_2039_JJA, DecSeasTasMin2030_2039_SON,
#                                DecSeasTasMin2040_2049_DJF, DecSeasTasMin2040_2049_MAM, DecSeasTasMin2040_2049_JJA, DecSeasTasMin2040_2049_SON,
#                                DecSeasTasMin2050_2059_DJF, DecSeasTasMin2050_2059_MAM, DecSeasTasMin2050_2059_JJA, DecSeasTasMin2050_2059_SON,
#                                DecSeasTasMin2060_2069_DJF, DecSeasTasMin2060_2069_MAM, DecSeasTasMin2060_2069_JJA, DecSeasTasMin2060_2069_SON,
#                                DecSeasTasMin2070_2079_DJF, DecSeasTasMin2070_2079_MAM, DecSeasTasMin2070_2079_JJA, DecSeasTasMin2070_2079_SON)
#
#
## TasMax
#DecSeasTasMax2000_2009_DJF<- calc(stack(SeasTasMax2000_2079[[1]], SeasTasMax2000_2079[[5]], SeasTasMax2000_2079[[9]], SeasTasMax2000_2079[[13]], SeasTasMax2000_2079[[17]], 
#                                        SeasTasMax2000_2079[[21]], SeasTasMax2000_2079[[25]], SeasTasMax2000_2079[[29]], SeasTasMax2000_2079[[33]], SeasTasMax2000_2079[[37]]), mean)
#DecSeasTasMax2010_2019_DJF<- calc(stack(SeasTasMax2000_2079[[41]], SeasTasMax2000_2079[[45]], SeasTasMax2000_2079[[49]], SeasTasMax2000_2079[[53]], SeasTasMax2000_2079[[57]], 
#                                        SeasTasMax2000_2079[[61]], SeasTasMax2000_2079[[65]], SeasTasMax2000_2079[[69]], SeasTasMax2000_2079[[73]], SeasTasMax2000_2079[[77]]), mean)
#DecSeasTasMax2020_2029_DJF<- calc(stack(SeasTasMax2000_2079[[81]], SeasTasMax2000_2079[[85]], SeasTasMax2000_2079[[89]], SeasTasMax2000_2079[[93]], SeasTasMax2000_2079[[97]], 
#                                        SeasTasMax2000_2079[[101]], SeasTasMax2000_2079[[105]], SeasTasMax2000_2079[[109]], SeasTasMax2000_2079[[113]], SeasTasMax2000_2079[[117]]), mean)
#DecSeasTasMax2030_2039_DJF<- calc(stack(SeasTasMax2000_2079[[121]], SeasTasMax2000_2079[[125]], SeasTasMax2000_2079[[129]], SeasTasMax2000_2079[[133]], SeasTasMax2000_2079[[137]], 
#                                        SeasTasMax2000_2079[[141]], SeasTasMax2000_2079[[145]], SeasTasMax2000_2079[[149]], SeasTasMax2000_2079[[153]], SeasTasMax2000_2079[[157]]), mean)
#DecSeasTasMax2040_2049_DJF<- calc(stack(SeasTasMax2000_2079[[161]], SeasTasMax2000_2079[[165]], SeasTasMax2000_2079[[169]], SeasTasMax2000_2079[[173]], SeasTasMax2000_2079[[177]], 
#                                        SeasTasMax2000_2079[[181]], SeasTasMax2000_2079[[185]], SeasTasMax2000_2079[[189]], SeasTasMax2000_2079[[193]], SeasTasMax2000_2079[[197]]), mean)
#DecSeasTasMax2050_2059_DJF<- calc(stack(SeasTasMax2000_2079[[201]], SeasTasMax2000_2079[[205]], SeasTasMax2000_2079[[209]], SeasTasMax2000_2079[[213]], SeasTasMax2000_2079[[217]], 
#                                        SeasTasMax2000_2079[[221]], SeasTasMax2000_2079[[225]], SeasTasMax2000_2079[[229]], SeasTasMax2000_2079[[233]], SeasTasMax2000_2079[[237]]), mean)
#DecSeasTasMax2060_2069_DJF<- calc(stack(SeasTasMax2000_2079[[241]], SeasTasMax2000_2079[[245]], SeasTasMax2000_2079[[249]], SeasTasMax2000_2079[[253]], SeasTasMax2000_2079[[257]], 
#                                        SeasTasMax2000_2079[[261]], SeasTasMax2000_2079[[265]], SeasTasMax2000_2079[[269]], SeasTasMax2000_2079[[273]], SeasTasMax2000_2079[[277]]), mean)
#DecSeasTasMax2070_2079_DJF<- calc(stack(SeasTasMax2000_2079[[281]], SeasTasMax2000_2079[[285]], SeasTasMax2000_2079[[289]], SeasTasMax2000_2079[[293]], SeasTasMax2000_2079[[297]], 
#                                        SeasTasMax2000_2079[[301]], SeasTasMax2000_2079[[305]], SeasTasMax2000_2079[[309]], SeasTasMax2000_2079[[313]], SeasTasMax2000_2079[[317]]), mean)
#
#DecSeasTasMax2000_2009_MAM<- calc(stack(SeasTasMax2000_2079[[2]], SeasTasMax2000_2079[[6]], SeasTasMax2000_2079[[10]], SeasTasMax2000_2079[[14]], SeasTasMax2000_2079[[18]], 
#                                        SeasTasMax2000_2079[[22]], SeasTasMax2000_2079[[26]], SeasTasMax2000_2079[[30]], SeasTasMax2000_2079[[34]], SeasTasMax2000_2079[[38]]), mean)
#DecSeasTasMax2010_2019_MAM<- calc(stack(SeasTasMax2000_2079[[42]], SeasTasMax2000_2079[[46]], SeasTasMax2000_2079[[50]], SeasTasMax2000_2079[[54]], SeasTasMax2000_2079[[58]], 
#                                        SeasTasMax2000_2079[[62]], SeasTasMax2000_2079[[66]], SeasTasMax2000_2079[[70]], SeasTasMax2000_2079[[74]], SeasTasMax2000_2079[[78]]), mean)
#DecSeasTasMax2020_2029_MAM<- calc(stack(SeasTasMax2000_2079[[82]], SeasTasMax2000_2079[[86]], SeasTasMax2000_2079[[90]], SeasTasMax2000_2079[[94]], SeasTasMax2000_2079[[98]], 
#                                        SeasTasMax2000_2079[[102]], SeasTasMax2000_2079[[106]], SeasTasMax2000_2079[[110]], SeasTasMax2000_2079[[114]], SeasTasMax2000_2079[[118]]), mean)
#DecSeasTasMax2030_2039_MAM<- calc(stack(SeasTasMax2000_2079[[122]], SeasTasMax2000_2079[[126]], SeasTasMax2000_2079[[130]], SeasTasMax2000_2079[[134]], SeasTasMax2000_2079[[138]], 
#                                        SeasTasMax2000_2079[[142]], SeasTasMax2000_2079[[146]], SeasTasMax2000_2079[[150]], SeasTasMax2000_2079[[154]], SeasTasMax2000_2079[[158]]), mean)
#DecSeasTasMax2040_2049_MAM<- calc(stack(SeasTasMax2000_2079[[162]], SeasTasMax2000_2079[[166]], SeasTasMax2000_2079[[170]], SeasTasMax2000_2079[[174]], SeasTasMax2000_2079[[178]], 
#                                        SeasTasMax2000_2079[[182]], SeasTasMax2000_2079[[186]], SeasTasMax2000_2079[[190]], SeasTasMax2000_2079[[194]], SeasTasMax2000_2079[[198]]), mean)
#DecSeasTasMax2050_2059_MAM<- calc(stack(SeasTasMax2000_2079[[202]], SeasTasMax2000_2079[[206]], SeasTasMax2000_2079[[210]], SeasTasMax2000_2079[[214]], SeasTasMax2000_2079[[218]], 
#                                        SeasTasMax2000_2079[[222]], SeasTasMax2000_2079[[226]], SeasTasMax2000_2079[[230]], SeasTasMax2000_2079[[234]], SeasTasMax2000_2079[[238]]), mean)
#DecSeasTasMax2060_2069_MAM<- calc(stack(SeasTasMax2000_2079[[242]], SeasTasMax2000_2079[[246]], SeasTasMax2000_2079[[250]], SeasTasMax2000_2079[[254]], SeasTasMax2000_2079[[258]], 
#                                        SeasTasMax2000_2079[[262]], SeasTasMax2000_2079[[266]], SeasTasMax2000_2079[[270]], SeasTasMax2000_2079[[274]], SeasTasMax2000_2079[[278]]), mean)
#DecSeasTasMax2070_2079_MAM<- calc(stack(SeasTasMax2000_2079[[282]], SeasTasMax2000_2079[[286]], SeasTasMax2000_2079[[290]], SeasTasMax2000_2079[[294]], SeasTasMax2000_2079[[298]], 
#                                        SeasTasMax2000_2079[[302]], SeasTasMax2000_2079[[306]], SeasTasMax2000_2079[[310]], SeasTasMax2000_2079[[314]], SeasTasMax2000_2079[[318]]), mean)
#
#DecSeasTasMax2000_2009_JJA<- calc(stack(SeasTasMax2000_2079[[3]], SeasTasMax2000_2079[[7]], SeasTasMax2000_2079[[11]], SeasTasMax2000_2079[[15]], SeasTasMax2000_2079[[19]], 
#                                        SeasTasMax2000_2079[[23]], SeasTasMax2000_2079[[27]], SeasTasMax2000_2079[[31]], SeasTasMax2000_2079[[35]], SeasTasMax2000_2079[[39]]), mean)
#DecSeasTasMax2010_2019_JJA<- calc(stack(SeasTasMax2000_2079[[43]], SeasTasMax2000_2079[[47]], SeasTasMax2000_2079[[51]], SeasTasMax2000_2079[[55]], SeasTasMax2000_2079[[59]], 
#                                        SeasTasMax2000_2079[[63]], SeasTasMax2000_2079[[67]], SeasTasMax2000_2079[[71]], SeasTasMax2000_2079[[75]], SeasTasMax2000_2079[[79]]), mean)
#DecSeasTasMax2020_2029_JJA<- calc(stack(SeasTasMax2000_2079[[83]], SeasTasMax2000_2079[[87]], SeasTasMax2000_2079[[91]], SeasTasMax2000_2079[[95]], SeasTasMax2000_2079[[99]], 
#                                        SeasTasMax2000_2079[[103]], SeasTasMax2000_2079[[107]], SeasTasMax2000_2079[[111]], SeasTasMax2000_2079[[115]], SeasTasMax2000_2079[[119]]), mean)
#DecSeasTasMax2030_2039_JJA<- calc(stack(SeasTasMax2000_2079[[123]], SeasTasMax2000_2079[[127]], SeasTasMax2000_2079[[131]], SeasTasMax2000_2079[[135]], SeasTasMax2000_2079[[139]], 
#                                        SeasTasMax2000_2079[[143]], SeasTasMax2000_2079[[147]], SeasTasMax2000_2079[[151]], SeasTasMax2000_2079[[155]], SeasTasMax2000_2079[[159]]), mean)
#DecSeasTasMax2040_2049_JJA<- calc(stack(SeasTasMax2000_2079[[163]], SeasTasMax2000_2079[[167]], SeasTasMax2000_2079[[171]], SeasTasMax2000_2079[[175]], SeasTasMax2000_2079[[179]], 
#                                        SeasTasMax2000_2079[[183]], SeasTasMax2000_2079[[187]], SeasTasMax2000_2079[[191]], SeasTasMax2000_2079[[195]], SeasTasMax2000_2079[[199]]), mean)
#DecSeasTasMax2050_2059_JJA<- calc(stack(SeasTasMax2000_2079[[203]], SeasTasMax2000_2079[[207]], SeasTasMax2000_2079[[211]], SeasTasMax2000_2079[[215]], SeasTasMax2000_2079[[219]], 
#                                        SeasTasMax2000_2079[[223]], SeasTasMax2000_2079[[227]], SeasTasMax2000_2079[[231]], SeasTasMax2000_2079[[235]], SeasTasMax2000_2079[[239]]), mean)
#DecSeasTasMax2060_2069_JJA<- calc(stack(SeasTasMax2000_2079[[243]], SeasTasMax2000_2079[[247]], SeasTasMax2000_2079[[251]], SeasTasMax2000_2079[[255]], SeasTasMax2000_2079[[259]], 
#                                        SeasTasMax2000_2079[[263]], SeasTasMax2000_2079[[267]], SeasTasMax2000_2079[[271]], SeasTasMax2000_2079[[275]], SeasTasMax2000_2079[[279]]), mean)
#DecSeasTasMax2070_2079_JJA<- calc(stack(SeasTasMax2000_2079[[283]], SeasTasMax2000_2079[[287]], SeasTasMax2000_2079[[291]], SeasTasMax2000_2079[[295]], SeasTasMax2000_2079[[299]], 
#                                        SeasTasMax2000_2079[[303]], SeasTasMax2000_2079[[307]], SeasTasMax2000_2079[[311]], SeasTasMax2000_2079[[315]], SeasTasMax2000_2079[[319]]), mean)
#
#DecSeasTasMax2000_2009_SON<- calc(stack(SeasTasMax2000_2079[[4]], SeasTasMax2000_2079[[8]], SeasTasMax2000_2079[[12]], SeasTasMax2000_2079[[16]], SeasTasMax2000_2079[[20]], 
#                                        SeasTasMax2000_2079[[24]], SeasTasMax2000_2079[[28]], SeasTasMax2000_2079[[32]], SeasTasMax2000_2079[[36]], SeasTasMax2000_2079[[40]]), mean)
#DecSeasTasMax2010_2019_SON<- calc(stack(SeasTasMax2000_2079[[44]], SeasTasMax2000_2079[[48]], SeasTasMax2000_2079[[52]], SeasTasMax2000_2079[[56]], SeasTasMax2000_2079[[60]], 
#                                        SeasTasMax2000_2079[[64]], SeasTasMax2000_2079[[68]], SeasTasMax2000_2079[[72]], SeasTasMax2000_2079[[76]], SeasTasMax2000_2079[[80]]), mean)
#DecSeasTasMax2020_2029_SON<- calc(stack(SeasTasMax2000_2079[[84]], SeasTasMax2000_2079[[88]], SeasTasMax2000_2079[[92]], SeasTasMax2000_2079[[96]], SeasTasMax2000_2079[[100]], 
#                                        SeasTasMax2000_2079[[104]], SeasTasMax2000_2079[[108]], SeasTasMax2000_2079[[112]], SeasTasMax2000_2079[[116]], SeasTasMax2000_2079[[120]]), mean)
#DecSeasTasMax2030_2039_SON<- calc(stack(SeasTasMax2000_2079[[124]], SeasTasMax2000_2079[[128]], SeasTasMax2000_2079[[132]], SeasTasMax2000_2079[[136]], SeasTasMax2000_2079[[140]], 
#                                        SeasTasMax2000_2079[[144]], SeasTasMax2000_2079[[148]], SeasTasMax2000_2079[[152]], SeasTasMax2000_2079[[156]], SeasTasMax2000_2079[[160]]), mean)
#DecSeasTasMax2040_2049_SON<- calc(stack(SeasTasMax2000_2079[[164]], SeasTasMax2000_2079[[168]], SeasTasMax2000_2079[[172]], SeasTasMax2000_2079[[176]], SeasTasMax2000_2079[[180]], 
#                                        SeasTasMax2000_2079[[184]], SeasTasMax2000_2079[[188]], SeasTasMax2000_2079[[192]], SeasTasMax2000_2079[[196]], SeasTasMax2000_2079[[200]]), mean)
#DecSeasTasMax2050_2059_SON<- calc(stack(SeasTasMax2000_2079[[204]], SeasTasMax2000_2079[[208]], SeasTasMax2000_2079[[212]], SeasTasMax2000_2079[[216]], SeasTasMax2000_2079[[220]], 
#                                        SeasTasMax2000_2079[[224]], SeasTasMax2000_2079[[228]], SeasTasMax2000_2079[[232]], SeasTasMax2000_2079[[236]], SeasTasMax2000_2079[[240]]), mean)
#DecSeasTasMax2060_2069_SON<- calc(stack(SeasTasMax2000_2079[[244]], SeasTasMax2000_2079[[248]], SeasTasMax2000_2079[[252]], SeasTasMax2000_2079[[256]], SeasTasMax2000_2079[[260]], 
#                                        SeasTasMax2000_2079[[264]], SeasTasMax2000_2079[[268]], SeasTasMax2000_2079[[272]], SeasTasMax2000_2079[[276]], SeasTasMax2000_2079[[280]]), mean)
#DecSeasTasMax2070_2079_SON<- calc(stack(SeasTasMax2000_2079[[284]], SeasTasMax2000_2079[[288]], SeasTasMax2000_2079[[292]], SeasTasMax2000_2079[[296]], SeasTasMax2000_2079[[300]], 
#                                        SeasTasMax2000_2079[[304]], SeasTasMax2000_2079[[308]], SeasTasMax2000_2079[[312]], SeasTasMax2000_2079[[316]], SeasTasMax2000_2079[[320]]), mean)
#
#DecSeasTasMax2000_2079 <- stack(DecSeasTasMax2000_2009_DJF, DecSeasTasMax2000_2009_MAM, DecSeasTasMax2000_2009_JJA, DecSeasTasMax2000_2009_SON,
#                                DecSeasTasMax2010_2019_DJF, DecSeasTasMax2010_2019_MAM, DecSeasTasMax2010_2019_JJA, DecSeasTasMax2010_2019_SON,
#                                DecSeasTasMax2020_2029_DJF, DecSeasTasMax2020_2029_MAM, DecSeasTasMax2020_2029_JJA, DecSeasTasMax2020_2029_SON,
#                                DecSeasTasMax2030_2039_DJF, DecSeasTasMax2030_2039_MAM, DecSeasTasMax2030_2039_JJA, DecSeasTasMax2030_2039_SON,
#                                DecSeasTasMax2040_2049_DJF, DecSeasTasMax2040_2049_MAM, DecSeasTasMax2040_2049_JJA, DecSeasTasMax2040_2049_SON,
#                                DecSeasTasMax2050_2059_DJF, DecSeasTasMax2050_2059_MAM, DecSeasTasMax2050_2059_JJA, DecSeasTasMax2050_2059_SON,
#                                DecSeasTasMax2060_2069_DJF, DecSeasTasMax2060_2069_MAM, DecSeasTasMax2060_2069_JJA, DecSeasTasMax2060_2069_SON,
#                                DecSeasTasMax2070_2079_DJF, DecSeasTasMax2070_2079_MAM, DecSeasTasMax2070_2079_JJA, DecSeasTasMax2070_2079_SON)
#
## Tas
#DecSeasTas2000_2009_DJF<- calc(stack(SeasTas2000_2079[[1]], SeasTas2000_2079[[5]], SeasTas2000_2079[[9]], SeasTas2000_2079[[13]], SeasTas2000_2079[[17]], 
#                                     SeasTas2000_2079[[21]], SeasTas2000_2079[[25]], SeasTas2000_2079[[29]], SeasTas2000_2079[[33]], SeasTas2000_2079[[37]]), mean)
#DecSeasTas2010_2019_DJF<- calc(stack(SeasTas2000_2079[[41]], SeasTas2000_2079[[45]], SeasTas2000_2079[[49]], SeasTas2000_2079[[53]], SeasTas2000_2079[[57]], 
#                                     SeasTas2000_2079[[61]], SeasTas2000_2079[[65]], SeasTas2000_2079[[69]], SeasTas2000_2079[[73]], SeasTas2000_2079[[77]]), mean)
#DecSeasTas2020_2029_DJF<- calc(stack(SeasTas2000_2079[[81]], SeasTas2000_2079[[85]], SeasTas2000_2079[[89]], SeasTas2000_2079[[93]], SeasTas2000_2079[[97]], 
#                                     SeasTas2000_2079[[101]], SeasTas2000_2079[[105]], SeasTas2000_2079[[109]], SeasTas2000_2079[[113]], SeasTas2000_2079[[117]]), mean)
#DecSeasTas2030_2039_DJF<- calc(stack(SeasTas2000_2079[[121]], SeasTas2000_2079[[125]], SeasTas2000_2079[[129]], SeasTas2000_2079[[133]], SeasTas2000_2079[[137]], 
#                                     SeasTas2000_2079[[141]], SeasTas2000_2079[[145]], SeasTas2000_2079[[149]], SeasTas2000_2079[[153]], SeasTas2000_2079[[157]]), mean)
#DecSeasTas2040_2049_DJF<- calc(stack(SeasTas2000_2079[[161]], SeasTas2000_2079[[165]], SeasTas2000_2079[[169]], SeasTas2000_2079[[173]], SeasTas2000_2079[[177]], 
#                                     SeasTas2000_2079[[181]], SeasTas2000_2079[[185]], SeasTas2000_2079[[189]], SeasTas2000_2079[[193]], SeasTas2000_2079[[197]]), mean)
#DecSeasTas2050_2059_DJF<- calc(stack(SeasTas2000_2079[[201]], SeasTas2000_2079[[205]], SeasTas2000_2079[[209]], SeasTas2000_2079[[213]], SeasTas2000_2079[[217]], 
#                                     SeasTas2000_2079[[221]], SeasTas2000_2079[[225]], SeasTas2000_2079[[229]], SeasTas2000_2079[[233]], SeasTas2000_2079[[237]]), mean)
#DecSeasTas2060_2069_DJF<- calc(stack(SeasTas2000_2079[[241]], SeasTas2000_2079[[245]], SeasTas2000_2079[[249]], SeasTas2000_2079[[253]], SeasTas2000_2079[[257]], 
#                                     SeasTas2000_2079[[261]], SeasTas2000_2079[[265]], SeasTas2000_2079[[269]], SeasTas2000_2079[[273]], SeasTas2000_2079[[277]]), mean)
#DecSeasTas2070_2079_DJF<- calc(stack(SeasTas2000_2079[[281]], SeasTas2000_2079[[285]], SeasTas2000_2079[[289]], SeasTas2000_2079[[293]], SeasTas2000_2079[[297]], 
#                                     SeasTas2000_2079[[301]], SeasTas2000_2079[[305]], SeasTas2000_2079[[309]], SeasTas2000_2079[[313]], SeasTas2000_2079[[317]]), mean)
#
#DecSeasTas2000_2009_MAM<- calc(stack(SeasTas2000_2079[[2]], SeasTas2000_2079[[6]], SeasTas2000_2079[[10]], SeasTas2000_2079[[14]], SeasTas2000_2079[[18]], 
#                                     SeasTas2000_2079[[22]], SeasTas2000_2079[[26]], SeasTas2000_2079[[30]], SeasTas2000_2079[[34]], SeasTas2000_2079[[38]]), mean)
#DecSeasTas2010_2019_MAM<- calc(stack(SeasTas2000_2079[[42]], SeasTas2000_2079[[46]], SeasTas2000_2079[[50]], SeasTas2000_2079[[54]], SeasTas2000_2079[[58]], 
#                                     SeasTas2000_2079[[62]], SeasTas2000_2079[[66]], SeasTas2000_2079[[70]], SeasTas2000_2079[[74]], SeasTas2000_2079[[78]]), mean)
#DecSeasTas2020_2029_MAM<- calc(stack(SeasTas2000_2079[[82]], SeasTas2000_2079[[86]], SeasTas2000_2079[[90]], SeasTas2000_2079[[94]], SeasTas2000_2079[[98]], 
#                                     SeasTas2000_2079[[102]], SeasTas2000_2079[[106]], SeasTas2000_2079[[110]], SeasTas2000_2079[[114]], SeasTas2000_2079[[118]]), mean)
#DecSeasTas2030_2039_MAM<- calc(stack(SeasTas2000_2079[[122]], SeasTas2000_2079[[126]], SeasTas2000_2079[[130]], SeasTas2000_2079[[134]], SeasTas2000_2079[[138]], 
#                                     SeasTas2000_2079[[142]], SeasTas2000_2079[[146]], SeasTas2000_2079[[150]], SeasTas2000_2079[[154]], SeasTas2000_2079[[158]]), mean)
#DecSeasTas2040_2049_MAM<- calc(stack(SeasTas2000_2079[[162]], SeasTas2000_2079[[166]], SeasTas2000_2079[[170]], SeasTas2000_2079[[174]], SeasTas2000_2079[[178]], 
#                                     SeasTas2000_2079[[182]], SeasTas2000_2079[[186]], SeasTas2000_2079[[190]], SeasTas2000_2079[[194]], SeasTas2000_2079[[198]]), mean)
#DecSeasTas2050_2059_MAM<- calc(stack(SeasTas2000_2079[[202]], SeasTas2000_2079[[206]], SeasTas2000_2079[[210]], SeasTas2000_2079[[214]], SeasTas2000_2079[[218]], 
#                                     SeasTas2000_2079[[222]], SeasTas2000_2079[[226]], SeasTas2000_2079[[230]], SeasTas2000_2079[[234]], SeasTas2000_2079[[238]]), mean)
#DecSeasTas2060_2069_MAM<- calc(stack(SeasTas2000_2079[[242]], SeasTas2000_2079[[246]], SeasTas2000_2079[[250]], SeasTas2000_2079[[254]], SeasTas2000_2079[[258]], 
#                                     SeasTas2000_2079[[262]], SeasTas2000_2079[[266]], SeasTas2000_2079[[270]], SeasTas2000_2079[[274]], SeasTas2000_2079[[278]]), mean)
#DecSeasTas2070_2079_MAM<- calc(stack(SeasTas2000_2079[[282]], SeasTas2000_2079[[286]], SeasTas2000_2079[[290]], SeasTas2000_2079[[294]], SeasTas2000_2079[[298]], 
#                                     SeasTas2000_2079[[302]], SeasTas2000_2079[[306]], SeasTas2000_2079[[310]], SeasTas2000_2079[[314]], SeasTas2000_2079[[318]]), mean)
#
#DecSeasTas2000_2009_JJA<- calc(stack(SeasTas2000_2079[[3]], SeasTas2000_2079[[7]], SeasTas2000_2079[[11]], SeasTas2000_2079[[15]], SeasTas2000_2079[[19]], 
#                                     SeasTas2000_2079[[23]], SeasTas2000_2079[[27]], SeasTas2000_2079[[31]], SeasTas2000_2079[[35]], SeasTas2000_2079[[39]]), mean)
#DecSeasTas2010_2019_JJA<- calc(stack(SeasTas2000_2079[[43]], SeasTas2000_2079[[47]], SeasTas2000_2079[[51]], SeasTas2000_2079[[55]], SeasTas2000_2079[[59]], 
#                                     SeasTas2000_2079[[63]], SeasTas2000_2079[[67]], SeasTas2000_2079[[71]], SeasTas2000_2079[[75]], SeasTas2000_2079[[79]]), mean)
#DecSeasTas2020_2029_JJA<- calc(stack(SeasTas2000_2079[[83]], SeasTas2000_2079[[87]], SeasTas2000_2079[[91]], SeasTas2000_2079[[95]], SeasTas2000_2079[[99]], 
#                                     SeasTas2000_2079[[103]], SeasTas2000_2079[[107]], SeasTas2000_2079[[111]], SeasTas2000_2079[[115]], SeasTas2000_2079[[119]]), mean)
#DecSeasTas2030_2039_JJA<- calc(stack(SeasTas2000_2079[[123]], SeasTas2000_2079[[127]], SeasTas2000_2079[[131]], SeasTas2000_2079[[135]], SeasTas2000_2079[[139]], 
#                                     SeasTas2000_2079[[143]], SeasTas2000_2079[[147]], SeasTas2000_2079[[151]], SeasTas2000_2079[[155]], SeasTas2000_2079[[159]]), mean)
#DecSeasTas2040_2049_JJA<- calc(stack(SeasTas2000_2079[[163]], SeasTas2000_2079[[167]], SeasTas2000_2079[[171]], SeasTas2000_2079[[175]], SeasTas2000_2079[[179]], 
#                                     SeasTas2000_2079[[183]], SeasTas2000_2079[[187]], SeasTas2000_2079[[191]], SeasTas2000_2079[[195]], SeasTas2000_2079[[199]]), mean)
#DecSeasTas2050_2059_JJA<- calc(stack(SeasTas2000_2079[[203]], SeasTas2000_2079[[207]], SeasTas2000_2079[[211]], SeasTas2000_2079[[215]], SeasTas2000_2079[[219]], 
#                                     SeasTas2000_2079[[223]], SeasTas2000_2079[[227]], SeasTas2000_2079[[231]], SeasTas2000_2079[[235]], SeasTas2000_2079[[239]]), mean)
#DecSeasTas2060_2069_JJA<- calc(stack(SeasTas2000_2079[[243]], SeasTas2000_2079[[247]], SeasTas2000_2079[[251]], SeasTas2000_2079[[255]], SeasTas2000_2079[[259]], 
#                                     SeasTas2000_2079[[263]], SeasTas2000_2079[[267]], SeasTas2000_2079[[271]], SeasTas2000_2079[[275]], SeasTas2000_2079[[279]]), mean)
#DecSeasTas2070_2079_JJA<- calc(stack(SeasTas2000_2079[[283]], SeasTas2000_2079[[287]], SeasTas2000_2079[[291]], SeasTas2000_2079[[295]], SeasTas2000_2079[[299]], 
#                                     SeasTas2000_2079[[303]], SeasTas2000_2079[[307]], SeasTas2000_2079[[311]], SeasTas2000_2079[[315]], SeasTas2000_2079[[319]]), mean)
#
#DecSeasTas2000_2009_SON<- calc(stack(SeasTas2000_2079[[4]], SeasTas2000_2079[[8]], SeasTas2000_2079[[12]], SeasTas2000_2079[[16]], SeasTas2000_2079[[20]], 
#                                     SeasTas2000_2079[[24]], SeasTas2000_2079[[28]], SeasTas2000_2079[[32]], SeasTas2000_2079[[36]], SeasTas2000_2079[[40]]), mean)
#DecSeasTas2010_2019_SON<- calc(stack(SeasTas2000_2079[[44]], SeasTas2000_2079[[48]], SeasTas2000_2079[[52]], SeasTas2000_2079[[56]], SeasTas2000_2079[[60]], 
#                                     SeasTas2000_2079[[64]], SeasTas2000_2079[[68]], SeasTas2000_2079[[72]], SeasTas2000_2079[[76]], SeasTas2000_2079[[80]]), mean)
#DecSeasTas2020_2029_SON<- calc(stack(SeasTas2000_2079[[84]], SeasTas2000_2079[[88]], SeasTas2000_2079[[92]], SeasTas2000_2079[[96]], SeasTas2000_2079[[100]], 
#                                     SeasTas2000_2079[[104]], SeasTas2000_2079[[108]], SeasTas2000_2079[[112]], SeasTas2000_2079[[116]], SeasTas2000_2079[[120]]), mean)
#DecSeasTas2030_2039_SON<- calc(stack(SeasTas2000_2079[[124]], SeasTas2000_2079[[128]], SeasTas2000_2079[[132]], SeasTas2000_2079[[136]], SeasTas2000_2079[[140]], 
#                                     SeasTas2000_2079[[144]], SeasTas2000_2079[[148]], SeasTas2000_2079[[152]], SeasTas2000_2079[[156]], SeasTas2000_2079[[160]]), mean)
#DecSeasTas2040_2049_SON<- calc(stack(SeasTas2000_2079[[164]], SeasTas2000_2079[[168]], SeasTas2000_2079[[172]], SeasTas2000_2079[[176]], SeasTas2000_2079[[180]], 
#                                     SeasTas2000_2079[[184]], SeasTas2000_2079[[188]], SeasTas2000_2079[[192]], SeasTas2000_2079[[196]], SeasTas2000_2079[[200]]), mean)
#DecSeasTas2050_2059_SON<- calc(stack(SeasTas2000_2079[[204]], SeasTas2000_2079[[208]], SeasTas2000_2079[[212]], SeasTas2000_2079[[216]], SeasTas2000_2079[[220]], 
#                                     SeasTas2000_2079[[224]], SeasTas2000_2079[[228]], SeasTas2000_2079[[232]], SeasTas2000_2079[[236]], SeasTas2000_2079[[240]]), mean)
#DecSeasTas2060_2069_SON<- calc(stack(SeasTas2000_2079[[244]], SeasTas2000_2079[[248]], SeasTas2000_2079[[252]], SeasTas2000_2079[[256]], SeasTas2000_2079[[260]], 
#                                     SeasTas2000_2079[[264]], SeasTas2000_2079[[268]], SeasTas2000_2079[[272]], SeasTas2000_2079[[276]], SeasTas2000_2079[[280]]), mean)
#DecSeasTas2070_2079_SON<- calc(stack(SeasTas2000_2079[[284]], SeasTas2000_2079[[288]], SeasTas2000_2079[[292]], SeasTas2000_2079[[296]], SeasTas2000_2079[[300]], 
#                                     SeasTas2000_2079[[304]], SeasTas2000_2079[[308]], SeasTas2000_2079[[312]], SeasTas2000_2079[[316]], SeasTas2000_2079[[320]]), mean)
#
#DecSeasTas2000_2079 <- stack(DecSeasTas2000_2009_DJF, DecSeasTas2000_2009_MAM, DecSeasTas2000_2009_JJA, DecSeasTas2000_2009_SON,
#                             DecSeasTas2010_2019_DJF, DecSeasTas2010_2019_MAM, DecSeasTas2010_2019_JJA, DecSeasTas2010_2019_SON,
#                             DecSeasTas2020_2029_DJF, DecSeasTas2020_2029_MAM, DecSeasTas2020_2029_JJA, DecSeasTas2020_2029_SON,
#                             DecSeasTas2030_2039_DJF, DecSeasTas2030_2039_MAM, DecSeasTas2030_2039_JJA, DecSeasTas2030_2039_SON,
#                             DecSeasTas2040_2049_DJF, DecSeasTas2040_2049_MAM, DecSeasTas2040_2049_JJA, DecSeasTas2040_2049_SON,
#                             DecSeasTas2050_2059_DJF, DecSeasTas2050_2059_MAM, DecSeasTas2050_2059_JJA, DecSeasTas2050_2059_SON,
#                             DecSeasTas2060_2069_DJF, DecSeasTas2060_2069_MAM, DecSeasTas2060_2069_JJA, DecSeasTas2060_2069_SON,
#                             DecSeasTas2070_2079_DJF, DecSeasTas2070_2079_MAM, DecSeasTas2070_2079_JJA, DecSeasTas2070_2079_SON)
#



#DecSeasTasMin2000_2079 <- stackApply(SeasTasMin2000_2079, indices=yt, fun=mean, na.rm=TRUE) 
#DecSeasTasMax2000_2079 <- stackApply(SeasTasMax2000_2079, indices=yt, fun=mean, na.rm=TRUE) 
#DecSeasTas2000_2079 <- stackApply(SeasTas2000_2079 , indices=yt, fun=mean, na.rm=TRUE) 

