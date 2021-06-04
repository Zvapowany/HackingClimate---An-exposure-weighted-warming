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




