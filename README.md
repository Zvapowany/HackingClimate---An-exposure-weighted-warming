# HackingClimate-An-exposure-weighted-warming

An exposure weighted warming - UCL Hackaton

Group 2 members: Marcin Gajewski, Salma Al-Zadjali, Yasna Palmeiro Silva

The aim of the project was to calculate population weighted exposure projection for years 2000-2080 for United Kingdom and see how it compares to simple avarage temperature projection. 
Population-weighted averages assure that degree-day averages reflect conditions in more densely populated areas of the country. Hence, it is a powerfull tool in asessing exposure to climate change anomalies. 

In order to calculate population weighted temperature anomalies we used UK Climate Projection 2018 (UKCP18) data, supplied by MetOffice, which was analysed for: 
1. Mean air temperature at 1.5m (°C)
2. Maximum air temperature at 1.5m (°C)
3. Minimum air temperature at 1.5m (°C)
       
In order to obtain meaningful results we used monthly data to calculate temperature anomalies compared to historical baseline (1980-2000). 
In order to calculate population-weighted average anomalies, we used population estimations from 2000 onward.

As a result we have obtained linear plots showing population weighted and area weighted temperature projections for UK, as well as individual countries of UK (England, Northern Ireland, Scotland, Wales):

a) [Mean temperature projections](Mean_temperature_projection.md)

b) [Mean maximum temperature projections](Max_temperature_projection.md)

c) [Mean minimum temperature projections](Min_temperature_projection.md)

'''
###############################################################################
#   TEMPERATURE PROJECTIONS
###############################################################################

tasmin.nc <- nc_open("/nfs/cfs/home4/rejb/rejbypa/UKCP18_R/mon/tasmin_rcp85_land-rcm_uk_12km_01_mon_198012-208011.nc")
print(tasmin.nc)
names(tasmin.nc$var) #"tasmin"
tasmin <- brick("/nfs/cfs/home4/rejb/rejbypa/UKCP18_R/mon/tasmin_rcp85_land-rcm_uk_12km_01_mon_198012-208011.nc")
'''
