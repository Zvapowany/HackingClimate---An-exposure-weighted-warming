##Loading libraries
library(ncdf4)                  # needed to read netcdf data
library(ncdf4.helpers)          # additional support functions for netcdf data
library(fields)                 # provides image.plot()
library(rgdal)                  # load shapefiles
library(abind)                  # join arrays together

##Reading data
fnm <- "~/UKCP18/day/tas_rcp85_land-rcm_uk_12km_01_day_20201201-20301130.nc"
nc <- nc_open(fnm)  # open a connection to the netcdf file
tas <- ncvar_get(nc, "tas") # read the data into R
dim(tas)  # check the dimensions of your variable. Daily UKCP18 data is stored in arrays of 82 · 112 · 3600, where the dimensions are eastings · northings · date.

##Labelling dimensions

dimnames(tas) <- list("E" = ncvar_get(nc, "projection_x_coordinate"),
                      "N" = ncvar_get(nc, "projection_y_coordinate"),
                      "date" = substr(nc.get.time.series(nc, "tas"),1,10))  # label the dimensions of your array with the dimensions of the netcdf - this means that you can refer to them by E/N/date

head(dimnames(tas)$date); tail(dimnames(tas)$date) # dates are stored as a date-time string: we're only really interested in the date, so we use 'substring' to extract the date only
