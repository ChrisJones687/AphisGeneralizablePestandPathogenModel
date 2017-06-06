#---------------------------------------------------------------------------------------------------------------
# Name:         GeneralizedPestandPathogen.r
# Purpose:      Generalizable spread model of pests and pathogens over forest or agriculutural ecosystems.
# Author:       Chris Jones
# Email:        cjones1688@gmail.com
# Created:      06/06/2017
# Copyright:    (c) 2017 by Chris Jones
# License:      GNU General Public License (GPL)
# Software:     In developement using R version 3.3.1 (http://www.r-project.org/)
#-----------------------------------------------------------------------------------------------------------------------

#install packages
#install.packages(c("rgdal","raster","lubridate","CircStats","Rcpp", "rgrass7", "optparse", "plotrix", "ncdf4", "dismo", "sp"))

#load packages:
suppressPackageStartupMessages(library(raster))    #Raster operation and I/O. Depends R (≥ 2.15.0)
suppressPackageStartupMessages(library(rgdal))     #Geospatial data abstraction library. Depends R (≥ 2.14.0)
suppressPackageStartupMessages(library(lubridate)) #Make dealing with dates a little easier. Depends R (≥ 3.0.0)
suppressPackageStartupMessages(library(CircStats)) #Circular Statistics - Von Mises distribution
suppressPackageStartupMessages(library(Rcpp))      #Seamless R and C++ Integration. Depends R (≥ 3.0.0)
suppressPackageStartupMessages(library(plotrix))   #Add text annotations to plot
suppressPackageStartupMessages(library(ncdf4))     #work with NetCDF datasets
suppressPackageStartupMessages(library(dismo))     #work regressio
suppressPackageStartupMessages(library(sp))        #work with NetCDF datasets