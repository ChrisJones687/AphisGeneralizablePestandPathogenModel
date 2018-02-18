#---------------------------------------------------------------------------------------------------------------
# Name:         GeneralizedPestandPathogen.r
# Purpose:      Generalizable spread model of pests and pathogens over forest or agriculutural ecosystems.
# Author:       Chris Jones
# Email:        cjones1688@gmail.com
# Created:      06/06/2017
# Copyright:    (c) 2017 by Chris Jones and Qiang (Jack) Zhang
# License:      GNU General Public License (GPL)
# Software:     In developement using R version 3.4 (http://www.r-project.org/)
#-----------------------------------------------------------------------------------------------------------------------

## install packages (only needed on first run of the model)
#install.packages(c("rgdal","raster","lubridate","CircStats","Rcpp", "rgrass7", "optparse", "plotrix", "ncdf4", "dismo", "sp"))

## load packages:
suppressPackageStartupMessages(library(raster))    #Raster operation and I/O. Depends R (≥ 2.15.0)
suppressPackageStartupMessages(library(rgdal))     #Geospatial data abstraction library. Depends R (≥ 2.14.0)
suppressPackageStartupMessages(library(lubridate)) #Make dealing with dates a little easier. Depends R (≥ 3.0.0)
suppressPackageStartupMessages(library(CircStats)) #Circular Statistics - Von Mises distribution
suppressPackageStartupMessages(library(Rcpp))      #Seamless R and C++ Integration. Depends R (≥ 3.0.0)
suppressPackageStartupMessages(library(plotrix))   #Add text annotations to plot
suppressPackageStartupMessages(library(ncdf4))     #work with NetCDF datasets
suppressPackageStartupMessages(library(dismo))     #Regression for ecological datasets
suppressPackageStartupMessages(library(sp))        #Classes and methods for spatial data

## Define data paths for climate and tamarisk leaf beetle data 
data_path = "C:\\Users\\cmjone25\\Desktop\\Daymet Data\\"
data_path2 = "C:\\Users\\cmjone25\\Dropbox\\Projects\\APHIS\\TamariskLeafBeetleData"

## Set data variables for calculations (noting assumption)
# Days in a year
maxdays <- 365

# Reproduction length (Average over the course of the reproductive season) and temperature range
reproductive_period <- 21
min_temperature <- 0
max_temperature <- 35

# Host density (tamarisk)
tamarisk_rate = 25
N_max=100 # maximum density

# Set critical daylength
critical_daylength <- 43500

# Mortality rate (will make this dependant on critical day length (i.e. more mortality if diapause starts early at a location) will need to be an estimate for many of the species)
mortality = 10

## Read in data into model (This is the section that you will change for different species)
# Climate maximum and minimum temperature and daylength (2007 to 2011).
tmax07 = stack(paste(data_path,"tmax2007small.tif",sep="/"))
tmax08 = stack(paste(data_path,"tmax2008small.tif",sep="/"))
tmax09 = stack(paste(data_path,"tmax2009small.tif",sep="/"))
tmax10 = stack(paste(data_path,"tmax2010small.tif",sep="/"))
tmax11 = stack(paste(data_path,"tmax2011small.tif",sep="/"))
tmax <- addLayer(tmax07,tmax08,tmax09,tmax10,tmax11)

tmin07 = stack(paste(data_path,"tmin2007small.tif",sep="/"))
tmin08 = stack(paste(data_path,"tmin2008small.tif",sep="/"))
tmin09 = stack(paste(data_path,"tmin2009small.tif",sep="/"))
tmin10 = stack(paste(data_path,"tmin2010small.tif",sep="/"))
tmin11 = stack(paste(data_path,"tmin2011small.tif",sep="/"))
tmin <-addLayer(tmin07,tmin08,tmin09,tmin10,tmin11)

# Read in day length data and prepare data for simulations
daylength = stack(paste(data_path,"daylengthsmall.tif",sep="/"))
daylengthData = addLayer(daylength,daylength,daylength,daylength,daylength)
dlmax = cellStats(daylength, 'max') # Get maximum day length for each day of the year for our location
dlfirst = min(which((dlmax >= critical_daylength)==TRUE)) # find the first day the year when critical day length is reached
dllast = max(which((dlmax >= critical_daylength)==TRUE)) # find the last day of the year when critical day length is reached

# Read in host species data (in this case we are testing different precentages of tamarisk within 1 km of the major rivers)
riverAssumption = raster(paste(data_path2,"TamariskRiverAssumption1km.tif",sep="/"))

# Read in TLB data for model runs and validation
p07 = raster(paste(data_path2,"2007TLB.tif",sep="/"))
p08 = raster(paste(data_path2,"2008TLB.tif",sep="/"))
p09 = raster(paste(data_path2,"2009TLB.tif",sep="/"))
p10 = raster(paste(data_path2,"2010TLB.tif",sep="/"))
p11 = raster(paste(data_path2,"2011TLB.tif",sep="/"))

## Create data frames for calculations
dp07 <- as.data.frame(p07)
dp08 <- as.data.frame(p08)
dp09 <- as.data.frame(p09)
dp10 <- as.data.frame(p10)
dp11 <- as.data.frame(p11)
host <- as.data.frame(riverAssumption)

# Create neighborhood functions
neighbor4_find <- function(cell_id,nrow,ncol){
  neighbors <- vector()
  cell_w <- cell_id - 1
  if((cell_id-1)%%ncol != 0) neighbors <- c(neighbors,cell_w)
  cell_n <- cell_id - ncol
  if(cell_id > ncol) neighbors <- c(neighbors,cell_n)
  cell_e <- cell_id + 1
  if(cell_id%%ncol != 0) neighbors <- c(neighbors,cell_e)
  cell_s <- cell_id + ncol
  if(cell_id <=ncol*(nrow-1)) neighbors <- c(neighbors,cell_s)
  neighbors <- c(neighbors,cell_id)
  return(neighbors)
}

neighbor8_find <- function(cell_id,nrow,ncol){
  neighbors <- neighbor4_find(cell_id,nrow,ncol)
  cell_nw <- cell_id - ncol - 1
  if(cell_id > ncol & (cell_id-1)%%ncol != 0) neighbors <- c(neighbors,cell_nw)
  cell_ne <- cell_id - ncol + 1
  if(cell_id > ncol & cell_id%%ncol != 0) neighbors <- c(neighbors,cell_ne)
  cell_sw <- cell_id + ncol - 1 
  if(cell_id <=ncol*(nrow-1) & (cell_id-1)%%ncol != 0) neighbors <- c(neighbors,cell_sw)
  cell_se <- cell_id + ncol + 1 
  if(cell_id <=ncol*(nrow-1) & cell_id%%ncol != 0) neighbors <- c(neighbors,cell_se)
  neighbors <- c(neighbors,cell_id)
  return(neighbors)
}

drawPoints <- function(cellid){
  x <-cellid %/% n_column
  y <-cellid %% n_column
  points(x,y)
}

####
# Φ_infulence = b_k * I_it * (S_jt/N_max) * K_k(d_ij;α_k)/d_ij
# i: local cell
# j: simulated cell
# t: time step
# k: short or long distance, ignore long distance dispersal in assumption
# b: amount of spores, determined by spore_kernal function
# I_it: Number of infected host units in cell i at time t
# S_jt: Number of susceptible host units in cell i at time t
# N_max: Maximum number of host units in any cell
# K: Dispersal kernel between cells i and j
# d_ij: distance between i and j, consider the distance between adjacent cells equals to 1
# α: scale parameters (short distance or long distance), abandoned in our simulation
####
#Revised function:
# Φ_infulence = b * I_it * (S_jt/N_max) * K(d_ij)

influence_calculate <- function(cell_id, n_row, n_column, start_time, beetle_num, tamariskmap){
  
  neighbors <- neighbor8_find(cell_id, n_row, n_column)
  spore <- spore_kernal(cell_id, start_time)
  I_it <- beetle_num
  neighbor_probabality <- vector()
  influence_sum <- 0
  influencemap <- data.frame(matrix(0, n_row*n_column))
  outcomemap <- data.frame(matrix(0, n_row*n_column))
  
  #calculate the influence of the trees in the neighborhood
  for(n in neighbors){
    #Going to add the distance to water
    #Going to replace the fake tamarisk data
    influencemap[n,] <- tamariskmap[n,] / N_max
  }
  
  for(n in neighbors){
    influence_sum <- influence_sum + influencemap[n,]
  }
  
  for(n in neighbors){
    neighbor_probabality <-c(neighbor_probabality,influencemap[n,]/influence_sum)
  }
  
  spreads <- sample(neighbors, size= I_it * spore, replace=TRUE, prob=neighbor_probabality)
  #print(neighbor_probabality)
  #print(spreads)
  
  for (i in spreads){
    outcomemap[i,] <- outcomemap[i,] + 1
  }
  return(outcomemap)
}



#dispersal kernal is a function of distance
#dispersal kernal is a function to generate probability to disperse the spore
#dispersal kernal could determine if the beetles survive locally or not

dispersal_kernal <- function(){
  
}

#Spore kernal is a function of temperature

spore_kernal <- function(cell_id,start_time){
  for (round in 0:reproductive_period){
    print(tmax[[start_time+round]][cell_id])
    print(tmin[[start_time+round]][cell_id])
    if(tmax[[start_time+round]][cell_id]<max_temperature & tmin[[start_time+round]][cell_id]>min_temperature)
    {}
    else{
      return(1)
    }
  }
  return(400)
}


#Main Simulation

n_row <- 335
n_column <- 239
plot(0:335, 0:335, type = "n")

#Initializing influence data
dinfluence <- data.frame(matrix(0, n_row*n_column))
colnames(dinfluence) <- c('influence')

#Initializing Tamarisk data
dtamarisk <- host * tamarisk_rate
dtamarisk[is.na(dtamarisk)] <- 0

#Initializing basemap data
start_present <- dp07
start_present[is.na(start_present)] <- 0
basemap <- start_present

#Initialize Simulation Stack
#Simulation Stack stores the inforamtion of cells to simulate
day = 1
new_present <- basemap
newstart <- which(!new_present$X2007TLB == 0)
simuStack <- data.frame(CellID = integer(),
                        BeetleNum = integer(),
                        Day = integer(),
                        stringsAsFactors=FALSE)
simuHeading = c("CellID","BeetleNum","Day")
names(simuStack) <- simuHeading

for (p in newstart){
  simuStack <- rbind(simuStack,data.frame("CellID"=p, "BeetleNum"=basemap[p,],"Day" =1))
  drawPoints(p)
}


#Simulation
while(is.data.frame(simuStack) && nrow(simuStack)!=0){
  currentCell <- simuStack[1,]
  simuStack <- simuStack[-c(1),]
  cell_id <- currentCell$CellID
  beetle_number <- currentCell$BeetleNum
  start_day <- currentCell$Day
  if(daylengthData[[start_day]][cell_id] > critical_daylength){
    cat(sprintf("Day: %i\t Cell: %i\n",start_day,cell_id))
    newmap <- influence_calculate(cell_id, n_row, n_column, start_time=start_day,beetle_number,dtamarisk)
    start_day <- start_day + 21
    if(start_day < maxdays){
      new_points <- which(newmap != 0)
      for (p in new_points){
        simuStack <- rbind(simuStack,data.frame("CellID"=p, "BeetleNum"=newmap[p,],"Day" = start_day))
      }
      simuStack <- simuStack[with(simuStack, order(simuStack$Day)), ]
    }
    new_points <- which(newmap != 0)
    for (p in new_points){
      simuStack <- rbind(simuStack,data.frame("CellID"=p, "BeetleNum"=newmap[p,],"Day" = start_day))
      drawPoints(p)
    }
    simuStack <- simuStack[with(simuStack, order(simuStack$Day)), ]
  }
  else{
    cat(sprintf("Day %i on Cell %i is not long enough\n",start_day,cell_id))
    start_day <- start_day + 1
    if(start_day < maxdays){
      simuStack <- rbind(simuStack,data.frame("CellID"=cell_id, "BeetleNum"=beetle_number,"Day" = start_day))
      simuStack <- simuStack[with(simuStack, order(simuStack$Day)), ]
    }
  }
}


