library(ncdf)

#years of interest:
start <- 2000
end <- 2014

#list of ALL weather layers
lst <- dir('./layers/weather', pattern='\\.img$', full.names=T)

#sublist of weather coefficients
Mlst <- lst[grep("_m", lst)] #M = moisture; 
Clst <- lst[grep("_c", lst)] #C = temperature;

Mlst_sub <- grep(paste(as.character(seq(start,end)), collapse="|"), Mlst, value=TRUE)  #use only the raster files matching the years of interest
Clst_sub <- grep(paste(as.character(seq(start,end)), collapse="|"), Clst, value=TRUE)  #use only the raster files matching the years of interest

Mstack <- stack(Mlst_sub) 
Cstack <- stack(Clst_sub) 

#----------------
# Make dimensions
#----------------
xvals <- xFromCol(Mstack, col=1:ncol(Mstack))
yvals <- yFromRow(Mstack, row=1:nrow(Mstack))
 
nx <- length(xvals)
ny <- length(yvals)
nt <- length(Mstack@layers) #MUST be same as Cstack@layers

xdim <- dim.def.ncdf( 'X', 'meters_easting', xvals )
ydim <- dim.def.ncdf( 'Y', 'meters_northing', yvals )
tdim <- dim.def.ncdf('time', units='weeks since 2000-01-01', nt, unlim=TRUE)

#---------
# Make var
#---------
mv <- -9999     # missing value
var_mcf <- var.def.ncdf( "Mcoef", "units", list(xdim,ydim,tdim), mv )
var_ccf <- var.def.ncdf( "Ccoef", "units", list(xdim,ydim,tdim), mv )

#---------------------
# Make new output file
#---------------------
output_fname <- paste('weatherCoeff_', start, '_', end, '.nc', sep='')
ncid_new <- create.ncdf( output_fname, list(var_mcf, var_ccf))

#-------------------------------
# Put data into the file
#-------------------------------
for(i in 1:nt) {
	data2d_m <- as.matrix(Mstack[[i]])   #one should take the transpose of the matrix if you want to match Mstack[[i]] at this stage
	data2d_c <- as.matrix(Cstack[[i]])
	put.var.ncdf(ncid_new, var_mcf, data2d_m, start=c(1,1,i), count=c(-1,-1,1))
	put.var.ncdf(ncid_new, var_ccf, data2d_c, start=c(1,1,i), count=c(-1,-1,1))
} 
#--------------------------
# Close our new output file
#--------------------------
close.ncdf( ncid_new )


#check if things are right
nc <- open.ncdf(output_fname)
#Var <- get.var.ncdf(nc, "var_name_here")  #get variable of interest into an R array
Var <- get.var.ncdf(nc, "Mcoef")
#plot original raster and netcdf slice...they must match!!
plot(Mstack[[10]])
plot(raster(Var[,,10]))  #if you used the transpose in the LOOP above, then remember to also use t(Var[,,]) here!
close.ncdf(nc)