##WEATHER RANKING:
lst <- dir('D:\\SOD-modeling\\layers\\weather\\weatherCoeff_19902008_20082014PRISM', pattern='\\.img$', full.names=T)
Mlst <- lst[grep("_m", lst)]
Clst <- lst[grep("_c", lst)]

years <- 1990:2014

prod_mc_TOT <- rep(0, length(years))
cnt <- 0

for (yr in years){
	
	cat(paste('year', yr), '\n')
	cnt <- cnt + 1

	Mlst_sel <- grep(paste0(as.character(yr),"_"), Mlst, value=TRUE)
	Clst_sel <- grep(paste0(as.character(yr),"_"), Clst, value=TRUE)
	
	Mstack <- stack(Mlst_sel) #M = moisture; 
	Cstack <- stack(Clst_sel) #C = temperature;
	
	prod_mc <- rep(0, nlayers(Mstack))

	for (wk in 1:nlayers(Mstack)){
		
		cat(paste('week', wk), '\n')
		prod_mc_ras <- Mstack[wk] * Cstack[wk]

		#sum over all cells of the raster for current week
		prod_mc[wk] <- sum(prod_mc_ras, na.rm=T)
	}
	
	#sum over all weeks of the current year
	prod_mc_TOT[cnt] <- mean(prod_mc, na.rm=T)

}