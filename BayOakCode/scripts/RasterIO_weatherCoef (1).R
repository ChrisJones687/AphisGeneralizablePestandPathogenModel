library(rgdal)
library(raster)

setwd("D:\\SCENIC\\PRISM")
fldCoef <- "WeatherCoeff"
if (!file.exists(file.path(fldCoef))) dir.create(file.path(getwd(),fldCoef))
#reference raster used to clip, project, and resample 
r <- raster("D:\\SOD-modeling\\layers\\UMCA_den_100m.img")

#list of folders
lstFld <- list.files(pattern='\\.zip')
lstName <- sapply(lstFld, function(x) unlist(strsplit(x, '\\.'))[1])

for (dir in lstName){
        
        fldPrj <- "NAD27_aea"
        if (!file.exists(file.path(dir,"NAD27_aea"))) dir.create(file.path(getwd(),dir,fldPrj))
        
        #lstRast <- list.files(dir, pattern='\\.tiff') 
        lstPrec <- grep("pcpn_gt2.5", list.files(dir), value=TRUE)
        message(paste0("Starting computation for precipitation data: ", substr(lstPrec[1],12,15)))
        cat('\n')
        for (lyr in lstPrec){

                rst <- raster(file.path(dir, lyr))
                message(paste0("Projecting raster ", lyr, "..."))
                outName <- paste0(sub(".tiff$", "", lyr), "_prj")
                rst_prj <- projectRaster(from = rst, to = r,
                                         method = "ngb",
                                         filename = file.path(dir,fldPrj,outName),
                                         format = "HFA",
                                         dataType = 'INT1U',
                                         overwrite=TRUE)
                message(paste0("...Saved as ", file.path(dir,fldPrj,outName)))
                message("Computing precipitation coefficients")
                m_coef <- rst_prj / 7
                #save coef raster to common folder
                cnt <- na.omit(as.numeric(unlist(strsplit(unlist(lyr), "[^0-9]+"))))[4]
                m_coef_outName <- ifelse(cnt < 10, paste0(substr(lstPrec[1],12,15),"_m0", cnt),
                                                   paste0(substr(lstPrec[1],12,15),"_m", cnt))
                writeRaster(m_coef, filename = file.path(fldCoef,m_coef_outName),
                                    format = "HFA",
                                    dataType = 'FLT4S',
                                    overwrite=TRUE)
                message(paste0("...Saved as ", file.path(fldCoef,m_coef_outName)))
        }
        lstTmean <- grep("avgt_mean", list.files(dir), value=TRUE)
        cat('\n\n') 
        message(paste0("Starting computation for avg temperature data: ", substr(lstTmean[1],11,14)))
        cat('\n')     
        for (lyr in lstTmean){
                
                rst <- raster(file.path(dir, lyr))
                message(paste0("Projecting raster ", lyr, "..."))
                outName <- paste0(sub(".tiff$", "", lyr), "_prj")
                rst_prj <- projectRaster(from = rst, to = r,
                                         method = "ngb",
                                         filename = file.path(dir,fldPrj,outName),
                                         format = "HFA",
                                         dataType = 'INT1U',
                                         overwrite=TRUE)
                message(paste0("...Saved as ", file.path(dir,fldPrj,outName)))
                message("Computing avg temperature coefficients")
                c_coef <- -0.066 + 0.056 * rst_prj - 0.0036 * (rst_prj - 15)**2 - 0.0003 * (rst_prj - 15)**3
                c_coef[c_coef < 0] <- 0 #restrain condition
                #save coef raster to common folder
                cnt <- na.omit(as.numeric(unlist(strsplit(unlist(lyr), "[^0-9]+"))))[2]
                c_coef_outName <- ifelse(cnt < 10, paste0(substr(lstTmean[1],11,14),"_c0", cnt),
                                         paste0(substr(lstPrec[1],12,15),"_c", cnt))
                writeRaster(c_coef, filename = file.path(fldCoef,c_coef_outName),
                            format = "HFA",
                            dataType = 'FLT4S',
                            overwrite=TRUE)
                message(paste0("...Saved as ", file.path(fldCoef,c_coef_outName)))
        }
}
message('DONE')
