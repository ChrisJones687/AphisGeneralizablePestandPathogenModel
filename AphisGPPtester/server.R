library(leaflet)
library(RColorBrewer)
library(scales)
library(lattice)
library(dplyr)
library(shiny)
library(raster)
library(rgdal)

# writeZip <- function(x, file, filename, format, ...) {
#   if (format=="ESRI Shapefile") {
#     writeOGR(x, "./", filename, format, overwrite_layer=T, check_exists=T)
#   } else {
#     writeRaster(x, filename, format, bylayer=F, overwrite=T, ...)
#   }
#   f <- list.files(pattern=paste0(strsplit(filename, ".", fixed=T)[[1]][1], ".*"))
#   zip(paste0(filename, ".zip"), f, flags="-9Xjm", zip="zip")
#   file.copy(paste0(filename, ".zip"), file)
#   file.remove(paste0(filename, ".zip"))
# }

server <- function(input, output) {
  options(shiny.maxRequestSize=65000*1024^2) 
  
  # Creates the text file that is downloaded upon model completion
  #output$model <- renderPrint({c("Model inputs are:",input$wind, input$windData, input$temp, input$tempData, input$precip, input$precipData)})
  
  # Used to set the initial zoom of the map and color of the rasters
  r <- raster("C:\\Users\\chris\\Dropbox\\Projects\\APHIS\\BayOakCode/layers/UMCA_den_100m.img")
  pal <- colorNumeric(c("#0C2C84","#41B6C4","#FFFFCC"), values(r), na.color = "transparent")
  
  #Creates the text saying the model is running when the action button is pressed
  modeltext <- eventReactive(input$run, {"Model is now running"})
  output$modelText <- renderText({modeltext()})
  
  #inTotalSpeciesData <- reactive(input$totalSpeciesData)
  #inTSD <- raster(inTotalSpeciesData())
  
  modelRun <- observeEvent(input$run, 
                           {pest("./layers/UMCA_den_100m.img","./layers/OAKS_den_100m.img","./layers/TPH_den_100m.img", "./layers/init_2000_cnt.img",
                                 input$start, input$end, input$seasonQ, input$seasonMonths[1],input$seasonMonths[2], 
                                 input$sporeRate, input$windQ, input$windDir, './layers/weather/weatherCoeff_2000_2014.nc')
    })
  
  inTempData <- renderText(input$tempData$datapath)
  # Plot the data
  output$plotData <- renderLeaflet({
    if (is.null(input$windData)==TRUE && is.null(input$tempData)==TRUE && is.null(input$totalSpeciesData)==TRUE){
      leaflet(height = "300px") %>%
        addProviderTiles("Esri.WorldImagery", group="background 1") %>%
        addRasterImage(r, opacity=0.0) %>%
        addLegend("bottomright", pal = pal, values = values(r),
                title = "Host species",
                #labFormat = labelFormat(prefix = "$"),
                opacity = 1)
    } else if (is.null(input$windData)==TRUE && is.null(input$tempData)==TRUE && is.null(input$totalSpeciesData)==FALSE) {
      inTotalSpeciesData <- input$totalSpeciesData
      leaflet(height = "300px") %>%
        addProviderTiles("Esri.WorldImagery", group="background 1") %>%
        addRasterImage({raster(inTotalSpeciesData$datapath)}, opacity=0.6) %>%
        addLegend("bottomright", pal = pal, values = values(r),
                title = "Host species",
                #labFormat = labelFormat(prefix = "$"),
                opacity = 1)
    } else if (is.null(input$windData)==FALSE && is.null(input$tempData)==TRUE) {
      inWind <- input$windData
      leaflet(height = "300px") %>%
        addProviderTiles("Esri.WorldImagery", group="background 1") %>%
        addRasterImage({raster(inWind$datapath)}, opacity=0.6)
    } else if (is.null(input$windData)==TRUE && is.null(input$tempData)==FALSE){
      inTemp <- input$tempData
      leaflet(height = "300px") %>%
        addProviderTiles("Esri.WorldImagery", group="background 1") 
        addRasterImage({raster(inTemp$datapath)}, opacity=0.6)
    } else if (is.null(input$windData)==FALSE && is.null(input$tempData)==FALSE) {
      inTemp <- input$tempData
      inWind <- input$windData
      leaflet(height = "300px") %>%
        addProviderTiles("Esri.WorldImagery", group="background 1") %>%
        addRasterImage({raster(inTemp$datapath)}, opacity=0.6) %>%
        addRasterImage({raster(inWind$datapath)}, opacity=0.4)
    }

  })
  
  output$extentMatch <- renderText({
    if (is.null(input$windData)==FALSE && is.null(input$tempData)==FALSE) {
      inTemp <- input$tempData
      inWind <- input$windData
      exMatch <- extent(inTemp)==extent(inWind)
      if (exMatch==TRUE){
        print("Your extents match")
      } else if (exMatch == FALSE) {
        print("Your extents don't match would you like to resample")
      }
    }
  })
  
  # Allows for the downloading of the user manual when the download link is pressed
  output$pdf <- downloadHandler("generalizablepestandpathogenmodel.pdf", content = function(file){
    file.copy("C:\\Users\\Chris\\Desktop\\Generalizable Pest and Pathogen Model User.pdf",file)
  })
  
  # output$saveData <- downloadHandler(
  #   filename, function(file) {
  #     x <- "some data to write"
  #     path <- "path to write files to"
  #     format <- "argument passed to writeRaster or writeOGR (needed in my case)"
  #     writeZip(x, file, path, format)
  #   })
  # 
}
