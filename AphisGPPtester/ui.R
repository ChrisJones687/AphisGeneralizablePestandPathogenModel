library(shiny)
#library(shinyjs)
library(raster)
library(rgdal)

ui <- fluidPage(title = "Generalizable Pest and Pathogen Model",#theme = "bootstrap.css",
  style = "background-color: black",

  
  # Add Title to App
  h1("Generalizable Pest and Pathogen Model", style = "color: green"),
  # Create a sidebar for variable inputs that react to user changes
  fluidRow(
    column(width = 4,
           # Create panel for Species paramaters
           wellPanel(h3("Species Parameters", icon("bug")),
                     style = "background-color: #54ACC1; border: #ADBD60; color: black",
                     textInput("pest", "Enter the name of the species being simulated (currently not used exept to keep track of model run information)"),
                     
                     # Create date range box for simulation
                     #dateRangeInput("date", "Start date to end date"),
                     numericInput("start", "Year to start simulation", value = 2000, min = 1960, max = 2020),
                     numericInput("end", "Year to end simulation", value = 2010, min = 1960, max = 2020),
                     # Create seasonality box
                     selectInput(inputId = "seasonQ", label = "Does spread occur during specific seasons (Currently set up to take in specific months but will be set to work with daylength as well)?", choices = c("YES","NO")),
                     conditionalPanel(
                       condition = "input.seasonQ == 'YES'", sliderInput("seasonMonths", "Months that contribute to the spread of the pest or pathogen?", value = c(1,9), min =1, max =12, step =1)
                     ),
                     numericInput(inputId ="sporeRate", label = "Enter the rate of reproduction or spore production.", value = "4.4", min=0, max = 100, step = 0.1)
           ),
           
           # Create panel for Host variables
           wellPanel(h3("Hosts", icon("tree")),
                     style = "background-color: #ADBD60; border: #ADBD60; color: black",
           # Create text box for pest or pathogen being simulated
           
           # Create input for host data
           selectInput(inputId = "hostQ", label = "Is the pest or pathogen a generalist or specialist?", choices = c("specialist","generalist")),
           conditionalPanel(
             condition = "input.hostQ == 'generalist'", numericInput("hostMulti", "How many hosts affect the spread of the pest or pathogen?", value = 2, min =1, max =10, step =1)
           ),
           conditionalPanel(
             condition = "input.hostQ == 'generalist' && input.hostMulti == 2",
             fileInput("hostDataM1", "Select your raster file for your 1st host!", accept = c(".tif")),
             fileInput("hostDataM2", "Select your raster file for your 2nd host!", accept = c(".tif"))
           ),
           conditionalPanel(
             condition = "input.hostQ == 'specialist'", fileInput("hostDataSingle", "Select your raster file for your host data", accept = c(".tif"))
           ),
           fileInput("totalSpeciesData", "Select your raster file for total species data", accept = c(".tif"))
           ),
           
           # Create panel for Environmental variables
           wellPanel(h3("Environmental Effects", icon("sun-o"), style = "color: black"), style = "color: black;background-color: #E7B15F; border: #E7B15F",
           # Create box asking if wind affects spread and if wind data is available
           selectInput(inputId = "windQ", label = "Does wind affect spread and do you have wind data?", choices = c("NO","YES")),
           conditionalPanel(
             condition = "input.windQ == 'YES'", selectInput("windType", "Is your wind data raster or direction form!", choices = c("Raster","Direction"))
           ),
           
           # add a data input box if yes is selected
           conditionalPanel(
             condition = "input.windQ == 'YES' && input.windType == 'Raster'", fileInput("windData", "Select your raster file for your wind data!", accept = c(".tif"))
           ),
           
           # add a data type option for wind in case detailed wind data isn't available but predominate wind direction
           conditionalPanel(
             condition = "input.windQ == 'YES' && input.windType == 'Direction'", selectInput("windDir", "What is the predominate wind direction!", choices = c("NE","E","SE","S","SW","W","NW","N"))
           ),
           
           # Create box asking if temperature affects spread and if Temperature data is available
           selectInput(inputId = "temp", label = "Does temperature affect spread and do you have temperature data?", choices = c("NO","YES")),
           # add a data input box if yes is selected
           conditionalPanel(
             condition = "input.temp == 'YES'", fileInput("tempData", "Select your raster file for your wind data!", accept = c(".tif", ".nc"))
           ),
           # Create box asking if precipitation affects spread and if precipitation data is available
           selectInput(inputId = "precip", label = "Does precipitation affect spread and do you have temperature data?", choices = c("NO","YES")),
           # add a data input box if yes is selected
           conditionalPanel(
             condition = "input.precip == 'YES'", fileInput("precipData", "Select your raster file for your wind data!", accept = c(".tif"))
           )

           # Create an Action button that will run the model when pressed
           
           ),
           actionButton("run", "Run Model")
    ),
    
    column(width = 8,
      # Create a text file with all input values 
      #textOutput("model"),
      
      # Let the user know if their extents match
      #textOutput("extentMatch")
      # Create a text box saying the the model is running while model is processes data
      verbatimTextOutput("modelText"),
      # Create a download link for the user Manual
      downloadLink("pdf", "Download User Manual ", icon("cloud-download"), style = "color: green"),

      
      leafletOutput("plotData", height = "600px")
    )
))

