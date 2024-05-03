require(shinythemes)
require(shinyWidgets)
require(plyr)
require(data.table)
require(tidyverse)
require(ggplot2)
require(gridExtra)
require(waiter)
require(pspline)
require(ggplate)

ui <- fluidPage(
  
  use_waiter(),
  useAttendant(),

  
  theme = shinytheme("slate"),
  titlePanel("Shiny Arrhythmia Detection Tool"),
  
  sidebarLayout(
    sidebarPanel(
      tags$h3("CSV Timepoint Imports"),
      fileInput("CSV1", "Choose CSV File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv"),
                width = "100%"),
      br(),
      
      tags$h3("Select Detection Parameters"),
      if (exists("CSV1")) {
        selectInput("search", label = "Select Well to Analyse...", choices = df[[1]])
      } else {
        selectInput("search", label = "Select Well to Analyse...", choices = "No CSV")
      },
      fluidRow(column(12, tabsetPanel(
        tabPanel("Trace & Peaks",
                 fluidRow(sliderInput("smooth_slider", "Smoothing Value:",
                                      min = 0.001, max = 0.02, value = 0.002, step = 0.001),
                          sliderInput("peak_slider", "Points to detect peak:",
                                      min = 1, max = 20, value = 6, step = 1),
                          sliderInput("trough_slider", "Points to detect trough:",
                                      min = 1, max = 20, value = 5, step = 1),
                          sliderInput("pref_slider", "Percentage Peak Cut-Off:",
                                      min = 0, max = 1, value = 0.6, step = 0.01),
                          numericInput("Trim_value", "Number of tail points to remove",
                                       value = 100))),
        tabPanel("Noise & Fibrillation", 
                 fluidRow(sliderInput("noise_slider", "Noise Amplitude Cutoff:",
                                               min = 0, max = 5000, value = 900, step = 10),
                          sliderInput("fib_slider", "Fibrillation Peak to Peak detection:",
                                      min = 0, max = 10, value = 1.2, step = 0.1))),
                        
        tabPanel("Arrhythmia",
                 fluidRow(sliderInput("arrh_slider", "Flag Troughs Above This Threshold:",
                                      min = 0, max = 1, value = 0.15, step = 0.01),
                          sliderInput("arrh_slider2", "Trough SD Percent Fluctuation:",
                                      min = 0, max = 100, value = 10, step = 0.5))),
        tabPanel("Warnings",
                 fluidRow(sliderInput("warning_slider1", "Rhythm Warning: SD distance between peaks:",
                                      min = 0, max = 10, value = 1.2, step = 0.1),
                          sliderInput("warning_slider2", "Trace Wandering Warning:",
                                      min = 0, max = 1, value = 0.1, step = 0.01)))
      ))),
      actionButton("run_analysis", "Apply Parameters"),
      
     br(),
      
      conditionalPanel(
        condition = "output.csvSelected",
        shinyWidgets::progressBar(id = "progress", value = 0, display_pct = TRUE)
      )
    ),
    mainPanel(
      uiOutput("tb"),
    )
  )
)
