server <- function(input, output, session) {
  

  
  source("Find_peaks_V2.0.R")
  
  options(shiny.maxRequestSize = 1024*1024*1024)
  
  well <- reactiveVal(NULL)
  peaksReactive <- reactiveVal(NULL)
  pstatsReactive <- reactiveVal(NULL)
  xminReac <- reactiveVal(NULL)
  xmaxReac <- reactiveVal(NULL)
  
  read_batch_with_progress <- function(file_path,nrows,no_batches){
    progress <- Progress$new(session, min = 1,max = no_batches)
    progress$set(message = "Importing Datasets...")
    seq_length <- ceiling(seq.int(from = 2, to = nrows-2,length.out = no_batches+1))
    seq_length <- seq_length[-length(seq_length)]
    
    #read the first line
    df <- read.csv(file_path,nrows = 1, header = T, check.names = F)
    col_names <- colnames(df)
    
    for(i in seq_along(seq_length)){
      progress$set(value = i)
      if(i == no_batches) chunk_size = -1 else chunk_size = seq_length[i+1] - seq_length[i]
      df_temp <- fread(file_path, skip = seq_length[i], nrows = chunk_size,header = F, check.names = F)
      colnames(df_temp) <- col_names
      df <- rbind(df,df_temp)
    }
    
    progress$close()
    return(df)
  }
  
  df <- reactive({
    req(input$CSV1)
    n_rows <- length(count.fields(input$CSV1$datapath))
    df_out <- read_batch_with_progress(input$CSV1$datapath,n_rows,10)
    return(df_out)
  })
  
  observe({
    output$sum <- renderPrint({
      print(head(df(), 10))
    }) 
  }) 
  
  observe({
    if (!is.null(df())) {
      updateSelectInput(session, "search", choices = df()[[1]])
    }
  })
  
  ## Generate Table for Data Set~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #output$table <- renderDataTable({
  #  df()
  #})
  
  output$tb <- renderUI({
    if (is.null(input$CSV1)) {
      return()
    } else {
      tabsetPanel(
        tabPanel("Peaks", tableOutput("peaks_table")),
        tabPanel("P.Stats", tableOutput("pstats_table")),
        tabPanel("Trace", fluidRow(plotOutput("graph"),
                                   fluidRow(column(4, tableOutput("pstats_table_short")),
                                            column(4, uiOutput("xslider")),
                                            column(4, uiOutput("graphcheckbox")))
                                  )
                 ),
         tabPanel("Full Analysis",
             attendantBar("progress-bar"),
             actionButton("perform_analysis", "Perform Full Analysis of CSV"),
             tableOutput("count_results"),
             fileInput("platePlan", "Import Plate Plan", multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")),
             textInput("result_name", "Enter Name for Results", value = "Results"),
             actionButton("export_count_results", "Export Results")),
        tabPanel("Plate View",
                 plotOutput("plateLayout"),
                 fluidRow(column(4, p(strong("Score Legend")),
                                    p(code("Noise - Black")),
                                    p(code("No Arrhythmia - Blue")),
                                    p(code("Fibrilation - Yellow")),
                                    p(code("Arrhythmia - Red")),
                                    p(code("Other - Pink"))),
                          column(4, uiOutput("plateRadio")),
                          column(4, selectInput("plateSize", "Select Plate Size",
                                                choices = c("96", "384"), selected = "384"))),
                 fluidRow(column(4, p(strong("Warning Legend")),
                                    p(code("No Warning - White")),
                                    p(code("Warning - Pink")),
                                    p(code("Rhythm Warning - Gold")),
                                    p(code("Wander Warning - Blue")),
                                    p(code("Multiple Warnings - Red"))))
      
      ))
    }
  })

  ## Select Well Math LONG....~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  observeEvent(input$run_analysis, {
    if (!is.null(input$search)) {
    if (input$search %in% df()[, 1]) {  
      well_data <- df()[df()[, 1] == input$search, ]
      well_data <- well_data %>% pivot_longer(2:length(well_data), names_to = "Time", values_to = ("OD"))
      
      TrimVal <- input$Trim_value
      well_data <- well_data %>% group_by(Well) %>%
        slice(1:(n() - TrimVal)) %>%
        ungroup()
      
      well_data$OD <- as.numeric(well_data$OD)
      well_data$Time <- as.numeric(well_data$Time)
      well_data$Smooth <- NA
      
      ## Smoothing~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      smooth_val <- input$smooth_slider
      well_data <- well_data %>%
          mutate(Smooth = suppressWarnings(predict(loess(OD ~ Time, span = smooth_val))))
      well(well_data)
      
      ##Slider values inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      peak_val <- input$peak_slider
      trough_val <- input$trough_slider
      noise_val <- input$noise_slider
      arrh_val <- input$arrh_slider
      arrh_val2 <- input$arrh_slider2
      pref_val <- input$pref_slider
      fib_val <- input$fib_slider
      warning_val <- input$warning_slider1
      wander_val <- input$warning_slider2
      
      ## Finding Peaks~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      t0 <- min(well_data$Time, na.rm = TRUE)
      tn <- max(well_data$Time, na.rm = TRUE)
      
      xminReac(t0)
      xmaxReac(tn)
      
      pks <- find_peaks(well_data, "Smooth", top = peak_val, bot = trough_val, pref = pref_val)
      pks <- pks %>% group_by(Well, shape) %>% slice(-c(1,n()))## removes first and last peak/troth
      peaksReactive(pks)
      
      ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ## Peak Math and Foundation for arrhtyhmia detection
      tTot <- (tn - t0)
      Peaks <- pks %>% filter(shape == "max") %>% group_by(Well) %>% arrange(Well, Time) %>% ## generate poincaré plot stats
        mutate(RR = (Time - lag(Time, default = NA))/1000, 
               RRn = (lead(Time, default = NA) - Time)/1000, 
               dRR = RR - RRn,
               Freq = (n() / (tTot/1000)) * 60)
      Peaks <- Peaks %>% ungroup() %>% 
        mutate(Rhythm = ifelse(RR > RRn*1.1, 0,
                               ifelse(RR < RRn*0.9, 0, 1)), 
               RRnmin = RRn*0.9, 
               Rnmax = 1.1*RRn)
      
      PeaksCSV1 <- pks
      
      p.stats <- Peaks %>%  #generate ellipse statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        group_by(Well) %>% summarise(meanRR = mean(RR, na.rm = TRUE), 
                                     medianRR = median(RR, na.rm = TRUE), 
                                     meandRR = mean(dRR, na.rm = TRUE),
                                     SD_length = sd(RRn, na.rm = TRUE), 
                                     SD_width = sd(dRR, na.rm = TRUE)/sqrt(2),
                                     Freq = mean(Freq, na.rm = TRUE),
                                     BVR = sum(abs(dRR), na.rm = TRUE)/length(dRR)*sqrt(2))
      p.stats$SD_length <- as.numeric(p.stats$SD_length)
      p.stats$SD_width <- as.numeric(p.stats$SD_width)
      p.stats$meanRR <- as.numeric(p.stats$meanRR)
      
      p.stats_ratio_ellipse <- p.stats %>%
        na.omit() %>%
        mutate(
          ratio = mean(c(SD_length, SD_width), na.rm = TRUE) / meanRR,
          ellipse = SD_length / SD_width
        )
      
      pks_abnormality <- pks %>%
        group_by(Well) %>%
        filter(shape == "min") %>%
        summarise(
          minTroughAmp = min(Smooth),
          maxTroughAmp = max(Smooth),
          diffTroughAmp = maxTroughAmp - minTroughAmp,
          sdTroughAmp = sd(Smooth)
        )
      
      pks_fib <- pks %>%
        group_by(Well) %>%
        arrange(desc(Smooth)) %>%
        summarise(
          lowSmooth = mean(tail(Smooth, 10), na.rm = T),
          hiSmooth = mean(head(Smooth, 10), na.rm = T),
          Amplitude = hiSmooth - lowSmooth
        )
      
      pks_max <- pks %>%
        group_by(Well) %>%
        filter(shape == "max") %>%
        arrange(desc(Time)) %>%
        summarise(
          Wander = mean(head(Smooth, 5), na.rm = T) - mean(tail(Smooth, 5), na.rm = T)
        )
      
      p.stats_essemble <- p.stats_ratio_ellipse %>%
        left_join(pks_abnormality, by = "Well") %>%
        left_join(pks_fib, by = "Well") %>%
        left_join(pks_max, by = "Well")
      
      p.stats_userP <- p.stats_essemble %>%
        group_by(Well) %>%
        summarise(
          userAmpDetectP = Amplitude*arrh_val + lowSmooth,
          FibArrhythDetect = Amplitude*(arrh_val*2) + lowSmooth,
          sdPercentDetect = sdTroughAmp/Amplitude*100
        )
      
      p.stats_final <- p.stats_essemble %>%
        left_join(p.stats_userP)
      
      ArrhythmiaDetection_NEW <- p.stats_final %>% 
        mutate(Noise = ifelse(Amplitude < noise_val, TRUE, FALSE),
                                                          
        Arrhythmia = ifelse((maxTroughAmp > userAmpDetectP &
                               Noise == FALSE) |
                              
                            (sdPercentDetect > arrh_val2 &
                               Noise == FALSE), TRUE, FALSE),
        
        Fib = ifelse((Amplitude > (noise_val*1.5) &
                        meanRR < fib_val &
                        Arrhythmia == FALSE &
                        Noise == FALSE) |
                       
                       (Amplitude > (noise_val*1.5) &
                           meanRR < fib_val &
                           maxTroughAmp < FibArrhythDetect &
                           Arrhythmia == TRUE &
                           Noise == FALSE), TRUE, FALSE),
        
        Arrhythmia = ifelse((maxTroughAmp > FibArrhythDetect & Fib == TRUE & Noise == FALSE) |
                              maxTroughAmp > userAmpDetectP & Noise == FALSE |
                              sdPercentDetect > arrh_val2 & Noise == FALSE, TRUE, FALSE),
        
        Warning = ifelse((Amplitude < (noise_val*1.5) &
                            Noise == FALSE) |
                           
                           (meanRR < 1 & Amplitude < (noise_val*2) &
                              Noise == FALSE), TRUE, FALSE),
        
        Rhythm.Warning = ifelse((SD_length > warning_val &
                                   Noise == FALSE &
                                   Arrhythmia == FALSE &
                                   Fib == FALSE), TRUE, FALSE),
        
        Wander.Warning = ifelse((abs(Wander) > (Amplitude*wander_val) &
                                   Noise == FALSE), TRUE, FALSE)
        )
      
      p.stats_t <- suppressWarnings(as_tibble(cbind(nms = names(ArrhythmiaDetection_NEW), t(ArrhythmiaDetection_NEW))))
      pstatsReactive(p.stats_t)
      
    } else {
      showNotification("Value not found. Please enter a valid search value.", type = "warning")
    }
    }
    
  })
  
  output$search_input <- renderUI({
    textInput("search", label = "Search by Value in the First Column")
  })
  
  
  ## Target Well Tab Data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$well_table <- renderTable({
    if (!is.null(well())) {
      well()
    }
  })
  
  ## Peaks Tab Data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$peaks_table <- renderTable({
    if (!is.null(peaksReactive())) {
      peaksReactive()
    }
  })
  
  ## Pstats Tab Data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$pstats_table <- renderTable({
    if (!is.null(pstatsReactive())) {
      pstatsReactive()
    }
  })
  
  output$pstats_table_short <- renderTable({
    if (!is.null(pstatsReactive())) {
      pstatsReactive()
      pstats_data <- pstatsReactive()
      extracted_rows <- tail(pstats_data, 6)
      extracted_rows
    }
  })
  
  observe({
    if (!is.null(well())) {
      
      output$graphcheckbox <- renderUI({
      checkboxGroupInput("show_noise_val", "Show Parameters", c("% Peak Cut Off", "Noise", "Arrhythmia Threshold"))
      })
    }
  })
  
  
  ## Graph X axis slider
  
  observe({
    if (!is.null(well())) {
      
      t0 <- as.numeric(xminReac())/1000
      tn <- as.numeric(xmaxReac())/1000
      
      output$xslider <- renderUI({
        sliderInput("xslider", "X Axis Range", 
                    min = t0, max = tn, value = c(t0, tn))
      })
    }
  })
  
  ## Graph Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$graph <- renderPlot({
    if (!is.null(well())) {
      

        baseline <- as.numeric(pstatsReactive() %>% slice(15) %>% pull(2))
        amp_value <- baseline + as.numeric(pstatsReactive() %>% slice(17) %>% pull(2))
        noise_val <- baseline + as.numeric(input$noise_slider)
        peakcut <- baseline + ((amp_value - baseline) * as.numeric(input$pref_slider))
        arrhythmiaCut <- baseline + ((amp_value - baseline) * as.numeric(input$arrh_slider))
        
        well_adjusted <- well() %>% mutate(Seconds = Time/1000)
        peaks_adjusted <- peaksReactive() %>% mutate(Seconds = Time/1000)
        Time0 <- input$xslider[1]
        Timemax <- input$xslider[2]
      
      
      
     p <- ggplot(well_adjusted, aes(x = Seconds, y = Smooth)) + 
        theme_classic() +
        geom_point(size = 1, colour = "black") +
        geom_line(size = 0.7, colour = "black") + 
        geom_point(data = peaks_adjusted, aes(x = Seconds, y = Smooth, colour = shape), size = 3) +
        geom_hline(yintercept = baseline, color = "red1", size = 1) +
        geom_hline(yintercept = amp_value, color = "blue", size = 1) +
        xlim(Time0, Timemax)
      
     if ("% Peak Cut Off" %in% input$show_noise_val) {
       p <- p + geom_hline(yintercept = peakcut, linetype = "dashed", color = "brown1", size = 1)
     }
     
     if ("Noise" %in% input$show_noise_val) {
       p <- p + geom_hline(yintercept = noise_val, linetype = "dashed", color = "deeppink1", size = 1)
     }
     
     if ("Arrhythmia Threshold" %in% input$show_noise_val) {
       p <- p + geom_hline(yintercept = arrhythmiaCut, linetype = "dashed", color = "skyblue", size = 1)
     }
     
     print(p)
     
    }
  })
  
  
  ## Full CSV Analysis~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  observeEvent(input$perform_analysis, {
    
    att <- Attendant$new("progress-bar")
    

    
        All_data <- df()
        counter <- 100 / (nrow(All_data) *2)
        
        All_data <- All_data %>% pivot_longer(2:length(All_data), names_to = "Time", values_to = ("OD"))
        
        TrimVal <- input$Trim_value
        All_data <- All_data %>% group_by(Well) %>%
          slice(1:(n() - TrimVal)) %>%
          ungroup()
        
        All_data$OD <- as.numeric(All_data$OD)
        All_data$Time <- as.numeric(All_data$Time)
        All_data$Smooth <- NA
        
        ## Smoothing~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        smooth_val <- input$smooth_slider
        
        for (well_name in unique(All_data$Well)) {
          
          subset_data1 <- All_data[All_data$Well == well_name, ]
          subset_data1 <- subset_data1 %>%
            mutate(Smooth = suppressWarnings(predict(loess(OD ~ Time, span = smooth_val))))
          
          All_data[All_data$Well == well_name, ]$Smooth <- subset_data1$Smooth
          
          att$inc(counter, text = "Applying Smoothing...")
        }

        ##Slider values inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        peak_val <- input$peak_slider
        trough_val <- input$trough_slider
        noise_val <- input$noise_slider
        arrh_val <- input$arrh_slider
        arrh_val2 <- input$arrh_slider2
        pref_val <- input$pref_slider
        fib_val <- input$fib_slider
        warning_val <- input$warning_slider1
        wander_val <- input$warning_slider2
        
        ## Finding Peaks~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        t0 <- min(All_data$Time, na.rm = TRUE)
        tn <- max(All_data$Time, na.rm = TRUE)
        pks1 <- data.frame()
        subset_data1 <- data.frame()
        
        for (well_name in unique(All_data$Well)) {
          
          subset_data1 <- All_data[All_data$Well == well_name, ]
          subset_data1 <- find_peaks(subset_data1, "Smooth", top = peak_val, bot = trough_val, pref = pref_val)
          
          pks1 <- rbind(pks1, subset_data1)
          att$inc(counter, text = "Finding Peaks & Searching for Arrhythmias...")
        }

        pks1 <- pks1 %>% group_by(Well, shape) %>% slice(-c(1,n()))## removes first and last peak/troth

        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ## Peak Math and Foundation for arrhtyhmia detection
        
        tTot <- (tn - t0)
        
        Peaks <- pks1 %>% filter(shape == "max") %>% group_by(Well) %>% arrange(Well, Time) %>% ## generate poincaré plot stats
          mutate(RR = (Time - lag(Time, default = NA))/1000, 
                 RRn = (lead(Time, default = NA) - Time)/1000, 
                 dRR = RR - RRn,
                 Freq = (n() / (tTot/1000)) * 60)
        Peaks <- Peaks %>% ungroup() %>% 
          mutate(Rhythm = ifelse(RR > RRn*1.1, 0,
                                 ifelse(RR < RRn*0.9, 0, 1)), 
                 RRnmin = RRn*0.9, 
                 Rnmax = 1.1*RRn)
        
        p.stats <- Peaks %>%  #generate ellipse statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          group_by(Well) %>% summarise(meanRR = mean(RR, na.rm = TRUE), 
                                       medianRR = median(RR, na.rm = TRUE), 
                                       meandRR = mean(dRR, na.rm = TRUE),
                                       SD_length = sd(RRn, na.rm = TRUE), 
                                       SD_width = sd(dRR, na.rm = TRUE)/sqrt(2),
                                       Freq = mean(Freq, na.rm = TRUE),
                                       BVR = sum(abs(dRR), na.rm = TRUE)/length(dRR)*sqrt(2))
        p.stats$SD_length <- as.numeric(p.stats$SD_length)
        p.stats$SD_width <- as.numeric(p.stats$SD_width)
        p.stats$meanRR <- as.numeric(p.stats$meanRR)
        
        p.stats_ratio_ellipse <- p.stats %>%
          na.omit() %>%
          mutate(
            ratio = mean(c(SD_length, SD_width), na.rm = TRUE) / meanRR,
            ellipse = SD_length / SD_width
          )
        
        pks_abnormality <- pks1 %>%
          group_by(Well) %>%
          filter(shape == "min") %>%
          summarise(
            minTroughAmp = min(Smooth),
            maxTroughAmp = max(Smooth),
            diffTroughAmp = maxTroughAmp - minTroughAmp,
            sdTroughAmp = sd(Smooth)
          )
        
        pks_fib <- pks1 %>%
          group_by(Well) %>%
          arrange(desc(Smooth)) %>%
          summarise(
            lowSmooth = mean(tail(Smooth, 10), na.rm = T),
            hiSmooth = mean(head(Smooth, 10), na.rm = T),
            Amplitude = hiSmooth - lowSmooth
          )
        
        pks_max <- pks1 %>%
          group_by(Well) %>%
          filter(shape == "max") %>%
          arrange(desc(Time)) %>%
          summarise(
            Wander = mean(head(Smooth, 5), na.rm = T) - mean(tail(Smooth, 5), na.rm = T)
          )
        
        p.stats_essemble <- p.stats_ratio_ellipse %>%
          left_join(pks_abnormality, by = "Well") %>%
          left_join(pks_fib, by = "Well") %>%
          left_join(pks_max, by = "Well")
        
        p.stats_userP <- p.stats_essemble %>%
          group_by(Well) %>%
          summarise(
            userAmpDetectP = Amplitude*arrh_val + lowSmooth,
            FibArrhythDetect = Amplitude*(arrh_val*2) + lowSmooth,
            sdPercentDetect = sdTroughAmp/Amplitude*100
          )
        
        p.stats_final <- p.stats_essemble %>%
          left_join(p.stats_userP)
        
        ArrhythmiaDetection_NEW <- p.stats_final %>% 
          mutate(Noise = ifelse(Amplitude < noise_val, TRUE, FALSE),
                 
                 Arrhythmia = ifelse((maxTroughAmp > userAmpDetectP &
                                        Noise == FALSE) |
                                       
                                       (sdPercentDetect > arrh_val2 &
                                          Noise == FALSE), TRUE, FALSE),
                 
                 Fib = ifelse((Amplitude > (noise_val*1.5) &
                                 meanRR < fib_val &
                                 Arrhythmia == FALSE &
                                 Noise == FALSE) |
                                
                                (Amplitude > (noise_val*1.5) &
                                   meanRR < fib_val &
                                   maxTroughAmp < FibArrhythDetect &
                                   Arrhythmia == TRUE &
                                   Noise == FALSE), TRUE, FALSE),
                 
                 Arrhythmia = ifelse((maxTroughAmp > FibArrhythDetect & Fib == TRUE & Noise == FALSE) |
                                       maxTroughAmp > userAmpDetectP & Noise == FALSE |
                                       sdPercentDetect > arrh_val2 & Noise == FALSE, TRUE, FALSE),
                 
                 Warning = ifelse((Amplitude < (noise_val*1.5) &
                                     Noise == FALSE) |
                                    
                                    (meanRR < 1 & Amplitude < (noise_val*2) &
                                       Noise == FALSE), TRUE, FALSE),
                 
                 Rhythm.Warning = ifelse((SD_length > warning_val &
                                            Noise == FALSE &
                                            Arrhythmia == FALSE &
                                            Fib == FALSE), TRUE, FALSE),
                 
                 Wander.Warning = ifelse((abs(Wander) > (Amplitude*wander_val) &
                                            Noise == FALSE), TRUE, FALSE)
          )
        
        att$set(100, text = "Analysis Complete!")
        
        count_results <- ArrhythmiaDetection_NEW %>% 
          summarise(
            Arrhythmias_Detected = sum(Arrhythmia == T),
            No_Arrhythmia = sum(Arrhythmia == F & Noise == F),
            Noise = sum(Noise == T),
            Fibrilation = sum(Fib == T),
            Warning = sum(Warning == T),
            Rhythm_Warning = sum(Rhythm.Warning == T),
            Wandering_Trace = sum(Wander.Warning == T))
        
        output$count_results <- renderTable({
          count_results
        })
        
        mergedData <- reactiveVal()
        resultsText <- reactiveVal()
        plateName <- reactiveVal()
        
        thedata <- reactive(ArrhythmiaDetection_NEW)

        # Importing Plate Plan
        observeEvent(input$platePlan, {
          req(input$platePlan)
          
          thePlate <- fread(input$platePlan$datapath)
          currentData <- thedata()
          mergedData <- merge(currentData, thePlate, by = "Well", all = TRUE)
          mergedData <- mergedData %>%
            select(names(thePlate), everything()) %>% bind_rows()
          mergedData(mergedData)
        })
        
        #Export Save file name
        observeEvent(input$result_name, {
          updateTextInput(session, "result_name", value = input$result_name)
          TextString <- paste0(input$result_name, ".csv")
          GraphSave <- paste0("_PlateLayout_", input$result_name)
          resultsText(TextString)
          plateName(GraphSave)
        })
        
        #Output CSV
        observeEvent(input$export_count_results, {
            req(!is.null(mergedData()) || !is.null(thedata()))
          
            if (!is.null(mergedData()) && nrow(mergedData()) > 0) {
              write.csv(mergedData(), resultsText())
            } else {
              write.csv(thedata(), resultsText())
            } 
          })
        
        # Plate Plan View with color
        output$plateRadio <- renderUI({
          
          if (!is.null(mergedData()) && nrow(mergedData()) > 0) {
            
            radioButtons("plateRadio", "Select Plate View",
                         choices = c("Score",
                                     "Warnings",
                                     "Amplitude",
                                     "Frequency",
                                     "Line",
                                     "Condition",
                                     "Treatment"))
            
          } else {
            radioButtons("plateRadio", "Select Plate View",
                         choices = c("Score",
                                     "Warnings",
                                     "Amplitude",
                                     "Frequency"))
          } 
          
        })
        
        output$plateLayout <- renderPlot({
          
          plateTitle <- as.character(plateName())
          plateTitle <- paste0(input$plateRadio, plateTitle)
          plateSize <- reactiveVal(input$plateSize)
          plateSize <- as.character(plateSize())
          
          if (!is.null(mergedData()) && nrow(mergedData()) > 0) {
            
            data_store <- data.frame(mergedData())
            data_store <- data_store %>% mutate(
              Score = ifelse(Noise == TRUE, 0,
                             ifelse(Noise == FALSE & Fib == FALSE & Arrhythmia == FALSE, 1,
                                    ifelse(Fib == TRUE, 2,
                                           ifelse(Arrhythmia == TRUE & Fib == FALSE, 3, 4)))),
              Warnings = ifelse(Warning == FALSE & Rhythm.Warning == FALSE & Wander.Warning == FALSE, 0,
                               ifelse(Warning == TRUE, 1,
                                      ifelse(Rhythm.Warning == TRUE, 2,
                                             ifelse(Wander.Warning == TRUE, 3, 4))))
            )
            if ("Score" %in% input$plateRadio) {
              
              plate_plot(
                data = data_store,
                position = Well,
                value = Score,
                plate_size = plateSize,
                title = plateTitle,
                show_legend = FALSE,
                limits = c(0,4),
                colour = c("gray0",
                           "dodgerblue",
                           "gold",
                           "red",
                           "deeppink")
              )
            } else {
              
              if ("Warnings" %in% input$plateRadio) {
                
                plate_plot(
                  data = data_store,
                  position = Well,
                  value = Warnings,
                  title = plateTitle,
                  plate_size = plateSize,
                  show_legend = FALSE,
                  limits = c(0,4),
                  colour = c("white",
                             "deeppink",
                             "gold",
                             "dodgerblue",
                             "red")
                )
                
              } else {
                
                if ("Amplitude" %in% input$plateRadio) {
                  
                  plate_plot(
                    data = data_store,
                    position = Well,
                    value = Amplitude,
                    title = plateTitle,
                    plate_size = plateSize
                  )
                  
                } else {
                  
                  if ("Frequency" %in% input$plateRadio) {
                    
                    plate_plot(
                      data = data_store,
                      position = Well,
                      value = Freq,
                      title = plateTitle,
                      plate_size = plateSize
                    )
                    
                  } else {
                    
                    if ("Line" %in% input$plateRadio) {
                      
                      plate_plot(
                        data = data_store,
                        position = Well,
                        value = Line,
                        title = plateTitle,
                        legend_n_row = 6,
                        plate_size = plateSize
                      )
                      
                    } else {
                      
                      if ("Condition" %in% input$plateRadio) {
                        
                        plate_plot(
                          data = data_store,
                          position = Well,
                          value = Condition,
                          title = plateTitle,
                          legend_n_row = 6,
                          plate_size = plateSize
                        )
                        
                      } else {
                        
                        if ("Treatment" %in% input$plateRadio) {
                          
                          plate_plot(
                            data = data_store,
                            position = Well,
                            value = Treatment,
                            title = plateTitle,
                            legend_n_row = 6,
                            plate_size = plateSize
                          )
                        } else {
                          
                          print("That Column Can Not Be Found")
                          
                          }
                  }
              } 
            }
                }}}
          
          } else {
            
            data_store <- data.frame(thedata())
            data_store <- data_store %>% mutate(
              Score = ifelse(Noise == TRUE, 0,
                             ifelse(Noise == FALSE & Fib == FALSE & Arrhythmia == FALSE, 1,
                                    ifelse(Fib == TRUE, 2,
                                           ifelse(Arrhythmia == TRUE & Fib == FALSE, 3, 4)))),
              Warnings = ifelse(Warning == FALSE & Rhythm.Warning == FALSE & Wander.Warning == FALSE, 0,
                               ifelse(Warning == TRUE, 1,
                                      ifelse(Rhythm.Warning == TRUE, 2,
                                             ifelse(Wander.Warning == TRUE, 3, 4))))
            )
            
            if ("Score" %in% input$plateRadio) {
              
              plate_plot(
                data = data_store,
                position = Well,
                value = Score,
                plate_size = plateSize,
                title = plateTitle,
                show_legend = FALSE,
                limits = c(0,4),
                colour = c("gray0",
                           "dodgerblue",
                           "gold",
                           "red",
                           "deeppink")
              )
            } else {
              
              if ("Warnings" %in% input$plateRadio) {
              
              plate_plot(
                data = data_store,
                position = Well,
                value = Warnings,
                title = plateTitle,
                show_legend = FALSE,
                plate_size = plateSize,
                limits = c(0,4),
                colour = c("white",
                           "deeppink",
                           "gold",
                           "dodgerblue",
                           "red")
              )
              
            } else {
              
              if ("Amplitude" %in% input$plateRadio) {
                
                plate_plot(
                  data = data_store,
                  position = Well,
                  value = Amplitude,
                  title = plateTitle,
                  plate_size = plateSize
                )
                
              } else {
                
              plate_plot(
                data = data_store,
                position = Well,
                value = Freq,
                title = plateTitle,
                plate_size = plateSize
              )}
          } 
          }
          }
     })
  })
}
  
