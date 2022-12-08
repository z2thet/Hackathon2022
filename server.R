# AE app 

# Load packages ####
library(shiny)
library(shinycssloaders)
library(tidyverse)
library(rio)
library(janitor)
library(DT)
library(skimr)
library(rlang)
library(ggpubr)
library(ggExtra)
library(ggsci)
library(moments)
library(gridExtra)
library(ggnewscale)
library(lubridate)
library(survival)
library(survminer)

 

# Define server logic required ####
shinyServer(function(input, output, session) {
  # ~~~~~ UPLOAD DATA ####
    # Upload AE data  
    AE_data_upload <- selected_data_upload <- eventReactive(input$data_file, {
        if (is.null(input$data_file)) return(NULL)
        rio::import(input$data_file$datapath) %>% 
            # if integer, conver it to numeric
            mutate_if(is.integer, as.numeric) %>% 
            # if any numeric variable has less than or equal to 5 unique categories,
            # convert it to factor
            mutate_if(~ is.numeric(.) & n_distinct(.) <= 5, as.factor) #%>% select(-any_of(Initials))
    })
    # Upload DEMOGRAPHICS data  
    demographics_upload <- eventReactive(input$demographics_file, {
      if (is.null(input$demographics_file)) return(NULL)
      rio::import(input$demographics_file$datapath) %>% 
        # if integer, conver it to numeric
        mutate_if(is.integer, as.numeric) %>% 
        # if any numeric variable has less than or equal to 5 unique categories,
        # convert it to factor
        mutate_if(~ is.numeric(.) & n_distinct(.) <= 5, as.factor)  #%>% select(-any_of(Initials))
    })
    # Upload follow up  data  
    fu_data_upload0  <- eventReactive(input$fu_file, {
      if (is.null(input$fu_file)) return(NULL)
      rio::import(input$fu_file$datapath) %>% 
        # if integer, conver it to numeric
        mutate_if(is.integer, as.numeric) %>% 
        # if any numeric variable has less than or equal to 5 unique categories,
        # convert it to factor
        mutate_if(~ is.numeric(.) & n_distinct(.) <= 5, as.factor) # %>% select(-any_of(Initials))
    })
    # Upload drug administration data  
    da_data_upload0  <- eventReactive(input$da_file, {
      if (is.null(input$da_file)) return(NULL)
      rio::import(input$da_file$datapath) %>% 
        # if integer, conver it to numeric
        mutate_if(is.integer, as.numeric) %>% 
        # if any numeric variable has less than or equal to 5 unique categories,
        # convert it to factor
        mutate_if(~ is.numeric(.) & n_distinct(.) <= 5, as.factor) # %>% select(-any_of(Initials))
    })
    
  # ~~~~~ IMPORT DEM0 DATA ~~~~~~ ####
    # data sets for displaying data before user upload ###
    demo_data <- reactive({
        #rio::import("demo_data.xls") %>% 
      rio::import("demo_ae_data.csv") %>% 
            # if integer, conver it to numeric
            mutate_if(is.integer, as.numeric) %>% 
            # if any numeric variable has less than or equal to 5 unique categories,
            # convert it to factor
            mutate_if(~ is.numeric(.) & n_distinct(.) <= 5, as.factor)
    })
    demo_data2 <- reactive({
      #rio::import("demo_data2.xls") %>% 
        rio::import("demo_demo_data.csv") %>%
        # if integer, conver it to numeric
        mutate_if(is.integer, as.numeric) %>% 
        # if any numeric variable has less than or equal to 5 unique categories,
        # convert it to factor
        mutate_if(~ is.numeric(.) & n_distinct(.) <= 5, as.factor)
    })
    demo_fu_data <- reactive({
      rio::import("demo_fu_data.csv") %>% 
        # if integer, conver it to numeric
        mutate_if(is.integer, as.numeric) %>% 
        # if any numeric variable has less than or equal to 5 unique categories,
        # convert it to factor
        mutate_if(~ is.numeric(.) & n_distinct(.) <= 5, as.factor)
    })
    demo_da_data <- reactive({
      rio::import("demo_da_data.csv") %>% 
        # if integer, conver it to numeric
        mutate_if(is.integer, as.numeric) %>% 
        # if any numeric variable has less than or equal to 5 unique categories,
        # convert it to factor
        mutate_if(~ is.numeric(.) & n_distinct(.) <= 5, as.factor)
    })
    # End demo data import ###

  # ~~~~~ USE DEMO DATA ~~~~~ ####
    # if demographics data not uploaded filter out missing sequence no ###
    AE_data <- selected_data <- reactive({ if (is.null(input$data_file)) {
      
      output <- demo_data()  %>% clean_names() 
      
      if ( "cdus_ctcae_toxicity_type_code" %in% names(output)){output <- output %>% 
          mutate(cdus_toxicity_type_code = cdus_ctcae_toxicity_type_code)}
      
      } else {
        output <- AE_data_upload() %>% clean_names() %>% filter(is.na(sequence_no) == FALSE)
      
      if ( "cdus_ctcae_toxicity_type_code" %in% names(output)){output <- output %>% 
          mutate(cdus_toxicity_type_code = cdus_ctcae_toxicity_type_code)}
      
      }
      
      output
      })
    
    demographics_data <- reactive({  if (is.null(input$demographics_file)) {demo_data2()  %>% clean_names() } else { demographics_upload() %>% clean_names() %>% 
                                                                                                                     filter(is.na(sequence_no) == FALSE) }  })
    fu_data_upload <- reactive({  if (is.null(input$fu_file)) {demo_fu_data()  %>% clean_names() } else { fu_data_upload0() %>% clean_names() %>% 
        filter(is.na(sequence_no) == FALSE) }  })
    da_data_upload <- reactive({  if (is.null(input$da_file)) {demo_da_data()  %>% clean_names() } else { da_data_upload0() %>% clean_names() %>% 
        filter(is.na(sequence_no) == FALSE) }  })
    

  # ~~~~~ Make toxicity data by merging AE data and Demographics #####
    toxicity_data <- reactive({  
      
      #if (is.null(input$fu_file)) { } else {
      AE_data <-  AE_data() %>% clean_names() %>%  filter(is.na(sequence_no) == FALSE)  
      demographics_data <-  demographics_data() %>% clean_names() %>%   filter(is.na(sequence_no) == FALSE) 
      
      # print("dims AE and Demos:")
      # print(dim(AE_data))
      # print(dim(demographics_data))
      # 
      if("onset_date" %in% names(AE_data)){AE_data$onset_date_of_ae <- AE_data$onset_date }
      if("cdus_toxicity_type_code" %in% names(AE_data)){AE_data$cdus_ctcae_toxicity_type_code <- AE_data$cdus_toxicity_type_code }
      if("toxicity" %in% names(AE_data)){AE_data$cdus_ctcae_toxicity_type_code <- AE_data$toxicity } 
      if(!"attribution_possible" %in% names(AE_data)){AE_data <- AE_data %>% mutate(attribution_possible = case_when(attribution == "Possible" ~ "Yes" , 
                                                                                                                     attribution != "Possible" ~ "Not  Applicable" ),
                                                                                    attribution_probable = case_when(attribution == "Probable" ~ "Yes", 
                                                                                                                     attribution != "Probable" ~ "Not  Applicable" ),
                                                                                    attribution_definite = case_when(attribution == "Definite" ~ "Yes", 
                                                                                                                     attribution != "Definite" ~ "Not  Applicable" ))
      }
      
      
      
     # AE_data$sequence_no <- as.numeric(AE_data$sequence_no)
     # demographics_data$sequence_no <- as.numeric(demographics_data$sequence_no)
      
  
      AE_data <- AE_data %>% select(-c( form, form_desc))
      # print("dims of AE and DEmos:")
        # print(names(AE_data))
        # print(names(demographics_data))   
      toxicity.data <-  merge(AE_data %>%  clean_names() %>% filter(is.na(onset_date_of_ae) == FALSE  & is.na(resolved_date) == FALSE ),
                              demographics_data %>% clean_names() %>% filter(is.na(on_treatment_date)==FALSE) %>% select(sequence_no, on_treatment_date),
                              by = "sequence_no")
      
      # print("dims tox data:")
      # print(dim(toxicity.data))
      # print("names tox data:")
      # names(toxicity.data)
      
      #---calculate AE time---
      if( "start_date_of_course" %in% names(toxicity.data)){
      toxicity.data$start_date_of_course_cycle <- toxicity.data$start_date_of_course} else {
      toxicity.data$start_date_of_course_cycle <- toxicity.data$start_date
        
      }
      
      print(toxicity.data %>% select(onset_date_of_ae, on_treatment_date))
      print(str(toxicity.data %>% select(onset_date_of_ae, on_treatment_date),1,1))
      
      toxicity.data <- toxicity.data %>% mutate(pid = sequence_no, time = start_date_of_course_cycle,
                                                AE.time = difftime(onset_date_of_ae, on_treatment_date, units = "days"))#
      toxicity.data <- toxicity.data[!is.na(toxicity.data$time), ]

      toxicity.data
      #}
      
      }) # end reactive toxicity data code chunk stuff thing 
    
  # ~~~~~ render datatables TO DISPLAY ~~~~~~~~#####
    
    output$selected_data <- DT::renderDataTable({
        datatable(AE_data(), extensions = "FixedHeader", options = list(
            scrollX = TRUE, pageLength = 15, fixedHeader = FALSE, autoWidth = TRUE, searchHighlight = TRUE, 
            initComplete = JS(
                "function(settings, json) {",
                "$(this.api().table().header()).css({'background-color': '#58b0e3', 'color': '#fff'});",
                "}")), 
            rownames = FALSE)
    })
    
    output$demographics_data <- DT::renderDataTable({
      datatable(demographics_data(), extensions = "FixedHeader", options = list(
        scrollX = TRUE, pageLength = 15, fixedHeader = FALSE, autoWidth = TRUE, searchHighlight = TRUE, 
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#58b0e3', 'color': '#fff'});",
          "}")), 
        rownames = FALSE)
    })
    
    output$followup_data <- DT::renderDataTable({
      datatable(fu_data_upload(), extensions = "FixedHeader", options = list(
        scrollX = TRUE, pageLength = 15, fixedHeader = FALSE, autoWidth = TRUE, searchHighlight = TRUE, 
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#58b0e3', 'color': '#fff'});",
          "}")), 
        rownames = FALSE)
    })
    
    output$tox_data <- DT::renderDataTable({
      datatable(toxicity_data(), extensions = "FixedHeader", options = list(
        scrollX = TRUE, pageLength = 15, fixedHeader = FALSE, autoWidth = TRUE, searchHighlight = TRUE, 
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#58b0e3', 'color': '#fff'});",
          "}")), 
        rownames = FALSE)
    })
    
    output$drug_data <- DT::renderDataTable({
      datatable(da_data_upload(), extensions = "FixedHeader", options = list(
        scrollX = TRUE, pageLength = 15, fixedHeader = FALSE, autoWidth = TRUE, searchHighlight = TRUE, 
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#58b0e3', 'color': '#fff'});",
          "}")), 
        rownames = FALSE)
    })   
    

  # MERGE AE AND on_treatment date from DEMOGRAPHICS changing names in uploaded data if necessary ####
    AEandDemoData <- AE_data_test2 <- reactive({ 
      AE_data <-  AE_data()  %>% clean_names()
         demographics_data  <-  demographics_data()   %>% clean_names()
         
        if("onset_date" %in% names(AE_data)){AE_data$onset_date_of_ae <- AE_data$onset_date }
        if("cdus_toxicity_type_code" %in% names(AE_data)){AE_data$cdus_ctcae_toxicity_type_code <- AE_data$cdus_toxicity_type_code }
        if("toxicity" %in% names(AE_data)){AE_data$cdus_ctcae_toxicity_type_code <- AE_data$toxicity } 
        if(!"attribution_possible" %in% names(AE_data)){AE_data <- AE_data %>% mutate(attribution_possible = case_when(attribution == "Possible" ~ "Yes" , 
                                                                                                                                         attribution != "Possible" ~ "Not  Applicable" ),
                                                                                                        attribution_probable = case_when(attribution == "Probable" ~ "Yes", 
                                                                                                                                         attribution != "Probable" ~ "Not  Applicable" ),
                                                                                                        attribution_definite = case_when(attribution == "Definite" ~ "Yes", 
                                                                                                                                         attribution != "Definite" ~ "Not  Applicable" ))
        }
        
 
       #AE_data$sequence_no <- as.numeric(AE_data$sequence_no)
       #demographics_data$sequence_no <- as.numeric(demographics_data$sequence_no)
 
        merge(AE_data %>%  clean_names() %>% filter(is.na(onset_date_of_ae) == FALSE  & is.na(resolved_date) == FALSE ),
              demographics_data %>% clean_names() %>% filter(is.na(on_treatment_date)==FALSE) %>% select(sequence_no, on_treatment_date),
                                     by = "sequence_no")
    })
 
  # start: AE days tab vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv  ####
  # Calculate number of days of AEs and number of unique Days of AEs to display with render table ### 
    output$selected_data_AEDAYS <- DT::renderDataTable({
 

      a2 <- AEandDemoData()  %>%  filter(is.na(onset_date_of_ae) == FALSE  & is.na(resolved_date ) == FALSE & is.na(on_treatment_date ) == FALSE) %>% 
        dplyr::mutate(attribution_possible = replace_na(attribution_possible, 'Not  Applicable'),
                      attribution_probable = replace_na(attribution_probable, 'Not  Applicable'),
                      attribution_definite = replace_na(attribution_definite, 'Not  Applicable')) %>% 
        group_by(sequence_no) %>% 
        mutate( treatment_related = !apply(as.matrix(vars(name2)),1,function(x) all(x=='Not  Applicable')),
                treatment_related = factor(treatment_related,level=c(F,T),label=c('No','Yes')),
                t1=as.numeric(difftime(onset_date_of_ae,on_treatment_date,units='days')),
                t2=as.numeric(difftime(resolved_date,on_treatment_date,units='days')),
                t12=as.numeric(difftime(resolved_date,onset_date_of_ae,units='days')),
                code=cdus_ctcae_toxicity_type_code,
                index = as.numeric(factor(code)), 
                listdaynumber = map2(t1,t2,function(.x, .y){.x:.y}),
                totalnumberofdayswithAEs = length(unlist( listdaynumber)),
                numberofdayswithAEs = length(unique(unlist( listdaynumber)))
        )
      
      
      # Select only what we need to show ###
      a3 <- as.data.frame(a2 %>% select(sequence_no, totalnumberofdayswithAEs,  numberofdayswithAEs) %>%  
                                 distinct(sequence_no, .keep_all = TRUE))
      
       
      
      # datatable(a3, extensions = "FixedHeader", options = list(
      #   scrollX = TRUE, pageLength = 15, fixedHeader = FALSE, autoWidth = TRUE, searchHighlight = TRUE, 
      #   initComplete = JS(
      #     "function(settings, json) {",
      #     "$(this.api().table().header()).css({'background-color': '#58b0e3', 'color': '#fff'});",
      #     "}")), 
      #   rownames = FALSE)
      
      datatable(a3 , 
        rownames = FALSE)
    })
    
    AEdays <- reactive({  
      a2 <- AEandDemoData()  %>%  filter(is.na(onset_date_of_ae) == FALSE  & is.na(resolved_date ) == FALSE & is.na(on_treatment_date ) == FALSE) %>% 
        dplyr::mutate(attribution_possible = replace_na(attribution_possible, 'Not  Applicable'),
                      attribution_probable = replace_na(attribution_probable, 'Not  Applicable'),
                      attribution_definite = replace_na(attribution_definite, 'Not  Applicable')) %>% 
        group_by(sequence_no) %>% 
        mutate( treatment_related = !apply(as.matrix(vars(name2)),1,function(x) all(x=='Not  Applicable')),
                treatment_related = factor(treatment_related,level=c(F,T),label=c('No','Yes')),
                t1=as.numeric(difftime(onset_date_of_ae,on_treatment_date,units='days')),
                t2=as.numeric(difftime(resolved_date,on_treatment_date,units='days')),
                t12=as.numeric(difftime(resolved_date,onset_date_of_ae,units='days')),
                code=cdus_ctcae_toxicity_type_code,
                index = as.numeric(factor(code)), 
                listdaynumber = map2(t1,t2,function(.x, .y){.x:.y}),
                totalnumberofdayswithAEs = length(unlist( listdaynumber)),
                numberofdayswithAEs = length(unique(unlist( listdaynumber)))
        )
      
      
      # Select only what we need to show ###
      a3 <- as.data.frame(a2 %>% select(sequence_no, totalnumberofdayswithAEs,  numberofdayswithAEs) %>%  
                            distinct(sequence_no, .keep_all = TRUE))
      
      a3
       
    })
    
    # Download handler for AE data
    {
    output$AEDAYSdownload <- downloadHandler( 
      filename = function(){"AEdays.csv"}, 
      content = function(fname){
        write.csv(AEdays(), fname)
      }
    )
    }
    
  # end: AE days tab ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  ####
  
     
  # start: AE plot tab vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv  ####
  # ~~~~~~AE swimmers plot ~~~~~~~~~~~~~~~~~~ 
    observe({ updateSelectizeInput(session, "single_var_input", choices = unique(AEandDemoData()$sequence_no),  
                                   #selected = character(0), 
                                   selected = unique(AEandDemoData()$sequence_no)[1],
                                   server = TRUE) })
  
 
    # Custom theme for ggplot
    my_theme = theme(axis.title = element_text(size = 20), axis.text = element_text(size = 15), legend.text = element_text(size = 15), strip.text = element_text(size = 15))
    
  # AE LINE-swimmers PLOT ######
    single_var_plot <- reactive({
      
      
      w1<-AEandDemoData() %>% filter(is.na(onset_date_of_ae) == FALSE  & is.na(resolved_date) == FALSE  & is.na(on_treatment_date) == FALSE)
      # print("dim of w1 for plot 1")
      # print(dim(w1))
      subject<-input$single_var_input  
      
      w1 <- w1 %>% filter(sequence_no == subject)
      name1<-c('onset_date_of_ae','cdus_ctcae_toxicity_type_code', 'resolved_date','on_treatment_date' ,'grade')
      name2<-c('attribution_possible','attribution_probable', 'attribution_definite')
      w2<-w1%>%select(c(name1,name2))
      w2[,name2][is.na(w2[,name2])] <- 'Not  Applicable'
      w2$treatment_related<-!apply(as.matrix(w2[,name2]),1,function(x) all(x=='Not  Applicable'))
      w2$treatment_related<-factor(w2$treatment_related,level=c(F,T),label=c('No','Yes'))
      w2=w2%>%mutate(t1=as.numeric(difftime(onset_date_of_ae,on_treatment_date,units='days')),
                     t2=as.numeric(difftime(resolved_date,on_treatment_date,units='days')),
                     t12=as.numeric(difftime(resolved_date,onset_date_of_ae,units='days')),
                     code=cdus_ctcae_toxicity_type_code)
      w2$index<-as.numeric(factor(w2$code))
      w2$grade<-as.factor(w2$grade)
      
       
      
      if("start_date_of_drug" %in% names( da_data_upload())) {
        da_data_subject <- da_data_upload() %>% filter(sequence_no == subject) 
        }else{
        
        da_data_subject <- da_data_upload() %>% mutate(start_date_of_drug = start_date) %>% filter(sequence_no == subject)
 
        
      }
      
      if("drug" %in% names( da_data_subject)) {
        da_data_subject <-da_data_subject %>% filter(sequence_no == subject) %>% 
          select(sequence_no, start_date_of_drug, cycle, drug) %>% 
          filter(!drug %in% c("Not  Applicable")) %>% 
          filter(is.na(drug) == FALSE)
        
        
      }else{
        
        da_data_subject <- da_data_subject %>% mutate(drug = level) %>% filter(sequence_no == subject) %>% 
          select(sequence_no, start_date_of_drug, cycle, drug) %>%
          filter(!drug %in% c("Not  Applicable")) %>%
          filter(is.na(drug) == FALSE)

        
      }
      
      
      # dim(da_data_subject)
      # print(names(da_data_subject))
      
      w2.2 <- w2 %>% select(on_treatment_date) %>% distinct(on_treatment_date) 
      # print("some AE data w2")
      #  print(w2.2)
      #  print("dim AE data w2")
      #  print(dim(w2.2))
      #  
     # w2.2 <- w2 %>% select(on_treatment_date) %>% distinct(on_treatment_date) 
       
      w3 <- merge(  da_data_subject, w2.2) %>% mutate(dadaynumber = as.numeric(difftime(start_date_of_drug ,on_treatment_date,units='days')) )
      # print("  vvvvvvvvvvvvvvvvvv DA data w3 for plot 1")
      # print(w3)
      # print("dim DA data w3 for plot 1")
      # print(dim(w3)) 
      # print("  ^^^^^^^^^^^^^^^^^^ DA data w3 for plot 1")
      # 
      # 
      # print("AE DATA FOR PLOT 1 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
      # print(w2)
      # print("AE DATA FOR PLOT 1 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
      UL <- max(w2$t2)+10
      plot1 <- ggplot(w2, aes(t1, index, color=grade, shape=treatment_related, label = code)) +
        geom_point(aes(t2, index),size=2)+
        geom_segment(aes(xend = t2, yend = index), size = 1, lineend = "butt")+
        xlab('days')+ylab(paste0("Seq #: ", subject))+
        xlim(c(-4,UL)) +
         scale_y_continuous(breaks=w2$index,labels=w2$code)+
        theme(axis.text=element_text(size=18),
              axis.title=element_text(size=18,face="bold"),
              legend.text = element_text(size=14),
              legend.title = element_text(size=18)) + 
        ggnewscale::new_scale_color() +
        geom_vline(data = w3, aes(xintercept = dadaynumber, color = drug))  
      
    })
    
    output$single_var_plot <- renderPlot({   print(single_var_plot())    })
    
    
  # NEW AE swimmers plot EARLY cut off ####
    single_var_plot_early <- reactive({
      
      
      w1<-AEandDemoData() %>% filter(is.na(onset_date_of_ae) == FALSE  & is.na(resolved_date) == FALSE  & is.na(on_treatment_date) == FALSE)
      # print("names in w1")
      # print(names(w1))
      subject<-input$single_var_input  
      
      w1 <- w1 %>% filter(sequence_no == subject)
      name1<-c('onset_date_of_ae','cdus_ctcae_toxicity_type_code', 'resolved_date','on_treatment_date' ,'grade')
      name2<-c('attribution_possible','attribution_probable', 'attribution_definite')
      w2<-w1%>%select(c(name1,name2))
      w2[,name2][is.na(w2[,name2])] <- 'Not  Applicable'
      w2$treatment_related<-!apply(as.matrix(w2[,name2]),1,function(x) all(x=='Not  Applicable'))
      w2$treatment_related<-factor(w2$treatment_related,level=c(F,T),label=c('No','Yes'))
      w2=w2%>%mutate(t1=as.numeric(difftime(onset_date_of_ae,on_treatment_date,units='days')),
                     t2=as.numeric(difftime(resolved_date,on_treatment_date,units='days')),
                     t12=as.numeric(difftime(resolved_date,onset_date_of_ae,units='days')),
                     code=cdus_ctcae_toxicity_type_code)
      w2$index<-as.numeric(factor(w2$code))
      w2$grade<-as.factor(w2$grade)
      
      
      
      if("start_date_of_drug" %in% names( da_data_upload())) {
        da_data_subject <- da_data_upload() %>% filter(sequence_no == subject) 
      }else{
        
        da_data_subject <- da_data_upload() %>% mutate(start_date_of_drug = start_date) %>% filter(sequence_no == subject)
        
        
      }
      
      if("drug" %in% names( da_data_subject)) {
        da_data_subject <-da_data_subject %>% filter(sequence_no == subject) %>% 
           select(sequence_no, start_date_of_drug, cycle, drug) %>% 
          filter(!drug %in% c("Not  Applicable")) %>% 
          filter(is.na(drug) == FALSE)
        
        
      }else{
        
        da_data_subject <- da_data_subject %>% mutate(drug = level) %>% filter(sequence_no == subject) %>% 
          select(sequence_no, start_date_of_drug, cycle, drug) %>%
          filter(!drug %in% c("Not  Applicable")) %>%
          filter(is.na(drug) == FALSE)
        
        
      }
      
      
      # print("names(da_data_subject) and DIM:")
      # dim(da_data_subject)
      # print(names(da_data_subject))
      # 
      # print("AE DATA FOR PLOT 2 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
      # print("names and DIM AE data w2 vvvvvv what ?????")
      # print(names(w2))
      # print(dim(w2))
      # print("AE DATA FOR PLOT 2 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
      #w2.2 <- w2 %>% select(on_treatment_date) %>% distinct(on_treatment_date) 
      w2.2early <- w2 %>% distinct(on_treatment_date, .keep_all = TRUE) 

      
      
      # print("AE DATA FOR PLOT 2 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
      # print("  AE data w2.2 early")
      # print(w2.2early)
      # print("dim AE data w2.2early")
      # print(dim(w2.2early))
      # print("AE DATA FOR PLOT 2 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
      # print("names AE data w2.2early")
      # print(names(w2.2early))
      # w2.2 <- w2 %>% select(on_treatment_date) %>% distinct(on_treatment_date) 
      
      w3 <- merge(  da_data_subject, w2.2early) %>% mutate(dadaynumber = as.numeric(difftime(start_date_of_drug ,on_treatment_date,units='days')) )
      # print("some DA data w3")
      # print(w3)
      # print("dim DA data w3")
      # print(dim(w3)) 
      # 
      # 
      # print("names w3")
      # print(names(w3)) 
      
      UL <- max(w2$t2)+10
      early_AE_plot_manuscript.fun<-function(AE.data= w2,early.time=365,k=1.1,k1=5,early.AE.status=FALSE)
      {
        vline.fun<-function(early.AE.status.tmp,early.time.tmp) if(early.AE.status.tmp)    geom_vline(xintercept = early.time.tmp,linetype="dotted")
        geom_label.no_early_AE.fun<-function(early.AE.status.tmp) if(!early.AE.status.tmp) geom_label(aes(x=t2-(t12/2),y=index,label = t12), inherit.aes = F,hjust=0.5,size = 7)
        geom_label.early_AE.fun<-function(early.AE.status.tmp) if(early.AE.status.tmp) geom_label(aes(x=t2-(t12/2),y=index,label = t12.early), inherit.aes = F,hjust=0.5,size = 7)
        
        # w1<-AE.data
        g.cols=c("pink","purple",  'red',"blue",'black' )
        w1<-AE.data %>% filter(is.na(resolved_date ) == FALSE  & is.na(on_treatment_date ) == FALSE)
        
        subject<-unique(w1$sequence_no)
        name1<-c('onset_date_of_ae','cdus_ctcae_toxicity_type_code', 'resolved_date','on_treatment_date',#'AE.time',
                 'grade')
        name2<-c('attribution_possible','attribution_probable', 'attribution_definite')
        w2<-w1%>%select(c(name1,name2))
        #  w2$treatment_related<-!apply(as.matrix(w2[,name2]),1,function(x) all(x=='Not  Applicable'))
        w2$treatment_related<-!apply(as.matrix(w2[,name2]),1,function(x) all(x%in%c('Not  Applicable',NA)))
        w2$treatment_related<-factor(w2$treatment_related,level=c(F,T),label=c('No','Yes'))
        w2=w2%>%mutate(t1=as.numeric(difftime(onset_date_of_ae,on_treatment_date,units='days')),
                       t2=as.numeric(difftime(resolved_date,on_treatment_date,units='days')),
                       t12=as.numeric(difftime(resolved_date,onset_date_of_ae,units='days'))+1,
                       code=cdus_ctcae_toxicity_type_code)
        w2$index<-as.numeric(factor(w2$code))
        
        if(early.AE.status){
          w2=w2%>%mutate(AE.early.indicator=(t1<=early.time)*(t2>=early.time))%>%
            mutate(t12.early=ifelse(AE.early.indicator%in%1,paste(early.time-t1+ifelse(early.time==0,1,0),'+',t2-early.time+1,sep=''),t12))
        }
        
        # print("vvvvvvvvvvv THIS IS THE DATA TO PLOT vvvvvvvvvvvvvvvv")
        # print(w2)
        # print("^^^^^^^^^^^ THIS IS THE DATA TO PLOT ^^^^^^^^^^^^^^^^")
        plot1 <- ggplot(w2, aes(t1, index,color=grade,shape=treatment_related, label = code)) +
          geom_point(aes(t2, index),size=2*k1)+
          geom_segment(aes(xend = t2, yend = index), size = 1*k1, lineend = "butt")+
          geom_label.no_early_AE.fun(early.AE.status.tmp=early.AE.status)+
          geom_label.early_AE.fun(early.AE.status.tmp=early.AE.status)+
          #    geom_label(aes(x=t2-(t12/2),y=index,label = t12), inherit.aes = F,hjust=1,size = 7)+
          xlab('days')+ylab('')+
          scale_y_continuous(breaks=w2$index,labels=w2$code)+
          scale_x_continuous(breaks=(0:floor(max(w2$t2)/early.time))*early.time,limits=c(0, max(w2$t2)))+
          #    xlim(0,NA)+
          scale_color_manual(values=g.cols, drop=FALSE) +
          vline.fun(early.AE.status.tmp=early.AE.status,early.time.tmp=early.time)+
          theme(legend.position="top",
                axis.text=element_text(size=20*k),
                axis.title=element_text(size=20*k,face="bold"),
                legend.text = element_text(size=16*k),
                legend.title = element_text(size=20*k))
        
        w2$listdaynumber <- map2(w2$t1,w2$t2,function(.x, .y){.x:.y})
        w2$totalnumberofdayswithAEs = length(unlist( w2$listdaynumber))
        w2$numberofdayswithAEs <- length(unique(unlist( w2$listdaynumber)))
        #list(data= w2, plot=plot1)
        plot1
      }
      early_AE_plot_manuscript.fun()
    })
    
    output$single_var_plot_early <- renderPlot({   print(single_var_plot_early())    })
    
    
    
  # allow users to download the AE swimmers plot ####
   { output$download_single_var_plot <- downloadHandler(
        filename = function() {
            paste0("plot_of_", as.character(input$single_var_input), "_", Sys.Date(), ".png")
        },
        content = function(file) { ggsave(file, plot = single_var_plot(),  width = 14,
                                          height = 8,
                                          units = c("in"),#, "cm", "mm", "px"),
                                          dpi = 300 )   }
    )
   }
    
    
  # allow users to download the AE swimmers plot ####
    { output$download_single_var_plot_early <- downloadHandler(
      filename = function() {
        paste0("plot_of_", as.character(input$single_var_plot_early), "_", Sys.Date(), ".png")
      },
      content = function(file) { ggsave(file, plot = single_var_plot_early(),  width = 14,
                                        height = 8,
                                        units = c("in"),#, "cm", "mm", "px"),
                                        dpi = 300 )   }
    )
    }
    
  # end: AE days tab ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  #### 
    
  # Calculate all AE measures!!!!!!!!!!!!. ####
    alldataoutput <- reactive({ 
      
      #---functions-----
      toxicity.out.fun <- function(toxicity.data = toxicity.data, AE.time.cutoff = NULL) {
        #---identify subject without AE at the AE.time.cutoff--
        id.tmp <- sort(unique(toxicity.data$pid))
        print("HOW MANY PATIENTS")
        print(length(id.tmp))
        print(paste("NROW original tox:", NROW(toxicity.data)))
        if (!is.null(AE.time.cutoff)) { print("WE HAVE A CUTOFF!!!")
          toxicity.data <- toxicity.data %>% filter(AE.time < AE.time.cutoff)  # TRUNCATE TIME
          print(paste("NROW after truncating tox:", NROW(toxicity.data)))
          print(paste("now how many patients?:",length(sort(unique(toxicity.data$pid)))))
        }
        
        # print(names(toxicity.data ))
        # print(names(toxicity.data[, (2:6)]))
       # toxicity.data.by.id <- by(toxicity.data[, -(2:6)], toxicity.data$pid, data.frame) # makes a data frame for each id.
        #LINE bELOW CHAINGE 20221207
        toxicity.data.by.id <- by(toxicity.data, toxicity.data$pid, data.frame) # makes a data frame for each id.
        
        # print("1......... str toxicity.data.by.id :")
        # print(class(toxicity.data.by.id ))
        # print(length(toxicity.data.by.id ))
        # print(names(toxicity.data.by.id[[1]]))
        
        null.status <- sapply(toxicity.data.by.id, is.null) # tests if null?
        # print("2........... str null.status :")
        # print(class(null.status ))
        # print(length(null.status ))
        # print(table(null.status ))
        
        toxicity.data.by.id <- toxicity.data.by.id[!null.status] # takes if null.status is FALSE
        # print("3.................str toxicity.data.by.id NOT NULL STAtUS :")
        # print(class(toxicity.data.by.id ))
        # print(length(toxicity.data.by.id ))
        # print(str(toxicity.data.by.id,1,1 ))
        
        id.in.data <- names(toxicity.data.by.id)
        id.no.AE <- id.tmp[!id.tmp %in% id.in.data] # gets patient numbers with out AE 
        # name.toxicity <- c( # WHERE IS THIS USED ? nowhere I think?
        #   "pid", "start_date_of_course_cycle", "onset_date_of_ae", "resolved_date", "cdus_ctcae_toxicity_type_code",
        #   "toxicity_category", "adverse_event_description", "grade", "attribution_total_therapy_regimen",
        #   "is_this_event_an_ir_ae", "dose_limiting_toxicity", "action", "serious_adverse_event", "sae_reported", "outcome"
        # )
        # name.toxicity<-c('pid','start_date_of_course_cycle','Onset.Date.of.AE',"Resolved.Date","CDUS.CTCAE.Toxicity.Type.Code", "Toxicity.Category","Adverse.Event.Description","Grade", "Attribution" ,"Dose.Limiting.Toxicity", "Action", "Serious.Adverse.Event" , "SAE.Reported", "Outcome")
        
        name.group <- c("cdus_ctcae_toxicity_type_code", "toxicity_category")
        name.time <- c("onset_date_of_ae", "resolved_date")
        name.grade <- "grade"
        name.treatment.related <- c("attribution_possible", "attribution_probable", "attribution_definite")
        
        toxicity.type.name <- names(table(as.vector(toxicity.data[, name.group[1]])))
        # print("5............ toxicity.type.names")
        # print(toxicity.type.name)
        toxicity.category.name <- names(table(as.vector(toxicity.data[, name.group[2]])))
        
        table1 <- table(toxicity.data$cdus_ctcae_toxicity_type_code, toxicity.data$toxicity_category)
        
        # THIS GETs each type of AE in each category by looking at each column of table1
        toxicity.type.within.category <- apply(table1, 2, function(x) dimnames(table1)[[1]][x != 0])
        
        #    name.tox.summary<-c("freq.all.grade","duration.all.grade", "freq.grade12","duration.grade12","freq.grade3","duration.grade3","freq.all.grade.treatment.related","duration.all.grade.treatment.related","freq.grade12.treatment.related","duration.grade12.treatment.related","freq.grade3.treatment.related","duration.grade3.treatment.related")
        
        name.tox.summary <-
          c(
            "all.grade.occurrence", "all.grade.fre", "all.grade.duration",
            "grade12.occurrence", "grade12.fre", "grade12.duration",
            "grade3.occurrence", "grade3.fre", "grade3.duration",
            "all.grade.treatment.related.occurrence", "all.grade.treatment.related.fre", "all.grade.treatment.related.duration",
            "grade12.treatment.related.occurrence", "grade12.treatment.related.fre", "grade12.treatment.related.duration",
            "grade3.treatment.related.occurrence", "grade3.treatment.related.fre", "grade3.treatment.related.duration"
          )
        
        
        duration.fun <- function(x, index.tmp, AE.time.cutoff.tmp) {
          x <- x[index.tmp, , drop = F]
          # print(names(x))
          # print(x %>% select(resolved_date,onset_date_of_ae))
          AE.whole.duration <- as.numeric(difftime(x$resolved_date, x$onset_date_of_ae, units = "days")) + 1 #--add 1 day to avoid 0 for same day of onset and resolved
          # print("AE.whole.duration")
          # print(AE.whole.duration)
          if (!is.null(AE.time.cutoff.tmp)) {
            max.AE <- AE.time.cutoff.tmp - x$AE.time + 1 #--add 1 day to avoid 0 for same day of onset and initial treatment day
            ans <- ifelse(AE.whole.duration > max.AE, max.AE, AE.whole.duration)
          } else {
            ans <- AE.whole.duration
          }
          
          # print(ans)
          ans
        }
        
        #---generate long format data--
        tmp1 <- numeric()
        for (i in 1:length(toxicity.type.name))
        {
          
          #--new---
          # print("toxicity.data.by.id")
          #  print(names(toxicity.data.by.id[[1]]))
          # print(name.group %in% names(toxicity.data.by.id[[1]]))
          tmp0 <- sapply(
            toxicity.data.by.id,
            function(x) {
              # print(dim(x))
              # print(names(x))
              # print(name.group[1])
              # print(toxicity.type.name[i])
              index.all <- x[, name.group[1]] == toxicity.type.name[i]
              index.grade12 <- (as.numeric(as.vector(x[, name.grade])) < 3)
              index.grade3 <- (as.numeric(as.vector(x[, name.grade])) >= 3)
              index.tretament.related <- apply(x[, name.treatment.related], 1, function(x) any(x != "Not  Applicable"))
              all.grade.fre <- sum(index.all, na.rm = T)
              all.grade.occurrence <- as.numeric(all.grade.fre > 0)
              all.grade.duration <- sum(duration.fun(x, index.tmp = index.all, AE.time.cutoff.tmp = AE.time.cutoff), na.rm = T)
              grade12.fre <- sum(index.all * index.grade12, na.rm = T)
              grade12.occurrence <- as.numeric(grade12.fre > 0)
              grade12.duration <- sum(duration.fun(x, index.tmp = index.all & index.grade12, AE.time.cutoff.tmp = AE.time.cutoff), na.rm = T)
              grade3.fre <- sum(index.all * index.grade3, na.rm = T)
              grade3.occurrence <- as.numeric(grade3.fre > 0)
              grade3.duration <- sum(duration.fun(x, index.tmp = index.all & index.grade3, AE.time.cutoff.tmp = AE.time.cutoff), na.rm = T)
              
              all.grade.treatment.related.fre <- sum(index.all * index.tretament.related, na.rm = T)
              all.grade.treatment.related.occurrence <- as.numeric(all.grade.treatment.related.fre > 0)
              all.grade.treatment.related.duration <- sum(duration.fun(x, index.tmp = index.all & index.tretament.related, AE.time.cutoff.tmp = AE.time.cutoff), na.rm = T)
              
              grade12.treatment.related.fre <- sum(index.all * index.grade12 * index.tretament.related, na.rm = T)
              grade12.treatment.related.occurrence <- as.numeric(grade12.treatment.related.fre > 0)
              grade12.treatment.related.duration <- sum(duration.fun(x, index.tmp = index.all & index.grade12 & index.tretament.related, AE.time.cutoff.tmp = AE.time.cutoff), na.rm = T)
              
              grade3.treatment.related.fre <- sum(index.all * index.grade3 * index.tretament.related, na.rm = T)
              grade3.treatment.related.occurrence <- as.numeric(grade3.treatment.related.fre > 0)
              grade3.treatment.related.duration <- sum(duration.fun(x, index.tmp = index.all & index.grade3 & index.tretament.related, AE.time.cutoff.tmp = AE.time.cutoff), na.rm = T)
              
              c(
                all.grade.occurrence, all.grade.fre, all.grade.duration,
                grade12.occurrence, grade12.fre, grade12.duration,
                grade3.occurrence, grade3.fre, grade3.duration,
                all.grade.treatment.related.occurrence, all.grade.treatment.related.fre, all.grade.treatment.related.duration,
                grade12.treatment.related.occurrence, grade12.treatment.related.fre, grade12.treatment.related.duration,
                grade3.treatment.related.occurrence, grade3.treatment.related.fre, grade3.treatment.related.duration
              )
            }
          )
          
            # print("VVVVVVVVVVVV str(tmp0,1,0)*******************")
          # print(name.tox.summary)
          # print(str(tmp0,1,0))
          rownames(tmp0) <- name.tox.summary # not working now ? 
          tmp0 <- as.data.frame(tmp0) %>%
            add_rownames(var = "measurement") %>%
            relocate(measurement)
          tmp0 <- tmp0 %>% pivot_longer(cols = -1, names_to = "pid", values_to = "value")
          tmp0$AE <- toxicity.type.name[i]
          tmp0$AE.category <- dimnames(table1)[[2]][table1[toxicity.type.name[i], ] != 0]
          
          tmp1 <- rbind(tmp1, tmp0)
        }
        
        #---summary function---
        
        
        AE.summary.fun <- function(data, var1, summary.status = T, id.no.AE.tmp = id.no.AE) {
          # data is a long format matrix with pid, AE, AE.catergory, measurement type, and the value
          # var1 = NULL for overall summary over all AEs
          if (summary.status) {
            data.summary.long <- data %>%
              group_by(pid, measurement) %>%
              dplyr::summarise(sum = sum(value))
            data.summary.wide <- data.summary.long %>% pivot_wider(names_from = measurement, values_from = sum)
            data.tmp <- data.summary.wide
            if (length(id.no.AE.tmp) > 0) {
              data.tmp.no.AE <- data.tmp[1:length(id.no.AE.tmp), , drop = F]
              data.tmp.no.AE$pid <- id.no.AE.tmp
              data.tmp.no.AE[, 2:dim(data.tmp.no.AE)[2]] <- 0
              data.tmp <- rbind(data.tmp, data.tmp.no.AE)
            }
            AE.data.list <- data.tmp
          } else {
            # this step is to sum the value within pid and AE or AE.category
            data.summary.long <- data %>%
              group_by(pid, {{ var1 }}, measurement) %>%
              dplyr::summarise(sum = sum(value))
            data.summary.wide <- data.summary.long %>% pivot_wider(names_from = measurement, values_from = sum)
            data.for.name <- data.summary.wide %>% select({{ var1 }})
            name1 <- names(table(data.for.name[, 2]))
            
            AE.data.list <- map(
              as.list(name1),
              function(x) {
                # under each AE or AE.category
                data.tmp <- data.summary.wide %>% filter({{ var1 }} %in% x)
                if (length(id.no.AE.tmp) > 0) {
                  data.tmp.no.AE <- data.tmp[1:length(id.no.AE.tmp), , drop = F]
                  data.tmp.no.AE$pid <- id.no.AE.tmp
                  data.tmp.no.AE[, 3:dim(data.tmp.no.AE)[2]] <- 0
                  data.tmp <- rbind(data.tmp, data.tmp.no.AE)
                }
                data.tmp
              }
            )
            names(AE.data.list) <- name1
          }
          AE.data.list
        } # END AE.summary.fun
        
        data.raw.long <- tmp1
        toxicity.whole.summary.data <- AE.summary.fun(data = tmp1, var1 = "", summary.status = T, id.no.AE.tmp = id.no.AE)
        toxicity.category.summary.data <- AE.summary.fun(data = tmp1, var1 = AE.category, summary.status = F, id.no.AE.tmp = id.no.AE)
        toxicity.type.summary.data <- AE.summary.fun(data = tmp1, var1 = AE, summary.status = F, id.no.AE.tmp = id.no.AE)
        #  
        
        list(
          id.no.AE = id.no.AE
          ,data.raw.long = data.raw.long
           ,toxicity.type.within.category = toxicity.type.within.category
          ,toxicity.category.summary.data = toxicity.category.summary.data
          ,toxicity.type.summary.data = toxicity.type.summary.data
          ,toxicity.whole.summary.data = toxicity.whole.summary.data
        )
      } # end toxicity.out.fun

      
      alldataoutput <- toxicity.out.fun(toxicity.data = toxicity_data(), AE.time.cutoff = NULL)    
      alldataoutput
      
      })
   
  # MERGING AE measures with follow up and demographics to make TOXDATA with all the measures ####
    toxicity.whole.summary.data <- reactive({  
   
      toxdata00 <- merge(alldataoutput()$toxicity.whole.summary.data,fu_data_upload(),by.x = "pid", by.y = "sequence_no")
      toxdata0 <- merge(toxdata00,demographics_data(),by.x = "pid", by.y = "sequence_no")
 

      # print("names(toxdata0)")
      # print(names(toxdata0))
      
      #LINEs BELOW ADDED 20221207
      
      # print( "I(off_treatment_date.y %in% names(toxdata0))")
      # print( I("off_treatment_date.y" %in% names(toxdata0)))
      # 
      if(I("off_treatment_date.y" %in% names(toxdata0)) == FALSE){
        print("off_treatment_date.y NOT IN DATA REPLACe WITH off_treatment_date")
      toxdata0 <- toxdata0  %>% mutate(off_treatment_date.y = off_treatment_date )}
     
      if(I("last_visit_date" %in% names(toxdata0)) == FALSE){
        print("LAST VISIT NOT IN DATA REPLACe WITH LAST_FOLLOWUP_DATE")
        toxdata0 <- toxdata0  %>% mutate(last_visit_date = last_followup_date )}
      
      # print("names(toxdata0)")
      # print(names(toxdata0))
      # print(toxdata0 %>% select(pid,off_treatment_date.y,on_treatment_date,last_followup_date,last_visit_date))
      # 
     #above ^^^added 20221207 
      toxdata1 <- toxdata0  %>% mutate(os_censor  = ifelse(is.na(expired_date) == FALSE, 1, 0)) %>%
        mutate(sequence_no = pid,
               #treatment.time  = off_treatment_date.y - on_treatment_date,
               treatment.time  = difftime(off_treatment_date.y, on_treatment_date, units = "days"),
               #lastdate = pmax( last_followup_date, last_visit_date, na.rm = TRUE),
                lastdate = as.Date(pmax( last_followup_date, last_visit_date, na.rm = TRUE)),
               
               os_time  = ifelse(os_censor  == 1, interval(on_treatment_date, expired_date) / months(1),
                                 interval(on_treatment_date, lastdate) / months(1)))
      # print("toxdata1....:")
      #print(toxdata1 %>% select(on_treatment_date, date_of_progression,last_followup_date,last_visit_date,lastdate))
      # print(str(toxdata1 %>% select(on_treatment_date,date_of_progression,  last_followup_date,last_visit_date,lastdate)))
      
      toxdata <- toxdata1 %>%   mutate(   
               #lastdate = pmax( last_followup_date, last_visit_date, na.rm = TRUE),
               lastdate = as.Date(pmax( last_followup_date, last_visit_date, na.rm = TRUE)),
               pfs_censor  = ifelse(is.na(date_of_progression) == FALSE, 1, 0),
               pfs_time   =    ifelse(pfs_censor  == 1, 
                                      interval(on_treatment_date, date_of_progression) / months(1),
                                      interval(on_treatment_date, lastdate) / months(1)) )  %>% 
        select(sequence_no, pid, everything())
      # print("toxdata....:")
      #  print(toxdata %>% select(on_treatment_date, date_of_progression,last_followup_date,last_visit_date,lastdate))
      # print(str(toxdata %>% select(on_treatment_date,date_of_progression,  last_followup_date,last_visit_date,lastdate)))
      # 
      
      toxdata
      
      })
    
    
    
    
  # start: AE measures tab VVVVVVVVVVVVVVVVVVVVVVVVV#####
    # Display toxicity.whole.summary.data() 
    output$toxicitytableoutput <- DT::renderDataTable({
      #todisplay <-alldataoutput()$toxicity.whole.summary.data
      todisplay <- toxicity.whole.summary.data()
      todisplay <- todisplay[,c(1,3:NCOL(todisplay))]
      datatable(todisplay, rownames = FALSE)
    })
    
    
    
    # Down load toxicity.whole.summary.data()       
    output$toxicitymeasuresdownload <- downloadHandler(
      filename = function(){"toxicitymeasures.csv"}, 
      content = function(fname){
        #write.csv(alldataoutput()$toxicity.whole.summary.data, fname)
        write.csv(toxicity.whole.summary.data(), fname)
      }
    )
    
  # end: AE measures tab ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^####  
    
  # start: COXph AE measures tab VVVVVVVVVVVVVVVVVVVVVVV####
    # KM plots and coxph models for AE measures ###  
    
    kmandcoxphinfofromtoxdata <- reactive({  
      
      
      kmandboxplotsfunction <- function(df){
        png.status<-F
        AE.km.and.boxplot.list<-coef.list<-list()
        k<-0
        #dir.tmp<-?"
        #howmany<-1:18 # temp for now  
        plotthesevars <-   c("all.grade.duration","all.grade.fre",  "all.grade.occurrence","all.grade.treatment.related.duration",
                             "all.grade.treatment.related.fre","all.grade.treatment.related.occurrence",  "grade12.duration","grade12.fre",
                             "grade12.occurrence","grade12.treatment.related.duration",  "grade12.treatment.related.fre","grade12.treatment.related.occurrence",
                             "grade3.duration","grade3.fre",  "grade3.occurrence","grade3.treatment.related.duration",  "grade3.treatment.related.fre","grade3.treatment.related.occurrence") 
        #for(h in 3:20)
        for (h in match(plotthesevars,names(df)))#[howmany])
        { print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ next !!!') 
          a4<-df;
          k<-k+1
          # print("The outcome: ")
          # print(names(a4)[h])
          # Make outcomes ###
          a4$y<-a4[,h] # original continuous 
          a4$AE.bin<-I(a4[,h]>0) # indicator TRUE if >0
          # print("outcomes:")
          # print(table(a4$AE.bin))
          # print(sort(a4$y))
          # print(summary(a4$y))
          # best responses ###
          # print("table of best response:")
          # print(table(a4$best_response))
          name.BOR<-c("Complete Response", "Partial Response","Stable Disease","Progressive Disease")
          a40<-a4%>%filter(best_response%in%name.BOR)
          name.BOR1<-names(table(a40$best_response)[table(a40$best_response)>0])
          # print("name.BOR1")
          # print(name.BOR1)
          # Make a factor 
          # print("make a factor with these levels:")
          # print(name.BOR[name.BOR%in%name.BOR1])
          # print("table best responses before as factor")
          # print(class(a40$best_response))
          # print(table(a40$best_response))
          a40$best_response<-factor(a40$best_response,level=name.BOR[name.BOR%in%name.BOR1])
          # print("table best responses after as factor")
          # print(table(a40$best_response))
          # print(class(a40$best_response))
          # Summarize N by best response.
          nn<-a40%>%group_by(best_response)%>%summarise(n=n())
          # print(nn)
          # print('~~~~~~~~~~~~~~~~~~~~~~~~~~end Best response stuff')
          # Plot boxplots ###
          plot1<-ggplot(a40, aes(y=y,x=best_response,fill=best_response)) +
            geom_boxplot(outlier.shape=NA) +
            geom_jitter(position = position_jitter(w = 0.2, h = 0))+
            theme(legend.position="none")+
            labs(title=names(a40)[h])+ylab(names(a40)[h])+
            ggpubr::stat_compare_means( aes(group=best_response),label = "p.signif", method="t.test", comparisons = combn(1:length(name.BOR1), 2, FUN = list))+
            geom_text(data=nn, aes(best_response,Inf,label = paste('(n=',n,')',sep='')), vjust = 1)+
            stat_summary(fun=mean, geom="point", shape=8, size=7, color="red", fill="red")
          #if(png.status) png(here::here(dir.tmp,paste(trial.id[i],'_',names(a4)[h],'_best_response.png',sep='')),width = 1024, height = 768,res = 120)
          # x11();
          # print(plot1)
          #if(png.status) dev.off()
          #plot(factor(a4$best_response,level=name.level),a3[,i],main=names(a3)[i])
          
          # print("try this:  --- Disease Control (DC: CR/PR/SD) vs PD----")
          #--- Disease Control (DC: CR/PR/SD) vs PD----
          #if(MCC=='MCC18494')  a40<-a40%>%mutate(disease_control=factor(best_response,label=c('DC','DC','PD')))
          #if(MCC=='MCC18321')  a40<-a40%>%mutate(disease_control=factor(best_response,label=c('DC','DC','DC','PD')))
          
          a40<-a40%>%mutate(disease_control=as.factor(case_when(best_response %in% c("Complete Response","Partial Response","Stable Disease") ~ "DC",
                                                                best_response %in% c("Progressive Disease") ~ "PD")))
          # print(a40%>%select(best_response,disease_control))
          # print(str(a40%>%select(best_response,disease_control)))                  
          #                      
          
          # print("table disease_control")
          # print(table(a40$disease_control))
          # 
          # print("length(table(a4$AE.bin))")
          # print(length(table(a4$AE.bin)))
          # print(table(a4$AE.bin))
          # print(summary(a4$y))
          
          name.BOR2<-names(table(a40$disease_control))
          nn<-a40%>%group_by(disease_control)%>%summarise(n=n())
          plot11<-ggplot(a40, aes(y=y,x=disease_control,fill=disease_control)) +
            geom_boxplot(outlier.shape=NA) +
            geom_jitter(position = position_jitter(w = 0.2, h = 0))+
            theme(legend.position="none")+
            labs(title=names(a40)[h])+ylab(names(a40)[h])+
            ggpubr::stat_compare_means( aes(group=disease_control),method="t.test", comparisons = combn(1:length(name.BOR2), 2, FUN = list))+
            geom_text(data=nn, aes(disease_control,Inf,label = paste('(n=',n,')',sep='')), vjust = 1)+
            stat_summary(fun=mean, geom="point", shape=8, size=7, color="red", fill="red")+xlab('')
          #if(png.status) png(here::here(dir.tmp,paste(trial.id[i],'_',names(a4)[h],'_disease_control.png',sep='')),width = 1024, height = 768,res = 120)
          # x11();
          # print(plot11)
          #if(png.status) dev.off()
          
          if(length(table(a4$AE.bin))>1)
          { #print("OS ------ OS -----------------OS ______")
            #print("The outcome: ")
            #print(names(a4)[h])
            #print("what is the length of length(table(a4$AE.bin))")
            # print(length(table(a4$AE.bin)))
            # print("OK RUN COX MODELS for OS")
            cox1<-coxph(Surv(os_time,os_censor)~y,data=a4)
            # print('OS')
            # print(summary(cox1))
            tmp99.os<-coef(summary(cox1))
            
            if(exists("cox1")==TRUE){rm(cox1)} 
            cox1<-coxph(Surv(os_time,os_censor)~AE.bin,data=a4)
            # print('OS.bin-------------------------------------------')
            # print(summary(cox1))
            tmp99.os<-rbind(tmp99.os,coef(summary(cox1)))
            
            if(exists("cox1")==TRUE){rm(cox1)} 
            # print("END OS ------ END OS -----------------END OS ______")
          } else {
            # print("OS ------ OS -----------------OS ______")
            # print("The outcome: ")
            # print(names(a4)[h])
            # print("PUT NAs in BIN OUTPUT because length(table(a4$AE.bin)) < 2")
            cox1<-coxph(Surv(os_time,os_censor)~y,data=a4)
            #print('OS')
            # print(summary(cox1))
            tmp99.os<-coef(summary(cox1))
            #print(tmp99.os)
            #cox1<-coxph(Surv(os_time,os_censor)~AE.bin,data=a4)
            tmp.os <- coef(summary(cox1))
            
            # print('OS.bin tmp.os')
            tmp.os<-rep(NA,5)
            # print(tmp.os)
            tmp99.os<-rbind(tmp99.os,tmp.os)
            # print('BOTH ')
            # print(tmp99.os)
            
            if(exists("cox1")==TRUE){rm(cox1)} 
            # print("END OS ------ END OS -----------------END OS ______")
          }
          #---PFS--
          if(length(table(a4$AE.bin))>1)
          { 
            # print("PFS ------ PFS -----------------PFS ______")
            # print("what is the length of length(table(a4$AE.bin))")
            # print(length(table(a4$AE.bin)))
            # print("table(a4$AE.bin)")
            # print(table(a4$AE.bin))
            # print("OK RUN COX MODELS for PFS")
            cox1<-coxph(Surv(pfs_time,pfs_censor)~y,data=a4)
            # print('PFS')
            # print(summary(cox1))
            tmp99.pfs<-coef(summary(cox1))
            cox1<-coxph(Surv(pfs_time,pfs_censor)~AE.bin,data=a4)
            # print('PFS.bin')
            # print(summary(cox1))
            tmp99.pfs<-rbind(tmp99.pfs,coef(summary(cox1)))
            # print("and plot KM plots for OS and PFS")
            fit3 <- survfit( Surv(os_time,os_censor) ~ AE.bin, data = a4 )
            ggsurv <- ggsurvplot(fit3, data = a4, pval = TRUE,surv.median.line = 'hv',break.time.by = 4,risk.table = TRUE,
                                 xlab='months',
                                 tables.height = 0.2,
                                 legend.title = "AE",
                                 legend.labs = c("No", "Yes"),
                                 palette = c("cyan","red"),
                                 title=paste('OS: ',names(a4)[h]))
            plot2<-ggsurv
            #plot2<-ggsurv$plot +theme_bw()+labs(title=paste('OS: ',names(a4)[h]))+xlab('months')+theme(legend.position="top")
            #if(png.status) png(here::here(dir.tmp,paste(trial.id[i],'_',names(a4)[h],'_OS.png',sep='')),width = 1024, height = 768,res = 120)
            # x11()
            # print(plot2)
            #if(png.status) dev.off()
            
            fit3 <- survfit( Surv(pfs_time,pfs_censor)~ AE.bin,data = a4 )
            ggsurv <- ggsurvplot(fit3,data = a4, pval = TRUE,surv.median.line ='hv',break.time.by = 4,risk.table = TRUE,
                                 xlab='months',
                                 tables.height = 0.2,
                                 legend.title = "AE",
                                 legend.labs = c("No", "Yes"),
                                 palette = c("cyan","red"),
                                 title=paste('PFS: ',names(a4)[h]))
            plot3<-ggsurv
            #    plot3<- ggsurv$plot +theme_bw()+labs(title=paste('PFS: ',names(a4)[h]))+xlab('months')+theme(legend.position="top")
            #if(png.status) png(here::here(dir.tmp,paste(trial.id[i],'_',names(a4)[h],'_PFS.png',sep='')),width = 1024, height = 768,res = 120)
            # x11()
            # print(plot3)
            #if(png.status) dev.off()
            
            if(exists("cox1")==TRUE){rm(cox1)} 
            # print("END PFS ------ END PFS -----------------END PFS ______")
          } else {
            
            # print("PFS ------ PFS -----------------PFS ______")
            # print("The outcome: ")
            # print(names(a4)[h])
            # print("PUT NAs in OUTPUT because length(table(a4$AE.bin)) < 2")
            # print("what is the length of length(table(a4$AE.bin))")
            # print(length(table(a4$AE.bin)))
            # print("the table...^^^^^^^^^^^^PFS bin")
            # print(table(a4$AE.bin))
            # print("OK RUN COX MODELS for PFS")
            cox1<-coxph(Surv(pfs_time,pfs_censor)~y,data=a4)
            # print('PFS')
            # print(summary(cox1))
            tmp99.pfs<-coef(summary(cox1))
            # print(tmp99.pfs)
            
            
            #cox1<-coxph(Surv(pfs_time,pfs_censor)~AE.bin,data=a4)
            # print('PFS.bin')
            tmp.pfs <- coef(summary(cox1))
            
            if(exists("cox1")==TRUE){rm(cox1)} 
            # print('pfs.bin tmp.pfs')
            tmp.pfs<-rep(NA,5)
            # print(tmp.pfs) 
            tmp99.pfs<-rbind(tmp99.pfs,tmp.pfs) 
            # print('BOTH PFS ')
            # print(tmp99.pfs)
            
            if(exists("cox1")==TRUE){rm(cox1)} 
            
            
            # print("and plot KM plots for OS and PFS")
            # fit3 <- survfit( Surv(os_time,os_censor) ~ AE.bin, data = a4 )
            # ggsurv <- ggsurvplot(fit3, data = a4, pval = TRUE,surv.median.line = 'hv',break.time.by = 4,risk.table = TRUE,
            #                      xlab='months',
            #                      tables.height = 0.2,
            #                      legend.title = "AE",
            #                      legend.labs = c("No", "Yes"),
            #                      palette = c("cyan","red"),
            #                      title=paste('OS: ',names(a4)[h]))
            plot2<-NA
            #plot2<-ggsurv$plot +theme_bw()+labs(title=paste('OS: ',names(a4)[h]))+xlab('months')+theme(legend.position="top")
            #if(png.status) png(here::here(dir.tmp,paste(trial.id[i],'_',names(a4)[h],'_OS.png',sep='')),width = 1024, height = 768,res = 120)
            # x11()
            # print(plot2)
            #if(png.status) dev.off()
            
            # fit3 <- survfit( Surv(pfs_time,pfs_censor)~ AE.bin,data = a4 )
            # ggsurv <- ggsurvplot(fit3,data = a4, pval = TRUE,surv.median.line ='hv',break.time.by = 4,risk.table = TRUE,
            #                      xlab='months',
            #                      tables.height = 0.2,
            #                      legend.title = "AE",
            #                      legend.labs = c("No", "Yes"),
            #                      palette = c("cyan","red"),
            #                      title=paste('PFS: ',names(a4)[h]))
            plot3<-NA
            
            if(exists("cox1")==TRUE){rm(cox1)}
            #print("END PFS ------ END PFS -----------------END PFS ______")
          }
          AE.km.and.boxplot.list[[k]]<-list(bp1=plot1, bp2=plot11, kmos = plot2, kmpfs=plot3)
          coef.list[[k]]<-list(os=tmp99.os,pfs=tmp99.pfs)
          # print("print coef.list")
          # print(coef.list)
        }
        
        names(coef.list)<-plotthesevars
        
        names(AE.km.and.boxplot.list)<-plotthesevars
        
        # make the coef.list into a data frame maybe for output? 
        
        step1<-lapply(coef.list, function(x){   as.data.frame(do.call(rbind,x), row.names = NULL) %>% 
            mutate(outcome =c("OS","OS","PFS","PFS"),
                   predictor = c("y","AE.bin","y","AE.bin")) })
        #step1names <- names(step1)
        step2<-do.call(rbind,step1)
        coef.data <- as_tibble(as.data.frame(step2) %>% mutate(AEmeasure  = rep( c("all.grade.duration","all.grade.fre",  "all.grade.occurrence","all.grade.treatment.related.duration",
                                                                                   "all.grade.treatment.related.fre","all.grade.treatment.related.occurrence",  "grade12.duration","grade12.fre",
                                                                                   "grade12.occurrence","grade12.treatment.related.duration",  "grade12.treatment.related.fre","grade12.treatment.related.occurrence",
                                                                                   "grade3.duration","grade3.fre",  "grade3.occurrence","grade3.treatment.related.duration",  "grade3.treatment.related.fre","grade3.treatment.related.occurrence") 
                                                                                 ,each = 4)))   %>% mutate(AEmeasure = paste(AEmeasure,".",outcome,".",predictor, sep='')) %>%
          select(AEmeasure, coef, HR=`exp(coef)`, `se(coef)`,      z, `Pr(>|z|)`)
        
        
        # end o'Loop de loop
        return(list(coef.data = coef.data, plots = AE.km.and.boxplot.list))
      }
      
      
      testing <- kmandboxplotsfunction( toxicity.whole.summary.data())
      testing
      
    })
    
  # Display Coxph output for AE measures ####
    output$kmcophinfo <- DT::renderDataTable({
      todisplay <-  kmandcoxphinfofromtoxdata()$coef.data   %>% 
        mutate(across(where(is.numeric), round, 4)) %>% 
        filter(is.na(coef)==FALSE)
      datatable(todisplay, rownames = FALSE)
      
    }, options = list(pageLength = 72))
    
  # KM and boxplots of responses####
    draw_plot <- function() {  kmandcoxphinfofromtoxdata()$plots  }
    
  # Download km and box plots ####
    output$plotall <- renderPlot({  draw_plot() })
    
  # Download KM plots pdf ####
    output$downloadplotspdf <- downloadHandler(
      filename = "Boxplots_and_KM_plots_report.pdf",
      content = function(file) {
        res <- rmarkdown::render(
          "template.Rmd",
          params = list(
            draw_plot = draw_plot
          )
        )
        file.rename(res, file)
      }
    )
  # end: Coxph AE measures tab ^^^^^^^^^^^^^^^^^^^^^^^####   
    
    
  # start: Forest plots OS and PF tab VVVVVVVVVVVVVVVV####
    # FOREST PLOTS FOR P-values ### 
    responsePvaluelistforestplots <- reactive({
      
      toxdataNOW <- toxicity.whole.summary.data()
      
      
      a400<-list(whole=alldataoutput()$toxicity.whole.summary.data)
 
      
      a400<-c(a400,
              alldataoutput()$toxicity.category.summary.data,
              alldataoutput()$toxicity.type.summary.data)
      a41now<-a400
      
      
      forestplots<-function(a41=a41now,a4.whole.summary=toxdataNOW){ # input will be a a list of 3 data sets and a data set with survival information, 
        N <-  length(a41)
        #--calcuate cox coef and p value for OS and PFS----
        coef.list.group<-list()
        for(j in 1:N)#length(a41))
          # for(j in 1:1)
        {
          # print("~~~~~~~~~LIST  NAME either AE or AE category ~~~~~~~~~~:")
          # print(names(a41)[j])
          AETYPEorCAT <-names(a41)[j]
          
          
          
          a4<-left_join(a41[[j]],a4.whole.summary[,c(2,21:dim(a4.whole.summary)[2])])
          #print(names(a4))
          # print(head(a4,5))
          png.status<-F
          coef.list<-list()
          k<-0
          
          measures<-c("all.grade.duration","all.grade.fre",
                      "all.grade.occurrence","all.grade.treatment.related.duration","all.grade.treatment.related.fre",       
                      "all.grade.treatment.related.occurrence", "grade12.duration","grade12.fre",                     
                      "grade12.occurrence","grade12.treatment.related.duration","grade12.treatment.related.fre",         
                      "grade12.treatment.related.occurrence","grade3.duration","grade3.fre",                            
                      "grade3.occurrence","grade3.treatment.related.duration","grade3.treatment.related.fre",          
                      "grade3.treatment.related.occurrence" )
          
          #for(h in 3:20)
          for(h in  match(measures,names(a4)) )
          {
            k<-k+1
            # print("~~~~~~~~~a4 names ~~~~~~~~~~:")
            # print(names(a41[[j]]))
            # print(names(a4.whole.summary[,c(2,21:dim(a4.whole.summary)[2])]))
            # 
            # print( names(a4))
            # print("~~~~~~~~~VAR NAME ~~~~~~~~~~:")
            # print(names(a4)[h])
            a4=a4%>%mutate(y=.data[[names(a4)[h]]], AE.bin=I(y>0))
            name.BOR<-c("Complete Response",   "Partial Response","Stable Disease","Progressive Disease")
            a40<-a4%>%filter(best_response%in%name.BOR)
            name.BOR1<-names(table(a40$best_response)[table(a40$best_response)>0])
            a40$best_response<-factor(a40$best_response,level=name.BOR[name.BOR%in%name.BOR1])
            #print("the data----------------------->")
            # print(a4 %>% select(os_time,os_censor,y))
            if(length(table(a4$AE.bin))>1)
            {
              if(length(table(a4$AE.bin))>1)
              {
                cox1<-coxph(Surv(os_time,os_censor)~y,data=a4)
                # print('OS')
                # print(summary(cox1))
                temp<-coef(summary(cox1)) 
                # print("what are the names of temp:")
                # print(temp)
                # print(names(temp))
                tmp99.os<- c(coef(summary(cox1)),
                             LCL= round(exp(temp[1]  - 1.96*temp[3]),3),
                             UCL= round(exp(temp[1]  + 1.96*temp[3]),3))
                # print( names(coef(summary(cox1))))
                names(tmp99.os)[1:5] <- c( "coef", "exp(coef)" , "se(coef)", "z", " Pr(>|z|)" ) #names(coef(summary(cox1)))
                
                if(exists("cox1")==TRUE){rm(cox1)} 
                cox1<-coxph(Surv(os_time,os_censor)~AE.bin,data=a4)
                temp<-coef(summary(cox1)) 
                # print('OS.bin')
                # print(summary(cox1))
                tmp99.os.bin<- c(coef(summary(cox1)),
                                 LCL= round(exp(temp[1] - 1.96*temp[3]),3),
                                 UCL= round(exp(temp[1] + 1.96*temp[3]),3))
                #names(tmp99.os.bin)[1:5] <- names(coef(summary(cox1)))
                names(tmp99.os.bin)[1:5] <- c( "coef", "exp(coef)" , "se(coef)", "z", " Pr(>|z|)" )
                #tmp99.os<-rbind(tmp99.os,coef(summary(cox1))) tmp99.os.bin
                tmp99.os<-rbind(tmp99.os,tmp99.os.bin)  
                rownames(tmp99.os)<- c("y","AE.binTRUE")
                
                if(exists("cox1")==TRUE){rm(cox1)} 
              }
              #---PFS-----
              if(length(table(a4$AE.bin))>1)
              {
                cox1<-coxph(Surv(pfs_time,pfs_censor)~y,data=a4)
                # print('PFS')
                # print(summary(cox1))        
                temp<-coef(summary(cox1)) 
                tmp99.pfs<- c(coef(summary(cox1)),
                              LCL= round(exp(temp[1] - 1.96*temp[3]),3),
                              UCL= round(exp(temp[1] + 1.96*temp[3]),3))
                #names(tmp99.pfs)[1:5] <- names(coef(summary(cox1)))
                names(tmp99.pfs)[1:5] <- c( "coef", "exp(coef)" , "se(coef)", "z", " Pr(>|z|)" )
                
                if(exists("cox1")==TRUE){rm(cox1)} 
                #tmp99.pfs<-coef(summary(cox1))
                cox1<-coxph(Surv(pfs_time,pfs_censor)~AE.bin,data=a4)
                temp<-coef(summary(cox1)) 
                # print('PFS.bin')
                # print(summary(cox1))
                tmp99.pfs.bin<- c(coef(summary(cox1)),
                                  LCL= round(exp(temp[1] - 1.96*temp[3]),3),
                                  UCL= round(exp(temp[1] + 1.96*temp[3]),3))
                #names(tmp99.pfs.bin)[1:5] <- names(coef(summary(cox1)))
                names(tmp99.pfs.bin)[1:5] <- c( "coef", "exp(coef)" , "se(coef)", "z", " Pr(>|z|)" )
                
                if(exists("cox1")==TRUE){rm(cox1)} 
                
                tmp99.pfs<-rbind(tmp99.pfs,tmp99.pfs.bin) 
                # print("rownames(tmp99,pfs)")
                rownames(tmp99.pfs)<- c("y","AE.binTRUE")
                # print(rownames(tmp99.pfs))
                
                if(exists("cox1")==TRUE){rm(cox1)} 
                # fit3 <- survfit( Surv(os_time,os_censor) ~ AE.bin, data = a4 )
                # fit3 <- survfit( Surv(pfs_time,pfs_censor)~ AE.bin,data = a4 )
                
              }
              coef.list<-c(coef.list,list(list(os=tmp99.os,pfs=tmp99.pfs)))
              names(coef.list)[length(coef.list)]<-names(a4)[h]
            }
            
          }  #--for h loop--
         # print("~~~~~~~do it again !!!!!!!!!!!!!!!!!!")
          #names(coef.list)<-names(a4)[3:20]
          coef.list.group[[j]]<-coef.list
        }  #--for j loop--
        
        names(coef.list.group)<-names(a41)[1:N]#length(a41)
        #names(coef.list.group)
        
        # print("*******************************************lets look at coef.list.group:")
        # print(coef.list.group[[1]])
        # print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^lets look at coef.list.group")
        # 
        # save(coef.list.group,file="MCC18786.coef.list.group.RData")
        
        
        # NOW TAKE THE OUTPUT OF THAT LOOP DE LOOP and use it to plot the p-values and HRs....--- 
        
        #---get occurrence (yes/no) based on dichotomized AE----
        fun1<-function(my.data.list=coef.list) # Input is coef.list.group a list of coxph results... coeffiences and pvalues etc... 
        {
          surv.name<-c('os','pfs')
          e2.comb<-numeric()
          for(i in 1:length(surv.name)) #loop over os and pfs vars
          {
            e1=sapply(my.data.list,function(x) x[[surv.name[i]]][,c(2,5:7)],simplify=F) # get the HR and pvalue
            # print("e1")
            # print(e1)
            e2<-data.frame(do.call("rbind",e1))    #combine them into a data frame
            # print("e2")
            # print(e2)
            dim.index<-dim(e2[!is.na(e2[,2]),])[1] #--ensure not all NA--
            if(dim.index>0)
            {
              #What the heck is going on here... mutate making a var to indicate OS or PFS and a var for AE type and measure type either continuous or Occurrence
              e2<-e2 %>% tibble::rownames_to_column(var='data.type') %>% #add_rownames(var='data.type') %>% 
                mutate(survival=surv.name[i], 
                       AE.type=rep(names(my.data.list),each=2),
                       data.type=sub('y.*','Continuous',sub('AE.*','Occurrence',data.type))) %>%
                relocate(AE.type)
              names(e2)[3:4]<-c('HR','p') # rename
              # print("e2 again ....")
              # print(e2)
              # print(dim(e2))
              e21<-e2%>%filter(data.type %in%'Occurrence')%>%slice(ends_with('occurrence',vars=AE.type))
              # print("e21, e22, e23")
              # print(head(e21,7))
              # print(dim(e21))
              e22<-e2%>%filter(data.type %in%'Continuous')%>%slice(ends_with('occurrence',vars=AE.type))%>%mutate(AE.type=sub('occurrence','sum_unique_AE',AE.type))
              # print(head(e22,7))
              # print(dim(e22))
              e23<-e2%>%filter(data.type %in%'Continuous')%>%slice(-ends_with('occurrence',vars=AE.type))
              # print(head(e23,12))
              # print(dim(e23))
              e2.comb<-rbind(e2.comb,rbind(e21,e22,e23))
              # print("ecombined")
              #print(e2.comb)
              # print(dim(e2.comb))
            }
          }
          
          e2.comb
        }
        
        #---generate AE survival p value plot for each AE----
        #---generate AE survival p value plot for each AE---
        AE.survival.p.plot.list<-list() 
        
        AE.survival.p.forestplot.list<-plot.data.AE.survival.p.forestplot.list<-list()
        k<-0
        e3<-numeric()
        for(i in 1:length(coef.list.group))
          # for(i in 1:N)
        {
          name.tmp<-names(coef.list.group)[i] # AE type ( or all = whole)
          e2<-fun1(my.data.list=coef.list.group[[i]]) # HERE WE call the function above. 
          
          if(length(e2)>0)
          {
            k<-k+1
            e2<-e2%>%mutate(AE=name.tmp) # make a variable for the TYPE
            # print("*******************************************lets look at e2 ")  
            # print(e2)
            # print(names(e3))
            # print("number of cols???")
            # print(NCOL(e3))
            # print(NCOL(e2))
            e3<-rbind(e3,e2) #rowbind something with e2 
            
            # NOW PLOTS saved in a list
            # AE.survival.p.plot.list[[k]]<-e2%>%ggplot(aes(y=AE.type,x=HR,fill=p<0.05))+geom_bar(stat = 'identity')+scale_x_log10()+facet_wrap(vars(survival))+labs(title=name.tmp)+
            #   labs(fill='P value')+scale_fill_manual(values=c("cyan","red"),breaks=c(F,T),labels=c('>0.05', '<0.05'))
            # names(AE.survival.p.plot.list)[k]<-name.tmp
            e3<-e3%>%relocate(AE)
            #as.data.frame(e3) -> e4
            # print("vvvvvvvAETYPEorCATvvvvvvvv")
            # print(AETYPEorCAT)
            
            # e4 is the plot data 
            
            e4<- as.data.frame(e2 %>% group_by(survival) %>% arrange(survival, AE.type) %>% mutate(index =  n():1))  %>% 
              # filter( is.finite(HR)==TRUE)  %>% 
              # filter(HR < 10 ) %>% 
              # filter(HR > .1) %>%
              # filter(is.finite(UCL)==TRUE) %>%
              mutate(
                HR = case_when(HR > 10 ~ NA_real_,
                               HR >= 0.01 & HR <= 10 ~ HR,
                               HR < 0.01 ~ NA_real_,
                               !is.finite(HR)~ NA_real_) ,
                UCL = case_when(HR > 10 ~ NA_real_,
                                HR >= 0.01 & HR <= 10 ~ UCL,
                                HR < 0.01 ~ NA_real_,
                                !is.finite(HR)~ NA_real_ ),
                LCL = case_when(HR > 10 ~ NA_real_,
                                HR >= 0.01 & HR <= 10 ~ LCL,
                                HR < 0.01 ~ NA_real_,
                                !is.finite(HR)~ NA_real_ ) ,
                UCL = case_when(!is.infinite(UCL)  ~ UCL,
                                is.infinite(UCL) ~ NA_real_ ),
                LCL = case_when(!is.infinite(LCL)  ~ LCL,
                                is.infinite(LCL) ~ NA_real_  ) 
                
                
                
              )
            
            
            e5<- as.data.frame(e3   %>% group_by(survival) %>% arrange(survival, AE.type) %>% mutate(index = n():1)) %>% 
              mutate(
                HR = case_when(HR > 10 ~ NA_real_,
                               HR >= 0.01 & HR <= 10 ~ HR,
                               HR < 0.01 ~ NA_real_,
                               !is.finite(HR)~ NA_real_) ,
                UCL = case_when(HR > 10 ~ NA_real_,
                                HR >= 0.01 & HR <= 10 ~ UCL,
                                HR < 0.01 ~ NA_real_,
                                !is.finite(HR)~ NA_real_ ),
                LCL = case_when(HR > 10 ~ NA_real_,
                                HR >= 0.01 & HR <= 10 ~ LCL,
                                HR < 0.01 ~ NA_real_,
                                !is.finite(HR)~ NA_real_ ) ,
                UCL = case_when(!is.infinite(UCL)  ~ UCL,
                                is.infinite(UCL) ~ NA_real_ ),
                LCL = case_when(!is.infinite(LCL)  ~ LCL,
                                is.infinite(LCL) ~ NA_real_  ) 
                
              )
            
            #VVVVVVVVVVVVVVVVGGPLOT VVVVVVVVVVV~~~~####
            #AE.survival.p.plot.list[[k]] <- 
            # print("length( e4$AE.type)")
            # print(length( e4$AE.type))
            # print( e4$AE.type)
            # print("e4$index")
            # print(e4$index)
            # print("VVVVVVVVVe4vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
            # print(e4 %>% select(AE.type,survival,index))
            plot.data.AE.survival.p.forestplot.list[[k]] <- e4
            AE.survival.p.forestplot.list[[k]]<- ggplot(e4, aes(y = index, x = HR, col = p<0.05)) + #
              scale_colour_manual(values=c("red","cyan"),breaks=c(T,F),labels=c('<0.05', '>0.05')) +
              facet_wrap(vars(survival)) + 
              geom_point(shape = 18, size = 3) +
              geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.25) +
              geom_vline(xintercept = 1, color = "black", linetype = "dashed", cex = .8, alpha = 0.5) +
              scale_y_continuous(name = "", breaks=e4$index, labels = e4$AE.type)+#, trans = "reverse") +  
              xlim(  c( -1,max(e4$UCL)+1)) +
              xlab("Hazard Ratio (95% CI)") +
              ylab(" ")+
              labs(title=name.tmp)
            names(AE.survival.p.forestplot.list)[k]<-name.tmp
            names(plot.data.AE.survival.p.forestplot.list)[k]<-name.tmp
            
          }
          
          
        }
        
        #---save survival p value data----
        responsePvaluelist <-list(coef=coef.list.group,
                                  coef.long=e3,
                                  coef.long4=plot.data.AE.survival.p.forestplot.list,
                                  filtered.coef.long=e5
                                  ,plot=AE.survival.p.forestplot.list
                                  
        )
      }
      
      forestplots(a41=a41now,a4.whole.summary=toxdataNOW)
    })
    
 
    output$forestplotcophinfo <- DT::renderDataTable({
      todisplay <-  responsePvaluelistforestplots()$coef.long %>% mutate(across(where(is.numeric), round, 4))
      #todisplay
      datatable(todisplay, rownames = FALSE)
    }, options = list(pageLength = 100))
    
     
    # draw_forestplot of responses ###  
    draw_forestplot <- function() {  responsePvaluelistforestplots()$plot   }

    # render draw_forestplot ###  
    output$plotallforest <- renderPlot({   draw_forestplot() })
    
  # Download FOREST PLOTS pdf ####
    output$downloadforestplotspdf <- downloadHandler(
      filename = "rendered_forest_plots_report.pdf",
      content = function(file) {
        res <- rmarkdown::render(
          "template_forestplots.Rmd",
          params = list(
            draw_forestplot = draw_forestplot 
          )
        )
        file.rename(res, file)
      }
    )
    
  # end: Forest plots OS and PF tab   ^^^^^^^^^^^^^^^^####  

    
  # start: Response testz plots tab VVVVVVVVVVVVVVVV#### 
  # testing_BOR_ans() reactive data set plots and t.test results... ####
      
      testing_BOR_ans   <- reactive({  
        toxdataNOW <- toxicity.whole.summary.data()
        
        
        a400<-list(whole=alldataoutput()$toxicity.whole.summary.data)
        
        
        a400<-c(a400,
                alldataoutput()$toxicity.category.summary.data,
                alldataoutput()$toxicity.type.summary.data)
        a41now<-a400
        
        
        responseplots<- function(a41=a41,adf = a1$toxicity.ans.summary_sub){
          
          a4.whole.summary<-adf;#a1$toxicity.ans.summary_sub
          
          #trial.id<-'MCC18494'
          coef.list.group<-AE.BOR.p.plot.list<-list()
          for(j in 1:length(a41))
          {
            a41.tmp<-ungroup(a41[[j]])
            var.drop<-c('AE.category','AE')[c('AE.category','AE')%in%names(a41.tmp)]
            #---get AE occurrence from sum.unique.AE
            a41.tmp<-select(a41.tmp,-all_of(var.drop))
            a41.tmp.sum.unique.AE<-a41.tmp%>%select(c(pid,ends_with('occurrence')))
            a41.tmp.occurrence<-as.list(a41.tmp.sum.unique.AE)[-1]%>%map_dfr(function(x) as.numeric(x>0))
            names(a41.tmp.sum.unique.AE)<-sub('occurrence','sum.unique',names(a41.tmp.sum.unique.AE))
            a41.tmp.fre.duraiton<-a41.tmp%>%select(c(pid,ends_with(c('fre','duration'))))
            a41.tmp.new<-cbind(a41.tmp.fre.duraiton,a41.tmp.sum.unique.AE%>%select(-pid),a41.tmp.occurrence)
            a4<-suppressMessages(left_join(a41.tmp.new,a4.whole.summary[,c(2,21:dim(a4.whole.summary)[2])]))
            name.BOR<-c("Complete Response",   "Partial Response","Stable Disease","Progressive Disease")
            name.BOR.short<-c("CR",   "PR","SD","PD")
            a40<-a4%>%filter(best_response%in%name.BOR)
            name.BOR1<-names(table(a40$best_response)[table(a40$best_response)>0])
            a40$best_response<-factor(a40$best_response,level=name.BOR[name.BOR%in%name.BOR1],label=name.BOR.short[name.BOR%in%name.BOR1])
            AE.var<-names(a40)[2:25]
            data.tmp<-a40%>%select(c(pid,best_response,all_of(AE.var)))%>%pivot_longer(cols=-(1:2),names_to ='type',values_to = 'value' )
            #---comparison of DC (CR/PR/SD) vs PD---
            tmp.com<-list('DC_vs_PD'=name.BOR.short,'PR_vs_PD'=name.BOR.short[c(1,2,4)],'SD_vs_PD'=name.BOR.short[c(3,4)])
            var1<-expand.grid(tmp.com,AE.var)
            var1$BOR<-rep(names(tmp.com),length(AE.var))
            
   
            if("AE.category" %in% names(a41[[j]])) {var1$AEtype <- unique(a41[[j]]$AE.category)} else {
              if("AE" %in% names(a41[[j]])) {var1$AEtype <- unique(a41[[j]]$AE)} else {var1$AEtype <- "whole"
              }
            }#end if statements
            
            # print("*****************************")
            # print( head(var1))
            options(warn=-1)
            tissue.test.ans<-  pmap(var1,~data.tmp%>%filter((best_response%in%..1)&(type%in%..2)))%>%
              map(function(x)
              {
                 
                df <- x %>% 
                  mutate(g1 = as.character(factor(x$best_response=='PD',level=c(T,F),label=c('PD','DC')))) %>%
                  select(value,g1)
                
                if( inherits(
                  try(tmp1<-suppressMessages(t.test(df$value~df$g1))
                      
                  )
                  , "try-error",T)
                  
                ){  
                  ans<-c(NA, NA) } else{
                    tmp1<-suppressMessages(t.test(df$value~df$g1))
                    
                    #print(tmp1);
                    ans<-c(tmp1$p.value,diff(tmp1$estimate))
                  }
                
              })
            
            
            tissue.test.ans.long<-cbind(var1,t(sapply(tissue.test.ans,c)))
            
           
            
            tissue.test.ans.long<-cbind(var1,t(sapply(tissue.test.ans,c)))%>%rename_all(~c('Best_response','Measurement.type','BOR','AE','p.value','difference'))%>%
              mutate(Measurement.type=factor(Measurement.type,level=sort(names(table(Measurement.type))))) 
            
          
            plot1<-tissue.test.ans.long%>%ggplot(aes(y=Measurement.type,x=difference,fill=p.value<0.05))+geom_bar(stat = 'identity')+
              #labs(title=paste(trial.id,names(a41)[j],sep='_'),fill='P value')+
              labs(title=paste("",names(a41)[j],sep=''),fill='P value')+
              scale_fill_manual(values=c("cyan","red"),breaks=c(F,T),labels=c('>0.05', '<0.05'))+
              ylab('')+xlab('Difference (PD as reference)')+
              facet_wrap(vars(BOR))
            plot2<-tissue.test.ans.long%>%filter(BOR%in%'DC_vs_PD')%>%ggplot(aes(y=Measurement.type,x=difference,fill=p.value<0.05))+geom_bar(stat = 'identity')+
              #labs(title=paste(trial.id,names(a41)[j],sep='_'),fill='P value')+
              labs(title=paste("",names(a41)[j],sep=''),fill='P value')+
              scale_fill_manual(values=c("cyan","red"),breaks=c(F,T),labels=c('>0.05', '<0.05'))+
              ylab('')+xlab('Difference (DC-PD )')
            
            AE.BOR.p.plot.list[[j]]<-list(all_comparison=plot1,DC_PD=plot2)
            coef.list.group[[j]]<-tissue.test.ans.long
          }
          
          names(AE.BOR.p.plot.list)<-names(coef.list.group)<-names(a41)
          
          
          twolists<-list(coef=coef.list.group,plot=AE.BOR.p.plot.list)
          
        }
        
        #testing_BOR_ans<- responseplots(a41=a41now,adf=toxdataNOW)
        
        #as_tibble(do.call(rbind,testing_BOR_ans$coef) %>% select(AE,everything()))
     
        testing_BOR_ans<- responseplots(a41=a41now,adf=toxdataNOW)
 
        testing_BOR_ans
        
        
        
         })
  # Display t-test pvalues and differences ####
      output$responsettestoutput <- DT::renderDataTable({
          todisplay <-  as_tibble(do.call(rbind,testing_BOR_ans() $coef) %>% 
                                    select(AE,everything(),-Best_response)) %>%
            mutate(comparison = case_when(BOR == "DC_vs_PD" ~ "(CR, PR, SD) vs PD",
                                          BOR == "PR_vs_PD" ~ "(CR, PR) vs PD",
                                          BOR == "SD_vs_PD" ~ "SD vs PD"),
                   p.value=round(p.value,6),
                   difference = round(difference,2))
          
          datatable(todisplay, rownames = FALSE)
          
        }, options = list(pageLength = 72))
        

    
    
    
    
    
    #  plots of responses ###  
    draw_responseplots <- function() {  testing_BOR_ans()$plot   }
    
    # render draw_responseplots ###  
    output$plotall_responseplots <- renderPlot({   draw_responseplots() })
    
    # Download FOREST PLOTS pdf ####
    output$download_responseplots <- downloadHandler(
      filename = "Response_plots_report.pdf",
      content = function(file) {
        res <- rmarkdown::render(
          "template_responseplots.Rmd",
          params = list(
            draw_responseplots = draw_responseplots 
          )
        )
        file.rename(res, file)
      }
    )
    
  # end: Response testz plots       ^^^^^^^^^^^^^^^^####     
      
    
  # start: correlation tab VVVVVVVVVVVVVVVV####    
    # TABLE OF SIGNIFICANT FINDINGS  ? or DURATION? ### 
      durationanalysis   <- reactive({  


        toxdataNOW <- toxicity.whole.summary.data() 
        # print("names(toxdataNOW)")
        # print(names(toxdataNOW))
        # 
        a0<-toxdataNOW
        
        cor.AE.treatment.time.ans<-cor.data<-list()
        
        #a1<-a0.list$`30`
        a1<-alldataoutput() 
        #a0<-a0.list$`30`$toxicity.ans.summary_sub
        a0<-a0%>%mutate(treatment.time=as.numeric(treatment.time))
        AE.var=names(a0)[3:20]
        a0=a0%>%select(-c(3:20))
        # print("names(a0)")
        # print(names(a0))
        
        # a400<-list(whole=a1$toxicity.ans.raw$toxicity.whole.summary.data)
        # a400<-c(a400,
        #         a1$toxicity.ans.raw$toxicity.category.summary.data,
        #         a1$toxicity.ans.raw$toxicity.type.summary.data)
        # a41<-a400
        
        a400<-list(whole=alldataoutput()$toxicity.whole.summary.data) #FROM 616
        
        a400<-c(a400,
                alldataoutput()$toxicity.category.summary.data,
                alldataoutput()$toxicity.type.summary.data)
        a41<-a400
        
        
        
        
        cor.plot.AE.treatment.time<-list()
        for(i in 1:length(a41))
        {
          tmp1<-left_join(a0,a41[[i]])%>%select(pid,treatment.time,AE.var)
          options(warn=-1)
          # print("~~~~~~~~DATA FOR COR()~~~~~~~~~~~~~")
          # print(tmp1[,c('treatment.time',AE.var)])
          cor1<-cor(tmp1[,c('treatment.time',AE.var)],method = 'pearson', use = "complete.obs")[1,]
          cor1<-data.frame(AE.measurement.type=names(cor1)[-1],r=round(cor1[-1],3))
          print(head(cor1))
          treatment.time<-tmp1$treatment.time
          cor1$pvalue<-apply(tmp1[,AE.var],2,function(x){ 
            
            # cor.test(x,treatment.time,method='pearson')$p.value
            
            if( inherits(
              try(cor.test(x,treatment.time,method='pearson')$p.value
                  
              )
              , "try-error",T)
              
            ){ NA } else{
              round(cor.test(x,treatment.time,method='pearson')$p.value,5)
              }
            }
            )#end apply
          # print(cor1)
          # print(dim(cor1))
          plot1<-cor1%>%ggplot(aes(y=AE.measurement.type,x=r,fill=pvalue<0.05))+geom_bar(stat = 'identity')+scale_fill_manual(values=c("cyan","red"),breaks=c(F,T),labels=c('>0.05', '<0.05'))+
            xlab('correlation coefficient (r)')+ylab('')+labs(title=names(a41)[i],fill='P value')
          cor.plot.AE.treatment.time[[i]]<-plot1
          cor.data[[i]]<-cor1
        }
        names(cor.plot.AE.treatment.time)<-names(cor.data)<-names(a41)
        cor.AE.treatment.time.ans<-list(cor=cor.data,plot=cor.plot.AE.treatment.time)
        #cor.AE.treatment.time.ans$trial.id=trial.id
    
        a2=cor.AE.treatment.time.ans$cor
        a3=sapply(a2[sapply(a2,length)>0],function(x) x%>%rename(p=pvalue, AE.type=AE.measurement.type)%>%filter((p<0.05)),simplify = F)
        
        a4=sapply(a3,dim)[1,]
        a31=a3[a4!=0]
        
        a5=sapply(a31,function(x)
        {
          x=x%>%mutate(group=paste(factor(p<0.05,level=c(T,F),label=c('Sig','NS')),factor(r>0,level=c(T,F),label=c('Pos_Cor','Neg_Cor')),sep='_'))
          tapply(as.vector(x$AE.type),x$group,sort)
        },simplify = F
        )
        
        fun.AE.summary<-function(x)
        {
          paste(gsub('NA','',gsub('grade3','High-Grade',gsub('grade12','Low-Grade',gsub('treatment\\.related','Trt',x)))),collapse='/')
        }
        
        fun.AE.paste<-function(x)
        {
          index1<-grep('treatment',x);
          if(length(index1)>0) x<-x[-index1]
          x
        }
        
        
        a6<-sapply(a5,function(y) as.vector(sapply(y,function(x) names(table(sub('\\.sum\\.unique','',sub('\\.sum_unique_AE','',sub('\\.occurrence','',sub('\\.fre','',sub('\\.duration','',x))))))),simplify = F)),simplify = F)
        
        a6.trt<-sapply(a6,function(y) unlist(sapply(y,function(x) {index1<-grep('treatment',x); if(length(index1)==0) NULL else
        {
          x<-x[index1];paste(x,collapse='/')
        }})))
        a6.trt<-a6.trt[!sapply(a6.trt,is.null)]
        a6.trt<-sapply(a6.trt,function(y) sapply(y,fun.AE.summary,simplify = F))
        
        a6.no_trt<-sapply(a6,function(y) unlist(sapply(y,function(x) {index1<-grep('treatment',x); if(length(index1)>0) x<-x[-index1];paste(x,collapse='/')})))
        a6.no_trt<-sapply(a6.no_trt,function(y) sapply(y,fun.AE.summary,simplify = F))
        
        a7<-sapply(a6,function(y) sapply(y,function(x) paste(sub('grade3','High-Grade',sub('grade12','Low-Grade',sub('treatment\\.related','Trt',x))),collapse='/')),simplify = F)
        
        name1<-c("Sig_Neg_Cor", "Sig_Pos_Cor")
        name2<-c("Negative Correlation","Positive Correlation")
        
        fun.group<-function(x,name.ref1=name1,name.ref2=name2)
        {
          name.tmp1<-names(x)
          index1<-name.ref1%in%name.tmp1
          if(sum(index1)<length(name.ref1)) {
            name.null<-rep('',length(name.ref1[!index1]))
            names(name.null)<-name.ref1[!index1]
            x<-c(x,name.null)
          }
          x<-x[match(name.ref1,names(x))]
          names(x)<-name.ref2
          x
        }
        
        table.AE.Response.sig.summary<-t(sapply(a7,function(x) fun.group(unlist(x))))
        
        q1=table.AE.Response.sig.summary
        q11=data.frame(AE=rownames(q1),q1)
        
        #---table of ind AE by AE category
        # print(names(AE_data()))
        # print(dim(AE_data()))
        #q2=table(AE_data()$cdus_ctcae_toxicity_type_code,AE_data()$toxicity_category)
        q2=table(AE_data()$cdus_toxicity_type_code,AE_data()$toxicity_category)
        
        #---get ind AE under AE category
        q3=q2[rownames(q2)%in%rownames(q1)[!rownames(q1)%in%c('whole',colnames(q2))],,drop=F]
        name1=rownames(q3)
        q31<-q3[,apply(q3,2,sum)!=0,drop=F]
        
        name2=colnames(q31)
        q5=apply(q31,1,function(x) name2[x>0])
        q6=data.frame(cate.AE=q5,AE=names(q5))
        q7=full_join(q6,q11)
        q7$cate.AE[(1:dim(q7)[1])[is.na(q7$cate.AE)]]<-q7$AE[(1:dim(q7)[1])[is.na(q7$cate.AE)]]
        q8=q7[order(q7$cate.AE),]
        q8$AE[q8$AE==q8$cate.AE]<-NA
        q8=q8[order(q8$cate.AE,q8$AE,na.last = F),]
        index1<-grep('whole',q8$cate.AE)
        if(length(index1)>0) q8<-q8[c(index1,(1:dim(q8)[1])[-index1]),]
        q8$AE[is.na(q8$AE)]<-''
        #write.csv(q8, file = "../manuscript/MCC18494_AE_treatment_duration_summary.csv", row.names = FALSE)
        #write.csv(q8, file = "../manuscript/MCC18321_AE_treatment_duration_summary.csv", row.names = FALSE)
        
        #q8
    list(q8=q8,plot=cor.AE.treatment.time.ans)
        
      })
  # Displaydurationanalysis table ####
      output$durationanalysistableoutput <- DT::renderDataTable({
        todisplay <- durationanalysis()$q8
        datatable(todisplay, rownames = FALSE)
        
      }, options = list(pageLength = 72))
      
     
      
      # Cor plots ###
      draw_plot_duration <- function() {  durationanalysis()$plot  }
      
      # render correlation plots ###
      output$plotduration <- renderPlot({  draw_plot_duration() })
      
      # Download correlation plots pdf ###
      output$downloaddurationplotspdf <- downloadHandler(
        filename = "duration_rendered_report.pdf",
        content = function(file) {
          res <- rmarkdown::render(
            "templateduration.Rmd",
            params = list(
              draw_plot_duration = draw_plot_duration
            )
          )
          file.rename(res, file)
        }
      )
  # end: duration plots        ^^^^^^^^^####
   
      
      
}) # end server ####
