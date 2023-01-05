library(shiny)
library(shinycssloaders)
library(ggplot2)
library(magrittr)
library(DT)
library(shinythemes)

# Increase uploaded file size to 30MB
options(shiny.maxRequestSize = 50*1024^2)

function(request){
  fluidPage(
    title = "AE plotting and calculation of days with AEs",
    # theme = shinytheme("spacelab"),
    # shinythemes::themeSelector(), 
    shinyjs::useShinyjs(),
    tags$head(
      tags$link(href = "style.css", rel = "stylesheet")
    ),
    div(id = "header",
        div(id = "title",
            "Analysis of Adverse Events"
        )
        # ,
        # div(id = "subtitle",
        #     "Adverse Events"),
        # div(id = "subsubtitle",
        #     "By:",
        #     tags$a(href = "mailto:zacha@moffitt.org", "Ram Thapa")
        # )
    ),
    fluidRow(
      column(
        width = 2, 
        class = "col-settings",
        
        fileInput("data_file", "Select AE file"),
        helpText("File size limited to 50 MB"),
        fileInput("demographics_file", "Select demographics file"),
        helpText("File size limited to 50 MB"),
        fileInput("fu_file", "Select followup file"),
        helpText("File size limited to 50 MB"),
        fileInput("da_file", "Select Drug administration file"),
        helpText("File size limited to 50 MB"),
        
        bookmarkButton()
      ),
      
      # tabs to show raw datasets ####
      column(
        width = 9, 
        class = "col-tabs",
        
   
        tabsetPanel( id = "panel", type = "tabs",
                     
                     tabPanel("Data", tabsetPanel(  
          tabPanel(
            div(icon("database"), "AE data"),
            DT::dataTableOutput("selected_data")
          ),  # tab to show raw data
          tabPanel(
            div(icon("database"), "Demographic data"),
            DT::dataTableOutput("demographics_data")
          ),  # tab to show raw data
          tabPanel(
            div(icon("database"), "Follow up data"),
            DT::dataTableOutput("followup_data")
          ), # tab to show raw data
          
          tabPanel(
            div(icon("database"), "Drug administration data"),
            DT::dataTableOutput("drug_data")
          ),
          tabPanel(
            div(icon("database"), "Toxicity data"),
            DT::dataTableOutput("tox_data")
          )
          
          )) #END tabPanel("Data", tabsetPanel
          ,  
     
  #  ,
           
        
        
  # AE measures and plots
  
  tabPanel("AE Plots and Measures", tabsetPanel(        
          # SWIMMERS PLOTS TAB ####
          tabPanel(
            div(icon("chart-line"), "AE plot"),
            div(
              style = "position:relative",
              fluidRow(
                column(2,
                       selectizeInput("single_var_input", "Patient number:", choices = NULL,
                                      options = list(placeholder = "Select a patient",
                                                     onInitialize = I('function() { this.setValue(""); }')))
                ),
                column(2,
                       numericInput("AEplot_EarlyAECut","Early AE Time Point:", value = NULL, step = 1, min = 1)
                ),
                column(3,
                       checkboxInput("AEplot_ShowTime","Display Time Annotation", value = T)
                )
              ),
 
              br(),
              downloadButton('download_single_var_plot_early',"Save plot"),
              br(),
              withSpinner(plotOutput("single_var_plot_early", width = "100%", height = "500px"))
              
            )
          ),
          
          
          # AE days ####S
          tabPanel(
            div(icon("database"), "AE table"),
            div(#downloadButton('AEDAYSdownload',"Download the data"),
                br(),
                br(),
                br(),
                withSpinner(  DT::dataTableOutput("AEcounttables"))
            )
          ),
          
          
          
          # AE days ####S
          tabPanel(
            div(icon("database"), "AE days data"),
            div(downloadButton('AEDAYSdownload',"Download the data"),
                br(),
                br(),
                br(),
                withSpinner(  DT::dataTableOutput("selected_data_AEDAYS"))
            )
          ),
          
          
          
 
          
          # AE measures tab ####
           tabPanel(
             div(icon("database"), "AE measures"),
             fluidRow(
               column(2,
                      div(downloadButton('toxicitymeasuresdownload',"Download the data"))
               ),
               column(2,
                      numericInput("AEmeasuresEarlyAEcut","Early AE Time Point:", value = NULL, step = 1, min = 1)
               ),
               column(3,
                      uiOutput("rendAEmeasuresAEcatselect")
               ),
               column(3,
                      uiOutput("rendAEmeasuresAEtypeselect")
               )
             ),
             #div(downloadButton('toxicitymeasuresdownload',"Download the data"),
             br(),
             br(),
             br(),
             withSpinner(DT::dataTableOutput("toxicitytableoutput"))
             #)
           )
          
          )) # END NEXT  tabPanel("AE Plots and Measures", tabsetPanel( 
  ,
        
          
  tabPanel("Survival Analysis", tabsetPanel( 
          tabPanel( #https://stackoverflow.com/questions/46471756/download-pdf-report-in-shiny
            div(icon( "fa-solid fa-code-branch" ), "Coxph AE measures"),  
            div(
              actionButton("gogocoxmodels", "Run Coxph models" ),
              br()),
            fluidRow(
              column(2,
                     div(downloadButton('downloadplotspdf',"Download plots")),
                     ),
              column(4,
                     uiOutput("rendBPKPplotSelection")
                     )
            ),
            withSpinner(plotOutput("BPKP_PlotOutput", width = "100%", height = "450px"), type = 6),
            #div(downloadButton('downloadplotspdf',"Download plots"),
                br(),
                br(),
                br(),
                withSpinner( DT::dataTableOutput("kmcophinfo"))
            #)
          ),
          #Forest plots OS and PFS tab ####
          tabPanel( #https://stackoverflow.com/questions/46471756/download-pdf-report-in-shiny
            div(icon("database"), "Forest plots OS and PFS"),
            div(
              actionButton("goforstplotstests", "Run tests" ),
              br()),
            fluidRow(
              column(2,
                     div(downloadButton('downloadforestplotspdf',"Download Forest plots")),
              ),
              column(4,
                     uiOutput("rendForestplotSelection")
              )
            ),
            withSpinner(plotOutput("Forest_PlotOutput", width = "100%", height = "450px"), type = 6),
            #div(downloadButton('downloadforestplotspdf',"Download Forest plots"),
                br(),
                br(),
                br() 
                ,                        
                withSpinner(DT::dataTableOutput("forestplotcophinfo"))
            #)
          )
        
  )) #End   tabPanel("Survival Analysis", tabsetPanel( 
          ,
        
  tabPanel("Response and Correlations", tabsetPanel( 
        
          # Response tests tab ####
          tabPanel( #https://stackoverflow.com/questions/46471756/download-pdf-report-in-shiny
            div(icon("database"), "Response tests"),
            div(
              actionButton("goresponsetests", "Run tests" ),
              br(),
              downloadButton('download_responseplots',"Download Response plots"),
              br(),
              br(),
              br() 
              ,                        
              withSpinner(DT::dataTableOutput("responsettestoutput"))
            )
          ),
          
          # Correlation tab ####
          tabPanel(
            div(icon("database"), "Correlation"),
            div(
              actionButton("goCorrelationtests", "Run tests" ),
              br(),
              downloadButton('downloaddurationplotspdf',"Download duration plots"),
              br(),
              br(),
              br(),
              withSpinner(DT::dataTableOutput("durationanalysistableoutput"))
            )
          )
  )) # END   tabPanel("Response and Correlations", tabsetPanel( 
        ,
  tabPanel("Tables and Reports", tabsetPanel(  
         # Survival summary tab #### 
          tabPanel(
            div(icon("database"), "Survival summary"),
            div(
              p("Please run the 'Coxph AE measures' and 'Forest plots OS and PFS tabs' before downloading the survival summary report." ),
              br(),
              downloadButton('survival_summarydownload',"Download Survival summary"),
              br(),
              br(),
              br(),
              withSpinner(DT::dataTableOutput("survivalsummarytable"))
            )
          ),
         # Response summary tab #### 
          tabPanel(
            div(icon("database"), "Response summary"),
            div( p("Please run the 'Response tests tab' before attempting to download the 'Response summary report'." ),
                 br(),
               downloadButton('response_summarydownload',"Download Response summary"),
              br(),
              br(),
              br(),
              withSpinner(DT::dataTableOutput("responsesummarytable"))
            )
          ),
          
         # Summary Report download tab ####
          tabPanel(
            div(icon("database"), "Summary Report"),
            div(
              downloadButton('downloadsumaryReport',"Download summary report"),
              br()  
              
            )
          )
  ))
         #End  tabPanel("Tables and Reports", tabsetPanel(  
         ,
          # tab to documentation/feedback
          tabPanel(
            div(icon("info-circle"), "Documentation"),
            # p("To use this app upload a file that has adverse event data stored in table format with first row as a header.
            #   Adverse event data from Encore must be merged with the on treatment date typically found in the demographics data set in Encore.
            #   The file can be of any format i.e. excel (.xlsx, .xls), csv, text (.txt), SAS (.sas7bdat), SPSS (.sav), Stata (.dta) etc. 
            #   and it should contain only one table (for example, a table in the first sheet of an Excel file)."),
            # hr(),
            # # p("Numerical variable with unique values less than or equal to 5 is treated as factor variable. Data type conversion is not available in the app for now."),
            # # hr(),
            # p(" Variables needed (please ensure variables are named exactly as below):    "),
            # tags$div(tags$ul(
            #     tags$li(tags$span("sequence_no")),
            #     tags$li(tags$span("onset_date_of_ae")),
            #     tags$li(tags$span("cdus_ctcae_toxicity_type_code")),
            #     tags$li(tags$span("resolved_date")),
            #     tags$li(tags$span("on_treatment_date")), 
            #     tags$li(tags$span("grade")),
            #     tags$li(tags$span("attribution_possible")),
            #     tags$li(tags$span("attribution_probable")),
            #     tags$li(tags$span("attribution_definite")) 
            #                 ) 
            #          ),
            # hr(),
            
            p("This shiny application was developed by xxx yyy at the Biostatistics Core. Please let them know if you find any errors or have any questions", HTML(paste0(a(href = "mailto:ram.thapa@moffitt.org", "by email."))))
          ) 
        
        ),#END ORigiNAL TABSETPANEL
        hr(),
        p(HTML('<a href = "https://intranet.moffitt.org/display/Biostats/Biostatistics+and+Bioinformatics+Home" target="_blank"> <img src="biostatistics.jpg" align = "left" height = "70px" width = "448"/></a>')),
        p(HTML('<a href = "https://www.moffitt.org/" target="_blank"> <img src = "MOFFITT_2c_RGB.jpg" align = "right" height = "80px" width = "320px"/></a>')),
        tags$br(),
        tags$br(),
        tags$br(),
        tags$br(),
        tags$br(),
        tags$br()
        
      
    )
  )
  ) 
}

