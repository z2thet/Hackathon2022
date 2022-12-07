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
                tabsetPanel(
                    id = "mainnav", 
                    # tab to show raw data
                    tabPanel(
                        div(icon("database"), "Raw AE data"),
                         DT::dataTableOutput("selected_data")
                    ),  # tab to show raw data
                    tabPanel(
                        div(icon("database"), "Raw demographic data"),
                        DT::dataTableOutput("demographics_data")
                    ),  # tab to show raw data
                    tabPanel(
                      div(icon("database"), "Follow up data"),
                      DT::dataTableOutput("followup_data")
                    ), # tab to show raw data
                    tabPanel(
                      div(icon("database"), "Toxicity data"),
                      DT::dataTableOutput("tox_data")
                    ), # tab to show raw data
                    tabPanel(
                      div(icon("database"), "Drug administration data"),
                      DT::dataTableOutput("drug_data")
                    ),
                    
                     
                    tabPanel(
                        div(icon("chart-line"), "AE plot"),
                        div(
                            style = "position:relative",
                            selectizeInput("single_var_input", "Patient number:", choices = NULL,
                                           options = list(placeholder = "Select a patient",
                                                          onInitialize = I('function() { this.setValue(""); }'))),
                            br(),
                            downloadButton("download_single_var_plot", "Save plot"),
                            br(),
                            withSpinner(plotOutput("single_var_plot", width = "100%", height = "500px")),
                            br(),
                            downloadButton('download_single_var_plot_early',"Save plot"),
                            br(),
                            withSpinner(plotOutput("single_var_plot_early", width = "100%", height = "500px"))
                            
                            # br(), br(),
                            # verbatimTextOutput("single_var_summary")
                        )
                    ),
                   tabPanel(
                       div(icon("database"), "AE days data"),
                       div(downloadButton('AEDAYSdownload',"Download the data"),
                           br(),
                           br(),
                           br(),
                           withSpinner(  DT::dataTableOutput("selected_data_AEDAYS"))
                       )
                   ),
                   tabPanel(
                     div(icon("database"), "AE measures"),
                     div(downloadButton('toxicitymeasuresdownload',"Download the data"),
                         br(),
                         br(),
                         br(),
                         withSpinner(DT::dataTableOutput("toxicitytableoutput"))
                     )
                   ),
 
                   tabPanel( #https://stackoverflow.com/questions/46471756/download-pdf-report-in-shiny
                     div(icon("database"), "Coxph AE measures"),
                     div(downloadButton('downloadplotspdf',"Download plots"),
                         br(),
                         br(),
                         br(),
                         withSpinner( DT::dataTableOutput("kmcophinfo"))
                   )
                  ),
                  tabPanel( #https://stackoverflow.com/questions/46471756/download-pdf-report-in-shiny
                    div(icon("database"), "Forest plots OS and PFS"),
                    div(downloadButton('downloadforestplotspdf',"Download Forest plots"),
                        br(),
                        br(),
                        br() 
                        ,                        
                        withSpinner(DT::dataTableOutput("forestplotcophinfo"))
                    )
                  ),
                  tabPanel( #https://stackoverflow.com/questions/46471756/download-pdf-report-in-shiny
                    div(icon("database"), "Response tests"),
                     div(
                       #downloadButton('downloadforestplotspdf',"Download Response plots"),
                        br(),
                        br(),
                        br() 
                        ,                        
                        withSpinner(DT::dataTableOutput("responsettestoutput"))
                    )
                  ),
                  
                   
                  tabPanel(
                    div(icon("database"), "Correlation"),
                     div(
                        downloadButton('downloaddurationplotspdf',"Download duration plots"),
                         br(),
                        br(),
                        br(),
                        withSpinner(DT::dataTableOutput("durationanalysistableoutput"))
                    )
                  ),
                  
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
                ),
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

