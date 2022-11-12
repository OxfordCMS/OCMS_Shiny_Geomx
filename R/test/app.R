
library(shiny)
library(plotly)
library(ggplot2)
library(GeomxTools)
library(shinythemes)
library(shinycssloaders)
library(tidyr)
library(formattable)

# TEMPORARY DATA UNTIL READING IN FILE WORKING
geomx_dat <- readRDS("/well/powrie/users/utu691/devel/shiny_geomx/data/hca_geomx.rds")
dat <- as.data.frame(sData(geomx_dat))
dat$Sample_ID <- rownames(dat)

# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("cosmo"),
                
                # Application title
                titlePanel("shiny geomx"),
                
                tabsetPanel(
                  tabPanel("Overview", fluid = TRUE,
                           sidebarLayout(
                             sidebarPanel("Import geomx .rds file", position="left",
                                          fileInput("rds", "", multiple=FALSE)
                                          
                             ),
                             mainPanel(h6("The purpose of this app is to provide an interface in order to assess\n
                         the quality of geomx spatial transcriptomic data. It allows you to visualise\n
                         the data in terms of segments that PASS/FAIL QC given the user-specified\n
                         parameters. The separate tabs represent different processing/analysis\n
                         types."),
                                       tableOutput("nsegments") %>% withSpinner(),
                                       textOutput("read_counts_title"),
                                       plotlyOutput("read_counts") %>% withSpinner()
                             )
                           )
                  ),
                  tabPanel("Select samples", fluid = TRUE,
                           sidebarLayout(
                             sidebarPanel("Sample selection", position="left",
                             ),
                             mainPanel(h6("Select samples from table below"),
                                       DT::dataTableOutput('selected_samples') %>% withSpinner(),
                             ),
                           )
                  ),
                  tabPanel("Segment QC", fluid = TRUE,
                           sidebarLayout(
                             sidebarPanel("Segment QC paramters", position="left",
                                          numericInput("min_segment_reads", "Min. segment reads", 1000),
                                          numericInput("percent_trimmed", "Percent trimmed", 80),
                                          numericInput("percent_stitched", "Percent stitched", 80),
                                          numericInput("percent_aligned", "Percent aligned", 75),
                                          numericInput("percent_saturation", "Percent saturation", 50),
                                          numericInput("min_negative_count", "Min. negative count", 1),
                                          numericInput("max_ntcc_count", "Max. NTC count", 9000),
                                          numericInput("min_nuclei", "Min. nuclei", 20),
                                          numericInput("min_area", "Min. area", 1000),
                             ),
                             mainPanel(h6("Samples that are selected are displayed below:"),
                                       DT::dataTableOutput("samples_selected"),
                                       DT::dataTableOutput("segment_qc_table"),
                                       formattableOutput("segment_pass_table")
                             ),
                           )
                  ),
                  tabPanel("Probe QC", fluid = TRUE,
                           sidebarLayout(
                             sidebarPanel("Probe QC", position="left",
                                          radioButtons("local_removal", "Perform local removal", choices = c("Yes", "No"))
                             ),
                             mainPanel(h6("Probe QC is performed on negative control probes. They are used to determine
                                the background in subsequent processing and so at this point if they are outliers
                                they are removed. You can perform local removal i.e. at the level of segment if you
                                wish and can select in the sidebar. By default only global removal is performed."),
                                       h6("Samples that are selected are displayed below:"),
                                       #DT::dataTableOutput("samples_selected"),
                                       DT::dataTableOutput("probe_qc_table")
                             ),
                           )
                  )
                )
)

# Define server logic required to draw a histogram
options(shiny.maxRequestSize=50*1024^2)
server <- function(input, output) {
  
  ###############################################
  # Summary stats  
  ###############################################
  # Summary table
  summary_df <- reactive({
    dat <- as.data.frame(sData(geomx_dat))
    summary_data <- data.frame("nsegments" = nrow(dat),
                               "ave_read_count" = round(mean(dat$Raw), digits = 0),
                               "ave_deduplicated_reads" = round(mean(dat$DeduplicatedReads), digits=0))
    summary_data
  })
  
  output$nsegments <- renderTable({
    summary_df()
  })
}

# Run the application 
shinyApp(ui = ui, server = server)