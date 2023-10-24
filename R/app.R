library(shiny)
library(plotly)
library(ggplot2)
library(GeomxTools)
library(shinythemes)
library(shinycssloaders)
library(shinyBS)
library(tidyr)
library(formattable)
library(ComplexHeatmap)
library(preprocessCore)
library(ggrepel)
library(shinyBS)

# TEMPORARY DATA UNTIL READING IN FILE WORKING
#geomx_dat <- readRDS("~/devel/shiny_geomx/data/hca_geomx.rds")
#geomx_dat@analyte <- "RNA"
#dat <- as.data.frame(sData(geomx_dat))
#dat$Sample_ID <- rownames(dat)

# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("cosmo"),

        
#    tags$style(type="text/css",
#               ".shiny-output-error { visibility: hidden; }",
#               ".shiny-output-error:before { visibility: hidden; }"
#    ),
    
    # Application title
    titlePanel("shiny geomx"),
                
    tabsetPanel(
        tabPanel("Overview", fluid = TRUE,
            sidebarLayout(
                sidebarPanel("Import geomx .rds file", position="left",
                              div(id="input_div",
                              fileInput("rds", "", multiple=FALSE, accept=c(".rds", ".Rds", ".RDS")),
                              bsTooltip("input_div", "Input file that is .Rds file of stored NanoStringGeoMxSet", placement = "bottom", trigger = "hover",
                                       options = NULL)
                              )
                              
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
                   h6("To include all samples do not select anything. Otherwise select required samples")
                   ),
                   mainPanel(h6("Select samples from table below"),
                             DT::dataTableOutput('selected_samples') %>% withSpinner(),
                             plotlyOutput("scatter_count_area") %>% withSpinner(),
                             plotlyOutput("scatter_count_area_log") %>% withSpinner(),
                   ),
                )
        ),
        tabPanel("Segment QC", fluid = TRUE,
                 sidebarLayout(
                   sidebarPanel("Segment QC paramters", position="left",
                                div(id="segment_reads_div",
                                numericInput("min_segment_reads", "Min. segment reads", 1000),
                                bsTooltip("segment_reads_div",
                                          "Segments will pass if they have greater than this number of reads",
                                          placement = "bottom",
                                          trigger = "hover",
                                          options = NULL)
                                ),
                                div(id="percent_trimmed_div",
                                numericInput("percent_trimmed", "Percent trimmed", 80),
                                bsTooltip("percent_trimmed_div",
                                          "Segments will pass if more than this % passed post-trimming",
                                          placement = "bottom",
                                          trigger = "hover",
                                          options = NULL)
                                ),
                                div(id="percent_stitched_div",
                                numericInput("percent_stitched", "Percent stitched", 80),
                                bsTooltip("percent_stitched_div",
                                          "Segments will pass if more than this % passed post-stitching",
                                          placement = "bottom",
                                          trigger = "hover",
                                          options = NULL)
                                ),
                                div(id="percent_aligned_div",                                
                                numericInput("percent_aligned", "Percent aligned", 75),
                                bsTooltip("percent_aligned_div",
                                          "Segments will pass if more than this % are successfully aligned",
                                          placement = "bottom",
                                          trigger = "hover",
                                          options = NULL)
                                ),
                                div(id="percent_saturation_div",
                                numericInput("percent_saturation", "Percent saturation", 50),
                                bsTooltip("percent_saturation_div",
                                          "% Sequencing saturation ([1-deduplicated reads/aligned reads]%). <50% suggests more sequencing required",
                                          placement = "bottom",
                                          trigger = "hover",
                                          options = NULL)
                                ),
                                div(id="min_neg_count_div",
                                numericInput("min_negative_count", "Min. negative count", 1),
                                bsTooltip("min_neg_count_div",
                                          "This is the geometric mean of the several unique negative probes in the GeoMx panel that do not target mRNA and establish the background count level per segment",
                                          placement = "bottom",
                                          trigger = "hover",
                                          options = NULL)
                                ),
                                div(id="max_ntcc_count_div",
                                numericInput("max_ntcc_count", "Max. NTC count", 9000),
                                bsTooltip("max_ntcc_count_div",
                                          "No Template Control (NTC) count: values >1,000 could indicate contamination for the segments associated with this NTC; however, in cases where the NTC count is between 1,000- 10,000, the segments may be used if the NTC data is uniformly low (e.g. 0-2 counts for all probes)",
                                          placement = "bottom",
                                          trigger = "hover",
                                          options = NULL)
                                ),
                                div(id="min_nuclei_div",
                                numericInput("min_nuclei", "Min. nuclei", 20),
                                bsTooltip("min_nuclei_div",
                                          ">100 nuclei per segment is generally recommended; however, this cutoff is highly study/tissue dependent and may need to be reduced; what is most important is consistency in the nuclei distribution for segments within the study",
                                          placement = "bottom",
                                          trigger = "hover",
                                          options = NULL)
                                ),                               
                                div(id="min_area_div",
                                numericInput("min_area", "Min. area", 1000),
                                bsTooltip("min_area_div",
                                          "Area generally correlates with number of nuclei. Not neccessary to have a strict cutoff",
                                          placement = "bottom",
                                          trigger = "hover",
                                          options = NULL)
                                ),    
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
                                div(id="probe_qc_div",
                                radioButtons("local_removal", "Perform local removal",
                                             choices = c("Yes", "No"),
                                             selected="No"),
                                bsTooltip("probe_qc_div",
                                          "Local removal will remove probes within segments in which they perform poorly and leave others untouched. WARNING: Not implemented",
                                           placement = "bottom",
                                           trigger = "hover",
                                           options = NULL)
                                ),
                   ),
                   mainPanel(h6("Probe QC is performed on negative control probes. They are used to determine
                                the background in subsequent processing and so at this point if they are outliers
                                they are removed. You can perform local removal i.e. at the level of segment if you
                                wish and can select in the sidebar. By default only global removal is performed."),
                             h6("Segments that are selected are displayed below:"),
                             DT::dataTableOutput("segments_for_qc"),
                             DT::dataTableOutput("probe_qc_table"),
                             DT::dataTableOutput("filtered_geomx_summary")
                   ),
                )
        ),
        tabPanel("Gene-level QC", fluid = TRUE,
                 sidebarLayout(
                   sidebarPanel("Limit of quantification", position="left",
                                
                                div(id="std_div",
                                numericInput("n_std_deviations", "Number of SDs above neg mean", 2),
                                bsTooltip("std_div",
                                          "Number of standard deviations above the negative control probes that amn endogenous probe has to be in order to be called as expressed",
                                          placement = "bottom",
                                          trigger = "hover",
                                          options = NULL)
                                ),
                                div(id="min_cutoff_div",
                                numericInput("min_cutoff", "Minimum cutoff", 2),
                                bsTooltip("min_cutoff_div",
                                          "This is a bit of a catch for probes that pass the Std cutoff but are still too low to be considered useful. Leave as 2",
                                          placement = "bottom",
                                          trigger = "hover",
                                          options = NULL)
                                ),
                                div(id="filter_threshold_div",
                                numericInput("filter_threshold", "% gene detection cut-off to remove segment", 10),
                                bsTooltip("filter_threshold_div",
                                          "Based on the above parameters what % of genes should represent a high quality segment? Segments below this threshold are removed",
                                          placement = "bottom",
                                          trigger = "hover",
                                          options = NULL)
                                ),      
                                actionButton("filter_segments", "Filter segments"),
                                
                                div(id="goi_div",
                                textAreaInput("goi_list", "Gene detection rates (paste gene list)"),
                                bsTooltip("goi_div",
                                          "Type or paste a list of genes (symbols) to retrieve the % of segments they are detected in. Must be a sinlge column list with no seperators",
                                          placement = "bottom",
                                          trigger = "hover",
                                          options = NULL)
                                ),
                                actionButton("search_genes", "Search"),
                                div(id="gene_filter_div",
                                numericInput("gene_filter_threshold", "Keep gene if detected in more than this % of segments", 10),
                                bsTooltip("gene_filter_div",
                                          "Filter genes out of the dataset if they are not present in at least this % of segments",
                                          placement = "bottom",
                                          trigger = "hover",
                                          options = NULL)
                                ),
                                actionButton("filter_genes", "Filter genes")
                   ),
                   mainPanel(h6("The geometric mean is taken for probes that belong to the same gene and a single
                                counts matrix is produced"),
                             h6("Gene Detection. The proportion of genes detected per segment is displayed"),
                             plotOutput("gene_detection") %>% withSpinner(),
                             h6("Explore genes for detection rates across segments"),
                             DT::dataTableOutput("filtered_stats") %>% withSpinner(),
                             DT::dataTableOutput("gois") %>% withSpinner(),
                             plotOutput("gois_heatmap"),
                             DT::dataTableOutput("filtered_segments_genes_summary")
                             
                   ),
                 )
        ),
        tabPanel("Normalization", fluid = TRUE,
                 sidebarLayout(
                   sidebarPanel("Normalization", position="left",
                                div(id="normalisation_div",
                                radioButtons("norm_method",
                                             "Choose method",
                                             choices=c("Quantile"),
                                             selected="Quantile"),
                                bsTooltip("normalisation_div",
                                          "Choose normalisation method",
                                          placement = "bottom",
                                          trigger = "hover",
                                          options = NULL)
                                ),
                                actionButton("normalize", "Normalize"),
                                textAreaInput("genes_to_view", "Gene expression"),
                                actionButton("search_genes_to_view", "Search"),
                                h6("Downloads"),
                                downloadButton("download_rds", "Download .Rds file"),
                                h6(""),
                                downloadButton("download_normalized", "Download normalized (.tsv) file"),
                                h6(""),
                                downloadButton("download_metadata", "Download metadata (.tsv), file")
                                
                   ),
                   mainPanel(h6("Normalization across segments is performed according to user-specified method"),
                             h6("Pre-normalization gene expression distribution across samples"),
                             plotOutput("prenorm_distribution") %>% withSpinner(),
                             h6("After hitting normalize, gene expression distribution appears below"),
                             plotOutput("postnorm_distribution") %>% withSpinner(),
                             plotOutput("neg_vs_pos_plot") %>% withSpinner(),
                             h6("PCA coloured by segment - post-normalization"),
                             plotlyOutput("postnorm_pca") %>% withSpinner(),
                             h6("Expression of selected genes of interest are shown in a heatmap below"),
                             bsModal("plotModal1", "Heatmap genes of interest", "search_genes_to_view", size = "large", plotOutput("genes_to_view_heatmap", inline = TRUE))
                                                          
                   ),
                   )
        )
    )
)

# Define server logic required to draw a histogram
options(shiny.maxRequestSize=50*1024^2)
server <- function(input, output) {
  
    # input .rds file
    geomx_dat <- reactive({
      filename <- input$rds
      if (is.null(filename)){
        return()
      }
      geomx_dat <- readRDS(filename$datapath)
      geomx_dat@analyte <- "RNA"
      colnames(pData(geomx_dat)) <- tolower(colnames(pData(geomx_dat)))
      pData(geomx_dat)$file <- rownames(pData(geomx_dat))
      geomx_dat
    })

      
    ###############################################
    # Summary stats  
    ###############################################
    # Summary table
    summary_df <- reactive({
        if (is.null(input$rds)){
          return()
        }
        dat <- as.data.frame(sData(geomx_dat()))
        dat$Sample_ID <- rownames(dat)
        dat$Sample_ID <- rownames(dat)
        summary_data <- data.frame("nsegments" = nrow(dat),
                                   "ave_read_count" = round(mean(dat$Raw), digits = 0),
                                   "ave_deduplicated_reads" = round(mean(dat$DeduplicatedReads), digits=0))
        summary_data
    })

    output$nsegments <- renderTable({
      summary_df()
    })

    output$read_counts_title <- renderText({
      "Summary of read counts before and after de-duplication"
    })
    
    # read summary plot
    summary_plot <- reactive({
        if (is.null(input$rds)){
          return()
        }
        dat <- as.data.frame(sData(geomx_dat()))
        dat$Sample_ID <- rownames(dat)
        dat.m <- dat[,c("Sample_ID", "Raw", "DeduplicatedReads")]
        dat.m <- tidyr::pivot_longer(!Sample_ID,
                                     values_to="count",
                                     data = dat.m)
                
        dat.m <- dat.m[order(dat.m$count, decreasing = TRUE),]
        dat.m$Sample_ID <- factor(dat.m$Sample_ID, levels=unique(dat.m$Sample_ID))        
        
        plt <- ggplot(dat.m, aes(x=Sample_ID, y=count)) +
          geom_point()  +
          coord_flip() +
          facet_wrap(~name) +
          theme_classic() +
          geom_hline(yintercept=10000, colour="red") +
          theme(axis.text.y = element_blank())
        plotly::ggplotly(plt)
    })
    
    output$read_counts <- renderPlotly({
      summary_plot()
    })

    ###############################################
    # Sample selection
    ###############################################
    sample_table_select <- reactive({
      pData(geomx_dat()) %>% select(file, slide.name, roi, segment)
    })
    
    output$selected_samples <- DT::renderDataTable(server = FALSE, {
        DT::datatable(sample_table_select(), filter = 'top',
                    rownames = FALSE,
                    options = list(
                    pageLength = 15,
                    scrollX = TRUE,
                    dom = 'lfrtip',
                    filter = 'top'))
    })    

    scatter_area_counts <- function(dat, logT=FALSE){
      dat <- sData(dat)
      dat$label <- paste0(dat$file, ":", dat$slide, ":", dat$roi, ":", dat$segment)
      if (logT == TRUE){
        dat$count <- log10(dat$DeduplicatedReads)
        dat$area_ <- log10(dat$area)
        labx <- "log10(Area)"
        laby <- "log10(Count)"
        }
      else{
          dat$count <- dat$DeduplicatedReads
          dat$area_ <- dat$area
          labx <- "Area"
          laby <- "Count"
      }
      p <- ggplot(dat, aes(x=area_, y=count, color=label)) +
        geom_point() +
        theme_classic() +
        theme(legend.position="none") +
        xlab(labx) +
        ylab(laby)
      return(plotly::ggplotly(p))
    }
    
    scatterplot_area_counts_raw <- reactive({
      scatter_area_counts(geomx_dat())
    })
    scatterplot_area_counts_log <- reactive({
      scatter_area_counts(geomx_dat(), logT=TRUE)
    })
    
    output$scatter_count_area <- renderPlotly({
      scatterplot_area_counts_raw()
    })
    
    output$scatter_count_area_log <- renderPlotly({
      scatterplot_area_counts_log()
    })
    
    
    
    ###############################################
    # Segment QC
    ###############################################
    output$samples_selected = DT::renderDataTable({
      s = input$selected_samples_rows_selected
      if (length(s)) {
        sample_table_select()[s,]}
      else{
        sample_table_select()
      }
    })
    
    # filter the data - if nothing is selected in the table then take all
      geomx_dat_filtered <- reactive({
      to_keep <- input$selected_samples_rows_selected
      if (is.null(to_keep)){
        geomx_dat_filtered <- geomx_dat()}
      else{
          geomx_dat_filtered <- geomx_dat()[,to_keep]
      }
      geomx_dat_filtered
    })
    
    run_segment_qc <- reactive({
      dat <- geomx_dat_filtered()
      qc_params <- list(
                        minSegmentReads = input$min_segment_reads,
                        percentTrimmed = input$percent_trimmed,
                        percentStitched = input$percent_stitched,
                        percentAligned = input$percent_aligned,
                        percentSaturation = input$percent_saturation,
                        minNegativeCount = input$min_negative_count,
                        maxNTCCount = input$max_ntcc_count,
                        minNuclei = input$min_nuclei,
                        minArea = input$min_area)
      
      qc_dat <- setSegmentQCFlags(dat, qcCutoffs = qc_params)
      
      
            
      # assign no template control counts to each segment based on the run NTC count
      qc_sdata <- sData(qc_dat)
      ntcs <- rownames(pData(qc_dat)[grep("No Template Control", pData(qc_dat)$slide.name),])
      ntc_counts <- sData(qc_dat)[ntcs,]$Aligned
      ntc_counts <- data.frame(NTC = ntcs,
                               Aligned = ntc_counts)
      ntc_counts$Run <- gsub("-[A-Z]-[A-Z].*.dcc", "", ntcs)
      rownames(ntc_counts) <- ntc_counts$Run
      
      # flag function
      flag_ntc <- function(x, ntc_counts, cutoff=input$max_ntcc_count){
        run <- gsub("-[A-Z]-[A-Z].*.dcc", "", x)
        count <- ntc_counts[run,]$Aligned
        if (is.na(count)){
          count <- 0
        }
        if (count >= cutoff){
          return(TRUE)}
        else{
          return(FALSE)
        }
      }
      protocolData(qc_dat)[["QCFlags"]]$HighNTC <- unlist(lapply(rownames(protocolData(qc_dat)[["QCFlags"]]),
                                                                 function(x) flag_ntc(x, ntc_counts)))
      qc_results <- protocolData(qc_dat)[["QCFlags"]]
      flag_columns <- colnames(qc_results)
      
      qc_summary <- data.frame(Pass = colSums(!qc_results[, flag_columns]),
                               Warning = colSums(qc_results[, flag_columns]))
      qc_results$QC_Status <- apply(qc_results, 1L, function(x) {
        ifelse(sum(x) == 0L, "PASS", "WARNING")
        })
            
      qc_summary["TOTAL FLAGS", ] <-
        c(sum(qc_results[, "QC_Status"] == "PASS", na.rm=TRUE),
          sum(qc_results[, "QC_Status"] == "WARNING", na.rm=TRUE))
      return(list(qc_summary, qc_results))
                  
    })
    
    output$segment_qc_table <- DT::renderDataTable({
      run_segment_qc()[[1]]
    })
    
    output$segment_pass_table <- renderFormattable({
      pass_warn <- formatter("span",
                              style = x~style("background-color" = ifelse(x == "PASS",
                                                            "lightgreen",
                                                            "orange")))
      qc_results <- run_segment_qc()[[2]]
      all_dat <- sData(geomx_dat_filtered())
      pdata <- pData(geomx_dat_filtered())
      qc_results <- data.frame(Sample_ID = rownames(qc_results),
                               Slide = pdata[rownames(qc_results),]$slide.name,
                               ROI = pdata[rownames(qc_results),]$roi,
                               QC_Status = qc_results$QC_Status,
                               Raw = all_dat[rownames(qc_results),]$Raw,
                               Trimmed = all_dat[rownames(qc_results),]$Trimmed,
                               Stitched = all_dat[rownames(qc_results),]$Stitched,
                               Aligned = all_dat[rownames(qc_results),]$Aligned,
                               Area = all_dat[rownames(qc_results),]$area)
      
      formattable(qc_results, list(QC_Status = pass_warn))
      
    })
    ###############################################
    # Probe QC
    ###############################################

    output$segments_for_qc <- DT::renderDataTable({
      s = input$selected_samples_rows_selected
      if (length(s)) {
        sample_table_select()[s,]}
      else{
        sample_table_select()
      }
    })
    
    # probe QC parameters are hardcoded
    # Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
    # FALSE if you do not want to remove local outliers
    
    probe_qc_summary <- reactive({
    
      if (input$local_removal == "Yes"){
        local_removal = TRUE
      }
      else{
        local_removal = FALSE  
      }

      geomx_dat_probe <- setBioProbeQCFlags(geomx_dat_filtered(), 
                                   qcCutoffs = list(minProbeRatio = 0.1,
                                                    percentFailGrubbs = 20), 
                                   removeLocalOutliers = local_removal)
    
      ProbeQCResults <- fData(geomx_dat_probe)[["QCFlags"]]
    
      # Define QC table for Probe QC
      qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                          Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                          Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                    & !ProbeQCResults$GlobalGrubbsOutlier))

      # show what is left after removing based on global
      # and local removal
      if (local_removal == FALSE){
        ProbeQCPassed <- subset(geomx_dat_probe,
                                fData(geomx_dat_probe)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
                                fData(geomx_dat_probe)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
      }
      else{
        ProbeQCPassed <- geomx_dat_probe
      }
      geomx_filtered_segments_probes <- ProbeQCPassed
      qc_removed <- as.data.frame(dim(ProbeQCPassed))
      colnames(qc_removed) <- c("Total")
      
      return(list(qc_df, qc_removed, geomx_filtered_segments_probes))
    })  

    output$probe_qc_table <- DT::renderDataTable({
      probe_qc_summary()[[1]]
    })
    
    output$filtered_geomx_summary <- DT::renderDataTable({
      probe_qc_summary()[[2]]
    })

###############################################
# Limits of quantification and aggregation    
###############################################
    
# first perform the aggregation
aggregate_genes <- reactive({
  
  target_data <- probe_qc_summary()[[3]]
  target_data <- aggregateCounts(target_data)
  target_data
})

# limit of quantification    
loq <- reactive({

  # Define LOQ SD threshold and minimum value
  cutoff <- input$n_std_deviations
  min_loq <- input$min_cutoff

  # aggregage
  dat_aggregated <- aggregate_genes()
  
  # Calculate LOQ per module tested
  LOQ <- data.frame(row.names = colnames(dat_aggregated))
  module <- fData(dat_aggregated)$Module
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"), module)
  if(all(vars[1:2] %in% colnames(pData(dat_aggregated)))) {
      LOQ[, module] <- pmax(min_loq, pData(dat_aggregated)[, vars[1]] * pData(dat_aggregated)[, vars[2]] ^ cutoff)
  }
  pData(dat_aggregated)$LOQ <- LOQ
  dat_aggregated
})    
    
# Filtering based on LOQ  
filter_loq <- reactive({
  
  dat_loq <- loq()
    
  module <- fData(dat_loq)$Module
  loq_mat <- t(esApply(dat_loq, MARGIN = 1,
                       FUN = function(x) { x > pData(dat_loq)$LOQ[,1]}))
  
  # ensure ordering since this is stored outside of the geomxSet
  loq_mat <- loq_mat[fData(dat_loq)$TargetName, ]
  
  # Save detection rate information to pheno data
  pData(dat_loq)$GenesDetected <- colSums(loq_mat, na.rm = TRUE)
  pData(dat_loq)$GeneDetectionRate <- pData(dat_loq)$GenesDetected / nrow(dat_loq)
  
  # Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
  pData(dat_loq)$DetectionThreshold <- 
    cut(pData(dat_loq)$GeneDetectionRate,
        breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
        labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))
  dat_loq
})
  
filtered_summary <- reactive({
  dat <- filter_loq()
  table(pData(dat)$DetectionThreshold,
        pData(dat)$segment)
})

gene_detection_plot <- reactive({
  
  dat <- filter_loq()
  
  # stacked bar plot of different cut points (1%, 5%, 10%, 15%)
  ggplot(pData(dat), aes(x = DetectionThreshold)) +
    geom_bar(aes(fill = segment)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_classic() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Gene Detection Rate",
         y = "Segments, #",
         fill = "Segment Type")
})

output$gene_detection <- renderPlot({
  gene_detection_plot()
})

output$filtered_summary <- DT::renderDataTable({
  filtered_summary()
})
  
filter_dataset <-  eventReactive(input$filter_segments, {
  dat <- filter_loq()
  dat <- dat[, pData(dat)$GeneDetectionRate*100 >= input$filter_threshold]
  dat
})

output$filtered_stats <- DT::renderDataTable({
  dat <- filter_dataset()
  as.data.frame(dim(dat))
})

gene_detection_rates <- eventReactive(input$search_genes, {
  genes <- unlist(strsplit(input$goi_list, "\n"))
  dat <- filter_dataset()
  
  # Calculate detection rate:
  loq_mat <- t(esApply(dat, MARGIN = 1,
                       FUN = function(x) { x > pData(dat)$LOQ[,1]}))
  
  # ensure ordering since this is stored outside of the geomxSet
  loq_mat <- loq_mat[fData(dat)$TargetName, ]
  loq_mat <- loq_mat[, colnames(dat)]
  fData(dat)$DetectedSegments <- rowSums(loq_mat, na.rm = TRUE)
  fData(dat)$DetectionRate <-
    fData(dat)$DetectedSegments / nrow(pData(dat))
  
  goi_df <- data.frame(
    Gene = genes,
    Number = fData(dat)[genes, "DetectedSegments"],
    DetectionRate = scales::percent(fData(dat)[genes, "DetectionRate"]))
  return(list(goi_df, dat, loq_mat))
})

output$gois <- DT::renderDataTable({
  gene_detection_rates()[[1]]
  
  
})
    
build_detected_gois <- eventReactive(input$search_genes, {
  
  loq_mat <- gene_detection_rates()[[3]]
  loq_mat[isTRUE(loq_mat)] <- 1
  loq_mat[isFALSE(loq_mat)] <- 0
  loq_mat
})

output$gois_heatmap <- renderPlot({
  genes <- unlist(strsplit(input$goi_list, "\n"))
  to_heat <- build_detected_gois()
  genes <- intersect(genes, rownames(to_heat))
  to_heat <- to_heat[genes,]
  
  # segment annotations
  dat <- gene_detection_rates()[[2]]
  col_annotation <- HeatmapAnnotation(segment = pData(dat)[colnames(to_heat),]$segment)
  
  h <- Heatmap(to_heat, top_annotation = col_annotation)
  draw(h)
})

#####################################################
#####################################################
#####################################################

filter_genes <- eventReactive(input$filter_genes, {
  dat <- gene_detection_rates()
  dat <- dat[[2]]
  negativeProbefData <- subset(fData(dat), CodeClass == "Negative")
  neg_probes <- unique(negativeProbefData$TargetName)

  threshold <- input$gene_filter_threshold
  threshold <- threshold/100
  dat <- dat[fData(dat)$DetectionRate >= threshold | fData(dat)$TargetName %in% neg_probes, ]
  dat
})

output$filtered_segments_genes_summary <- DT::renderDataTable({
  filtered_all <- filter_genes()
  as.data.frame(dim(filtered_all))
})

##########################################################################
##########################################################################
# Normalization
##########################################################################
##########################################################################

expression_distribution_plot <- reactive({
  dat <- filter_genes()
  dat_r <- reshape2::melt(exprs(dat))
  ggplot(dat_r, aes(x=Var2, y=log2(value))) +
    geom_violin() +
    theme_classic() +
    xlab("Segment") +
    ylab("Expression") +
    theme(axis.text.x=element_blank())
})

output$prenorm_distribution <- renderPlot({
  expression_distribution_plot()
})

normalize_expression <- function(mat, method=c("quantile")){
  
  # normalize data based on user-defined method
  dat_norm <- log2(as.data.frame(normalize.quantiles(as.matrix(mat), keep.names=TRUE)))
  return(dat_norm)
}

normalized <- eventReactive(input$normalize, {
  dat <- filter_genes()
  dat_raw <- exprs(dat)
  dat_norm <- normalize_expression(dat_raw, method=input$norm_method)
  dat_norm
})

post_norm_expression_plot <- reactive({
  dat <- normalized()
  dat_r <- reshape2::melt(dat)
  ggplot(dat_r, aes(x=variable, y=value)) +
    geom_violin() +
    theme_classic() +
    xlab("Segment") +
    ylab("Expression") +
    theme(axis.text.x=element_blank())
})

output$postnorm_distribution <- renderPlot({
    post_norm_expression_plot()
})

# Negative probe vs. positive probe distributions
neg_vs_pos_plot <- eventReactive(input$normalize, {
  dat <- filter_genes()
  negativeProbefData <- subset(fData(dat), CodeClass == "Negative")
  neg_probes <- unique(negativeProbefData$TargetName)
  dat_norm <- normalized()
  dat_pos <- dat_norm[!(rownames(dat_norm)) %in% neg_probes,]
  dat_neg <- dat_norm[rownames(dat_norm) %in% neg_probes,]
  dat_pos <- reshape2::melt(dat_pos)
  dat_pos$type <- rep("Endogenous", nrow(dat_pos))
  dat_neg <- reshape2::melt(dat_neg)
  dat_neg$type <- rep("Negative", nrow(dat_neg))
  dat_pos_v_neg <- dplyr::bind_rows(dat_pos, dat_neg)
  print(head(dat_pos_v_neg))
  ggplot(dat_pos_v_neg, aes(x=type, y=value, fill=type, group=type, alpha=0.5)) +
    geom_violin() +
    theme_classic() +
    scale_fill_manual(values=c("lightblue", "grey"))
})

output$neg_vs_pos_plot <- renderPlot({
  neg_vs_pos_plot()
})

###########################################################
# PCA functions
###########################################################

runPCA <- function(df, scale=TRUE){
  
  pc <- prcomp(t(df), scale=scale)
  return (pc)
}

getPCA <- function(pc){
  
  pcs <- data.frame(pc$x)
  return(pcs)
}

getVE <- function(pc, component="PC1"){
  
  pve <- base::summary(pc)$importance[,component][[2]]
  return (pve)
}

plotPrincipalComponents <- function(pc, metadata, colourby="none", shapeby="none", group="none", continuous=FALSE,  pcs=c("PC1", "PC2")){
  
  # covariate must be in same order as pc rownames
  
  # get variance explained for each component
  ve1 <- getVE(pc, component=pcs[1])
  ve2 <- getVE(pc, component=pcs[2])
  
  ve1 <- round(ve1, 2)*100
  ve2 <- round(ve2, 2)*100
  
  # get data frame of components
  pca <- data.frame(pc$x)
  
  # add conditions
  if (colourby == "none"){
    pca$condition <- "none"}else{
      pca$condition <- metadata[,colourby]}
  
  # add shape
  if (shapeby == "none"){
    pca$shape <- "none"}else{
      pca$shape <- metadata[,shapeby]}#
  
  #if (group == "none"){
  #  pca$group <- "none"}else{
  #    pca$group <- metadata[,group]}
  
  if (continuous==FALSE){
    pca$condition <- factor(pca$condition, levels=unique(pca$condition))
  }
  
  # plot
  pc1 <- pcs[1]
  pc2 <- pcs[2]
  
  # labels
  xlabel <- paste(pc1, ve1, sep=" (")
  xlabel <- paste(xlabel, "%", sep="")
  xlabel <- paste(xlabel, ")", sep="")
  ylabel <- paste(pc2, ve2, sep=" (")
  ylabel <- paste(ylabel, "%", sep="")	
  ylabel <- paste(ylabel, ")", sep="")
  
  n <- length(unique(pca$condition))
  colours <- rainbow(n, s=0.7, v=0.6)
  
  plot1 <- ggplot(pca, aes_string(x=pc1, y=pc2, colour="condition", shape="shape", label="group"))
  plot2 <- plot1 + geom_point(size=3)
  plot3 <- plot2 + theme_bw() 
  plot4 <- plot3 + xlab(xlabel) + ylab(ylabel)
  if (continuous==TRUE){
    plot4 <- plot4 + scale_colour_gradient()}
  else{
    plot4 <- plot4 + scale_colour_manual(values=colours)
  }
  return(plot4 + theme(legend.position="none")) 
}

pca_plot <- eventReactive(input$normalize, {
  dat <- normalized()
  all_dat <- filter_genes()
  metadata <- pData(all_dat)
  metadata$label <- paste0(metadata$file, ":", metadata$slide.name, ":", metadata$roi)
  pc <- runPCA(dat, scale=FALSE)
  plotPrincipalComponents(pc, metadata, shapeby="label", colourby = "segment", group="roi")
    
})

output$postnorm_pca <- renderPlotly({
  pca_plot()
})

########################################################################
# heatmap of raw counts and normalized values
########################################################################

plot_normalized_heatmap <- eventReactive(input$search_genes_to_view, {
  
  dat <- filter_genes()
  
  # GOIs
  genes <- unlist(strsplit(input$genes_to_view, "\n"))

  # raw data
  raw <- exprs(dat)
  genes <- intersect(genes, rownames(raw))
  raw <- raw[genes,]

  # normalised matrix
  norm <- normalized()
  norm <- norm[genes,]
  
  # scaled data
  norm_scaled <- as.data.frame(t(apply(norm, 1, scale)))
  colnames(norm_scaled) <- colnames(norm)
  
  # make colnames slide.name and roi
  pd <- pData(dat)
  roi <- gsub("=", "", gsub('"', '', pd$roi))
  
  # column labels informative - if long :(
  col_labels <- roi
  col_labels <- paste0(colnames(norm_scaled), "-", roi)
    
  # segment annotations
  col_annotation <- HeatmapAnnotation(segment = pData(dat)[colnames(raw),]$segment)
    
  h1 <- Heatmap(raw, top_annotation = col_annotation, name="Raw data", column_labels=col_labels, column_names_max_height = unit(20, "cm"))
  h2 <- Heatmap(norm, top_annotation = col_annotation, name = "Normalized data", column_labels=col_labels)
  h3 <- Heatmap(norm_scaled, top_annotation = col_annotation, name = "Normalized + scaled data", column_labels=col_labels)
  draw(h1 + h2 + h3)
})

output$genes_to_view_heatmap <- renderPlot({
  plot_normalized_heatmap()
}, width=1000, height=1000)

########################################################################
# Download object or .tsv files
########################################################################

output$download_rds <- downloadHandler(
   filename = function() {
     paste('geomx_data-', Sys.Date(), '.Rds', sep='')
   },
   content = function(file) {
     dat <- filter_genes()
     
     # need to add normalized slot
     attr(dat, "normalized") <- normalized()
     saveRDS(dat, file)
   }
)

output$download_normalized <- downloadHandler(

  filename = function() {
    paste('geomx_normalized_counts-', Sys.Date(), '.tsv', sep='')
  },
  content = function(file) {
    dat <- normalized()
    dat <- data.frame(Gene = rownames(dat), dat)
    write.table(dat, file, sep="\t", row.names=FALSE, quote=FALSE)
  }
)

output$download_metadata <- downloadHandler(

    filename = function() {
    paste('geomx_metadata-', Sys.Date(), '.tsv', sep='')
  },
  content = function(file) {
    dat <- filter_genes()
    metadata <- pData(dat)
    metadata <- data.frame(segment_id = rownames(metadata), metadata)
    write.table(metadata, file, sep="\t", row.names=FALSE, quote=FALSE)
  }
)

}

# Run the application 
shinyApp(ui = ui, server = server)
