#
# Geneva HRD test -- graphical interface to oncoscanR package via R shiny
# author: Yann Christinat
# company: HUG
# date: june 3, 2025
#


library(shiny)
library(oncoscanR)
library(magrittr)
library(shinycssloaders)


# Example analysis function – replace with your real logic
analyze_data <- function(file) {
  segments <- load_chas(file, oncoscanR::oncoscan_na33.cov)
  
  # Clean the segments: resctricted to Oncoscan coverage, LOH not overlapping
  # with copy loss segments, smooth&merge segments within 300kb and prune
  # segments smaller than 300kb.
  segs.clean <- trim_to_coverage(segments, oncoscanR::oncoscan_na33.cov) %>%
    adjust_loh() %>%
    merge_segments() %>%
    prune_by_size()
  
  # Get the number of nLST and TDplus
  wgd <- score_estwgd(segs.clean, oncoscanR::oncoscan_na33.cov)
  hrd <- score_nlst(segs.clean, wgd["WGD"], oncoscanR::oncoscan_na33.cov)
  
  n.td <- score_td(segs.clean)
  
  mbalt <- score_mbalt(segs.clean, oncoscanR::oncoscan_na33.cov, loh.rm=TRUE)
  
  hrd.label <- hrd["HRD"]
  if (mbalt['sample'] / mbalt['kit'] < 0.01)
    hrd.label <- paste(hrd["HRD"], "(no tumor?)")
  
  # Get the alterations into a single list.
  return(list(
    HRD = paste0(hrd.label, ", nLST=", hrd["nLST"]),
    TDplus = n.td$TDplus,
    avgCN = paste0(substr(as.character(wgd["avgCN"]), 1, 4), " (", wgd["WGD"], 
                   " genome-doubling event", ifelse(wgd["WGD"]>1, "s)", ")"))
  ))
}

ui <- fluidPage(
  titlePanel("Geneva HRD test App"),
  
  # Informational text
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose a ChAS segment file",
                accept = c(".csv", ".tsv", ".txt")),
      actionButton("analyze", "Run Analysis", class = "btn-primary"),
      br(), br(),
      downloadButton("download", "Download Results")
    ),
    
    mainPanel(
      h4("Information to user"),
      p(tagList("The app relies on the OncoscanR R package that handles Copy Number Variation analyses originating from the Applied Biosystems(TM) OncoScan(TM) CNV Assay (",
              a("oncoscanR package link", href="https://www.bioconductor.org/packages/release/bioc/html/oncoscanR.html"), "). ",
                "It allows computation of the HRD score used by the Geneva HRD test[Christinat 2023] and the tandem duplication plus score (TDplus) to identify CDK12-mutated tumors [Popova 2016].")),
      p(paste0("The HRD test is positive if the nLST score is greater or equal to 15. Of note, its predictive power for PARP inhibitors response has only been validated in ",
               "high-grade ovarian cancers as part of a research project[Christinat 2025] and it is not a companion diagnostic test.")),
      p("There is no recommandation on how to interpret the TDplus score but a number above 100 is generally indicative of a CDK12 mutation."),
      
      
      
      br(), br(),
      h4("Analysis Results"),
      fluidRow(
        column(4, strong("HRD status: "), withSpinner(textOutput("result1"))),
        column(4, strong("TDplus (CDK12-mut): "), withSpinner(textOutput("result2"))),
        column(4, strong("Average copy number: "), withSpinner(textOutput("result3")))
      ),
      
      br(), hr(), br(),
      
      h4("References"),
      tags$ul(
        tags$li(
          "Christinat et al., ",
          em("Geneva HRD test is predictive of survival benefit from olaparib and bevacizumab maintenance in ovarian cancer"),
          ", JCO Precision Oncology, in press, 2025."
        ),
        tags$li(
          "Christinat et al., ",
          em("Normalized LST Is an Efficient Biomarker for Homologous Recombination Deficiency and Olaparib Response in Ovarian Carcinoma"),
          ", JCO Precision Oncology, volume 7, 2023."
        ),
        tags$li(
          "Popova et al., ",
          em("Ovarian Cancers Harboring Inactivating Mutations in CDK12 Display a Distinct Genomic Instability Pattern Characterized by Large Tandem Duplications"),
          ", Cancer Res volume 76(7), 2016."
        )
      ),
        
        
      br(), hr(), br(),
      
      # Disclaimer in small font
      tags$div(
        style = "font-size: 9px; color: #666;",
        HTML(paste0("Disclaimer: The use of Oncoscan(TM) platform in this research does not imply an endorsement or recommendation of Thermo Fisher Scientific Inc. ",
                    "and its subsidiaries for the utilization of any specific algorithm or methodology with the Oncoscan(TM) platformm for HRD analysis. Thermo Fisher ",
                    "Scientific Inc. and its subsidiaries make no claims regarding the suitability, performance or efficacy of any algorithms or methodologies used in ",
                    "conjunction with the Oncoscan(TM) platformm for HRD analysis. Furthermore, Thermo Fisher Scientific Inc. and its subsidiaries take no ",
                    "responsibility, and anyone using any specific algorithms or methodologies in conjunction with the Oncoscan(TM) platformm for HRD analysis is solely ",
                    "responsible for researching, identifying, and obtaining any necessary third-party rights and ensuring that such use is in compliance with applicable laws and regulations."))
      )
    )
  )
)

server <- function(input, output, session) {
  results <- eventReactive(input$analyze, {
    req(input$file)
    analyze_data(input$file$datapath)
  })
  
  output$result1 <- renderText({ req(results()); results()$HRD })
  output$result2 <- renderText({ req(results()); results()$TDplus })
  output$result3 <- renderText({ req(results()); results()$avgCN })
  
  output$download <- downloadHandler(
    filename = function() {
      paste0("results_HRD_", tools::file_path_sans_ext(input$file$name), ".txt")
    },
    content = function(file) {
      req(results())
      fileConn <- file(file)
      writeLines(c(
        paste("Uploaded file:", input$file$name),
        paste("oncoscanR version:", as.character(packageVersion("oncoscanR"))),
        "",
        "Analysis Results",
        paste("  HRD status: ", results()$HRD),
        paste("  TDplus (CDK12-mut): ", results()$TDplus),
        paste("  Average copy number: ", results()$avgCN)
      ), fileConn)
      close(fileConn)
    }
  )
}

shinyApp(ui, server)
