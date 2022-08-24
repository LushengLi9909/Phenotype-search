library(shiny)
library(shinythemes)
library(shinyWidgets)
library(ontologyIndex)
library(DT)
library(plotly)
require(visNetwork)
source("phenotype_search.R")
source("cell-phenos.R")


ui <- fluidPage(
navbarPage(theme = shinytheme("flatly"),title = h2("Cell-Phenotype"),
           tabPanel(h3("Phenotype"),fluidRow(
             column(3, wellPanel(
               textInput("keywords", "Enter search keywords or phenotype term",
                         value = "anxiety, gut"),
               numericInput("qvalue", "q-value threshold", value = 0.005,min = 0,max =1,step=0.0005),
               numericInput("fold", "Minimum fold change", value = 1),
               numericInput("sd", "Minimum standard deviations from mean", value = 0),
               actionBttn(
                 inputId = "do",label = "Plot",style = "jelly", color = "success"),
               br(),
               br(),
               sliderInput("search_width", "Download image width (px)",
                           min = 400, max = 2000, value = 600, step = 100),
               sliderInput("search_height", "Download image height (px)",
                           min = 400, max = 2000, value = 600, step = 100),
               downloadLink("download", "Download figure")
               
           )),
           column(9, plotlyOutput("keyword_plot",height='100%'),
                  DTOutput('keyword_df')),
           column(2),
           column(10,h3("Cell-Phenotype Network:"),
                  h5("(if you enter one phenotype term, the network will show the descendant phenotypes)")),
           visNetworkOutput("Visnetwork", width = "100%", height = "660px"),
           DTOutput('network_df')
           ))
           )
)
           

server <- function(input, output) {
  plot_pheno <- eventReactive(input$do,plot_phenotype_counts(keyword_search_df(process_search_terms(input$keywords),q_threshold = input$qvalue,
                                                                    fold_threshold = input$fold,
                                                                    min_sd_from_mean = input$sd), input$keywords))
  output$keyword_plot <- plotly::renderPlotly(plot_pheno())
  table_pheno <- eventReactive(input$do,keyword_search_df(process_search_terms(input$keywords),q_threshold = input$qvalue,
                                              fold_threshold = input$fold,
                                              min_sd_from_mean = input$sd))
  output$keyword_df <- DT::renderDT(table_pheno(),
                                    rownames = FALSE,extensions = 'Buttons', 
                                    options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))

  output$download <- downloadHandler(
    filename =paste0("Phenotype_search_",input$keywords,Sys.Date(),".png"),
    content = function(filename) {
      png(filename, width = input$search_width, height = input$search_height)
      print(plot_pheno())
      dev.off()
    },
    contentType = "image/png")
  Visnet <- eventReactive(input$do,visnetwork_plot_full(phenotype_to_genes, all_results_merged,mpo,branch = input$keywords,q_threshold = input$qvalue,
                                                 fold_threshold = input$fold))
  output$Visnetwork <- renderVisNetwork(Visnet())
  sub_pheno <- eventReactive(input$do,subset_phenos(phenotype_to_genes, all_results_merged,mpo,branch = input$keywords,q_threshold = input$qvalue,
                                          fold_threshold = input$fold))
  output$network_df <- DT::renderDT(sub_pheno(),
                                    rownames = FALSE,extensions = 'Buttons', 
                                    options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
}

# Run the application 
shinyApp(ui = ui, server = server)
