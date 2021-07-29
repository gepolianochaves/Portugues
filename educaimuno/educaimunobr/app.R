#execute o comando abaixo para encontrar as bibliotecas
#options(repos = BiocManager::repositories())


library(shiny)
library(shinythemes)
library(shiny.i18n)
i18n <- Translator$new(automatic = TRUE)
i18n$set_translation_language("pt")

library(msa)
library(ape)
library(seqinr)


#library(adegenet)
#library(ggtree)
#library(ggplot2)
#library(stats)
#library(ips)



 getdfmutations<- function(fastafile) {

   sequences <- readAAStringSet(fastafile)
   sarsAln <- msa(sequences)
   conMat <- consensusMatrix(sarsAln)
   subset <- conMat[, 1:240]
   
   return(subset)
   
 }


 plottree <- function(fastafile) {
   
   
 }



# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                navbarPage(
                  theme = "cerulean",  # <--- To use a theme, uncomment this
                  "Visualizações de análises filogenéticas a partir dos genomas de Sars-Cov-2",
                  tabPanel("Frequência das Mutações",
                           sidebarPanel(
                             fileInput("fastafile", "Escolha o arquivo fasta",
                                       multiple = TRUE,
                                       accept = c(".fasta", ".fa")),
                             # Horizontal line ----
                             tags$hr(),
                             #actionButton("button_covid", "Dados covid")
                           ), # sidebarPanel
                           mainPanel(
                             h4("Tabela de frequências"),
                             tableOutput("contents")
                           ) # mainPanel
                           
                  ), # Navbar 1, tabPanel
                  tabPanel("Árvores filogenéticas", 
                           mainPanel(
                             imageOutput("contents2")
                           ) #
                           )
                ) # navbarPage
) # fluidPage


# Define server function  
server <- function(input, output) {
  
    output$contents <- renderTable({
     req(input$fastafile)
     sequences <- readAAStringSet(input$fastafile$datapath)
     sarsAln <- msa(sequences)
     conMat <- consensusMatrix(sarsAln)
     df <- conMat[, 1:10]
    return(df)})
    
    output$contents2 <- renderPlot(height=2000, {
      req(input$fastafile)
      sequences <- readAAStringSet(input$fastafile$datapath)
      sarsAln <- msa(sequences)
      sarsAln2 <- msaConvert(sarsAln, type="seqinr::alignment")
      d <- dist.alignment(sarsAln2, "identity")
      sarsTree <- nj(d)
      g <- plot(sarsTree, main="Phylogenetic Tree of Sars-Cov-2")
      return(g)
    })
    
    
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)


#rsconnect::deployApp('/mnt/D/PycharmProjects/local_educaimuno/educaimuno/educaimunobr')
#https://github.com/MSQ-123/CovidMutations

