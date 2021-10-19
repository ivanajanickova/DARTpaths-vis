#v37

#########################
#load necessary packages#
#########################

library(shiny)
library(shinydashboard)
library(reshape2)
library(ggplot2)
library(shinyWidgets)
library(wesanderson)
library(splitstackshape)
library(DT)

######################
#the global functions#
######################

#Functions related to list of pathways
pathname <- c("R-HSA-69239", "R-HSA-8937144", "R-HSA-5358351","R-HSA-381119", "R-HSA-195721","R-HSA-194138", "R-HSA-157858", "R-HSA-8869496"
              , "R-HSA-8864260","R-HSA-71262", "R-HSA-5674499", "R-HSA-888568","R-HSA-380612", "R-HSA-166665", "R-HSA-191273", "R-HSA-399710", "R-HSA-205017", "R-HSA-187042", "R-HSA-977443")
names(pathname) <- c("Synthesis of DNA (R-HSA-69239)","Aryl hydrocarbon receptor signaling (R-HSA-8937144)",
                     "Signaling by Hedgehog (R-HSA-5358351)", "Unfolded Protein Response (UPR) (R-HSA-381119)", "Signaling by WNT (R-HSA-195721)",
                     "Signaling by VEGF (R-HSA-194138)", "Gap junction trafficking and regulation (R-HSA-157858)",
                     "TFAP2A acts as a transcriptional repressor during retinoic acid induced cell differentiation (R-HSA-8869496)",
                     "Transcriptional regulation by the AP-2 (TFAP2) family of transcription factors (R-HSA-8864260)",
                     "Carnitine synthesis (R-HSA-71262)", "Negative feedback regulation of MAPK pathway (R-HSA-5674499)","GABA synthesis (R-HSA-888568)",
                     "Metabolism of serotonin (R-HSA-380612)", "Terminal pathway of complement (R-HSA-166665)", "Cholesterol biosynthesis (R-HSA-191273)",
                     "Activation of AMPA receptors (R-HSA-399710)", "NFG and proNGF binds to p75NTR (R-HSA-205017)", "NTRKA activation by NGF (R-HSA-187042)", "GABA receptor activation (R-HSA-977443)"
) 

choiceVec <- c(pathname)

#Functions related to loading results
loadfiles <- function(xvar){
  dir <- paste0('./',xvar,'_Enrichment_Results/')
  temp = list.files(path = dir, pattern="*_Result.txt", full.names = TRUE)
  myfiles = lapply(temp, function(x) { read.csv(x, skip = 3)})
  names(myfiles) <- c("celegans","mouse", "slimemould", "zebrafish")
  return(myfiles)
}

#Functions related to linking to external sources for each organism
ToLink <- function(txt) {
  paste0('<a href=https://wormbase.org/species/all/phenotype/',txt,">",txt,'</a>')
}

ToLinkZFIN <- function(txt) {
  paste0('<a href=https://zfin.org/',txt,">",txt,'</a>')
}

ToLinkDict  <- function(val) {
  sprintf('<a href="https://www.ebi.ac.uk/ols/search?q=%s" target="_blank" class="btn btn-primary">more info</a>',val)
}

ToLinkMus <- function(txt) {
  paste0('<a href=https://www.mousephenotype.org/data/phenotypes/',txt,">",txt,'</a>')
}

##################
#Start Shiny here#
##################

ui <- dashboardPage(
  dashboardHeader(title = "Phenotype Enrichment v1.0", titleWidth = 300
                  
  ),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    tags$style(HTML("

                    
                    .box.box-solid.box-primary>.box-header {
                    color:#fff;
                    background:#f44378
                    }
                    
                    .box.box-solid.box-primary{
                    border-bottom-color:#f44378;
                    border-left-color:#f44378;
                    border-right-color:#f44378;
                    border-top-color:#f44378;
                    }

                    .box.box-solid.box-warning>.box-header {
                    color:#fff;
                    background:#ecad68
                    }
                    
                    .box.box-solid.box-warning{
                    border-bottom-color:#ecad68;
                    border-left-color:#ecad68;
                    border-right-color:#ecad68;
                    border-top-color:#ecad68;
                    }
                    
                    ")),
    fluidRow(
      box(width = 2, title="Select pathway", solidHeader = TRUE, status = "primary", selectInput("selectpath", label = "from list of human pathway below:", 
                      c(pathname), size = 5,selectize = FALSE),
          
          actionBttn(
            inputId = "go",
            label = "Go",
            color = "danger",
            block = TRUE,
            size = "sm",
            style = "fill"
          )),
      box(width = 2, title="Select species", solidHeader = TRUE, status = "primary",
                 actionBttn(
                   inputId = "zebrafish",
                   label = "D. rerio",
                   color = "danger",
                   block = TRUE,
                   size = "sm",
                   style = "fill"
                 ),br(),
                 actionBttn(
                   inputId = "slimemould",
                   label = "D. discoideum",
                   color = "danger",
                   block = TRUE,
                   size = "sm",
                   style = "fill"
                 ),br(),
                 actionBttn(
                   inputId = "mouse",
                   label = "M. musculus",
                   color = "danger",
                   block = TRUE,
                   size = "sm",
                   style = "fill"
                 ),br(),
                 actionBttn(
                   inputId = "celegans",
                   label = "C. elegans",
                   color = "danger",
                   block = TRUE,
                   size = "sm",
                   style = "fill"
                 )),
      
      box(width = 4,  title = "Welcome!", status = "warning", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, h4("to Phenotype Enrichment tool!"),
          p(
            p(h6(strong("follow this step for analysis:")),
            p(h6(strong("1.")," select the human pathway you want to analyze, then click", strong("Go") )),
            p(h6(strong("Summary conservation")," box will appear barplot showing the number of orthologs of this pathway in the each species " )),
            p(h6(strong("2.")," select the specific species to see the enriched phenotypes")),
            p(h6(strong("3.")," once the table appears below, you can", strong("Download"), "it as well"))
            
          )))
      
         )
             ,

    fluidRow(
             box(width = 4, p(h4("Summary Conservation", textOutput("pathtext"), style="display:inline")), plotOutput("plot",height = 150)),
             box(width = 4,  title = "More on enrichment analysis", status = "warning", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                 p(strong("Pvalue"),"is calculated with hypergeometric test, only the ones that passed the", strong("qvalue"), " (with Benjamini-Hochberg FDR correction) is presented here as", strong("Enriched Phenotype"),
                   p(h6(strong("n"),": the total of annotated genes in the geneset (genes in pathway) that appears in the phenotype library")),
                   p(h6(strong("m"),": all the genes in phenotype library minus n")),
                   p(h6(strong("i"),": the total of given phenotype in the geneset in the geneset"))
                   
                 ))
             ),
    fluidRow(box(width = 4, p(h3("Enriched Phenotypes"),
                              p(em("(Select the species above to enrich)", style = "font-size:11pt"))))),
    fluidRow(column(width = 12, textOutput("noenrich"))),
    tags$style(HTML('table.dataTable th {background-color: pink !important;}')),
    fluidRow(
             column(width = 12, dataTableOutput('table')))
    ,

    
    # Also add some custom CSS to make the title background area the same
    # color as the rest of the header.
    tags$head(tags$style(HTML('
                              
                              
                              /* logo */
                              .skin-blue .main-header .logo {
                              background-color: #b51645;
                              }
                              
                              /* logo when hovered */
                              .skin-blue .main-header .logo:hover {
                              background-color: #7a0729;
                              }
                              
                              /* navbar (rest of the header) */
                              .skin-blue .main-header .navbar {
                              background-color: #f44378;
                              }
                              
                              /* toggle button when hovered  */                    
                              .skin-blue .main-header .navbar .sidebar-toggle:hover{
                              background-color: #ff69b4;
                              }
                              
                              ')))
    ))


server <- function(input, output) { 
  outVar <- reactive({
    loadfiles(input$selectpath)
  })
  
  outSum <- eventReactive(input$go,{
    df <- read.csv(paste0('./',input$selectpath,'_Enrichment_Results/',input$selectpath, '_summary_conservation', '.txt'),header=TRUE,sep=',')
    dat.m2 <- melt(df, id.var = "Species")
    dat.m2$variable <- reorder(dat.m2$variable, dat.m2$value)
    dat.m2$variable <- factor(dat.m2$variable, levels=rev(levels(dat.m2$variable)))
    return(dat.m2)
  })
  pathselect <- eventReactive(input$go,
                              {paste(names(choiceVec)[choiceVec == input$selectpath], "pathway")
                                
                              })
  output$pathtext <- renderText({pathselect()})
  output$table <- renderDataTable({
    #reactable(reactpath, filterable = TRUE, minRows = 10)
    if (is.null(v$data)) return()
    if (dim(v$data)[1] == 0) {
      return()
    }
    names(v$data) <- gsub("\\.", " ", names(v$data))
    datatable(v$data[-1], escape = FALSE, extensions = 'Buttons',
              options = list(scrollX = TRUE,autoWidth = FALSE, 
                             columnDefs = list(list(width = '350px', targets = c(2)),
                                               list(width = '150px', targets = c(1)),
                                               list(width = '400px', targets = c(9))),
                             dom = 'Bfrtip',
                             buttons = c('copy', 'csv', 'excel','print')
                             ))})
    
  output$nothing <- renderText({
    if(input$selectpath == "pleb"){
      out <- paste("No result to show")
    }
  })
  output$noenrich <- renderText({
    if (is.null(v$data)) return()
    if(dim(v$data)[1] == 0){
      out <- paste("No result to show")
    }
  })
  v <- reactiveValues(data = NULL)

  
  observeEvent(input$celegans, {
    v$data <- data.frame(outVar()$celegans) 
    if (dim(v$data)[1] == 0) {
      return()
    } else{
      v$data$Enriched.Phenotype <- ToLink(v$data$Enriched.Phenotype)
      v$data <- concat.split.list(v$data, split.col="Overlap.Genes", sep=",", drop = TRUE)
      v$data$Overlap.Genes_list <-lapply(v$data$Overlap.Genes, function(x)  
        paste0('<a href=https://wormbase.org/species/c_elegans/gene/',x,">",x,'</a>'))
    }
  })
  observeEvent(input$mouse, {
    v$data <- data.frame(outVar()$mouse) 
    if (dim(v$data)[1] == 0) {
      return()
    } else{
      v$data$Enriched.Phenotype <- ToLinkMus(v$data$Enriched.Phenotype)
      v$data <- concat.split.list(v$data, split.col="Overlap.Genes", sep=",", drop = TRUE)
      v$data$Overlap.Genes_list <-lapply(v$data$Overlap.Genes, function(x)  
        paste0('<a href=http://www.informatics.jax.org/marker/',x,">",x,'</a>'))
    }
  })
  observeEvent(input$zebrafish, {
    v$data <- data.frame(outVar()$zebrafish) 
    if (dim(v$data)[1] == 0) {
      return()
    } else{
      v$data$Enriched.Phenotype <- ToLinkZFIN(v$data$Enriched.Phenotype)
      v$data <- concat.split.list(v$data, split.col="Overlap.Genes", sep=",", drop = TRUE)
      v$data$Overlap.Genes_list <-lapply(v$data$Overlap.Genes, function(x)  
        paste0('<a href=https://zfin.org/',x,">",x,'</a>'))
    }
  })
  observeEvent(input$slimemould, {
    v$data <- data.frame(outVar()$slimemould) 
    if (dim(v$data)[1] == 0) {
      return()
    } else{
      v$data$Link <- ToLinkDict(v$data$Enriched.Phenotype)
      v$data <- concat.split.list(v$data, split.col="Overlap.Genes", sep=",", drop = TRUE)
      v$data$Overlap.Genes_list <-lapply(v$data$Overlap.Genes, function(x)  
      paste0('<a href=http://dictybase.org/gene/',x,">",x,'</a>'))
    }
  })
  
 
  output$plot <- renderPlot(
    ggplot(outSum(), aes(x = Species, y=value, fill=variable, label=value)) + 
                              geom_bar(position="stack", stat="identity")+ geom_text(size = 4, position = position_stack(vjust = 0.5), color = "white") +
                            #wont do position fill  
                            theme(
                                panel.grid.major = element_blank(), 
                                panel.grid.minor = element_blank(),
                                axis.ticks = element_blank(),
                                panel.background = element_rect(fill = "transparent",colour = NA),
                                plot.background = element_rect(fill = "transparent",colour = NA),
                                legend.text=element_text(size=12)
                              ) + coord_flip() + labs(x="",y="") + scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1"))
  )
  
}

shinyApp(ui, server)