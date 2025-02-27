options(shiny.trace = TRUE)

library(tidyverse)

function(input, output, session) {
  
  # Load and keep in cache the raw and analyzed datasets when they are selected in the UI
  to_analyze <- reactive({readRDS(paste0(input$dataset, "_toana.rds"))}) # %>% bindCache(input$dataset, cache=session)
  analyzed <- reactive({readRDS(paste0(input$dataset, "_ana.rds"))}) # %>% bindCache(input$dataset, cache=session)
  
  
  
  ## Table tab
  
  source("ROC.R", local=T)
  output$ROC_table <- ROC_table
  output$ROC_plot <- ROC_plot
  output$density_plot <- density_plot
  
  
  
  ## Ranked-residues tab
  
  observeEvent(
  # eventReactive(
    input$show_H1,
  # observe(
    {
      source("Heatmap_rank_server.R", local=T)
      to_ana <- to_analyze()
      ana <- analyzed()
      
      ht <- get_heatmap(to_ana, ana,
                  input$splits, dplyr::setdiff(unique(c(input$splits, input$annos)), c("None")), input$sort,
                  as.numeric(input$numcols))#, as.numeric(input$topn))
      
      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
        input = input, 
        output = output, 
        session = session,
        ht_list = ht,
        heatmap_id = "H1")
    },
  suspended=TRUE,
  # ignoreNULL=FALSE,
  ignoreInit=TRUE
    
  ) # %>% bindEvent(input$show_H1)
  
  
  ## Residues tab
  
  # observeEvent(
  eventReactive(
    # input$show_H2,
  # observe(
    {
      source("Heatmap_aas_server.R", local=T)
      to_ana <- to_analyze()
      ana <- analyzed()
      
      analysis <- get_analysis(ana, input$splits_aa, input$sort_aa)
      mat <- order_matrix(get_matrix(to_ana, input$dataset), analysis)
      
      ht <- get_heatmap(mat, analysis,
                        input$splits_aa, input$sort_aa,
                        dplyr::setdiff(unique(c(input$splits_aa, input$annos_aa)), c("None")))
      
      show_scatterplot <- get_show_scatterplot(mat)
      
      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
        input = input, 
        output = output, 
        session = session,
        ht_list = ht,
        click_action = show_scatterplot,
        heatmap_id = "H2")
    },
    suspended=TRUE,
    # ignoreNULL=FALSE,
    ignoreInit=TRUE
  ) # %>% bindEvent(input$show_H2)
}