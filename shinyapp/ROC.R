ROC_table <- DT::renderDT({
  analyzed() %>%
    DT::datatable(rownames=F, extensions = 'FixedColumns', width = 1250, filter = 'top', selection="single",
                  options = list(scrollX = TRUE, fixedColumns = list(leftColumns = 6), pageLength = 5)) %>%
    DT::formatStyle(c(1,5), "word-break" = "break-all") %>%
    DT::formatRound(grep("(pAUC|EF)", names(analyzed())))
}) %>%
  bindEvent(input$dataset)






ROC_plot <- renderPlot({
  sel <- slice(analyzed(), input$ROC_table_rows_selected)
  
  data <- to_analyze() %>%
    filter(Package == sel$Package & Filter == sel$Filter & Analysis == sel$Analysis) %>%
    arrange(rank)
  
  pROC::roc(response = data$hit, predictor = data$rank,
            quiet=F, direction =">", plot=T,
            print.auc=T, auc.polygon=T,
            partial.auc=c(95, 100), partial.auc.focus="se", partial.auc.correct=T,
            legacy.axes=T, percent=T,
            xlab="False Positives %", ylab="True Positives %",
            main = glue::glue_collapse(sep = " / ",
                                       x = c(sel$Package,
                                             paste("AUC", sel$pAUC_0),
                                             paste("Filter", sel$Filter),
                                             paste("Analysis", sel$Analysis))))
  pROC::roc(response = data$hit, predictor = data$rank,
            quiet=F, direction =">",
            plot=T, add=T, print.auc=T)
  
}) %>%
  bindEvent(input$ROC_table_rows_selected)






density_plot <- renderPlot({
  sel <- slice(analyzed(), input$ROC_table_rows_selected)
  
  plot(main="", density(to_analyze() %>%
                          filter(Package == sel$Package & Filter == sel$Filter & Analysis == sel$Analysis) %>%
                          pull(value)
                        )
       )
  
  
}) %>%
  bindEvent(input$ROC_table_rows_selected)