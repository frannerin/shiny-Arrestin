fluidPage(
  fluidRow(
    br(),
    p("Each row of the heatmap is a combination of network, filter and analysis, and the rows can be grouped and/or ordered by the selected variables. For each combination/row, the centrality values are normalized in the range 0-1 and residues are shown in sequence order, with sequence annotations of conservation and binding regions shown at the top of the heatmap. Clicking on a row will show the correlation plot between the centrality values and the ConSurf scores on the right. Selecting a specific area of the heatmap will show a zoomed subheatmap below."),
    br()
  ),
  
   fluidRow(
    column(6, 
           selectInput(inputId = "splits_aa", label = "Network info to use for grouping and sorting:",
                       choices = c("None", "Type of network", "Atom or angle", "Correlation metric", "Filter", "Analysis"),
                       selected = "None", multiple = T, selectize = T, width = "100%"),
           p("Variables of information of the networks to use to group them together for sorting and organized visualization.")),
  
    column(6,
           selectInput(inputId = "sort_aa", label = "Sort network combinations by this EF or pAUC:",
                       choices = c("norm_EF_1",   "norm_EF_2.5", "norm_EF_5",   "norm_EF_10",  "norm_EF_15",
                                   "NEF_3",       "NEF_5",       "NEF_10",      "NEF_20",  "NEF_50",
                                   "AUC", "Spearman corr"),
                       selected = "AUC", multiple = T, selectize = T, width = "100%"),
           p("Enrichment metrics, AUC and Spearman correlation to use for sorting the networks and/or network groups; the second and following variables will be used to break ties."))
    ),
   fluidRow(
     column(6,
            selectInput(inputId = "annos_aa", label = "Network info to show next to the heatmap:",
                        choices = c("Type of network", "Atom or angle", "Correlation metric", "Filter", "Analysis"),
                        selected = "Type of network", multiple = T, selectize = T, width = "100%"),
            p("Additional variables of information of the networks to show next to the Heatmap, without using them for sorting and thus without affecting the order of the networks in the Heatmap.")),
  
     column(6, wellPanel(actionButton("show_H2", "Update heatmap")))
   ),
   splitLayout(
     InteractiveComplexHeatmap::originalHeatmapOutput(heatmap_id="H2", title = "Original heatmap",
                                                      width=1000, height=900),
     InteractiveComplexHeatmap::HeatmapInfoOutput(heatmap_id="H2", #title = "Output",
                                                  output_ui = plotly::plotlyOutput("scatterplot")), 
     cellWidths = c(validateCssUnit('auto'), validateCssUnit(400))
   ),
   InteractiveComplexHeatmap::subHeatmapOutput(heatmap_id="H2", title = "Sub-heatmap", width=1000)
)