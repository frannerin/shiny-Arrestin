fluidPage(
  fluidRow(
    br(),
    p("Each row of the heatmap is a combination of network, filter and analysis, and the rows can be grouped and/or ordered by the selected variables. For each combination/row, residues have been ranked by the centrality value and are depicted in order, colored according to whether they are conserved (ConSurf score of 8 or 9) or not, with the aim to detect the rows where top results are enriched in conserved residues. All ordered/ranked residues can be visualized, or a number among the top can be selected. Selecting a specific area of the heatmap will show a zoomed subheatmap below."),
    br()
  ),
  
  fluidRow(
    column(4, 
      selectInput(inputId = "splits", label = "Network info to use for grouping and sorting:",
                  choices = c("None", "Type of network", "Atom or angle", "Correlation metric", "Filter", "Analysis"),
                  selected = "None", multiple = T, selectize = T, width = "100%"),
      p("Variables of information of the networks to use to group them together for sorting and organized visualization.")),
  
    column(4, 
      selectInput(inputId = "sort", label = "Sort network combinations by this EF or pAUC:",
                  choices = c("norm_EF_1",   "norm_EF_2.5", "norm_EF_5",   "norm_EF_10",  "norm_EF_15",  
                              "NEF_3",       "NEF_5",       "NEF_10",      "NEF_20",  "NEF_50",
                              "AUC", "Spearman corr"),
                  selected = "AUC", multiple = T, selectize = T, width = "100%"),
      p("Enrichment metrics, AUC and Spearman correlation to use for sorting the networks and/or network groups; the second and following variables will be used to break ties.")),
  
    column(4, 
      selectInput(inputId = "annos", label = "Network info to show next to the heatmap:",
                  choices = c("Type of network", "Atom or angle", "Correlation metric", "Filter", "Analysis"),
                  selected = "Type of network", multiple = T, selectize = T, width = "100%"),
      p("Additional variables of information of the networks to show next to the Heatmap, without using them for sorting and thus without affecting the order of the networks in the Heatmap."))
  ),
  
  fluidRow(
    column(4, 
      numericInput(inputId = "numcols", label = "Number of top-ranked residues to show from the networks:", 
                   value = 340, min = 5, max = 356, step = 1),
      p("The protein has 356 residues (the Delta-network is only calculated for 340) and thus networks can have up to this number of data points.")),
  
   column(2, wellPanel(actionButton("show_H1", "Update heatmap")))
  ),
  
  InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(heatmap_id="H1", layout="1|2|3", width1=1000, height1=900, width2=1000)
)