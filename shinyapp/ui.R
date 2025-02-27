fluidPage(
  
  h1("Allosteric networks of \u03B2-arrestin 1 computed with AlloViz"), br(),
  
  navlistPanel(
    tabPanel("Allosteric networks",
  p("Proteins are highly dynamic biomolecules that carry out essential functions in organisms, and they can undergo regulating conformational changes driven by mutations, interactions, and various factors. Allostery, which involves long-range signal transmission, is essential for this. Molecular Dynamics (MD) simulations offer a way to study protein dynamics at an atomic scale, and network analysis using graph theory provides insights into residue interactions. AlloViz unifies various network construction methods and simplifies network analysis and visualization, including delta-networks."),
  h3("Computing networks"),
  HTML("<p>The purpose of computing an allosteric network from an MD trajectory is to transform the data into a single inter-residue value for each residue pair in the protein, which represents their relationship. The position coordinates of the residues’ atoms over time can be used to calculate correlation metrics such as Pearson’s correlation coefficient or other generalized correlation coefficients that quantify the relationship between their movements throughout the trajectory. The time-series of the different dihedral angles of the residues, both backbone (Phi and Psi) and side-chain (Chi angles 1 through 4, depending on the side chain), can also be extracted from the simulations, and the correlations between the different residues’ dihedrals calculated using the same metrics. In addition, there exist different ways of modelling residue contacts to calculate contact frequencies, and of computing their interaction energies as well. A complete list of methods available in AlloViz can be accessed <a href='https://alloviz.readthedocs.io/en/latest/table.html'>here</a>.</p>"),
  h3("Analyzing networks"),
  p('Once networks are created, they can undergo diverse processing, including filtering and graph theory-based analyses. AlloViz offers different analysis options, such as the complete network ("All" filter), contact-based filtering using GetContacts ("GetContacts_edges"), spatial distance filtering ("Spatially_distant"), sequence distance filtering ("No_Sequence_Neighbors"), and combinations thereof. Among common network analyses, centrality is key, namely betweenness centrality ("btw") and current-flow betweenness centrality ("cfb").'),
  h3("Delta-networks"),
  p('Delta-networks, calculated from the subtraction of two allosteric networks, allow to compare the differences in information flow among different structures, such as different activation states of a receptor. The global structure and communication of the two systems compared is expected to be alike (except for just a few conformational changes due to the perturbation), and thus most residue pair interactions will have similar weights. The calculation of the delta-network reduces this noise and helps to highlight the subtle differences between the two networks. Each of the network’s edge weights of one structure ("other") are subtracted from the corresponding edge weights of a "reference" one.')
  ),
  
  tabPanel("\u03B2-arrestin 1 allosteric networks and sequence conservation",
    p("As a demonstration of the data-intensive power of AlloViz, we computed all available combinations of network, filter (”All” (no filter), ”Spatially distant”, ”GetContacts_edges” and the combination of ”GetContacts_edges” and ”No_Sequence_Neighbors”) and node analysis (betweenness centrality and current-flow betweenness centrality) for β-arrestin 1 in the inactive (apo) and active (V2RppWT-bound) conformations, and the delta-network thereof. We then investigated which of these analysis combinations was able to highlight better the residues that are important for the systems’ function. In this case, we chose to use residue sequence conservation scores provided by ConSurf to quantify the residues’ functional relevance, and systematically measured how much the computed node centralities’ top results (intuitively, the most important nodes in the network) are enriched with conserved residues, and how they correlate with the scores."),
    h3("Enrichment factors"),
    p("Enrichment factors (EFs) help assess whether specific entities are more prevalent or concentrated within a subset than would be expected by chance alone. In this case, we measure how many conserved residues are in the top X% (norm_EF) or top # (NEF) results of each combination of network, analysis and filter (ordering residues by the absolute value of the centrality analysis applied). They are calculated with the formula (c/n)/(C/N) wherere c and C are the number of conserved residues and n and N are the total number of residues in the top X%/# results (c, n) and in the complete dataset (C, N), respectively. EFs of the top X% (norm_EF) are normalized by dividing by the best possible EF achievable calculated as if all top results were conserved residues."),
    h3("AUC"),
    p("The Area Under the Curve (AUC) of the Receiver Operating Characteristic (ROC) curve is a metric used to evaluate the performance of a classifier in terms of the true positive rate (sensitivity) and the false positive rate (1-specificity) that it achieves at different thresholds of the predicting variable. In this case, it measures how well the networks' centralities 'predict' whether the residues are conserved or not."),
    h3("Spearman correlation"),
    p("Spearman's correlation is a statistical measure that assesses the strength and direction of the relationship between two variables. Unlike the Pearson correlation coefficient, which measures linear relationships, Spearman's correlation focuses on capturing monotonic relationships that can be linear or nonlinear but still exhibit a consistent trend. For this case, it's used to investigate how well the relative values of centrality relate to conservation scores.")),
    
  widths=c(2,9)),
  
  br(), br(),
  
  radioButtons(inputId = "dataset", label = "Select a dataset:", inline = T,
               choices = setNames(c("apo", "wt", "delta"), #"mut"
                                  c("Apo bArr1", "bArr1-WTpp", "Delta-network"))), #bArr1-Mutpp
  
  br(),
  
  tabsetPanel(
    tabPanel("Table of networks' results",
                          DT::DTOutput("ROC_table"),
                           fluidRow(
                             column(6, plotOutput("ROC_plot")),
                             column(6, plotOutput("density_plot"))
                           )
             ),
    
    tabPanel("Networks' predictive power", source("Heatmap_rank_ui.R", local=T)
             ),
    
    tabPanel("Networks' predicted residues", source("Heatmap_aas_ui.R", local=T)
    )
              )
  
)