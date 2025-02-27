library(tidyverse)
library(ComplexHeatmap)
library(shiny)
library(circlize)

sort_by_splits <- function(data, sort, splitvars) {
  for (splitvar in rev(splitvars)) {
    data <- data %>%
      group_by(across(all_of(splitvar))) %>%
      arrange(across(all_of(sort), ~desc(mean(., na.rm=T)))) %>%
      ungroup()
  }
  data
}



get_analysis <- function(analyzed, splits, sort) {
  analyzed %>%
    select(!contains("retrieved")) %>%
    unite(col="rowname", c(Package, Filter, Analysis), sep = "_", remove = F) %>%
    mutate(`Atom or angle` = str_replace(`Atom or angle`, "\\(.*", ""),
           Analysis = str_replace(Analysis, "_weight", ""),
           None = 0
    ) %>%
    sort_by_splits(sort, splits) %>%
    group_by(across(all_of(splits))) %>%
    arrange(desc(across(all_of(sort))), .by_group = T) %>%
    ungroup() %>%
    column_to_rownames()
}





source("arrestin_data.R")

cons_split <- c()
for (row in 1:nrow(top_cons)) {
  cons_split <- c(cons_split, c(rep(row, length(top_cons[row,]$ini:top_cons[row,]$end))))
}

cons_anno <- as.numeric(all_cons >= 8)





get_matrix <- function(to_analyze, dataset) {
  missing_res <- setdiff(cons$Residue, unique(to_analyze$Residue))
  if (dataset == "delta") {
    a <- 2
    b <- 1
  } else {
    a <- 1
    b <- 0
  }
  
  to_analyze %>%
    select(!c(`rank`, hit)) %>%
    pivot_wider(names_from = Residue, values_from = value,) %>%
    
    add_column(!!! missing_res) %>% 
      rename_with(.fn = ~gsub('"', "", ., fixed = T), ) %>%
      mutate(across(all_of(missing_res), ~ NA)) %>%
    unite(col="rowname", c(Package, Filter, Analysis), sep = "_") %>%
    column_to_rownames() %>%
    select(4:ncol(.)) %>%
    t() %>% data.frame() %>%
      mutate(across(everything(), ~ a*( (. - min(., na.rm=T)) / (max(., na.rm=T) - min(., na.rm=T)))-b )) %>%
    t() %>%
    as.matrix()
}


order_matrix <- function(mat, analysis) {
  mat[rownames(analysis),order(map_int(colnames(mat), ~as.integer(str_split(., ":")[[1]][2])))]
}







get_heatmap = function(mat, analysis, splits, sort, annos) {
  
  split <- analysis %>%
    unite("splits", all_of(splits)) %>%
    pull(splits)
  
  ans <- annos %>% 
    purrr::set_names() %>% 
    map(~str_replace(pull(analysis, .), "\\(.*", ""))
  
  matcolmeans = colMeans(mat, na.rm = T)
  
  hits <- Heatmap(
    mat, name="hits",
    cluster_rows = F, cluster_columns = F, 
    cluster_row_slices = F, cluster_column_slices = F,
    show_row_dend = F, show_column_dend = F,
    show_row_names = F, show_column_names = F,
    row_split = split, row_gap = unit(0, "mm"), row_title=NULL,
    column_title = paste0("Protein residues in order"),
    right_annotation = exec(rowAnnotation, !!!c(ans, list(annotation_name_side = "top", annotation_name_rot = 45,
                                                          annotation_name_offset = unit(rep(0.75, length(ans)), "mm")
							  ))),
    top_annotation = columnAnnotation(Regions = unname(cbind(regions1v, regions2v)),
                                      `Secondary structure` = ssv,
                                      `ConSurf score (1-9)` = all_cons,
                                      `ConSurf > 8` = anno_simple(cons_anno),
                                      `Res. average` = as.array(matcolmeans),
                                      col = list(`ConSurf > 8` = setNames(c("white", "black"), c("0", "1")),
                                                 `Regions` = setNames(palette.colors()[c(2:4,1,5)],#list(1, "white", 2, 3, "black", 4),
                                                                      na.omit(unique(c(regions1v, regions2v)))),
                                                 `Res. average` = colorRamp2(c(min(matcolmeans, na.rm = T), mean(matcolmeans, na.rm = T), max(matcolmeans, na.rm = T)), 
                                                                             c("blue", "white", "red"))),
                                      annotation_name_offset = unit(c(0, -1.25, -1.25, -1.25, -1.25), "mm"),
                                      annotation_height = unit(c(1, 1, 1, 1, 1), c("null", "null", "null", "null", "null")),
                                      na_col = "white",
                                      annotation_name_rot = 45
                                      )
  )
  
  
  efs <- Heatmap(analysis %>% select(contains("norm_EF")) %>% as.matrix(), name="efs",
                 cluster_columns = F, show_row_dend = F, show_row_names = F, 
                 width=unit(40,'pt'), column_title = "Top %\nEFs", column_names_side = "top"
  )
  
  
  nefs <- Heatmap(analysis %>% select(contains("NEF")) %>% as.matrix(), name="nefs",
                  cluster_columns = F, show_row_dend = F, show_row_names = F,
                  width=unit(40,'pt'), column_title = "Top N\nEFs", column_names_side = "top"
  )
  
  
  auc <- Heatmap(analysis %>% select(contains("AUC")) %>% as.matrix(), name="auc",
                   cluster_columns = F, show_row_dend = F, show_row_names = F, column_title_rot = 90,
                   width=unit(10,'pt'), column_labels = c("AUC"), column_names_side = "top"
  )
  
  
  corr <- Heatmap(analysis %>% select(contains("Spearman corr")) %>% as.matrix(), name="corr",
                  cluster_columns = F, show_row_dend = F, show_row_names = F, column_title_rot = 90,
                  width=unit(10,'pt'), column_labels = c("Spearman corr."), column_names_side = "top"
  )
  
  
  
  hm <- draw(nefs + efs + auc + corr + hits, main_heatmap="hits",
             merge_legend = T,
             heatmap_legend_side = "bottom")
  
  return(hm)
  
  
}





get_show_scatterplot = function (mat) {
  return(
    function(df, output) {
      output$scatterplot = plotly::renderPlotly({
        if(is.null(df)) {
          plotly::plot_ly() %>%
            plotly::layout(annotations = list(x = 1, y = 1, text = "Click on a heatmap cell.", showarrow=F))
        } else {
          i = df$row_index
          
          x = all_cons
          y = mat[i,]
          
          plotly::plot_ly(type="scatter") %>%
            plotly::add_trace(x=x, y=y, hoverinfo="text", showlegend=F,
                              text=paste(paste("Res:", colnames(mat)), 
                                         paste("ConSurf:", x), 
                                         paste("Network value:", sprintf('%.3f', y)), 
                                         sep="\n")) %>%
            plotly::layout(title = list(text = paste0("Spearman correlation = ",
                                                      sprintf('%.3f', cor(x, y, method="spearman", 
                                                                          use="pairwise.complete.obs"))),
                                        yanchor = "bottom"),
                           xaxis = list(title = "ConSurf score"),
                           yaxis = list(title = rownames(mat)[i]),
                           margin = list(t = 100)
            )
          }
      })
    }
  )
}
