library(tidyverse)
library(ComplexHeatmap)
library(shiny)



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




get_matrix <- function(to_analyze) {
  to_analyze %>%
    select(!c(value, Residue)) %>%
    pivot_wider(names_from=`rank`, values_from = hit,) %>%
    unite(col="rowname", c(Package, Filter, Analysis), sep = "_") %>%
    column_to_rownames() %>%
    select(4:ncol(.)) %>%
    as.matrix()
}




get_heatmap = function(to_analyze, analyzed, splits, annos, sort, numcols) {
  
  mat <- get_matrix(to_analyze)
  
  analysis <- get_analysis(analyzed, splits, sort)
  
  
  split <- analysis %>%
    unite("splits", all_of(splits)) %>%
    pull(splits)
  
  ans <- annos %>% 
    purrr::set_names() %>% 
    map(~str_replace(pull(analysis, .), "\\(.*", ""))
  
  
  hits <- Heatmap(
    mat[rownames(analysis),1:numcols], name="hits",
    cluster_rows=F, cluster_columns = F, 
    cluster_row_slices = F, cluster_column_slices = F,
    show_row_dend = F, show_column_dend = F,
    show_row_names = F, show_column_names = F,
    col = structure(1:2, names = c("0", "1")),
    row_split = split, row_gap = unit(0, "mm"), row_title=NULL, 
    column_title = paste0("Top ", numcols, "-ranked residues in each network, colored by conservation"),
    right_annotation = exec(rowAnnotation, !!!c(ans, list(annotation_name_side = "top")))
  )
  
  
  
  
  efs <- Heatmap(analysis %>% select(contains("norm_EF")) %>% as.matrix(), name="efs",
                 cluster_columns = F, show_row_dend = F, show_row_names = F, #show_column_names = F,
                 width=unit(55,'pt'), column_title = "Top %\nEFs", column_names_side = "top"
  )
  
  
  nefs <- Heatmap(analysis %>% select(contains("NEF")) %>% as.matrix(), name="nefs",
                  cluster_columns = F, show_row_dend = F, show_row_names = F, #show_column_names = F,
                  width=unit(55,'pt'), column_title = "Top N\nEFs", column_names_side = "top"
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
             adjust_annotation_extension = F, merge_legend = T,
             heatmap_legend_side = "bottom", annotation_legend_side = "right")
  
  return(hm)
  
  
}
