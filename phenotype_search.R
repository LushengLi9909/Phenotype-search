require(stringr)
require(cowplot)
require(reshape2)
require(scales)
require(ggplot2)
library(RColorBrewer)
phenotype_to_genes <- read.delim("data/Phenotype_MPO.txt",header = TRUE,sep=",")
all_results_merged <- readRDS("data/mpo_results_merged.rds")
cell_order <- read.delim("data/cell_order.txt",header = TRUE,sep=",")
#' Create Phenotype keyword search pattern 
#' 
#' This creates a regex search pattern to find all phenotypes containing the keywords 
#' entered by the user. The search expression is not case sensitive and 
#' uses an "or" operator, so the phenotype only needs to contain one of the keywords.
#' The function also removes any spaces. 
#' 
#' @param keywords phenotype related keywords separated by "," <string>
#' @returns  A regex pattern used to search for keywords <string>
#' @example 
#' process_search_terms("muscle, fatigue, heart")
#' @export
process_search_terms <- function(keywords) {
  if(keywords !=""){
  keywords <- strsplit(keywords, ",")[[1]]
  output = c()
  for (i in seq(1, length(keywords))){
    cur = keywords[i]
    while (substr(cur, 1,1) == " ") {
      cur = substr(cur, 2, nchar(cur))
    }
    while (substr(cur, nchar(cur), nchar(cur)) == " ") {
      cur = substr(cur, 1, nchar(cur)-1)
    }
    output = append(output, paste0("(?i)", cur))
  }
  output <- paste(output, collapse = "|")
  return(output)
  }
}


#' Subset the RD EWCE results using phenotype keywords 
#' 
#' This finds all results related to phenotypes that contain the kewords entered 
#' into the \code{process_search_terms} function. It further subests by
#' fold change and q value thresholds. 
#' 
#' @param Terms the regex search pattern generated by \code{process_search_terms} <string>
#' @param q_threshod q value maximum threshold <numeric>
#' @param fold_threhold fold change minimum threshold <numberic>
#' @param min_sd_from_mean Z score threshold <numeric> - possibly remove?
#' 
#' @returns A data frame of subest of RD EWCE results 
#' @export 
keyword_search_df <- function(Terms, 
                              q_threshold = 0.05,
                              fold_threshold = 1,
                              min_sd_from_mean = 0) {
  if(Terms !=""){
  # remove as charcater below ?
  Phenos <- as.character(unique(phenotype_to_genes$Phenotype[stringr::str_detect(phenotype_to_genes$Phenotype, pattern = Terms)]))
  DF <- all_results_merged[all_results_merged$Phenotype %in% Phenos &
                             all_results_merged$q <= q_threshold & 
                             all_results_merged$fold_change >= fold_threshold &
                             all_results_merged$sd_from_mean >= min_sd_from_mean, ]
  if(length(DF$CellType) != 0){
  return(DF)}
  }
}


#' Plot bar chart of RD EWCE subest 
#' 
#' Should redo this figure. Possibly reuse a better one from the report
#' 
#' @param DF The RD EWCE results subset to be visualised <dataframe>
#' @param keywords the keyword string (used to create title)
#' 
#' @returns a bar chart <ggplot>
#' @export

plot_phenotype_counts <- function(DF, keywords) {
  if(keywords !="" ){
    if(length(DF) != 0){
      if(length(DF$CellType) != length(unique(DF$CellType))){
        DF$CellType = factor(DF$CellType,levels = cell_order$x)
        plot <- ggplot(DF, aes(x = CellType)) + 
    geom_bar(fill="#69b3a2", color="#e9ecef") + 
    # cowplot::theme_cowplot() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle=45,hjust=1,vjust = 1, size = 10,face = "bold"),axis.title.x = element_blank(),
          plot.title = element_text(margin=margin(0,0,30,0),size = 12)) + 
    labs(title =paste0("Phenotype enrichment counts associated with: ",keywords), 
          y = "Number of Enrichments") 
    #scale_y_continuous(breaks = scales::pretty_breaks(), 
    #                   expand = expansion(mult = c(0, .1))) +
    #theme(legend.position = "top", legend.key.size = unit(0.25, "cm"),
    #      legend.text = element_text(size = 10), 
    #      legend.title = element_blank(),
    #      )+
    #scale_fill_continuous(name = "Fold Change")
  return (plot)} else{
    plot <- ggplot(DF, aes(x = CellType,y = fold_change)) + 
      geom_bar(stat="identity",fill="#69b3a2", color="#e9ecef") + 
      # cowplot::theme_cowplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle=45,hjust=1,vjust = 1, size = 10,face = "bold"),axis.title.x = element_blank(),
            plot.title = element_text(margin=margin(0,0,30,0),size = 12)) + 
      labs(title =paste0("Phenotype enrichment associated with: ",keywords), 
            y = "Fold Change")
    return(plot)
  }}}
}