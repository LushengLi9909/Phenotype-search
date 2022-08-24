library(ontologyIndex)
library(dplyr)
require(visNetwork)
mpo <- get_ontology("data/MPheno_OBO.ontology")
phenotype_to_genes <- read.delim("data/Phenotype_MPO.txt",header = TRUE,sep=",")
all_results_merged <- readRDS("data/mpo_results_merged.rds")

disease_descriptions <- read.delim("data/Phenotype_descriptions.txt",header = FALSE)
colnames(disease_descriptions) = c("MPID","description")
rownames(disease_descriptions) = disease_descriptions$MPID

#' Get mpo term definition
#'
#' This gets the disease description from a data frame of disease descriptions.
#' The rows names of the data frame are the mpo ID and the description column is
#' called "description". It also adds new lines to the description so that the
#' hover box in the web app does not get too wide. This is done by calling the
#' \code{newlines_to_definition} function.
#'
#' @param ontologyId The mpo Id of the term (string)
#' @param disease_descriptions A data frame of disease descriptions corresponding to each mpo Id
#'
#' @return The disease description with new lines added.
#'
#' @examples
#' mpo_get_term_definition("HP:123456", disease_descriptions)
#'
#' @export
##
mpo_get_term_definition <- function(ontologyId, disease_descriptions) {
  definition = disease_descriptions[ontologyId,"description"]
  definition = newlines_to_definition(definition)
  return (definition)
}


#' Get mpo term description for a list of terms
#'
#' Need to redo this without a for loop (use \code{lapply} or something). This just
#' applys the \code{mpo_get_term_definition} function to a character vector of terms.
#'
#' @param ontologyId_list A character vector of mpo Ids
#' @param disease_descriptions A data frame of disease descriptions for all mpo Id
#'
#' @retun A named vector of disease descriptions, with mpo Id as names and descriptions
#' as values.
#'
#' @examples
#' mpo_term_definition_list(mpo_terms_char_vector, Disease_description_df)
#'
#' @export
##
mpo_term_definition_list <- function(ontologyId_list, disease_descriptions) {
  term_details <- c()
  for (term in ontologyId_list) {
    term_details[term] <- mpo_get_term_definition(term, disease_descriptions)
  }
  return (term_details)
}


#' Add new lines to disease description
#'
#' Adds new lines to the description so that hover boxes dont get too wide.
#'
#' @param definition A disease description string
#' @param line_length A integer representing the desired words per line.
#'
#' @returns The disease description with newline symbols added every nth word.
#'
#' @examples
#' newlines_to_definition(disease_description, 10)
#' @export
##
newlines_to_definition <- function(definition, line_length = 10) {
  definition = strsplit(definition, split = " ")[[1]]
  if (length(definition) > line_length) {
    remainder = length(definition) %% line_length
    n_new_lines = floor((length(definition)/line_length))
    new_line_index = seq(line_length,(n_new_lines*line_length),line_length)
    definition[new_line_index] = paste0("\n", definition[new_line_index])
  }
  definition = paste(definition,collapse = " ")
  return(definition)
}

#' Get absolute ontology level
#'
#' This gets the absolute ontology level of a term (without consideration for the
#' particular subset of the data you are looking at, as in the find_parent function).
#'
#' @param mpo The mpo ontology data object
#' @param term_id mpo term ID <string>
#' @example get_ont_level(mpo,"HP:0000003")
#' @return returns the ontology level <numeric>
#' @export
get_ont_level = function(mpo,term_id) {
  children = unique(setdiff(unlist(mpo$children[term_id]), term_id))
  if (length(children) == 0) {
    return(0)
  } else {
    return(1 + get_ont_level(mpo,children)) #<- recursion..
  }
}

#' Subset RD EWCE Results
#'
#' This subsets  the Rare disease EWCE results by cell type, q threshold and fold change.
#'
#' @param phenotype_to_genes The list of mpo terms with their assocaited gene lists taken from mpo website
#' @param all_results_merged The dataframe of RD EWCE Results
#' @param mpo The mpo ontology data object
#' @param cell_type A string representing the cell type of interest.
#' @param q_threshold The q threshold. The subset of results will have a q lower than this
#' @param fold_threshold The fold change threshold. The subest of results will have a fold change greater than this.
#'
#' @returns A data frame of results taken from the main data frame of results
#' @export
#
subset_phenos = function(phenotype_to_genes, all_results_merged,mpo,branch, q_threshold, fold_threshold) {
  if(branch != ""){
  if(!(branch %in% unique(all_results_merged$Phenotype))){
    Terms = process_search_terms(branch)
    Phenos <- as.character(unique(phenotype_to_genes$Phenotype[stringr::str_detect(phenotype_to_genes$Phenotype, pattern = Terms)]))
    phenos = all_results_merged[all_results_merged$Phenotype %in% Phenos & all_results_merged$q <= q_threshold & all_results_merged$fold_change >= fold_threshold,]
  } else {
  branch_id = paste(mpo$id[match(branch,mpo$name)])
  branch_descendants = paste(get_descendants(mpo,branch_id))
  branch_descendants_names = paste(mpo$name[branch_descendants])
  phenos = all_results_merged[all_results_merged$Phenotype %in% branch_descendants_names & all_results_merged$q <= q_threshold & all_results_merged$fold_change >= fold_threshold,]
}
   if(length(phenos$CellType) != 0){
  phenos = add_mpo_termid_col(phenos, phenotype_to_genes, mpo)
  phenos = phenos[!is.na(phenos$mpo_term_Id),]
  #phenos <- phenos %>% group_by(mpo_term_Id) %>% filter (! duplicated(mpo_term_Id))
  ontlevz <- c()
  for (p in unique(phenos$mpo_term_Id)) {
    ontlevz[p] <- get_ont_level(mpo,p)
  }
  phenos$ontlvl <- ontlevz[phenos$mpo_term_Id]
  mpo_term = mpo_term_definition_list(ontologyId_list = phenos$mpo_term_Id,
                                  disease_descriptions = disease_descriptions)
  phenos$descriptions = mpo_term[phenos$mpo_term_Id]
  return (phenos)}
}}



add_mpo_termid_col = function(cells, phenotype_to_genes, mpo) {
  mpotermID = c()
  for (p in cells$Phenotype){
    termid = phenotype_to_genes$MPID[phenotype_to_genes$Phenotype == p][1]
    mpotermID = append(mpotermID, termid)
  }
  cells$mpo_term_Id =mpotermID
  return(cells)
}
#' Make network object
#'
#' This uses the network package to coerce the  into a
#' network object. It also adds the fold change, label, and relative ontology level
#' parameters to each node in the network.
#'
#' @param phenos The subset of the results to be plotted
#' @param mpo The mpo ontology data object
#'
#' @returns A ggnetowrk graph/ network object of a subset of the  EWCE results.
#' @export
visnetwork_plot = function(phenos, mpo) {
  if(length(phenos) != 0){
  adjacency = data.frame()
  mpo_id = phenos$mpo_term_Id
  phenotype = phenos$Phenotype
  
  for( i in 1:length(mpo_id)){
    for( d in 1:length(mpo_id)){
      if(mpo_id[i] %in% (mpo$children[mpo_id[d]][[1]])){
        adjacency = rbind(adjacency,data.frame("from"=phenotype[d],"to"=phenotype[i]))
      }} 
  }
  for( c in unique(phenos$CellType)){
    for( p in phenos[phenos$CellType == c , "Phenotype"]){
      adjacency = rbind(adjacency,data.frame("from"=c,"to"=p))
    }
  }
  adjacency = adjacency[!duplicated(adjacency),]
  phenos_lev = unique(phenos[,c("Phenotype","ontlvl")])
  Npheno = length(unique(phenos$Phenotype))
  Ncell = length(unique(phenos$CellType))
  vertices <- data.frame(name=c(unique(phenos$Phenotype),unique(phenos$CellType)),
                         molecular=c(rep('Phenotype',Npheno),rep('Celltype',Ncell)),
                         heirarchy =c(as.numeric(phenos_lev$ontlvl),rep(0,Ncell)),
                         term_definitions = c(as.character(unique(phenos$descriptions)),rep('Celltype',Ncell)),
                         label = c(as.character(unique(phenos$Phenotype)),as.character(unique(phenos$CellType)))
  )
  nodes <- data.frame( id=c(unique(phenos$Phenotype),unique(phenos$CellType)),
                       label=c(as.character(unique(phenos$Phenotype)),as.character(unique(phenos$CellType))),
                       value =c(as.numeric(phenos_lev$ontlvl),rep(0,Ncell)),
                       shape=c(rep("box",Npheno),rep("ellipse",Ncell)),
                       color=c(rep("lightblue",Npheno),rep("orange",Ncell)))
  
  edges <- data.frame(from = c(adjacency$from),to = c(adjacency$to),arrows = c("to"))
  visNetwork(nodes, edges)
  }  
}

# PLOT USING GGNETWORK
#' Plot RD EWCE results subset as interactive network plot
#'
#' This coerces the network object into a plot with ggnetwork (which assigns
#' x, y coordinates to each node in the network). The hover box with results and
#' disease description for each node is also added here in the
#' \code{pheno_ggnetwork$hover} column. Once the x and y coordiantes have been
#' added, it can be plot using ggplot2.
#' @param phenoNet The network object created using \code{create_network_object}
#' @param phenos The subset of results to be plotted (data frame)
#' @param disease_descriptions The data frame of all disease descriptions, This is
#' where new lines are added using the mpo_term_definition_list function.
#' @returns A interactive plot of the network of phenotypes in the selected subset of results.
#' @export


#' Create interactive network plot start to finish (including subset data)
#'
#' This puts all the functions together from gettig the subest of results to creating
#' the final interactive plot.
#'
#'  @param penotype_to_genes The phenotype gene lists taken from the mpo
#'  @param all_results_merged The RD EWCE Results
#'  @param mpo The mpo ontology data object
#'  @param disease_descriptions The dataframe of all disease descriptions in the mpo
#'  @param cell_type The cell type of interest to be plotted
#'  @param q_threshold The q value threshold for the subset of results to be plotted
#'  @param fold_threshold The minimum fold change in specific expression for the subest of results to be plotted
#'
#'  @return A interactive network plot of the selected subset of results from RD EWCE analysis
#'  @export
visnetwork_plot_full = function(phenotype_to_genes, all_results_merged, mpo, disease_descriptions,branch, q_threshold, fold_threshold){
  if(branch != ""){
  phenos = subset_phenos(phenotype_to_genes, all_results_merged,mpo,branch=branch, q_threshold = q_threshold, fold_threshold = fold_threshold)
  network_plot = visnetwork_plot(phenos,mpo)
  return(network_plot)
  }
}


#visnetwork_plot_full(phenotype_to_genes, all_results_merged, mpo, disease_descriptions,branch= "abnormal consumption behavior", q_threshold= 0.05, fold_threshold=1)
