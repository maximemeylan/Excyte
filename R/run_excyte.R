#' Master function to run preprocessing,phenograph clustering and umap computation
#' @param fcs_dir directory or vector containing fcs files to be used
#' @param downsampling numeric indicating the number event to randomly select from each fcs, if the number of events request is bigger than the number of event in the  fcs, all event are selected
#' @param channels vector containing channels to select. Can be "all" to select all channels, "with_desc" to select channels with a marker description or a vector a channels.
#' @param k numeric indicating the number of neighbor for phenograph and umap computation
#' @export
run_excyte <- function(fcs_dir,
                       downsampling=3000,
                       channels="all",
                       k=30){
  #pre-process fcs
  processed_fcs_obj <- pre_process_fcs(fcs_dir = fcs_dir,downsampling = downsampling)
  #compute phenograph membership for each event
  pheno_obj <- compute_phenograph(processed_fcs_obj,channels = channels,k = k)
  #compute umap for each event
  umap_obj <- compute_umap(processed_fcs_obj,channels = channels,k=k)
  return(list("processed_fcs_obj"=processed_fcs_obj,"pheno_obj"=pheno_obj,"umap_obj"=umap_obj))
}
#' Rerun the excyte pipeline on selected phenograph clusters
#' @param excyte_obj list of object obtained from an initial run with the excyte pipeline
#' @param channels vector containing channels to select. Can be "all" to select all channels, "with_desc" to select channels with a marker description or a vector a channels.
#' @param k numeric indicating the number of neighbor for phenograph and umap computation
#' @param downsampling numeric indicating the number event to randomly select from each fcs, if the number of events request is bigger than the number of event in the  fcs, all event are selected
#' @param clusters_id vector of character containing the ID of the phenograph clusters to rerun the excyte pipeline on

#' @export
rerun_excyte <- function(excyte_obj,
                         clusters_id=NA,
                         downsampling=3000,
                         channels="all",
                         k=30){
  if(is.na(clusters_id)){
    stop("Please submit clusters ID")
  }
  message("Excyte re-running with selected clusters: ",paste0(clusters_id," "))
  event_to_select <- excyte_obj$pheno_obj$processed_fcs$Phenograph_membership %in% clusters_id
  excyte_obj$processed_fcs_obj$processed_fcs <- excyte_obj$processed_fcs_obj$processed_fcs[event_to_select,]
  #compute new phenograph membership for selected events
  pheno_obj <- compute_phenograph(processed_fcs_obj = excyte_obj$processed_fcs_obj,channels = channels,k = k)
  #compute umap for selected events
  umap_obj <- compute_umap(excyte_obj$processed_fcs_obj,channels = channels,k=k)
  return(list("processed_fcs_obj"= excyte_obj$processed_fcs_obj,"pheno_obj"=pheno_obj,"umap_obj"=umap_obj))
}
#' Compute phenograph membership for each event
#' @param processed_fcs_obj list containing a datraframe of processed intensities for each event and informations of channel used
#' @param channels vector containing channels to select. Can be "all" to select all channels, "with_desc" to select channels with a marker description or a vector a channels.
#' @param k numeric indicating the number of neighbor for phenograph and umap computation
#' @import Rphenograph
#' @export
compute_phenograph <- function(processed_fcs_obj,channels=c("all","with_desc")[1],k=30){
  processed_fcs <- query_extract(processed_fcs_obj,channels=channels)
  channels_to_use <- setdiff(colnames(processed_fcs),"sample_id")
  message("\nComputing Phenograph clustering with channels: \n",paste0(channels_to_use,collapse = "\t"))

  #compute phenograph
  phenograph_obj <- Rphenograph(processed_fcs[,channels_to_use],k = k)
  processed_fcs$Phenograph_membership <- factor(paste0("c",phenograph_obj[[2]]$membership))
  phenograph_perc <- t(sapply(unique(processed_fcs$sample_id),function(y){
    all_pop <- table(processed_fcs[processed_fcs$sample_id == y,"Phenograph_membership"])
    perc <- all_pop/sum(all_pop)
  }))
  return(list("phenograph"=phenograph_obj,"phenograph_percentage"=phenograph_perc,"processed_fcs"=processed_fcs))
}
#' Compute Umap coordinates for each event
#' @param processed_fcs_obj list containing a datraframe of processed intensities for each event and informations of channel used
#' @param channels vector containing channels to select. Can be "all" to select all channels, "with_desc" to select channels with a marker description or a vector a channels.
#' @param k numeric indicating the number of neighbor for phenograph and umap computation
#' @import umap
#' @importFrom stats setNames
#' @export

compute_umap <- function(processed_fcs_obj,channels=c("all","with_desc")[1],k=30){
  processed_fcs<- query_extract(processed_fcs_obj,channels=channels)
  channels_to_use <- setdiff(colnames(processed_fcs),"sample_id")
  message("\nComputing Umap with channels: \t",paste0(channels_to_use,collapse = "\t"))
  #compute umap
  umap_obj <- umap(processed_fcs[,channels_to_use],method = "umap-learn",k=k)
  umap_obj_2D <- setNames(data.frame(umap_obj$layout,check.names = F),c("X","Y"))
  return(list("umap_obj"=umap_obj,"umap_2D"=umap_obj_2D,"channels_used"=channels_to_use))
}

