
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
