#' @import Rphenograph

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
