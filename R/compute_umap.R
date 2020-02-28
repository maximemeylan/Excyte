#' @import umap
#' @importFrom stats setNames
#' @export

compute_umap <- function(processed_fcs_obj,channels=c("all","with_desc")[1],samples="all",k=30){
  processed_fcs<- query_extract(processed_fcs_obj,channels=channels)
  channels_to_use <- setdiff(colnames(processed_fcs),"sample_id")
  message("\nComputing Umap with channels: \t",paste0(channels_to_use,collapse = "\t"))
  #compute umap
  umap_obj <- umap(processed_fcs[,channels_to_use],method = "umap-learn",k=k)
  umap_obj_2D <- setNames(data.frame(umap_obj$layout,check.names = F),c("X","Y"))
  return(list("umap_obj"=umap_obj,"umap_2D"=umap_obj_2D,"channels_used"=channels_to_use))
}
