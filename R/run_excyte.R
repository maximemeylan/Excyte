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
