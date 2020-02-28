#' @export
plot_excyte <- function(excyte_obj,cut_top_99th=T,show_perc=T,alpha=0.5){
  #generate umaps for selected channels and save plots in a list
  umap_channels <- plot_umap(umap_2D = excyte_obj$umap_obj$umap_2D,
                             processed_fcs_df = excyte_obj$processed_fcs_obj,
                             channels =  excyte_obj$umap_obj$channels_used,
                             cut_top_99th = cut_top_99th,
                             alpha=alpha)
  #plot umap to visualize phenograph memberships
  umap_phenograph <- plot_phenograph(umap_2D = excyte_obj$umap_obj$umap_2D,
                                     phenograph_obj = excyte_obj$pheno_obj,
                                     alpha=alpha)

  #display a heatmap of channels intensities for each phenograph membership
  heatmap <- plot_cluster_profile_heatmap(pheno_obj = excyte_obj$pheno_obj,
                                          angle=45,
                                          cut_top_99th = T,
                                          show_perc = show_perc)
  #display cell repartition according to phenograph memberships for each channels and save plots in a list
  ridges <- plot_ridge(excyte_obj$pheno_obj)
  return(list("umap_channels"=umap_channels,"umap_phenograph"=umap_phenograph,"heatmap"=heatmap,"ridges"=ridges))
}
