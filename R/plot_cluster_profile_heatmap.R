#' @import pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @export
#'
plot_cluster_profile_heatmap <- function(pheno_obj,
                                         palette="default",
                                         channels_order=NULL,
                                         memberships_order=NULL,
                                         display_marker_names=T,
                                         show_perc=T,
                                         cut_top_99th=T,
                                         ...){
  all_channels <- attr(pheno_obj$processed_fcs,"all_channels")
  processed_fcs <- pheno_obj$processed_fcs
  channels_to_use <-  setdiff(colnames(processed_fcs),c("sample_id","Phenograph_membership"))
  clusters_id <- unique(processed_fcs$Phenograph_membership)

  mean_intensity_pheno <- sapply(clusters_id,function(x){
    colMeans(processed_fcs[processed_fcs$Phenograph_membership==x,channels_to_use])
  })
  if(show_perc){
    mean_perc <- colMeans(pheno_obj$phenograph_percentage)
    colnames(mean_intensity_pheno)<- paste0(clusters_id," ",round(mean_perc[clusters_id],digits = 3)*100,"%")
  }else{
    colnames(mean_intensity_pheno)<- clusters_id
  }
  if(cut_top_99th){
    #normalized_mean_intensity_pheno <- apply(mean_intensity_pheno,2,function(x) norm_range(x,range=c(min(x),quantile(x,probs = 0.98))))
    normalized_mean_intensity_pheno <- apply(mean_intensity_pheno,2,function(x) norm_range(x,c(range=quantile(x,probs = 0.01),quantile(x,probs = 0.99))))

  }else{
    normalized_mean_intensity_pheno <- apply(mean_intensity_pheno,2,function(x) norm_range(x,range=c(0,1)))
  }
  if(display_marker_names==T){
    rownames(normalized_mean_intensity_pheno) <- sapply(rownames(normalized_mean_intensity_pheno), function(x)paste(x,all_channels[which(all_channels[,1]==x),2],sep=" / "))
  }
  if(palette=="default"){
    palette <- colorRampPalette(brewer.pal(name = "PuBu",n=9))(100)
  }
  if(!is.null(channels_order)&!is.null(memberships_order)){
    pheatmap(t(normalized_mean_intensity_pheno[channels_order,memberships_order]),
             color=palette,
             angle_col = 0,
             cluster_cols = F,
             silent = T,
             cluster_rows = F,
             ...)
  }else if(!is.null(channels_order)){
    pheatmap(t(normalized_mean_intensity_pheno[channels_order,]),
             color=palette,
             cluster_cols = F,
             silent = T,
             angle_col = 0,
             ...)
  }else if(!is.null(memberships_order)){
    pheatmap(t(normalized_mean_intensity_pheno[,memberships_order]),
             color=palette,
             cluster_rows = F,
             silent = T,
             angle_col = 0,
             ...)
  }else{
    pheatmap(t(normalized_mean_intensity_pheno),
             color=palette,
             silent = T,
             ...)
  }

}
