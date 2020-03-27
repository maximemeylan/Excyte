#' Master function to call all available plot functions
#'
#' @param excyte_obj list of object obtained from an initial run with the excyte pipeline
#' @param cut_top_99th boolean value indicating if the top 1 percent should be remove to compute scaling
#' @param show_perc boolean indicating if percentages for each cluster should be displayed
#' @param alpha numeric indicating the transparency level when plotting the umap
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
                                     phenograph_obj = excyte_obj$phenograph_obj,
                                     alpha=alpha)

  #display a heatmap of channels intensities for each phenograph membership
  heatmap <- plot_cluster_profile_heatmap(phenograph_obj = excyte_obj$phenograph_obj,
                                          angle=45,
                                          cut_top_99th = T,
                                          show_perc = show_perc)
  #display cell repartition according to phenograph memberships for each channels and save plots in a list
  ridges <- plot_ridge(excyte_obj$phenograph_obj)
  return(list("umap_channels"=umap_channels,"umap_phenograph"=umap_phenograph,"heatmap"=heatmap,"ridges"=ridges))
}

#' Plot function that display phenograph memberships for each event on a umap representation
#' @param phenograph_obj list containing result of phenograph clustering and processed fcs
#' @param umap_2D matrix detailing coordinates of each event on the umap
#' @param alpha numeric indicating the transparency level when plotting the umap
#' @import ggplot2
#' @importFrom stats quantile
#' @importFrom ggrepel geom_text_repel
#' @export
plot_phenograph <- function(umap_2D,phenograph_obj,alpha=0.5){
  processed_fcs <- phenograph_obj$processed_fcs

  #subset processed_fcs based on umap corrdinates (safer when downsampling_umap != NULL)
  processed_fcs <- processed_fcs[rownames(umap_2D),]

  memberships <- factor(processed_fcs$Phenograph_membership)
  #define centers of each cluster
  cluster_pos <- as.data.frame(t(sapply(unique(memberships),function(x){
    ids <- which(memberships == x)
    x_cluster <- mean(umap_2D[ids,1])
    y_cluster <- mean(umap_2D[ids,2])
    return(c(x_cluster,y_cluster))
  })))
  colnames(cluster_pos) <- c("xpos","ypos")
  cluster_pos$cluster_id <- factor(unique(memberships))
  umap_2D$colour <- memberships
  p <- ggplot(umap_2D,aes(x = X,y=Y, colour =colour ))
  p <- p + theme_bw()
  p <- p + geom_point(alpha=alpha,size=0.01)
  p <- p + labs(x="",y="",color="phenograph membership")
  p <- p + guides(colour = guide_legend(override.aes = list(size=10)))
  suppressWarnings(p <- p + ggrepel::geom_text_repel(data = cluster_pos, aes(x =xpos,y=ypos, label = cluster_id),show_guide  = F,colour='black' ))
  return(p)
}

#' Plot function to display channels intensities according to event memberships
#' @import ggplot2
#' @import ggridges
#' @param phenograph_obj list containing result of phenograph clustering and processed fcs
#' @export
plot_ridge <- function(phenograph_obj){
  processed_fcs <- phenograph_obj$processed_fcs
  all_channels <- attr(processed_fcs,"all_channels")
  channels_to_plot <- intersect(colnames(processed_fcs),all_channels[,"name"])
  all_plots <- lapply(channels_to_plot,function(x){
    melted_df<- melt_df(processed_fcs[,c(x,"Phenograph_membership")],var_to_group ="Phenograph_membership")
    p <- ggplot(melted_df, aes(x = value, y = groups,fill = groups))
    p <- p + geom_density_ridges(scale = 4, rel_min_height = 0.045,alpha = 0.85)
    p <- p + theme_ridges()
    p <- p + theme(axis.title.y = element_blank(),axis.title.x = element_blank())
    p <- p + labs(title=paste(x,all_channels[which(all_channels[,1]==x),2],sep=" / "))
    p <- p + scale_fill_discrete(guide=FALSE)
    return(p)
  })
  return(all_plots)
}

#' Plot function that display a Umap of scaled intensities for event
#' @param umap_2D matrix detailing coordinates of each event on the umap
#' @param processed_fcs_df dataframe of processed intensities for each event and informations of channel used as attribute
#' @param channels vector of channels to select. Can be "all" to select all channels, "with_desc" to select channels with a marker description or a vector a channels
#' @param cut_top_99th boolean value indicating if the top 1 percent should be remove to compute scaling
#' @param alpha numeric indicating the transparency level when plotting the umap
#' @import ggplot2
#' @importFrom stats quantile
#' @importFrom RColorBrewer brewer.pal
#' @export
plot_umap <- function(umap_2D,processed_fcs_df,channels=c("channels_used","all","with_desc")[1],cut_top_99th=T,alpha=0.5){
  processed_fcs<- query_extract(processed_fcs_df,channels=channels)
  all_channels <- attr(processed_fcs,"all_channels")
  channels_to_use <-  setdiff(colnames(processed_fcs),"sample_id")

  #subset processed_fcs based on umap corrdinates (safer when downsampling_umap != NULL)
  processed_fcs <- processed_fcs[rownames(umap_2D),]

  all_plot <- lapply(channels_to_use, function(channel){
    title <- paste(channel,all_channels[which(all_channels[,1]==channel),2],sep=" / ")
    p <- ggplot(umap_2D,aes(x=X,y=Y))
    p <- p + theme_bw()
    if(cut_top_99th){
      p <- p + geom_point(size=0.05,
                          alpha=alpha,
                          aes(colour = norm_range(processed_fcs[,channel],c(quantile(processed_fcs[,channel],probs = 0.01),quantile(processed_fcs[,channel],probs = 0.99)))))
    }else{
      p <- p + geom_point(size=0.05,
                          alpha=alpha,
                          aes(colour = processed_fcs[,channel]))
    }
    p <- p + scale_color_gradientn(colours =rev(brewer.pal(name = "Spectral",n=11)))
    p <- p + ggtitle(title)
    p <- p + labs(x="",y="",color="intensity")
    return(p)
  })
}

#' Plot function that display a heatmap of channels intensities for each phenograph cluster
#'
#' @param phenograph_obj list containing result of phenograph clustering and processed fcs
#' @param palette vector of colors to be used
#' @param channels_order vector of ordered channels names
#' @param memberships_order vector of ordered memberships id
#' @param display_marker_names boolean indicating if markers names should be displayed after the channel name
#' @param cut_top_99th boolean value indicating if the top 1 percent should be remove to compute scaling
#' @param show_perc boolean indicating if percentages for each cluster should be displayed
#' @param ... arguments to be passed to pheatmap
#'
#' @import pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @export
plot_cluster_profile_heatmap <- function(phenograph_obj,
                                         palette="default",
                                         channels_order=NULL,
                                         memberships_order=NULL,
                                         display_marker_names=T,
                                         show_perc=T,
                                         cut_top_99th=T,
                                         ...){
  all_channels <- attr(phenograph_obj$processed_fcs,"all_channels")
  processed_fcs <- phenograph_obj$processed_fcs
  channels_to_use <-  setdiff(colnames(processed_fcs),c("sample_id","Phenograph_membership"))
  clusters_id <- unique(processed_fcs$Phenograph_membership)

  mean_intensity_pheno <- sapply(clusters_id,function(x){
    colMeans(processed_fcs[processed_fcs$Phenograph_membership==x,channels_to_use])
  })
  if(show_perc){
    mean_perc <- colMeans(phenograph_obj$phenograph_percentage[,clusters_id])
    colnames(mean_intensity_pheno)<- paste0(clusters_id,"  ",round(mean_perc,digits = 3)*100,"%")
  }else{
    colnames(mean_intensity_pheno)<- clusters_id
  }
  if(cut_top_99th){
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

