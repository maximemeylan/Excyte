#' Master function to call all available plot functions
#'
#' @param excyte_obj list of object obtained from an initial run with the excyte pipeline
#' @param cut_top_99th boolean value indicating if the top 1 percent should be remove to compute scaling
#' @param show_perc boolean indicating if percentages for each cluster should be displayed
#' @param alpha numeric indicating the transparency level when plotting the umap
#' @export
plot_excyte <- function(excyte_obj,cut_top_99th=T,show_perc=T,alpha=0.5,threshold="tertile",title='both'){
  #generate umaps for selected channels and save plots in a list
  umap_channels <- plot_umap(umap_2D = excyte_obj$umap_obj$umap_2D,
                             processed_fcs_df = excyte_obj$processed_fcs_obj,
                             channels =  excyte_obj$umap_obj$channels_used,
                             cut_top_99th = cut_top_99th,
                             alpha=alpha,
                             title=title)
  #plot umap to visualize phenograph memberships
  umap_phenograph <- plot_phenograph(umap_2D = excyte_obj$umap_obj$umap_2D,
                                     phenograph_obj = excyte_obj$phenograph_obj,
                                     alpha=alpha)

  #display a heatmap of channels intensities for each phenograph membership
  heatmap <- plot_cluster_profile_heatmap(phenograph_obj = excyte_obj$phenograph_obj,
                                          angle=45,
                                          cut_top_99th = T,
                                          show_perc = show_perc)
  #display cell repartition according to phenograph memberships for each channels and save plots in a list (clusters)
  ridges_clusters <- plot_ridge(excyte_obj$phenograph_obj,
                       type = "clusters",
                       threshold = threshold,
                       channel_names = title,
                       downsampling = 1000
                       )
  #display cell repartition according to phenograph memberships for each channels and save plots in a list (channels)
  ridges_channels <- plot_ridge(excyte_obj$phenograph_obj,
                                type = "channels",
                                threshold = threshold,
                                channel_names = title,
                                downsampling = 1000
  )
  return(list("umap_channels"=umap_channels,"umap_phenograph"=umap_phenograph,"heatmap"=heatmap,"ridges_clusters"=ridges_clusters,"ridges_channels"=ridges_channels))
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
  p <- p + theme_classic()
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
#' @param threshold character, draw the median, tertile or quartile
#' @param channels vector of channels to use, default uses all channels
#' @param cluster_to_use vector of cluster to use, default uses all clusters
#' @param type character, plot selected channels for each clusters (type="clusters") or selected clusters for each channel (type="channels")
#' @param downsampling numeric indicating the number of events to plot
#' @param channel_names character, edit channels names accordingly
#' @export
plot_ridge <- function(phenograph_obj,
                       channels="all",
                       cluster_to_use="all",
                       type=c("channels","clusters")[1],
                       threshold=c("median","tertile","quartile",NA)[1],
                       downsampling=NULL,
                       limits=c(-0.2,4.3),
                       channel_names=c("channel","marker","both")[3]){
  if(all(cluster_to_use!="all")){
    processed_fcs <- phenograph_obj$processed_fcs[phenograph_obj$processed_fcs$Phenograph_membership %in%  cluster_to_use ,]
  }else{
    processed_fcs <- phenograph_obj$processed_fcs
    cluster_to_use <- as.character(unique(processed_fcs$Phenograph_membership))
    cluster_to_use <- cluster_to_use[order(nchar(cluster_to_use), cluster_to_use)]
  }
  all_channels <- attr(processed_fcs,"all_channels")
  if(all(channels=="all")){
    channels <- intersect(colnames(processed_fcs),all_channels[,"name"])
  }
  if(!is.null(downsampling)){
    processed_fcs <- downsample(processed_fcs,downsampling)
  }
  if(channel_names == "both" | channel_names =="marker"){
    if(channel_names == "both"){
      edited_channels_name <- paste( channels,all_channels[match(channels,all_channels$name),"desc"],sep=" / ")
      edited_channels_name <- gsub(x = edited_channels_name,pattern = " / NA",replacement = "")
    }else if(channel_names =="marker"){
      edited_channels_name <- all_channels[match(channels,all_channels$name),"desc"]
      na_str <- is.na(edited_channels_name)
      edited_channels_name[na_str] <- channels[na_str]
    }
    order <- match(channels,colnames(processed_fcs))
    colnames(processed_fcs)[order] <- edited_channels_name
    channels <- edited_channels_name
  }
  #compute threshold values for each channels
  threshold_values <- apply(processed_fcs[,channels],2,function(x){
    if(threshold=="median"){
      return(median(x,na.rm=T))
    }
    if(threshold=="tertile"){
      return(quantile(x,probs = c(0.33,0.66),na.rm = T))
    }
    if(threshold=="quartile"){
      return(quantile(x,probs = c(0.25,0.5,0.75),na.rm = T))
    }
  })
  #if(threshold == "median") threshold_values <- t(threshold_values)
  #plot clusters for each channels
  #nasty nested ifelses but works...
  if(type == "channels"){
    all_plots <- lapply(channels,function(chan){
      melted_df<- melt_df(processed_fcs[,c(chan,"Phenograph_membership")],var_to_group ="Phenograph_membership")
      if(!is.na(threshold)){
        if(threshold == "median"){
          p <- ggplot(melted_df, aes(x = value, y = groups, fill = ifelse(..x..> threshold_values[chan], "higher than median", "lower than median")))
        }
        if(threshold == "tertile"){
          p <- ggplot(melted_df, aes(x = value, y = groups, fill = factor(ifelse(..x.. > threshold_values[2,chan],"third tertile",
                                                                            ifelse(..x.. > threshold_values[1,chan],"second tertile","first tertile")),
                                                                          levels=c("third tertile","second tertile","first tertile"))))
        }
        if(threshold == "quartile"){
          p <- ggplot(melted_df, aes(x = value, y = groups, fill = factor(ifelse(..x.. > threshold_values[3,chan],"fourth quartile",
                                                                                 ifelse(..x.. > threshold_values[2,chan],"third quartile",
                                                                                        ifelse(..x.. > threshold_values[1,chan],"second quartile","first quartile"))),
                                                                          levels=c("fourth quartile","third quartile","second quartile","first quartile"))))          }
      }
      p <- p + stat_density_ridges(geom = "density_ridges_gradient", scale = 3, rel_min_height = 0.045)
      p <- p + theme_ridges()
      p <- p + scale_fill_viridis_d(name=paste0(threshold," threshold"),direction = -1,option = "C",alpha = 0.6)
      p <- p + theme(axis.title.y = element_blank(),axis.title.x = element_blank())
      p <- p + labs(title=chan)
      p <- p + scale_x_continuous(limits = limits)
      #p <- p + theme(legend.position = "none")
      return(p)
    })
  }
  #plot channels for each cluster
  if(type == "clusters"){
    all_plots <- lapply(cluster_to_use,function(x){
      temp_df <- data.frame(t(processed_fcs[processed_fcs$Phenograph_membership == x,channels]))
      temp_df$chan <- rownames(temp_df)
      melted_df<- melt_df(temp_df,var_to_group ="chan")
      melted_df$groups <- factor(melted_df$groups,levels=channels)
      if(!is.na(threshold)){
        if(threshold == "median"){
          p <- ggplot(melted_df, aes(x = value, y = groups, fill = ifelse(..x.. > threshold_values[..y..],"higher than median","lower than median")))
        }
        if(threshold == "tertile"){
          p <- ggplot(melted_df, aes(x = value, y = groups, fill = factor(ifelse(..x.. > threshold_values[2,..y..],"third tertile",
                                                                          ifelse(..x.. > threshold_values[1,..y..],"second tertile","first tertile")),
                                                                          levels=c("third tertile","second tertile","first tertile"))))
        }
        if(threshold == "quartile"){
          p <- ggplot(melted_df, aes(x = value, y = groups, fill = factor(ifelse(..x.. > threshold_values[3,..y..],"fourth quartile",
                                                                          ifelse(..x.. > threshold_values[2,..y..],"third quartile",
                                                                                 ifelse(..x.. > threshold_values[1,..y..],"second quartile","first quartile"))),
                                                                          levels=c("fourth quartile","third quartile","second quartile","first quartile"))))
        }
      }
      p <- p + stat_density_ridges(geom = "density_ridges_gradient", scale = 3, rel_min_height = 0.045,alpha = 0.60)
      p <- p + theme_ridges()
      p <- p + scale_fill_viridis_d(name=paste0(threshold," threshold"),direction = -1,option = "C",alpha = 0.6)
      p <- p + theme(axis.title.y = element_blank(),axis.title.x = element_blank())
      p <- p + labs(title=x)
      p <- p + scale_x_continuous(limits = limits)
     # p <- p + theme(legend.position = "none")
      return(p)
    })
  }
  return(all_plots)
}

#' Plot function that display a Umap of scaled intensities for event
#' @param umap_2D matrix detailing coordinates of each event on the umap
#' @param processed_fcs_df dataframe of processed intensities for each event and informations of channel used as attribute
#' @param channels vector of channels to select. Can be "all" to select all channels, "with_desc" to select channels with a marker description or a vector a channels
#' @param cut_top_99th boolean value indicating if the top 1 percent should be removed to compute color scaling
#' @param alpha numeric indicating the transparency level when plotting the umap
#' @import ggplot2
#' @importFrom stats quantile
#' @importFrom RColorBrewer brewer.pal
#' @export
plot_umap <- function(umap_2D,
                      processed_fcs_df,
                      channels=c("channels_used","all","with_desc")[1],
                      cut_top_99th=T,
                      title=c("both","marker","channel")[1],
                      alpha=0.5){

  processed_fcs<- query_extract(processed_fcs_df,channels=channels)
  all_channels <- attr(processed_fcs,"all_channels")
  channels_to_use <-  setdiff(colnames(processed_fcs),"sample_id")

  #subset processed_fcs based on umap corrdinates (safer when downsampling_umap != NULL)
  processed_fcs <- processed_fcs[rownames(umap_2D),]

  all_plot <- lapply(channels_to_use, function(channel){
    if(title=="marker"){
      title_to_use <- all_channels[which(all_channels[,1]==channel),2]
      if(is.na(title_to_use)){
        title_to_use <-channel
      }
    }
    if(title=="both"){
      title_to_use <- paste(channel,all_channels[which(all_channels[,1]==channel),2],sep=" / ")
      title_to_use <- gsub(pattern = "/ NA",replacement = "",x = title_to_use)
    }
    if(title=='channel'){
      title_to_use <-channel
    }
    p <- ggplot(umap_2D,aes(x=X,y=Y))
    p <- p + theme_classic()
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
    p <- p + ggtitle(title_to_use)
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
    percs <- round(mean_perc,digits = 3)*100
    percs[percs<0.1] <- "< 0.1"
    colnames(mean_intensity_pheno)<- paste0(clusters_id,"  ",percs,"%")
  }else{
    colnames(mean_intensity_pheno)<- clusters_id
  }
  if(cut_top_99th){
    normalized_mean_intensity_pheno <- apply(mean_intensity_pheno,2,function(x) norm_range(x,c(range=quantile(x,probs = 0.01),quantile(x,probs = 0.99))))
  }else{
    normalized_mean_intensity_pheno <- apply(mean_intensity_pheno,2,function(x) norm_range(x,range=c(0,1)))
  }
  #change channels order
  if(!is.null(channels_order)){
    normalized_mean_intensity_pheno <- normalized_mean_intensity_pheno[channels_order,]
  }else{
    channels_order <- rownames(normalized_mean_intensity_pheno)
  }
  if(display_marker_names==T){
    new_rownames <- all_channels[match(channels_order,all_channels$name),"desc"]
    new_rownames[is.na(new_rownames)] <- channels_order[is.na(new_rownames)]
    rownames(normalized_mean_intensity_pheno) <- new_rownames
  }
  if(palette=="default"){
    palette <- colorRampPalette(brewer.pal(name = "PuBu",n=9))(100)
  }
  #change clusters order
  if(!is.null(memberships_order)){
    normalized_mean_intensity_pheno <- normalized_mean_intensity_pheno[,memberships_order]
  }
  if(!is.null(channels_order)&!is.null(memberships_order)){
    pheatmap(t(normalized_mean_intensity_pheno),
             color=palette,
             angle_col = 0,
             cluster_cols = F,
             silent = T,
             cluster_rows = F,
             ...)
  }else if(!is.null(channels_order)){
    pheatmap(t(normalized_mean_intensity_pheno),
             color=palette,
             cluster_cols = F,
             silent = T,
             angle_col = 0,
             ...)
  }else if(!is.null(memberships_order)){
    pheatmap(t(normalized_mean_intensity_pheno),
             color=palette,
             cluster_rows = F,
             silent = T,
             angle_col = 0,
             ...)
  }else{
    p <- pheatmap(t(normalized_mean_intensity_pheno),
             color=palette,
             silent = T,
             angle_col = 45,
             ...)
  }

}

