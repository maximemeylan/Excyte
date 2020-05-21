#' Master function to run preprocessing,phenograph clustering and umap computation
#' @param fcs_dir directory or vector containing fcs files to be used
#' @param downsampling numeric indicating the number event to randomly select from each fcs, if the number of events request is bigger than the number of event in the  fcs, all event are selected
#' @param channels vector containing channels to select. Can be "all" to select all channels, "with_desc" to select channels with a marker description or a vector a channels.
#' @param k numeric indicating the number of neighbor for phenograph and umap computation
#' @param downsampling_umap numeric indicating the number of events to sample to compute the umap
#' @param method define if the umap should be computed with the python package (umap-learn) or with the naive R implementation
#' @export
run_excyte <- function(fcs_dir,
                       downsampling=3000,
                       downsampling_umap=NULL,
                       channels="all",
                       method=c("umap-learn","naive")[1],
                       k=30){
  if(any(!is.null(downsampling_umap) & downsampling_umap > downsampling)){
    stop("the number of events to sample for umap should smaller than the number of event sampled for the pipeline")
  }
  #pre-process fcs
  processed_fcs_obj <- pre_process_fcs(fcs_dir = fcs_dir,downsampling = downsampling)
  #compute phenograph membership for each event
  phenograph_obj <- compute_phenograph(processed_fcs_obj,channels = channels,k = k)
  #compute umap coordinates for each event
  umap_obj <- compute_umap(processed_fcs_obj,channels = channels,k=k,downsampling_umap=downsampling_umap,method=method)
  #provide annoatation for each phenograph membership
  annotation <- annotate_clusters(phenograph_obj = phenograph_obj,channels=channels)
  return(list("processed_fcs_obj"=processed_fcs_obj,"phenograph_obj"=phenograph_obj,"umap_obj"=umap_obj,"annotation" =annotation))
}


#' Rerun the excyte pipeline on selected phenograph clusters
#' @param excyte_obj list of object obtained from an initial run with the excyte pipeline
#' @param channels vector containing channels to select. Can be "all" to select all channels, "with_desc" to select channels with a marker description or a vector a channels.
#' @param k numeric indicating the number of neighbor for phenograph and umap computation
#' @param downsampling numeric indicating the number event to randomly select from each fcs, if the number of events request is bigger than the number of event in the  fcs, all event are selected
#' @param clusters_id vector of character containing the ID of the phenograph clusters to rerun the excyte pipeline on
#' @param downsampling_umap numeric indicating the number of events to sample to compute the umap
#' @param method define if the umap should be computed with the python package (umap-learn) or with the naive R implementation
#' @export
#'
rerun_excyte <- function(excyte_obj,
                         clusters_id=NA,
                         downsampling=3000,
                         downsampling_umap=NULL,
                         channels="all",
                         method=c("umap-learn","naive")[1],
                         k=30){
  if(any(is.na(clusters_id))){
    stop("Please submit valid clusters ID")
  }
  if(any(!is.null(downsampling_umap) & downsampling_umap > downsampling)){
    stop("the number of events to sample for umap should smaller than the number of event sampled for the pipeline")
  }
  message("Excyte re-running with selected clusters: ",paste0(clusters_id," "))
  event_to_select <- excyte_obj$phenograph_obj$processed_fcs$Phenograph_membership %in% clusters_id
  excyte_obj$processed_fcs_obj$processed_fcs <- excyte_obj$processed_fcs_obj$processed_fcs[event_to_select,]
  #downsample if requested
  if(!is.null(downsampling)){
    excyte_obj$processed_fcs_obj$processed_fcs <- downsample(excyte_obj$processed_fcs_obj$processed_fcs,downsampling)
  }

  #compute new phenograph membership for selected events
  phenograph_obj <- compute_phenograph(processed_fcs_obj = excyte_obj$processed_fcs_obj,channels = channels,k = k)
  #compute umap for selected events
  umap_obj <- compute_umap(excyte_obj$processed_fcs_obj,channels = channels,k=k,downsampling_umap =downsampling_umap,method=method)
  #provide annoatation for each phenograph membership
  annotation <- annotate_clusters(phenograph_obj = phenograph_obj,channels=channels)
  return(list("processed_fcs_obj"= excyte_obj$processed_fcs_obj,"phenograph_obj"=phenograph_obj,"umap_obj"=umap_obj,"annotation" =annotation))
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
  spacer  <- ifelse(nchar(phenograph_obj[[2]]$membership) <2,"C_0","C_")
  processed_fcs$Phenograph_membership <- factor(paste0(spacer,phenograph_obj[[2]]$membership))

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
#' @param downsampling_umap numeric indicating the number of events to select to compute the umap
#' @param method define if the umap should be computed with the python package (umap-learn) or with the naive R implementation
#' @import umap
#' @importFrom stats setNames
#' @export

compute_umap <- function(processed_fcs_obj,channels=c("all","with_desc")[1],k=30,downsampling_umap=NULL,method=c("umap-learn","naive")[1]){
  processed_fcs<- query_extract(processed_fcs_obj,channels=channels)
  #downsample if requested
  if(!is.null(downsampling_umap)){
    processed_fcs <- downsample(processed_fcs,downsampling_umap)
  }
  channels_to_use <- setdiff(colnames(processed_fcs),"sample_id")
  message("\nComputing Umap with channels: \t",paste0(channels_to_use,collapse = "\t"))
  #compute umap
  umap_obj <- umap(processed_fcs[,channels_to_use],method = method,k=k)
  umap_obj_2D <- setNames(data.frame(umap_obj$layout,check.names = F),c("X","Y"))
  return(list("umap_obj"=umap_obj,"umap_2D"=umap_obj_2D,"channels_used"=channels_to_use))
}

#' Provide annotations for each cluster, based on intensity distribution of all events
#' @param phenograph_obj list containing result of phenograph clustering and processed fcs
#' @param thresholds character defining if threshold should be caracterized as the median, tertiles or quartiles
#' @param positivity_threshold numeric value between 0 and 1 defining the percentage of cells needed to call positivity to a threshold
#' @param channels vector of channels to use, default uses all channels
#' @param cluster_to_use vector of cluster to use, default uses all clusters
#' @param channel_names character, edit channels names accordingly
#' @export
annotate_clusters <- function(phenograph_obj,
                              channels="all",
                              cluster_to_use="all",
                              thresholds=c("median","tertile","quartile")[2],
                              positivity_threshold = 0.5,
                              channel_names=c("channel_only","marker_only","both")[3]
                              ){
  ## clean channels names and prepare data
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
  if(channel_names == "both" | channel_names =="marker_only"){
    if(channel_names == "both"){
      edited_channels_name <- paste( channels,all_channels[match(channels,all_channels$name),"desc"],sep=" / ")
      edited_channels_name <- gsub(x = edited_channels_name,pattern = " / NA",replacement = "")
    }else if(channel_names =="marker_only"){
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
    if(thresholds=="median"){
      return(median(x,na.rm=T))
    }
    if(thresholds=="tertile"){
      return(quantile(x,probs = c(0.33,0.66),na.rm = T))
    }
    if(thresholds=="quartile"){
      return(quantile(x,probs = c(0.25,0.5,0.75),na.rm = T))
    }
  })
  annot <- sapply(as.character(channels), function(chan){
    positive_cell_per_chan <- sapply(threshold_values[,chan],function(thresh){
      processed_fcs[,chan] > thresh
      })
    positive_cluster_per_chan <- sapply(cluster_to_use,function(cluster){
      ids <- processed_fcs$Phenograph_membership==cluster
      #compute percentage of positive cells for each thresold
      perc_sup_cell_per_chan <- colSums(positive_cell_per_chan[ids,]) / sum(ids)
      positivity_cell_per_chan <- as.logical(perc_sup_cell_per_chan > positivity_threshold)
      if(thresholds == "median"){
        if(identical(positivity_cell_per_chan,F)) return("low")
        if(identical(positivity_cell_per_chan,T)) return("high")
      }
      if(thresholds == "tertile"){
        if(identical(positivity_cell_per_chan,c(F,F))) return("low")
        if(identical(positivity_cell_per_chan, c(T,F))) return("medium")
        if(identical(positivity_cell_per_chan,c(T,T))) return("high")
      }
      if(thresholds == "quartile"){
        if(identical(positivity_cell_per_chan,c(F,F,F))) return("very_low")
        if(identical(positivity_cell_per_chan,c(T,F,F))) return("low")
        if(identical(positivity_cell_per_chan,c(T,T,F))) return("high")
        if(identical(positivity_cell_per_chan,c(T,T,T))) return("very_high")
      }
    })
  })
  annot <- data.frame(annot,check.names = F,check.rows = F)
  #aggregated vector to facilitate sorting
  annot$phenotype <- apply(annot,1,function(x) paste0(rbind(colnames(annot),x),collapse = "_"))
  return(annot)
}
