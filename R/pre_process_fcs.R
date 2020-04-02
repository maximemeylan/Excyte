#' Extract and aggreagte fluorescence intensity matrixx of FCS files.
#' Perfom logicle transformation and scaling, downsampling if necessary
#'
#' @importFrom flowCore read.flowSet
#' @importFrom flowCore pData
#' @importFrom flowCore fsApply
#' @importFrom flowCore estimateLogicle
#' @importFrom flowCore transform
#' @importFrom flowCore parameters
#' @importFrom flowCore exprs
#' @importFrom flowCore identifier
#' @importFrom flowCore logicleTransform
#' @importFrom flowCore transformList
#' @export
#' @param fcs_dir directory or vector containing fcs files to be used
#' @param downsampling number of event to randomly select from each fcs, if the number of events request is bigger than the number of event in the  fcs, all event are selected
#' @param rescale_all vector of two values indicating the range of the values to scale between
#'
#' @return a list containing the normalized aggregated dataframe and all_channels
#'

#Open fcs and put then in a flowset
pre_process_fcs <- function(fcs_dir,downsampling=NULL,rescale_all=c(0,4.5)){
  fs <- read.flowSet(fcs_dir,transformation = F,emptyValue = F)
  #get markers
  all_channels <- pData(parameters(fs[[1]]))[,c(1,2)]
  all_channels <- all_channels[as.vector(all_channels[,1] != "T- e" & all_channels[,1] != "Event"),]
  shape_marker <- grep('FSC|SSC',all_channels$name,value = T)
  channels_to_normalize <- setdiff(all_channels[,1],shape_marker)
  #transform values
  event_for_each_sample <- fsApply(fs,function(ff){
    param_T_list <- lapply(channels_to_normalize,function(chan){
      lgcl <- estimateLogicle(ff, channels = chan,type="data")
      param_T <- sapply(c("t","m", "a", "w"), function(param) as.numeric(format(as.vector(environment(lgcl@transforms[[chan]]@f)[[param]]), digits = 2)))
      if(param_T["w"] > 2 | param_T["w"] < 0 ){
        warning(paste0("estimateLogicle failed for channel: ", chan, " defaulting to standard parameters"))
        param_T["t"] <- 4000
        param_T["m"] <- 4.5
        param_T["a"] <- 0
        param_T["w"] <- 0.1
      }
      logicleTransform(transformationId = chan,
                       w = param_T["w"], t = param_T["t"], m = param_T["m"], a = param_T["a"] )
    })
    lgcl <- transformList(channels_to_normalize, param_T_list)
    ff <- transform(ff,lgcl)
    mat <- data.frame(exprs(ff),check.names = F)
   #linear scale for scatter values
   if(!is.null(shape_marker)){
    mat[,shape_marker]<- sapply(shape_marker,function(x) norm_range(mat[,x,drop=F]))
   }
   mat[,"sample_id"] <-identifier(ff)
   return(mat)
  })
  if(!is.null(downsampling)){
    event_for_each_sample <- lapply(event_for_each_sample,function(x){
      if(nrow(x) > downsampling){
        x[sample(nrow(x),downsampling),]
      }else{
        return(x)
      }
    })
  }
  processed_fcs_df <- data.frame(do.call("rbind", event_for_each_sample),check.names = F)
  if(!is.null(rescale_all)){
    processed_fcs_df[,all_channels[,1]]<- apply(processed_fcs_df[,all_channels[,1]],2,function(x) norm_range(x,rescale_all))
  }
  return(list("processed_fcs" =processed_fcs_df,"all_channels"=all_channels))
}
#' Query an aggregated dataframe to subset requested channels
#' @param processed_fcs_obj list containing a datraframe of processed intensities for each event and informations of channel used
#' @param channels vector containing channels to select. Can be "all" to select all channels, "with_desc" to select channels with a marker description or a vector a channels.
#' @export
query_extract <- function(processed_fcs_obj,channels=c("all","with_desc")[1]){
  all_channels <- processed_fcs_obj$all_channels
  processed_fcs_df <- processed_fcs_obj$processed_fcs
  #if channels are not specified
  if (all(channels == "with_desc")){
    channels_to_select <- all_channels[!is.na(all_channels[, "desc"]), 1]
    processed_fcs_df <- processed_fcs_df[, c(channels_to_select, "sample_id")]
    channels <- colnames(processed_fcs_df)[colnames(processed_fcs_df) != "sample_id"]
  }
  if (all(channels != "with_desc") & all(channels != "all")){
    processed_fcs_df <- processed_fcs_df[, c(channels, "sample_id")]
    channels <- colnames(processed_fcs_df)[colnames(processed_fcs_df) != "sample_id"]
  } else{
    channels <- all_channels[,1]
  }
  processed_fcs_df <- processed_fcs_df[,c(channels,"sample_id")]
  attr(processed_fcs_df,"all_channels") <- all_channels
  return(processed_fcs_df)
}
