#' @importFrom flowCore read.flowSet
#' @importFrom flowCore pData
#' @importFrom flowCore fsApply
#' @importFrom flowCore estimateLogicle
#' @importFrom flowCore transform
#' @export
#'

#Open fcs and put then in a flowset
pre_process_fcs <- function(fcs_dir,downsampling="none",rescale_all=c(0,4.5)){
  fs <- read.flowSet(path = fcs_dir,transformation = F,emptyValue = F)
  #get markers
  all_channels <- pData(parameters(fs[[1]]))[,c(1,2)]
  all_channels <- all_channels[as.vector(all_channels[,1] != "Time" & all_channels[,1] != "Event"),]
  shape_marker <- grep('FSC|SSC',all_channels$name,value = T)
  #transform values
  event_for_each_sample <- fsApply(fs,function(ff){
    lgcl <- estimateLogicle(ff, channels = setdiff(all_channels[,1],shape_marker),type="data")
    ff <- transform(ff,lgcl)
    mat <- data.frame(exprs(ff),check.names = F)
    #linear scale for scatter values
    if(!is.null(shape_marker)){
      mat[,shape_marker]<- sapply(shape_marker,function(x) norm_range(mat[,x,drop=F]))
    }
    mat[,"sample_id"] <-identifier(ff)
    return(mat)
  })
  if(downsampling !="none"){
    event_for_each_sample <- lapply(event_for_each_sample,function(x){
      if(nrow(x)>1000 & nrow(x) > downsampling){
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
