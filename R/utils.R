#' @importFrom reshape2 melt
norm_range <- function(x, range = c(0, 4.5)) {
  (x - min(x))/(max(x) - min(x)) * (range[2] - range[1]) + range[1]
}
melt_df <- function(df,var_to_group){
  df.m <- melt(df,id = var_to_group,varnames = c("variable","value"))
  colnames(df.m)[1] <- "groups"
  df.m$value <- as.numeric(as.character(df.m$value))
  return(df.m)
}
downsample <- function(df,number_to_sample){
  selected_events <- sapply(unique(df$sample_id),function(x){
    if(sum(df$sample_id == x)  >= number_to_sample){
      sample(x = rownames(df[df$sample_id == x,]),size = number_to_sample)
    }else{
      rownames(df[df$sample_id == x,])
    }
    })
  df <- df[unlist(selected_events),]
  return(df)
}


