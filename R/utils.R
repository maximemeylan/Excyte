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
