#' @import ggplot2
#' @importFrom stats quantile
#' @importFrom RColorBrewer brewer.pal
#' @export
plot_umap <- function(umap_2D,processed_fcs_df,channels=c("channels_used","all","with_desc")[1],cut_top_99th=T,alpha=0.5){
  processed_fcs<- query_extract(processed_fcs_df,channels=channels)
  all_channels <- attr(processed_fcs,"all_channels")
  channels_to_use <-  setdiff(colnames(processed_fcs),"sample_id")

  all_plot <- lapply(channels_to_use, function(channel){
    title <- paste(channel,all_channels[which(all_channels[,1]==channel),2],sep=" / ")
    p <- ggplot(umap_2D,aes(x=X,y=Y))
    p <- p + theme_bw()
    if(cut_top_99th){
      p <- p + geom_point(size=0.05,
                          alpha=alpha,
                          aes(colour = norm_range(processed_fcs[,channel],c(quantile(processed_fcs[,channel],probs = 0.01),quantile(processed_fcs[,channel],probs = 0.99)))))
      #aes(colour = norm_range(processed_fcs[,channel],c(min(processed_fcs[,channel]),quantile(processed_fcs[,channel],probs = 0.99)))))
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
