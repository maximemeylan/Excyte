#' @import ggplot2
#' @import ggridges

plot_ridge <- function(pheno_obj){
  processed_fcs <- pheno_obj$processed_fcs
  all_channels <- attr(processed_fcs,"all_channels")
  channels_to_plot <- intersect(colnames(processed_fcs),all_channels[,"name"])
  all_plots <- lapply(channels_to_plot,function(x){
    melted_df<- melt_df(processed_fcs[,c(x,"Phenograph_membership")],var_to_group ="Phenograph_membership")
    p <- ggplot(melted_df, aes(x = value, y = groups,fill = groups))
    p <- p + geom_density_ridges(scale = 4, rel_min_height = 0.045,alpha = 0.85)
    #p <- p + scale_y_discrete(expand = expand_scale(mult=c(0.01, 0.7)))
    p <- p + theme_ridges()
    p <- p + theme(axis.title.y = element_blank(),axis.title.x = element_blank())
    p <- p + labs(title=paste(x,all_channels[which(all_channels[,1]==x),2],sep=" / "))
    p <- p + scale_fill_discrete(guide=FALSE)
    return(p)
  })
  return(all_plots)
}
