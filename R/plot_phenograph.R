#' @import ggplot2
#' @importFrom stats quantile
#' @importFrom ggrepel geom_text_repel

plot_phenograph <- function(umap_2D,phenograph_obj,alpha=0.5){
  memberships <- factor(phenograph_obj$phenograph[[2]]$membership)
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
