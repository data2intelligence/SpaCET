#' @title Metrics
#' @description Calculate the ST metrics
#' @param ST An SpaCE object
#' @return An SpaCE object
#' @details This function omputes both UMI and gene counts across ST spots
#' @examples 
#' ST <- ST.metrics(ST)
#' @rdname ST.metrics
#' @export 
ST.metrics <- function(ST)
{
  visiualVector <- colSums(ST@input$counts)
  nUMI <- colSums(ST@input$counts)
  nGene <- colSums(ST@input$counts>0)
  
  ST@results$metrics <- rbind(nUMI,nGene)

  ST
}

#' @title ST metrics plot
#' @description Plot the quality control metrics for all ST spots
#' @param ST An SpaCE object
#' @param cols Colors, Default: c("lightblue", "blue", "darkblue")
#' @return A ggplot2 object
#' @details NULL
#' @examples ST.metrics.plot(ST, cols = c("lightblue", "blue", "darkblue"))
#' @rdname ST.metrics.plot
#' @export 
#' @importFrom cowplot plot_grid
ST.metrics.plot <- function( ST, cols = c("lightblue", "blue", "darkblue") )
{
  p1 <- visualSpatial(ST@results$metrics["nUMI",],ST@input$HEimage,cols,NULL,"UMI count","nUMI")
  p2 <- visualSpatial(ST@results$metrics["nGene",],ST@input$HEimage,cols,NULL,"Gene count","nGene")
  
  cowplot::plot_grid(p1, p2, nrow=1)
}
