#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ST PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ST PARAM_DESCRIPTION
#' @param cols PARAM_DESCRIPTION, Default: c("lightblue", "blue", "darkblue")
#' @param colsLimits PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[cowplot]{plot_grid}}
#' @rdname ST.metrics.plot
#' @export 
#' @importFrom cowplot plot_grid
ST.metrics.plot <- function( ST, cols = c("lightblue", "blue", "darkblue") )
{
  p1 <- visualSpatial(ST@results$metrics["nUMI",],ST@input$HEimage,cols,NULL,"UMI count","nUMI")
  p2 <- visualSpatial(ST@results$metrics["nGene",],ST@input$HEimage,cols,NULL,"Gene count","nGene")
  
  cowplot::plot_grid(p1, p2, nrow=1)
}
