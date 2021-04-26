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
#' @seealso 
#'  \code{\link[reshape]{melt}}
#' @rdname ST.deconv.res.violin
#' @export 
#' @importFrom reshape melt
ST.deconv.res.violin <- function(ST)
{
  library(Polychrome)
  mypal <- kelly.colors(22)
  
  colorVector <- c(
    mypal["purplishpink"],"B cell",
    mypal["purple"],"CAF",
    mypal["orange"],"Dendritic",
    mypal["olivegreen"],"Endothelial",
    mypal["greenishyellow"],"M1",
    mypal["reddishbrown"],"M2",
    mypal["yellowishbrown"],"Macrophage other",
    mypal["yellow"],"Malignant",
    mypal["orangeyellow"],"Mast",
    mypal["gray"],"Neutrophil",
    mypal["lightblue"],"NK",
    mypal["blue"],"T CD4 naive",
    mypal["buff"],"T CD4 other",
    mypal["yellowishpink"],"T CD8",
    mypal["green"],"T helper",
    mypal["violet"],"Treg",
    mypal["red"],"Macrophage",
    mypal["yellowgreen"],"T CD4",
    mypal["reddishorange"],"Unidentifiable"
  )
  
  colorMatrix <- matrix(colorVector,nrow=2)
  colnames(colorMatrix) <- colorMatrix[2,]
  myColors  <- colorMatrix[1,] # name : cell type
  
  deconv.res <- ST@results$fraction
  deconv.res.m <- reshape::melt(deconv.res)
  deconv.res.m[,1] <- factor(deconv.res.m[,1],levels=row.names(deconv.res))
  
  library(ggplot2)
  ggplot(deconv.res.m)+
    geom_violin(aes(x=X1,y=value,fill=X1), scale = "width")+
    scale_y_continuous(limits=c(0,1),breaks = seq(0, 1, by =0.2 ))+
    scale_fill_manual(values = myColors)+
    ylab("Cell fraction")+
    xlab("Cell lineage and sub-lineage")+
    theme_bw()+ 
    theme(
      panel.grid.minor = element_blank(),
      axis.text = element_text(size=12,colour = "black"),
      axis.title = element_text(size=14,colour = "black"),
      axis.text.x = element_text(angle = 30,hjust = 1),
      legend.position = "none"
    )
  
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ST PARAM_DESCRIPTION
#' @param cellTypes PARAM_DESCRIPTION
#' @param cols PARAM_DESCRIPTION, Default: c("blue", "yellow", "red")
#' @param colsLimits PARAM_DESCRIPTION, Default: c(0, 1)
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
#' @rdname ST.deconv.res.scatter
#' @export 
#' @importFrom cowplot plot_grid
ST.deconv.res.scatter <- function(
  ST,
  cellTypes,
  cols = c("blue", "yellow", "red"),
  colsLimits=c(0,1)
  )
{
  deconv.res <- ST@results$fraction
 
  for(cellType in cellTypes)
  {
    p <- visualSpatial(deconv.res[cellType,],ST@input$HEimage,cols,colsLimits,cellType,"Fraction")
    
    if(cellType!=cellTypes[1])
    {
      pp <- cowplot::plot_grid(pp, p, nrow=1, rel_widths=c(which(cellTypes==cellType)-1, 1))
    }else{
      pp <- p
    }
  }
  
  pp
}

