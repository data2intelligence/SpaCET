#' @title ST Visualization
#' @description FUNCTION_DESCRIPTION
#' @param visiualVector PARAM_DESCRIPTION
#' @param HEimage PARAM_DESCRIPTION
#' @param cols PARAM_DESCRIPTION
#' @param colsLimits PARAM_DESCRIPTION
#' @param titleName PARAM_DESCRIPTION
#' @param legendName PARAM_DESCRIPTION
#' @param legendPosition PARAM_DESCRIPTION
#' @param is.continuous PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[png]{readPNG}}
#'  \code{\link[grid]{grid.raster}},\code{\link[grid]{unit}}
#' @rdname visualSpatial
#' @export 
#' @importFrom png readPNG
#' @importFrom grid rasterGrob unit
visualSpatial <- function(visiualVector,HEimage,cols,colsLimits,titleName,legendName,legendPosition="right",is.continuous=TRUE)
{
  r <- png::readPNG(HEimage)
  rg <- grid::rasterGrob(r, width=grid::unit(1,"npc"), height=grid::unit(1,"npc"))
  nDiml <- dim(r)[1] # dim pixal
  
  coordi <- t(matrix(as.numeric(unlist(strsplit(names(visiualVector),"x"))),nrow=2))
  fig.df <- data.frame(x=nDiml-coordi[,1],y=coordi[,2])
  fig.df[["value"]] <- visiualVector
  rownames(fig.df) <- names(visiualVector)
  
  if(is.continuous==TRUE)
  {
    ggplot(fig.df,aes(x=x,y=y))+ 
      annotation_custom(rg) + # add background image
      geom_point(aes(colour=value))+
      scale_colour_gradientn(name=legendName, colours = cols, limits=colsLimits)+
      scale_x_continuous(limits = c(0, nDiml), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, nDiml), expand = c(0, 0)) +
      ggtitle(titleName)+
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = legendPosition
      )+coord_flip()
  }else{
    ggplot(fig.df,aes(x=x,y=y))+ 
      annotation_custom(rg) + # add background image
      geom_point(aes(colour=value))+
      scale_colour_manual(name=legendName,values=cols)+
      scale_x_continuous(limits = c(0, nDiml), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, nDiml), expand = c(0, 0)) +
      ggtitle(titleName)+
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = legendPosition
      )+coord_flip()
    
  }
  
}
