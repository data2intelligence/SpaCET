#' @title Deconvolution results violin plot
#' @description Visualize the deconvolution results in violin plot
#' @param ST An SpaCE object
#' @return A ggplot2 object
#' @details This function plot the deconvolution results across all ST spots
#' @examples 
#' 
#' @rdname ST.deconv.res.violin
#' @export 
#' @importFrom reshape2 melt
ST.deconv.res.violin <- function(ST)
{
  mypal <- c(
    "#f2f3f4",
    "#222222",
    "#f3c300",
    "#875692",
    "#f38400",
    "#a1caf1",
    "#be0032",
    "#c2b280",
    "#848482",
    "#008856",
    "#e68fac",
    "#0067a5",
    "#f99379",
    "#604e97",
    "#f6a600",
    "#b3446c",
    "#dcd300",
    "#882d17",
    "#8db600",
    "#654522",
    "#e25822",
    "#2b3d26"
  )
  
  names(mypal) <- c(
    "white",
    "black",
    "yellow",
    "purple",
    "orange",
    "lightblue",
    "red",
    "buff",
    "gray",
    "green",
    "purplishpink",
    "blue",
    "yellowishpink",
    "violet",
    "orangeyellow",
    "purplishred",
    "greenishyellow",
    "reddishbrown",
    "yellowgreen",
    "yellowishbrown",
    "reddishorange",
    "olivegreen"
  )
  
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
  deconv.res.m <- reshape2::melt(deconv.res)
  deconv.res.m[,1] <- factor(deconv.res.m[,1],levels=row.names(deconv.res))
  colnames(deconv.res.m)[1] <- "cellType"
  
  ggplot(deconv.res.m)+
    geom_violin(aes(x=cellType,y=value,fill=cellType), scale = "width")+
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

#' @title Deconvolution results scatter plot
#' @description Visualize the deconvolution results in scatter plot
#' @param ST An SpaCE object
#' @param cellTypes A vector of cell-types
#' @param cols Colors, Default: c("blue", "yellow", "red")
#' @param colsLimits Colors Limits of cell type fraction, Default: c(0, 1)
#' @return A ggplot2 object
#' @details DETAILS
#' @examples 
#' 
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

