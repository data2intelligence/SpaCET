#' @title Spatial feature visualization
#' @description Visualize multiple types of spatial features in ST data.
#' @param SpaCET_obj An SpaCET object.
#' @param spatialType Type of spatial features, i.e., "QualityControl", "GeneExpression", "CellFraction", and "LRNetworkScore". See ‘details’ for more information.
#' @param spatialFeatures A vector of spatial features.
#' @param sameScaleForFraction Indicate whether all cell types have the same scale for cell fraction.
#' @param nrow Row number of the combined panel for multiple spatial features.
#' @param imageBg logical: should the image be shown?
#' @return A ggplot2 object.
#' @details
#' `SpaCET.visualize.spatialFeature` is able to plot multiple types of spatial features, including "QC", "GeneExp", "CellFraction", and "LRNetworkScore".
#' "QualityControl" refers to quality control metrics. User can visualize the total counts of UMI and genes across all spots.
#' "GeneExpr" is set to visualize the expression level of genes.
#' "CellFraction" refers to the cell fraction of cell types.
#' "LRNetworkScore" is selected to show the Ligand-Receptor network score. See `SpaCET.CCI.LRNetworkScore` for how to calculate it.
#'
#' @examples
#' SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "QualityControl", spatialFeatures=c("UMI","Gene"))
#' SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "GeneExpression", spatialFeatures=c("EPCAM","MS4A1"))
#' SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "CellFraction", spatialFeatures=c("Malignant","Macrophage"))
#' SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "CellFraction", spatialFeatures="All", pointSize = 0.1, nrow=5)
#' SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "LRNetworkScore", spatialFeatures=c("Network_Score","Network_Score_pv"))
#'
#' @rdname SpaCET.visualize.spatialFeature
#' @export
SpaCET.visualize.spatialFeature <- function(
    SpaCET_obj,
    spatialType = c("QualityControl","GeneExpr","CellFraction","LRNetworkScore"),
    spatialFeatures = NULL,
    sameScaleForFraction = FALSE,
    pointSize = 1,
    nrow = 1,
    imageBg = TRUE
)
{
  if(spatialType == "QualityControl")
  {
    if(is.null(SpaCET_obj@results$metrics))
    {
      stop("Please run SpaCET.quality.control first.")
    }

    mat <- SpaCET_obj@results$metrics

    scaleType="color-continuous"
    colors = c("lightblue", "blue", "darkblue")
    legendName = "Count"
    limits = NULL
  }else if(spatialType == "GeneExpression"){
    mat <- SpaCET_obj@input$counts
    mat <- Matrix::t(Matrix::t(mat)*1e6/Matrix::colSums(mat))
    mat@x <- log2(mat@x+1)

    scaleType="color-continuous"
    colors = c("#4d9221", "yellow", "#c51b7d")
    legendName = "Expr"
    limits = NULL
  }else if(spatialType == "CellFraction"){
    if(is.null(SpaCET_obj@results$deconvolution$propMat))
    {
      stop("Please run cell type deconvolution first.")
    }

    mat <- SpaCET_obj@results$deconvolution$propMat

    if("All"%in%spatialFeatures) spatialFeatures <- rownames(SpaCET_obj@results$deconvolution$propMat)

    scaleType="color-continuous"
    colors = c("blue", "yellow", "red")
    legendName = "Fraction"

    if(sameScaleForFraction)
    {
      limits = c(0,1)
    }else{
      limits = NULL
    }
  }else if(spatialType == "LRNetworkScore"){
    if(is.null(SpaCET_obj@results$CCI$LRNetworkScore))
    {
      stop("Please run SpaCET.CCI.LRNetworkScore first.")
    }

    mat <- SpaCET_obj@results$CCI$LRNetworkScore
    mat[2,mat[2,]> 1.5] <- 1.5
    mat[2,mat[2,]< 0.5] <- 0.5
    mat[3,] <- -log10(mat[3,])

    scaleType="color-continuous"
    colors = c("blue","blue","blue","blue","cyan","cyan","yellow")
    limits = NULL
  }else{
    if(is.null(SpaCET_obj@results$CCI$interface))
    {
      stop("Please run SpaCET.identify.interface first.")
    }

    mat <- SpaCET_obj@results$CCI$interface

    spatialFeatures = c("interface")
    scaleType="color-discrete"
    colors = colors=c("black","darkgrey","#f3c300")
    legendName = "Spot"
    limits = NULL
  }

  for(spatialFeature in spatialFeatures)
  {
    if(spatialType == "LRNetworkScore"){
      if(spatialFeature == "Network_Score"){
        legendName = "Score"
      }else{
        legendName = "-log10pv"
      }
    }

    visiualVector <- mat[spatialFeature,]
    names(visiualVector) <- paste0(SpaCET_obj@input$spotCoordinates[,1],"x",SpaCET_obj@input$spotCoordinates[,2])

    p <- visualSpatial(
      visiualVector,
      image=SpaCET_obj@input$image,
      platform=SpaCET_obj@input$platform,
      scaleType=scaleType,
      colors=colors,
      pointSize=pointSize,
      pointAlpha=1,
      limits=limits,
      titleName=spatialFeature,
      legendName=legendName,
      imageBg=imageBg)

    if(exists("pp"))
    {
      pp <- pp + p
    }else{
      pp <- p
    }
  }

  pp <- pp + patchwork::plot_layout(nrow = nrow)

  if(sameScaleForFraction)
  {
    pp <- pp + patchwork::plot_layout(guides = "collect")
  }

  pp
}


visualSpatial <- function(
    visiualVector,
    image,
    platform,
    scaleType,
    colors,
    limits,
    pointSize,
    pointAlpha,
    titleName,
    legendName,
    imageBg
)
{
  library(ggplot2)

  if(grepl("visium", tolower(platform)))
  {
    coordi <- t(matrix(as.numeric(unlist(strsplit(names(visiualVector),"x"))),nrow=2))

    if(imageBg & !is.na(image$path))
    {
      xDiml <- dim(image$grob$raster)[1] # dim pixel
      yDiml <- dim(image$grob$raster)[2] # dim pixel
    }else{
      xDiml <- max(coordi[,1])
      yDiml <- max(coordi[,2])
    }

    fig.df <- data.frame(
      x=xDiml-coordi[,1],
      y=coordi[,2],
      value=visiualVector
    )
    rownames(fig.df) <- names(visiualVector)

    # 1. initiate plot
    p <- ggplot(fig.df,aes(x=x,y=y))

    # 2. add image
    if(imageBg & !is.na(image$path))
    {
      p <- p+
        annotation_custom(image$grob)
    }

    # 3. draw spot
    p <- p+
      geom_point(aes(colour=value), size=pointSize, alpha=pointAlpha)

    # 4. axis
    if(imageBg & !is.na(image$path))
    {
      p <- p+
        scale_x_continuous(limits = c(0, xDiml), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, yDiml), expand = c(0, 0))
    }

    # 5. set theme
    p <- p+
      ggtitle(titleName)+
      theme(
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()
      )+coord_flip()

    # 6. set color scale
    if(scaleType=="color-continuous")
    {
      p+scale_colour_gradientn(name=legendName, colours = colors, limits=limits)
    }else{
      p+scale_colour_manual(name=legendName,values=colors)
    }

  }else{

    if(imageBg& !is.na(image$path))
    {
      r <- jpeg::readJPEG(image)
      rg <- grid::rasterGrob(r, width=grid::unit(1,"npc"), height=grid::unit(1,"npc"))
      pix_x <- dim(r)[2]
      pix_y <- dim(r)[1]

      coordi <- t(matrix(as.numeric(unlist(strsplit(names(visiualVector),"x"))),nrow=2))
      fig.df <- data.frame(
        x=coordi[,1],
        y=pix_y-coordi[,2],
        value=visiualVector
      )

      ggplot(fig.df,aes(x=x,y=y))+
        annotation_custom(rg)+  # add background image
        geom_point(aes(colour=value), size=pointSize)+
        scale_colour_gradientn(name=legendName, colours = colors,limits=limits)+
        scale_x_continuous(limits = c(0, pix_x), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, pix_y), expand = c(0, 0)) +
        ggtitle(titleName)+
        theme_bw()+
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank()
        )
    }else{
      coordi <- t(matrix(as.numeric(unlist(strsplit(names(visiualVector),"x"))),nrow=2))
      fig.df <- data.frame(
        x=coordi[,1],
        y=coordi[,2],
        value=visiualVector
      )

      ggplot(fig.df,aes(x=x,y=y))+
        geom_point(aes(colour=value), size=pointSize)+
        scale_colour_gradientn(name=legendName, colours = colors, limits=limits)+
        ggtitle(titleName)+
        theme_bw()+
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank()
        )
    }


  }
}
