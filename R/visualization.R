#' @title Quality control metrics visualization
#' @description Visualize quality control metrics in ST dataset.
#' @param SpaCE_obj An SpaCE object.
#' @param itemQC Item for quality control metrics.i.e., "UMI" or "gene".
#' @param colors Legend color scale, Default: c("blue", "yellow", "red").
#' @return A ggplot2 object
#' @examples
#' SpaCE.visualize.metrics(SpaCE_obj,"EPCAM",imageBg = FALSE)
#' @rdname SpaCE.visualize.metrics
#' @export
SpaCE.visualize.metrics <- function(
    SpaCE_obj,
    itemQC=c("UMI","Gene"),
    colors = c("lightblue", "blue", "darkblue"),
    imageBg = TRUE
)
{
  visiualVector <- SpaCE_obj@results$metrics[itemQC,]
  names(visiualVector) <- paste0(SpaCE_obj@input$spotCoordinates[,1],"x",SpaCE_obj@input$spotCoordinates[,2])

  visualSpatial(
    visiualVector,
    image=SpaCE_obj@input$image,
    platform=SpaCE_obj@input$platform,
    scaleType="color-continuous",
    colors=colors,
    pointSize=1,
    pointAlpha=1,
    limits=NULL,
    titleName=itemQC,
    legendName="Count",
    imageBg=imageBg)
}

#' @title Gene expression visualization
#' @description Visualize gene expression level in ST dataset.
#' @param SpaCE_obj An SpaCE object.
#' @param geneSymbol Gene symbol.
#' @param colors Legend color scale, Default: c("blue", "yellow", "red").
#' @return A ggplot2 object
#' @examples
#' SpaCE.visualize.gene(SpaCE_obj,"EPCAM",imageBg = FALSE)
#' @rdname SpaCE.visualize.gene
#' @export
SpaCE.visualize.gene <- function(
    SpaCE_obj,
    gene,
    colors = c("blue", "yellow", "red"),
    imageBg = TRUE
)
{
  visiualVector <- log2(SpaCE_obj@input$counts[gene,]+1)
  names(visiualVector) <- paste0(SpaCE_obj@input$spotCoordinates[,1],"x",SpaCE_obj@input$spotCoordinates[,2])

  visualSpatial(
    visiualVector,
    SpaCE_obj@input$image,
    SpaCE_obj@input$platform,
    scaleType="color-continuous",
    colors=colors,
    pointSize=1,
    pointAlpha=1,
    limits=NULL,
    titleName=gene,
    legendName="Expr",
    imageBg=imageBg)
}

#' @title Cell type fraction visualization
#' @description Visualize cell type fraction in ST dataset.
#' @param SpaCE_obj An SpaCE object.
#' @param cellType Cell type name.
#' @param colors Legend color scale, Default: c("blue", "yellow", "red").
#' @param limits Value range, Default: c(0,1).
#' @return A ggplot2 object
#' @examples
#' SpaCE.visualize.deconvolution(SpaCE_obj,"Malignant")
#' @rdname SpaCE.visualize.deconvolution
#' @export
SpaCE.visualize.deconvolution <- function(
    SpaCE_obj,
    cellType,
    colors = c("blue", "yellow", "red"),
    limits = c(0,1),
    interactive=FALSE,
    imageBg = TRUE
)
{
  if(interactive)
  {
    visualSpatialInteractive(SpaCE_obj,cellType)
  }else{
    visiualVector <- SpaCE_obj@results$deconvolution[cellType,]
    names(visiualVector) <- paste0(SpaCE_obj@input$spotCoordinates[,1],"x",SpaCE_obj@input$spotCoordinates[,2])

    visualSpatial(
      visiualVector,
      SpaCE_obj@input$image,
      SpaCE_obj@input$platform,
      scaleType="color-continuous",
      colors=colors,
      pointSize=1,
      pointAlpha=1,
      limits=limits,
      titleName=cellType,
      legendName="Prop",
      imageBg=imageBg)
    }
}

#' @title Ligand-Receptor network score visualization
#' @description Visualize L-R network score in ST dataset.
#' @param SpaCE_obj An SpaCE object.
#' @param colors Legend color scale, Default: c("blue", "yellow", "red").
#' @param limits Value range, Default: c(0,2).
#' @return A ggplot2 object
#' @examples
#' SpaCE.visualize.LRNetworkScore(SpaCE_obj)
#' @rdname SpaCE.visualize.LRNetworkScore
#' @export
SpaCE.visualize.LRNetworkScore <- function(
    SpaCE_obj,
    colors = c("black","#081d58","#253494","blue","blue","blue","blue","cyan","cyan","yellow"),
    interactive=FALSE,
    imageBg = TRUE
)
{
    visiualVector <- SpaCE_obj@results$LRNetworkScore["Network_Score",]
    names(visiualVector) <- paste0(SpaCE_obj@input$spotCoordinates[,1],"x",SpaCE_obj@input$spotCoordinates[,2])
    visiualVector[visiualVector>1.5] <- 1.5

    visualSpatial(
      visiualVector,
      SpaCE_obj@input$image,
      SpaCE_obj@input$platform,
      scaleType="color-continuous",
      colors=colors,
      pointSize=1,
      pointAlpha=1,
      limits=range(visiualVector),
      titleName="Network_Score",
      legendName="Score",
      imageBg=imageBg)
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

  if(platform=="Visium")
  {
    coordi <- t(matrix(as.numeric(unlist(strsplit(names(visiualVector),"x"))),nrow=2))

    if(imageBg& !is.na(image))
    {
      r <- png::readPNG(image)
      rg <- grid::rasterGrob(r, width=grid::unit(1,"npc"), height=grid::unit(1,"npc"))
      xDiml <- dim(r)[1] # dim pixel
      yDiml <- dim(r)[2] # dim pixel
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

    p <- ggplot(fig.df,aes(x=x,y=y))

    if(imageBg& !is.na(image))
    {
      p <- p+annotation_custom(rg)+# add background image
        scale_x_continuous(limits = c(0, xDiml), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, yDiml), expand = c(0, 0))
    }

    p <- p+geom_point(aes(colour=value),size=pointSize,alpha=pointAlpha)+
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

    if(scaleType=="color-continuous")
    {
      p+scale_colour_gradientn(name=legendName, colours = colors, limits=limits)
    }else{
      p+scale_colour_manual(name=legendName,values=colors)
    }

  }else{

    if(imageBg& !is.na(image))
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
        geom_point(aes(colour=value))+
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
        geom_point(aes(colour=value))+
        scale_colour_gradientn(name=legendName, colours = colors,limits=limits)+
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

visualSpatialInteractive <- function(SpaCE_obj,gene)
{
  library(shiny)

  app <- list(
    ui=fluidPage(

      titlePanel("Interactive visualization"),

      sidebarLayout(

        sidebarPanel(
          p("Select a cell type",style="color:black; font-weight: bold; margin-bottom:33px"),
          selectInput("cellType", p("Cell type:",style="color:black; text-align:center"), choices = rownames(SpaCE_obj@results$deconvolution)),
          br(),
          sliderInput("pointSize", "Spot size", min=0, max=2, value=1, step=0.2),
          br(),
          sliderInput("pointAlpha", "Spot alpha", min=0, max=1, value=1, step=0.1),
          br(),
          style="background-color:papayawhip;border-left:8px solid orange"
        ),

        mainPanel(column(br(),plotOutput("overlayPlot"),br(), width = 10, style="border:1px solid black"))
      )

    ),

    server=function(input, output) {
      output$overlayPlot <- renderPlot({
        cellType <- input$cellType
        pointSize <- input$pointSize
        pointAlpha <- input$pointAlpha

        visiualVector <- SpaCE_obj@results$deconvolution[cellType,]
        names(visiualVector) <- paste0(SpaCE_obj@input$spotCoordinates[1,],"x",SpaCE_obj@input$spotCoordinates[2,])

        visualSpatial(
          visiualVector,
          SpaCE_obj@input$image,
          SpaCE_obj@input$platform,
          scaleType="color-continuous",
          colors=c("blue", "yellow", "red"),
          pointSize=pointSize,
          pointAlpha=pointAlpha,
          limits=c(0,1),
          titleName=cellType,
          legendName="Prop",
          imageBg=TRUE
        )

      })

    }

  )
  shiny::runApp(app)
}
