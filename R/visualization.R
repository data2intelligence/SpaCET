#' @title Spatial feature visualization
#' @description Visualize multiple types of spatial features in ST data.
#' @param SpaCET_obj A SpaCET object.
#' @param spatialType Type of spatial features, i.e., "QualityControl", "GeneExpression", "CellFraction", and "LRNetworkScore". See ‘details’ for more information.
#' @param spatialFeatures A vector of spatial features.
#' @param scaleTypeForGeneExpression Scale type of gene expression, i.e., "RawCounts","LogRawCounts","LogTPM/10", and "LogTPM".
#' @param sameScaleForFraction Indicate whether all cell types have the same scale for cell fraction.
#' @param colors A vector of colors for spots.
#' @param pointSize Size of spots.
#' @param pointAlpha Alpha transparency scales of spots. Must lie between 0 and 1.
#' @param nrow Row number of the combined panel for multiple spatial features.
#' @param imageBg Logical: should the image be shown as background?
#' @param imageSize Size of the image, i.e., "CompleteImage", "CaptureArea", "CustomizedArea".
#' @param CustomizedAreaScale A vector of four numbers (0~1) for scale of the Customized Area, i.e., x_left, x_right, y_bottom, y_top.
#' @param legend.position The position of the legend. Set it as "none" if you want to remove the legend.
#' @param interactive Logical: should the interactive app be activited?
#' @return A ggplot2 object.
#' @details
#' `SpaCET.visualize.spatialFeature` is able to plot multiple types of spatial features, including "QC", "GeneExp", "CellFraction", and "LRNetworkScore".
#' "QualityControl" refers to quality control metrics. User can visualize the total counts of UMI and genes across all spots.
#' "GeneExpression" is set to visualize the expression level of genes.
#' "CellFraction" refers to the cell fraction of cell types.
#' "LRNetworkScore" is selected to show the Ligand-Receptor network score. See `SpaCET.CCI.LRNetworkScore` for how to calculate it.
#'
#' @examples
#' SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "QualityControl", spatialFeatures=c("UMI","Gene"))
#' SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "GeneExpression", spatialFeatures=c("EPCAM","MS4A1"))
#' SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "CellFraction", spatialFeatures=c("Malignant","Macrophage"))
#' SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "LRNetworkScore", spatialFeatures=c("Network_Score","Network_Score_pv"))
#'
#' @rdname SpaCET.visualize.spatialFeature
#' @export
SpaCET.visualize.spatialFeature <- function(
    SpaCET_obj,
    spatialType,
    spatialFeatures = NULL,
    scaleTypeForGeneExpression = "LogTPM",
    sameScaleForFraction = FALSE,
    colors = NULL,
    pointSize = 1,
    pointAlpha = 1,
    nrow = 1,
    imageBg = TRUE,
    imageSize = "CaptureArea",
    CustomizedAreaScale = NULL,
    legend.position = "right",
    interactive = FALSE
)
{
  if(!imageSize%in%c("CompleteImage", "CaptureArea", "CustomizedArea"))
  {
    stop("Please set imageSize as one of them, CompleteImage, CaptureArea, CustomizedArea.")
  }

  if(interactive==FALSE)
  {

    if(spatialType == "QualityControl")
    {
      if(is.null(SpaCET_obj@results$metrics))
      {
        stop("Please run SpaCET.quality.control first.")
      }

      mat <- SpaCET_obj@results$metrics

      scaleType="color-continuous"
      if(is.null(colors)) colors = c("lightblue", "blue", "darkblue")
      legendName = "Count"
      limits = NULL
    }else if(spatialType == "GeneExpression"){
      mat <- SpaCET_obj@input$counts

      if(scaleTypeForGeneExpression=="RawCounts")
      {
        legendName = "Counts"
      }else if(scaleTypeForGeneExpression=="LogRawCounts"){
        mat@x <- log2(mat@x+1)
        legendName = "LogCounts"
      }else if(scaleTypeForGeneExpression=="LogTPM/10"){
        mat <- Matrix::t(Matrix::t(mat)*1e5/Matrix::colSums(mat))
        mat@x <- log2(mat@x+1)
        legendName = "LogTPM/10"
      }else if(scaleTypeForGeneExpression=="LogTPM"){
        mat <- Matrix::t(Matrix::t(mat)*1e6/Matrix::colSums(mat))
        mat@x <- log2(mat@x+1)
        legendName = "LogTPM"
      }else{
        stop("Please set scaleTypeForGeneExpression as one of four scale types, i.e., RawCounts, LogRawCounts, LogTPM/10, LogTPM.")
      }

      geneFlag <- spatialFeatures%in%rownames(mat)
      if(sum(geneFlag)!=length(geneFlag))
      {
        excluded <- spatialFeatures[!geneFlag]
        warning("The following genes are excluded because they are not official gene symbols.")
        message(excluded)
        spatialFeatures <- spatialFeatures[geneFlag]
      }

      scaleType="color-continuous"
      if(is.null(colors)) colors = c("#4d9221", "yellow", "#c51b7d")
      limits = NULL
    }else if(spatialType == "CellFraction"){
      if(is.null(SpaCET_obj@results$deconvolution$propMat))
      {
        stop("Please run cell type deconvolution first.")
      }

      mat <- SpaCET_obj@results$deconvolution$propMat

      if(!is.list(spatialFeatures))
      {
        if("All"%in%spatialFeatures) spatialFeatures <- rownames(mat)

      }else{
        if(is.null(names(spatialFeatures)) | ""%in%names(spatialFeatures))
        {
          stop("Please assign a name for each element of your cell-type list.")
        }

        for(i in names(spatialFeatures))
        {
          if( sum( spatialFeatures[[i]]%in%rownames(mat) ) != length(spatialFeatures[[i]]) )
          {
            wrongCellTypes <- paste0( spatialFeatures[[i]][!spatialFeatures[[i]]%in%rownames(mat)], collapse=", ")
            stop(paste0("The following cell-type names are not in the deconvolution results. Please check your input.\n", wrongCellTypes))
          }
          mat <- rbind(mat,i=colSums(mat[spatialFeatures[[i]],,drop=F]))
          rownames(mat)[nrow(mat)] <- i
        }
        spatialFeatures <- names(spatialFeatures)
      }

      scaleType="color-continuous"
      if(is.null(colors)) colors = c("blue", "yellow", "red")
      legendName = "Fraction"

      if(sameScaleForFraction)
      {
        limits = c(0,1)
      }else{
        limits = NULL
      }

    }else if(spatialType == "MostAbundantCellType"){
      if(is.null(SpaCET_obj@results$deconvolution$propMat))
      {
        stop("Please run cell type deconvolution first.")
      }

      mat <- SpaCET_obj@results$deconvolution$propMat

      scaleType="color-discrete"
      legendName = "Cell Type"
      limits = NULL

    }else if(spatialType == "CellTypeComposition"){
      if(is.null(SpaCET_obj@results$deconvolution$propMat))
      {
        stop("Please run cell type deconvolution first.")
      }

      mat <- SpaCET_obj@results$deconvolution$propMat

      scaleType="color-discrete"
      legendName = "Cell Type"
      limits = NULL

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
      if(is.null(colors)) colors = c("blue","blue","blue","blue","cyan","cyan","yellow")
      limits = NULL

    }else if(spatialType == "Interface"){
      if(is.null(SpaCET_obj@results$CCI$interface))
      {
        stop("Please run SpaCET.identify.interface first.")
      }

      mat <- SpaCET_obj@results$CCI$interface

      scaleType="color-discrete"
      legendName = "Spot"
      limits = NULL
    }else if(spatialType == "GeneSetScore"){
      if(is.null(SpaCET_obj@results$GeneSetScore))
      {
        stop("Please run SpaCET.GeneSetScore first.")
      }

      mat <- SpaCET_obj@results$GeneSetScore

      scaleType="color-continuous"
      if(is.null(colors)) colors = c("#91bfdb","#fee090","#d73027")
      legendName = "Score"
      limits = NULL
    }else if(spatialType == "SecretedProteinActivity"){
      if(is.null(SpaCET_obj@results$SecAct_output$SecretedProteinActivity))
      {
        stop("Please run SecAct.signaling.inference first.")
      }

      mat <- SpaCET_obj@results$SecAct_output$SecretedProteinActivity$zscore

      scaleType="color-continuous"
      if(is.null(colors)) colors = c("#b8e186","#b8e186","#b8e186","#de77ae","#c51b7d")
      legendName = "Activity"
      limits = NULL
    }else if(spatialType == "SignalingPattern"){
      if(is.null(SpaCET_obj @results $SecAct_output $pattern))
      {
        stop("Please run SecAct.signaling.pattern first.")
      }

      mat <- SpaCET_obj @results $SecAct_output $pattern $signal.H

      if("All"%in%spatialFeatures)
      {
        spatialFeatures <- rownames(mat)
      }else{
        spatialFeatures <- as.character(spatialFeatures)
      }

      scaleType="color-continuous"
      if(is.null(colors)) colors = c("#000004","#1A1042","#4A1079","#D9466B","#FCFDBF")
      legendName = "Signal"
      limits = NULL
    }else if(spatialType == "metaData"){

      mat <- t(SpaCET_obj@input$metaData)

      if(length(spatialFeatures) > 1)
      {
        stop("Please input one feature once.")
      }

      scaleType="color-discrete"
      legendName = spatialFeatures
      limits = NULL
    }else{
      stop("Please set spatialType as one of these spatial feature types, i.e., QualityControl, GeneExpression, CellFraction, MostAbundantCellType, CellTypeComposition, LRNetworkScore, Interface, GeneSetScore, SecretedProteinActivity, and SignalingPattern.")
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

      if(spatialType == "Interface"){
        if(is.null(spatialFeature))
        {
          stop("Please set spatialFeature.")
        }

        if(spatialFeature == "Interface"){
          if(is.null(colors)) colors=c("black","darkgrey","#f3c300")
        }else{
          if(is.null(colors)) colors=c("green","black","darkgrey","#f3c300")
        }
      }

      if(spatialType == "MostAbundantCellType")
      {
        if(!spatialFeature%in%c("MajorLineage","SubLineage"))
        {
          stop("Please set spatialFeatures as MajorLineage or SubLineage")
        }else{
          if(spatialFeature == "MajorLineage")
          {
            allCellTypes <- names(SpaCET_obj@results$deconvolution$Ref$lineageTree)
            if(!"Malignant"%in%allCellTypes & "Malignant"%in%rownames(SpaCET_obj@results$deconvolution$propMat))
            {
              allCellTypes <- c("Malignant",allCellTypes)
            }
          }else{
            allCellTypes <- unlist(SpaCET_obj@results$deconvolution$Ref$lineageTree)
            if(!"Malignant"%in%allCellTypes & "Malignant"%in%rownames(SpaCET_obj@results$deconvolution$propMat))
            {
              allCellTypes <- c("Malignant",allCellTypes)
            }
          }

          res_deconv_level <- mat[allCellTypes,,drop=F]

          Content <- sapply(1:dim(res_deconv_level)[2],function(x) names(sort(res_deconv_level[,x],decreasing=T))[1])
          names(Content) <- colnames(res_deconv_level)

          Content <- Content[order(match(Content,allCellTypes))]

          visiualVector <- Content
        }

      }else if(spatialType == "CellTypeComposition"){
        if(!spatialFeature%in%c("MajorLineage","SubLineage"))
        {
          stop("Please set spatialFeatures as MajorLineage or SubLineage")
        }else{
          if(spatialFeature == "MajorLineage")
          {
            allCellTypes <- names(SpaCET_obj@results$deconvolution$Ref$lineageTree)
            if(!"Malignant"%in%allCellTypes & "Malignant"%in%rownames(SpaCET_obj@results$deconvolution$propMat))
            {
              allCellTypes <- c("Malignant",allCellTypes)
            }
          }else{
            allCellTypes <- unlist(SpaCET_obj@results$deconvolution$Ref$lineageTree)
            if(!"Malignant"%in%allCellTypes & "Malignant"%in%rownames(SpaCET_obj@results$deconvolution$propMat))
            {
              allCellTypes <- c("Malignant",allCellTypes)
            }
          }

          res_deconv_level <- mat[allCellTypes,,drop=F]

          if("Unidentifiable"%in%rownames(mat))
          {
            res_deconv_level <- rbind(res_deconv_level, Unidentifiable=mat["Unidentifiable",] )
          }

          visiualVector <- as.data.frame(res_deconv_level)
        }

      }else{
        if(spatialType == "GeneExpression")
        {
          visiualVector <- sort(mat[spatialFeature,])
        }else{
          visiualVector <- mat[spatialFeature,]
        }
      }

      spotID <- names(visiualVector)
      names(visiualVector) <- paste0(
        SpaCET_obj@input$spotCoordinates[names(visiualVector),1],"x",
        SpaCET_obj@input$spotCoordinates[names(visiualVector),2])

      p <- visualSpatial(
        visiualVector,
        image=SpaCET_obj@input$image,
        platform=SpaCET_obj@input$platform,
        scaleType=scaleType,
        colors=colors,
        pointSize=pointSize,
        pointAlpha=pointAlpha,
        limits=limits,
        titleName=spatialFeature,
        legendName=legendName,
        legend.position=legend.position,
        imageBg=imageBg,
        imageSize=imageSize,
        CustomizedAreaScale=CustomizedAreaScale,
        spotID=spotID)

      library(patchwork)
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

  }else{
    if(!grepl("visium", tolower(SpaCET_obj@input$platform)))
    {
      stop("This function is only applicable to 10X Visium data.")
    }

    library(shiny)
    library(plotly)

    app <- list(
      ui=fluidPage(

        titlePanel("Interactive visualization from SpaCET"),

        sidebarLayout(

          sidebarPanel(
            p("1. Select a spatial feature.",style="color:black; font-weight: bold; font-size:18px;; margin-bottom:33px"),
            selectInput("spatialType", p("Spatial Feature Type:",style="color:black; text-align:center"), ""),
            selectInput("spatialFeature", p("Spatial Feature:",style="color:black; text-align:center"), ""),
            br(),
            p("2. Adjust the style of spots.",style="color:black; font-weight: bold; font-size:18px;; margin-bottom:33px"),
            sliderInput("pointSize", "Spot size", min=0, max=2, value=1, step=0.1),
            br(),
            sliderInput("pointAlpha", "Spot opacity", min=0, max=1, value=1, step=0.1),
            br(),
            style="background-color:papayawhip;border-left:8px solid orange",
            width = 3
          ),

          mainPanel(
            column(
              br(),
              plotlyOutput("overlayPlotly"),
              br(),
              width = 7,
              style="border:1px solid darkgrey; height:592px; border-right:18px"
            ),
            column(
              br(),
              DT::dataTableOutput("se"),
              br(),
              width = 5,
              style="border:1px solid darkgrey; height:592px"
            )
          )

        )

      ),

      server=function(input, output, session) {

        observe({
          spatialTypes <- c("CellFraction","QualityControl","SignalingPattern","LRNetworkScore","Interface")

          if(is.null(SpaCET_obj@results$deconvolution$propMat))
          {
            spatialTypes <- setdiff(spatialTypes,"CellFraction")
          }
          if(is.null(SpaCET_obj@results$metrics))
          {
            spatialTypes <- setdiff(spatialTypes,"QualityControl")
          }
          if(is.null(SpaCET_obj @results $SecAct_output $pattern $signal.H))
          {
            spatialTypes <- setdiff(spatialTypes,"SignalingPattern")
          }
          if(is.null(SpaCET_obj@results$CCI$LRNetworkScore))
          {
            spatialTypes <- setdiff(spatialTypes,"LRNetworkScore")
          }
          if(is.null(SpaCET_obj@results$CCI$interface))
          {
            spatialTypes <- setdiff(spatialTypes,"Interface")
          }
          updateSelectInput(session, "spatialType", choices = spatialTypes, selected=spatialTypes[1])
        })

        observe({
          if(input$spatialType=="QualityControl"){
            spatialFeatures <- rownames(SpaCET_obj@results$metrics)
          }else if(input$spatialType=="CellFraction"){
            spatialFeatures <- rownames(SpaCET_obj@results$deconvolution$propMat)
          }else if(input$spatialType=="SignalingPattern"){
            spatialFeatures <- rownames(SpaCET_obj @results $SecAct_output $pattern $signal.H)
          }else if(input$spatialType=="LRNetworkScore"){
            spatialFeatures <- rownames(SpaCET_obj@results$CCI$LRNetworkScore)[2:3]
          }else if(input$spatialType=="Interface"){
            spatialFeatures <- rownames(SpaCET_obj@results$CCI$interface)
          }else{
            spatialFeatures <- c()
          }
          updateSelectInput(session, "spatialFeature", choices = spatialFeatures, selected=spatialFeatures[1])
        })

        output$overlayPlotly <- renderPlotly({
          pointSize <- input$pointSize
          pointAlpha <- input$pointAlpha

          if(input$spatialType=="QualityControl"){
            mat <- SpaCET_obj@results$metrics
            scaleType="color-continuous"
            colors = c("lightblue", "blue", "darkblue")
            legendName = "Count"
            limits = NULL
          }else if(input$spatialType=="SignalingPattern"){
            mat <- SpaCET_obj @results $SecAct_output $pattern $signal.H
            scaleType="color-continuous"
            colors = c("#000004","#1A1042","#4A1079","#D9466B","#FCFDBF")
            legendName = "Signal"
            limits = NULL
          }else if(input$spatialType=="CellFraction"){
            mat <- SpaCET_obj@results$deconvolution$propMat
            scaleType="color-continuous"
            colors = c("blue", "yellow", "red")
            legendName = "Fraction"
            limits = NULL
          }else if(input$spatialType=="LRNetworkScore"){
            mat <- SpaCET_obj@results$CCI$LRNetworkScore
            mat[2,mat[2,]> 1.5] <- 1.5
            mat[2,mat[2,]< 0.5] <- 0.5
            mat[3,] <- -log10(mat[3,])

            scaleType="color-continuous"
            colors = c("blue","blue","blue","blue","cyan","cyan","yellow")

            if(input$spatialFeature == "Network_Score"){
              legendName = "Score"
            }else{
              legendName = "-log10pv"
            }

            limits = NULL
          }else{
            mat <- SpaCET_obj@results$CCI$interface

            scaleType="color-discrete"

            if(input$spatialFeature == "Interface"){
              colors=c("black","darkgrey","#f3c300")
            }else{
              colors=c("green","black","darkgrey","#f3c300")
            }

            legendName = "Spot"
            limits = NULL
          }

          if( input$spatialFeature%in%rownames(mat) )
          {
            visiualVector <- mat[input$spatialFeature,]
            spotID <- names(visiualVector)
            names(visiualVector) <- paste0(
              SpaCET_obj@input$spotCoordinates[names(visiualVector),1],"x",
              SpaCET_obj@input$spotCoordinates[names(visiualVector),2])

            g <- visualSpatial(
              visiualVector,
              image=SpaCET_obj@input$image,
              platform=SpaCET_obj@input$platform,
              scaleType=scaleType,
              colors=colors,
              pointSize=pointSize,
              pointAlpha=pointAlpha,
              limits=limits,
              titleName=input$spatialFeature,
              legendName=legendName,
              legend.position=legend.position,
              imageBg=imageBg,
              imageSize=imageSize,
              CustomizedAreaScale=CustomizedAreaScale,
              spotID=spotID
              )

            ggplotly(g,layerData=1,height=555) %>%
              plotly::layout(
                dragmode = "lasso",
                images=list(
                  list(
                    #source = base64enc::dataURI(file = SpaCET_obj@input$image$path),
                    xref = "paper",
                    yref = "paper",
                    x= 0,
                    y= 1,
                    sizex = 1,
                    sizey = 1,
                    opacity = 0.5,
                    sizing = "stretch",
                    layer = "bottom"
                  )
                )
              )
          }
        })


        output$se <- DT::renderDataTable(server = FALSE,{

          if(input$spatialType=="QualityControl")
          {
            mat <- SpaCET_obj@results$metrics
          }else if(input$spatialType=="CellFraction"){
            mat <- SpaCET_obj@results$deconvolution$propMat
          }else if(input$spatialType=="LRNetworkScore"){
            mat <- SpaCET_obj@results$CCI$LRNetworkScore
            mat[2,mat[2,]> 1.5] <- 1.5
            mat[2,mat[2,]< 0.5] <- 0.5
            mat[3,] <- -log10(mat[3,])
          }else{
            mat <- SpaCET_obj@results$CCI$interface
          }

          if( input$spatialFeature%in%rownames(mat) )
          {
            visiualVector <- mat[input$spatialFeature,]
            spotID <- names(visiualVector)
            names(visiualVector) <- paste0(
              SpaCET_obj@input$spotCoordinates[names(visiualVector),1],"x",
              SpaCET_obj@input$spotCoordinates[names(visiualVector),2])

            d <- event_data("plotly_selected")

            if(!is.null(d))
            {
              d[,3] <- round(d[,3],3)
              d[,4] <- round(d[,4],3)
              d <- cbind(d, xy="a")
              d <- cbind(d, SpotID="b")
              d <- cbind(d, Value="c")
              d[,"xy"] <- paste0(nrow(SpaCET_obj@input$image$grob$raster)-d[,4],"x",d[,3])
              d[,"SpotID"] <- spotID[names(visiualVector)%in%d[,"xy"]]
              d[,"Value"] <- visiualVector[names(visiualVector)%in%d[,"xy"]]
              d <- d[,c("SpotID","Value")]

              DT::datatable(d, caption = "Your selected spots in the middle panel.",
                            extensions = "Buttons",
                            options = list(paging = TRUE,
                                           scrollX=TRUE,
                                           searching = TRUE,
                                           ordering = TRUE,
                                           dom = 'Bfrtip',
                                           buttons = c('csv', 'excel'),
                                           pageLength=10,
                                           lengthMenu=c(10,20,50) ))
            }else{
              DT::datatable(data.frame("The selected spots in the middle panel will be listed here."),rownames = FALSE, colnames=NULL)
            }
          }

        })


      }

    )
    shiny::runApp(app)
  }
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
    legend.position,
    imageBg,
    imageSize,
    CustomizedAreaScale,
    spotID
)
{
  library(ggplot2)

  if(grepl("visium", tolower(platform)))
  {
    coordi <- t(matrix(as.numeric(unlist(strsplit(names(visiualVector),"x"))),nrow=2))

    if(imageBg & !is.na(image$path))
    {
      if(imageSize=="CaptureArea")
      {
        left_edge <- floor( min(coordi[,1])*0.95)
        right_edge <- ceiling( max(coordi[,1])*1.02 )
        bottom_edge <- floor( min(coordi[,2])*0.95 )
        top_edge <- ceiling( max(coordi[,2])*1.02 )

        image$grob$raster <- image$grob$raster[
          left_edge:right_edge,
          bottom_edge:top_edge
        ]

        coordi[,1] <- coordi[,1]-left_edge
        coordi[,2] <- coordi[,2]-bottom_edge
      }

      if(imageSize=="CustomizedArea")
      {
        if( length(CustomizedAreaScale)!=4 | sum(CustomizedAreaScale>=0 &CustomizedAreaScale<=1)!=4 )
        {
          stop("Please assign four numbers (0~1) to CustomizedAreaScale.")
        }
        range1 <- max(coordi[,1]) - min(coordi[,1])
        range2 <- max(coordi[,2]) - min(coordi[,2])
        CustomizedAreaScale_3 <- 1-CustomizedAreaScale[3]
        CustomizedAreaScale_4 <- 1-CustomizedAreaScale[4]

        left_edge <- floor( (min(coordi[,1]) + (range1 * CustomizedAreaScale_4)) * 0.95)
        right_edge <- ceiling( (min(coordi[,1]) + (range1 * CustomizedAreaScale_3)) * 1.02 )
        bottom_edge <- floor( (min(coordi[,2]) + (range2 * CustomizedAreaScale[1])) * 0.95 )
        top_edge <- ceiling( (min(coordi[,2])+ (range2 * CustomizedAreaScale[2])) * 1.02 )

        image$grob$raster <- image$grob$raster[
          left_edge:right_edge,
          bottom_edge:top_edge
        ]

        coordi[,1] <- coordi[,1]-left_edge
        coordi[,2] <- coordi[,2]-bottom_edge
      }

      xDiml <- dim(image$grob$raster)[1] # dim pixel
      yDiml <- dim(image$grob$raster)[2] # dim pixel

    }else{
      # no image
      xDiml <- max(coordi[,1])
      yDiml <- max(coordi[,2])
    }

    fig.df <- data.frame(
      x=xDiml-coordi[,1],
      y=coordi[,2],
      spotID=spotID
    )
    rownames(fig.df) <- names(visiualVector)

    if(is.vector(visiualVector))
    {
      fig.df <- cbind(fig.df, value=visiualVector)
    }else{
      fig.df <- cbind(fig.df, t(visiualVector) )
    }


    if(is.vector(visiualVector))
    {
      # 1. initiate plot
      p <- ggplot(fig.df,aes(x=x,y=y,text=spotID))

      # 2. add image
      if(imageBg & !is.na(image$path))
      {
        p <- p+
          annotation_custom(image$grob)
      }

      # 3. draw spot
      p <- p+
        geom_point(aes(colour=value), size=pointSize, alpha=pointAlpha)

    }else{

      # 1. initiate plot
      p <- ggplot()

      # 2. add image
      if(imageBg & !is.na(image$path))
      {
        p <- p+
          annotation_custom(image$grob)
      }

      # 3. draw spot
      p <- p+
        scatterpie::geom_scatterpie(
          aes(x=x, y=y, group=spotID, r=pointSize),
          data=fig.df,
          cols=colnames(fig.df)[4:ncol(fig.df)],
          color=NA)+
        scale_fill_manual(values=colors, name=legendName)
    }

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
        panel.border = element_blank(),
        legend.position = legend.position
      )+coord_flip()

    # 6. set color scale
    if(scaleType=="color-continuous")
    {
      p+scale_colour_gradientn(name=legendName, colours = colors, limits=limits)
    }else{
      if(is.vector(visiualVector))
      {
        p+scale_colour_manual(name=legendName, values=colors)
      }else{ # cell type composition no color
        p
      }
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
        geom_point(aes(colour=value), size=pointSize, alpha=pointAlpha)+
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
          panel.border = element_blank(),
          legend.position = legend.position
        )
    }else{
      coordi <- t(matrix(as.numeric(unlist(strsplit(names(visiualVector),"x"))),nrow=2))
      fig.df <- data.frame(
        x=coordi[,1],
        y=coordi[,2],
        value=visiualVector
      )

      p <- ggplot(fig.df,aes(x=x,y=y))+
        geom_point(aes(colour=value), size=pointSize, alpha=pointAlpha)+
        ggtitle(titleName)+
        theme_bw()+
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.position = legend.position
        )

      if(scaleType=="color-continuous")
      {
        p+scale_colour_gradientn(name=legendName, colours = colors, limits=limits)
      }else{
        p+scale_colour_manual(name=legendName, values=colors)
      }

    }


  }
}
