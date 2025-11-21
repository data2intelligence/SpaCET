#' @title Define the SpaCET object
#' @slot input The input data
#' @slot results The results
#'
setClass("SpaCET",
  slots = c(
    input = "list",
    results = "list"
  )
)


#' @title Create a SpaCET object from 10X Visium
#' @description Read an ST dataset to create a SpaCET object.
#' @param visiumPath Path to the Space Ranger output folder. See ‘details’ for more information.
#' @return A SpaCET object.
#' @details
#' If users are analyzing an ST data set from 10X Visium platform, they only need to input "visiumPath".
#' Please make sure that "visiumPath" points to the standard output folder of 10X Space Ranger,
#' which has both (1) sequencing data, i.e., `filtered_feature_bc_matrix.h5` file or `filtered_feature_bc_matrix` folder,
#'
#' The "filtered_feature_bc_matrix" folder includes \cr
#' "barcodes.tsv.gz": spot level barcodes; \cr
#' "features.tsv.gz": list of genes; \cr
#' "matrix.mtx.gz": (sparse) matrix of counts.
#'
#' and (2) image folder `spatial`.
#'
#' The "spatial" folder includes \cr
#' “tissue_positions_list.csv” : barcodes and spatial information; \cr
#' “tissue_lowres_image.png” : hematoxylin and eosin (H&E) image; \cr
#' “scalefactors_json.json” : scaling factors for adjusting the coordinates.
#'
#'
#' @examples
#' visiumPath <- file.path(system.file(package = "SpaCET"), "extdata/Visium_BC")
#' SpaCET_obj <- create.SpaCET.object.10X(visiumPath = visiumPath)
#'
#' @rdname create.SpaCET.object.10X
#' @export
#' @importFrom Matrix readMM
#' @importFrom jsonlite fromJSON
#'
create.SpaCET.object.10X <- function(visiumPath)
{
  if(!file.exists(visiumPath))
  {
    stop("The visiumPath does not exist. Please input the correct path.")
  }

  if(file.exists(paste0(visiumPath,"/filtered_feature_bc_matrix/matrix.mtx.gz")))
  {
    st.matrix.data <- Matrix::readMM(paste0(visiumPath,"/filtered_feature_bc_matrix/matrix.mtx.gz")) #dgT
    st.matrix.data <- methods::as(st.matrix.data, "CsparseMatrix")

    st.matrix.gene <- as.matrix(read.csv(paste0(visiumPath,"/filtered_feature_bc_matrix/features.tsv.gz"),as.is=T,header=F,sep="\t"))
    st.matrix.anno <- as.matrix(read.csv(paste0(visiumPath,"/filtered_feature_bc_matrix/barcodes.tsv.gz"),as.is=T,header=F,sep="\t"))

    if(ncol(st.matrix.gene)==1)
    {
      rownames(st.matrix.data) <- st.matrix.gene[,1]
    }else{
      rownames(st.matrix.data) <- st.matrix.gene[,2]
    }
    colnames(st.matrix.data) <- st.matrix.anno[,1]
  }else{
    st.matrix.data <- Seurat::Read10X_h5(filename = paste0(visiumPath,"/filtered_feature_bc_matrix.h5"))
  }

  jsonFile <- jsonlite::fromJSON(paste0(visiumPath,"/spatial/scalefactors_json.json"))

  if(file.exists(paste0(visiumPath,"/spatial/tissue_positions_list.csv")))
  {
    barcode <- read.csv(paste0(visiumPath,"/spatial/tissue_positions_list.csv"),as.is=T,header=FALSE)
    platform <- "Visium"
  }else if(file.exists(paste0(visiumPath,"/spatial/tissue_positions.csv"))){
    barcode <- read.csv(paste0(visiumPath,"/spatial/tissue_positions.csv"),as.is=T,header=TRUE)
    platform <- "Visium"
  }else{
    barcode <- as.data.frame(arrow::read_parquet(paste0(visiumPath,"/spatial/tissue_positions.parquet")))
    platform <- "VisiumHD"
  }

  colnames(barcode) <- c("barcode","in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres")
  rownames(barcode) <- barcode[,"barcode"]

  barcode <- barcode[barcode[,"in_tissue"]==1,]

  olp <- intersect(colnames(st.matrix.data),rownames(barcode))
  st.matrix.data <- st.matrix.data[,olp,drop=F]
  barcode <- barcode[olp,,drop=F]

  if(file.exists(paste0(visiumPath,"/spatial/tissue_lowres_image.png")))
  {
    imagePath <- paste0(visiumPath,"/spatial/tissue_lowres_image.png")
    barcode[["pixel_row"]] <- round(barcode[,"pxl_row_in_fullres"]*jsonFile$tissue_lowres_scalef,3)
    barcode[["pixel_col"]] <- round(barcode[,"pxl_col_in_fullres"]*jsonFile$tissue_lowres_scalef,3)
  }else{
    imagePath <- paste0(visiumPath,"/spatial/tissue_hires_image.png")
    barcode[["pixel_row"]] <- round(barcode[,"pxl_row_in_fullres"]*jsonFile$tissue_hires_scalef,3)
    barcode[["pixel_col"]] <- round(barcode[,"pxl_col_in_fullres"]*jsonFile$tissue_hires_scalef,3)
  }

  spotCoordinates <- barcode[,c("pixel_row","pixel_col","array_row","array_col")]
  rownames(spotCoordinates) <- paste0(barcode[,"array_row"],"x",barcode[,"array_col"])

  metaData <- barcode[,"barcode",drop=FALSE]
  rownames(metaData) <- paste0(barcode[,"array_row"],"x",barcode[,"array_col"])

  colnames(st.matrix.data) <- rownames(spotCoordinates)

  spotCoordinates[["coordinate_x_um"]] <- spotCoordinates[,"array_col"] * 0.5 * 100
  spotCoordinates[["coordinate_y_um"]] <- spotCoordinates[,"array_row"] * 0.5 * sqrt(3) * 100
  spotCoordinates[["coordinate_y_um"]] <- max(spotCoordinates[["coordinate_y_um"]]) - spotCoordinates[["coordinate_y_um"]]

  SpaCET_obj <- create.SpaCET.object(
    counts=st.matrix.data,
    spotCoordinates=spotCoordinates,
    metaData=metaData,
    imagePath=imagePath,
    platform=platform
  )

  SpaCET_obj
}


#' @title Create a SpaCET object
#' @description Read an ST dataset to create a SpaCET object.
#' @param counts Count matrix with gene name (row) x spot ID (column).
#' @param spotCoordinates Spot coordinate matrix with spot ID (row) x coordinates (column). This matrix should include two columns, i,e., X and Y coordinates, respectively, which represent the position of spots in H&E image.
#' @param imagePath Path to the H&E image file. Can be NA if it is not available.
#' @param platform A character string indicating the platform, i.e., "Visium", "OldST", or "Slide-Seq". "OldST" is the early in situ capturing method from which "Visium" was developed.
#' @return A SpaCET object.
#' @details
#' To create a SpaCET object, user need to input four parameters, i.e., "counts", "spotCoordinates", "imagePath", and "platform".
#' However, if analyzing the Visium data, it is more easy to use `create.SpaCET.object.10X` to read ST data.
#'
#' @examples
#' SpaCET_obj <- create.SpaCET.object(counts, spotCoordinates, imagePath, platform)
#'
#' @rdname create.SpaCET.object
#' @export
#'
create.SpaCET.object <- function(counts, spotCoordinates, metaData=NULL, imagePath=NA, platform)
{
  if(!identical(colnames(counts),rownames(spotCoordinates)))
  {
    stop("The Spot IDs in the counts and spotCoordinates matrices are not identical.")
  }

  if(!is.null(metaData))
  {
    if(!identical(colnames(counts),rownames(metaData)))
    {
      stop("The Spot IDs in the counts and metaData matrices are not identical.")
    }
  }

  if(is.na(imagePath))
  {
    rg <- NA
  }else{
    if(!file.exists(imagePath))
    {
      stop("The image under the imagePath does not exist. Please input the correct path. User can set imagePath=NA if the current ST dataset does not have a matched H&E image.")
    }else{
      if(grepl("visium", tolower(platform)))
      {
        r <- png::readPNG(imagePath)
        rg <- grid::rasterGrob(r, width=grid::unit(1,"npc"), height=grid::unit(1,"npc"))
      }
    }
  }

  st.matrix.data <- methods::as(counts, "CsparseMatrix")

  st.matrix.data <- rm_duplicates(st.matrix.data)

  SpaCET_obj <- methods::new("SpaCET",
    input=list(
      counts=st.matrix.data,
      spotCoordinates=spotCoordinates,
      metaData=metaData,
      image=list(path=imagePath,grob=rg),
      platform=platform
    )
  )

  SpaCET_obj
}


#' @title Filter spatial spots and calculate the QC metrics
#' @description Spots with less than `min.genes` expressed genes would be removed.
#' @param SpaCET_obj A SpaCET object.
#' @param min.genes Minimum number of expressed genes. Default: 1.
#' @return A SpaCET object.
#' @examples
#' SpaCET_obj <- SpaCET.quality.control(SpaCET_obj, min.genes=100)
#' SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "QualityControl", spatialFeatures=c("UMI","Gene"))
#'
#' @rdname SpaCET.quality.control
#' @export
#'
SpaCET.quality.control  <- function(SpaCET_obj, min.genes=1)
{
  message(paste0("Removing spots with less than ",min.genes," expressed genes."))

  st.matrix.data <- SpaCET_obj@input$counts
  expressed.genes<- Matrix::colSums(st.matrix.data>0)

  remaining.spots <- expressed.genes>=min.genes
  remaining.spots.num <- sum(remaining.spots)

  message(paste0(length(remaining.spots)-remaining.spots.num," spots are removed."))
  message(paste0(remaining.spots.num," spots are kept."))

  st.matrix.data <- st.matrix.data[,remaining.spots]

  SpaCET_obj@input$counts <- st.matrix.data

  SpaCET_obj@input$spotCoordinates <- SpaCET_obj@input$spotCoordinates[colnames(st.matrix.data),,drop=F]

  if(!is.null(SpaCET_obj@input$metaData))
  {
    SpaCET_obj@input$metaData <- SpaCET_obj@input$metaData[colnames(st.matrix.data),,drop=F]
  }

  SpaCET_obj@results$metrics <- rbind(
    UMI=Matrix::colSums(st.matrix.data),
    Gene=Matrix::colSums(st.matrix.data>0)
  )

  SpaCET_obj
}


#' @title Convert Seurat to SpaCET
#' @description Convert an Seurat object to a SpaCET object.
#' @param Seurat_obj An Seurat object.
#' @param platform A character string indicating the platform, i.e., "Visium", "OldST", or "Slide-Seq". "OldST" is the early in situ capturing method from which "Visium" was developed.
#' @param visiumPath Path to the Space Ranger output folder (Optional). If setting, this function will retrieve more information from the raw output.
#' @return A SpaCET object.
#' @examples
#' visiumPath <- file.path(system.file(package = "SpaCET"), "extdata/Visium_BC")
#' Seurat_obj <- Seurat::Load10X_Spatial(data.dir = visiumPath)
#' SpaCET_obj <- convert.Seurat(Seurat_obj, platform = "Visium")
#' SpaCET_obj <- SpaCET.quality.control(SpaCET_obj)
#' SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "QualityControl", spatialFeatures = c("UMI","Gene"))
#'
#' @rdname convert.Seurat
#' @export
#'
convert.Seurat <- function(Seurat_obj, platform, visiumPath=NULL)
{
  sliceNum <- length(Seurat_obj@images)

  if(sliceNum==1)
  {
    if(class(Seurat_obj@assays$Spatial)=="Assay5")
    {
      st.matrix.data <- Seurat_obj@assays$Spatial@layers$counts
      colnames(st.matrix.data) <- rownames(Seurat_obj@assays$Spatial@cells)
      rownames(st.matrix.data) <- rownames(Seurat_obj@assays$Spatial@features)
    }else{
      st.matrix.data <- Seurat_obj@assays$Spatial@counts
    }

    slice <- Seurat_obj@images[[1]]

    if(class(slice)=="VisiumV2")
    {
      st.matrix.data <- st.matrix.data[,slice@boundaries$centroids@cells]

      spotCoordinates <- data.frame(
        pixel_row=slice@boundaries$centroids@coords[,1] * slice@scale.factors$lowres,
        pixel_col=slice@boundaries$centroids@coords[,2] * slice@scale.factors$lowres
      )

      metaData <- data.frame(
        barcode=colnames(st.matrix.data)
      )

      if(is.null(visiumPath))
      {
        # VisiumV2 does not have spot ID
        colnames(st.matrix.data) <- paste0(
          spotCoordinates$pixel_row,"x",
          spotCoordinates$pixel_col
        )
      }else{
        if(file.exists(paste0(visiumPath,"/spatial/tissue_positions_list.csv")))
        {
          barcode <- read.csv(paste0(visiumPath,"/spatial/tissue_positions_list.csv"),as.is=T,header=FALSE)
        }else if(file.exists(paste0(visiumPath,"/spatial/tissue_positions.csv"))){
          barcode <- read.csv(paste0(visiumPath,"/spatial/tissue_positions.csv"),as.is=T,header=TRUE)
        }else{
          barcode <- as.data.frame(arrow::read_parquet(paste0(visiumPath,"/spatial/tissue_positions.parquet")))
        }
        colnames(barcode) <- c("barcode","in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres")
        rownames(barcode) <- barcode[,"barcode"]
        barcode <- barcode[colnames(st.matrix.data),]

        colnames(st.matrix.data) <- paste0(barcode[,"array_row"],"x",barcode[,"array_col"])

        spotCoordinates[["array_row"]] <- barcode[,"array_row"]
        spotCoordinates[["array_col"]] <- barcode[,"array_col"]
      }

    }else{ #V1
      spotCoordinates <- data.frame(
        pixel_row=slice@coordinates$imagerow * slice@scale.factors$lowres,
        pixel_col=slice@coordinates$imagecol * slice@scale.factors$lowres,
        array_row=slice@coordinates$row,
        array_col=slice@coordinates$col
      )

      metaData <- data.frame(
        barcode=colnames(st.matrix.data)
      )

      colnames(st.matrix.data) <- paste0(
        slice@coordinates$row,"x",
        slice@coordinates$col
      )

    }

    rownames(spotCoordinates) <- colnames(st.matrix.data)
    rownames(metaData) <- colnames(st.matrix.data)

    if("array_col"%in%colnames(spotCoordinates))
    {
      spotCoordinates[["coordinate_x_um"]] <- spotCoordinates[,"array_col"] * 0.5 * 100
      spotCoordinates[["coordinate_y_um"]] <- spotCoordinates[,"array_row"] * 0.5 * sqrt(3) * 100
      spotCoordinates[["coordinate_y_um"]] <- max(spotCoordinates[["coordinate_y_um"]]) - spotCoordinates[["coordinate_y_um"]]
    }

    st.matrix.data <- rm_duplicates(st.matrix.data)

    rg <- grid::rasterGrob(slice@image, width=grid::unit(1,"npc"), height=grid::unit(1,"npc"))

    SpaCET_obj <- methods::new("SpaCET",
      input=list(
        counts=st.matrix.data,
        spotCoordinates=spotCoordinates,
        metaData=metaData,
        image=list(path="FromSeurat",grob=rg),
        platform=platform
      )
    )

    SpaCET_obj
  }else{

    SpaCET_obj_list <- list()
    spotNum <- c(0)
    for(i in 1:sliceNum)
    {
      slice <- Seurat_obj@images[[i]]
      spotNum <- c(spotNum, nrow(slice@coordinates))

      spot_start <- sum(spotNum[1:i])+1
      spot_end <- spot_start+ nrow(slice@coordinates) - 1

      if(class(Seurat_obj@assays$Spatial)=="Assay5")
      {
        st.matrix.data <- Seurat_obj@assays$Spatial@layers$counts
        colnames(st.matrix.data) <- rownames(Seurat_obj@assays$Spatial@cells)
        rownames(st.matrix.data) <- rownames(Seurat_obj@assays$Spatial@features)
      }else{
        st.matrix.data <- Seurat_obj@assays$Spatial@counts
      }

      st.matrix.data <- st.matrix.data[,spot_start:spot_end]
      colnames(st.matrix.data) <- sapply(strsplit(colnames(st.matrix.data),"_",fixed=T),function(x) return(x[1])) #remove id

      if(class(slice)=="VisiumV2")
      {
        st.matrix.data <- st.matrix.data[,slice@boundaries$centroids@cells]

        spotCoordinates <- data.frame(
          pixel_row=slice@boundaries$centroids@coords[,1] * slice@scale.factors$lowres,
          pixel_col=slice@boundaries$centroids@coords[,2] * slice@scale.factors$lowres
        )

        metaData <- data.frame(
          barcode=colnames(st.matrix.data)
        )

        if(is.null(visiumPath))
        {
          # VisiumV2 does not have spot ID
          colnames(st.matrix.data) <- paste0(
            spotCoordinates$pixel_row,"x",
            spotCoordinates$pixel_col
          )
        }else{
          if(file.exists(paste0(visiumPath,"/spatial/tissue_positions_list.csv")))
          {
            barcode <- read.csv(paste0(visiumPath,"/spatial/tissue_positions_list.csv"),as.is=T,header=FALSE)
          }else if(file.exists(paste0(visiumPath,"/spatial/tissue_positions.csv"))){
            barcode <- read.csv(paste0(visiumPath,"/spatial/tissue_positions.csv"),as.is=T,header=TRUE)
          }else{
            barcode <- as.data.frame(arrow::read_parquet(paste0(visiumPath,"/spatial/tissue_positions.parquet")))
          }
          colnames(barcode) <- c("barcode","in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres")
          rownames(barcode) <- barcode[,"barcode"]
          barcode <- barcode[colnames(st.matrix.data),]

          colnames(st.matrix.data) <- paste0(barcode[,"array_row"],"x",barcode[,"array_col"])

          spotCoordinates[["array_row"]] <- barcode[,"array_row"]
          spotCoordinates[["array_col"]] <- barcode[,"array_col"]
        }

      }else{
        spotCoordinates <- data.frame(
          pixel_row=slice@coordinates$imagerow * slice@scale.factors$lowres,
          pixel_col=slice@coordinates$imagecol * slice@scale.factors$lowres,
          array_row=slice@coordinates$row,
          array_col=slice@coordinates$col
        )

        metaData <- data.frame(
          barcode=colnames(st.matrix.data)
        )

        colnames(st.matrix.data) <- paste0(
          slice@coordinates$row,"x",
          slice@coordinates$col
        )

      }

      rownames(spotCoordinates) <- colnames(st.matrix.data)
      rownames(metaData) <- colnames(st.matrix.data)

      if("array_col"%in%colnames(spotCoordinates))
      {
        spotCoordinates[["coordinate_x_um"]] <- spotCoordinates[,"array_col"] * 0.5 * 100
        spotCoordinates[["coordinate_y_um"]] <- spotCoordinates[,"array_row"] * 0.5 * sqrt(3) * 100
        spotCoordinates[["coordinate_y_um"]] <- max(spotCoordinates[["coordinate_y_um"]]) - spotCoordinates[["coordinate_y_um"]]
      }

      st.matrix.data <- rm_duplicates(st.matrix.data)

      rg <- grid::rasterGrob(slice@image, width=grid::unit(1,"npc"), height=grid::unit(1,"npc"))

      SpaCET_obj <- methods::new("SpaCET",
        input=list(
          counts=st.matrix.data,
          spotCoordinates=spotCoordinates,
          metaData=metaData,
          image=list(path="FromSeurat",grob=rg),
          platform=platform
        )
      )

      SpaCET_obj_list[[i]] <- SpaCET_obj
    }

    SpaCET_obj_list
  }
}


#' @title Add SpaCET to Seurat
#' @description Add deconvolution results from a SpaCET object to an Seurat object as a new assay.
#' @param SpaCET_obj A SpaCET object.
#' @param Seurat_obj A Seurat object.
#' @return An Seurat object.
#' @examples
#' visiumPath <- file.path(system.file(package = "SpaCET"), "extdata/Visium_BC")
#' Seurat_obj <- Seurat::Load10X_Spatial(data.dir = visiumPath)
#' SpaCET_obj <- convert.Seurat(Seurat_obj)
#' SpaCET_obj <- SpaCET.quality.control(SpaCET_obj)
#' SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType="BRCA", coreNo=6)
#' Seurat_obj <- addTo.Seurat(SpaCET_obj, Seurat_obj)
#' Seurat::DefaultAssay(Seurat_obj) <- "propMatFromSpaCET"
#' Seurat::SpatialFeaturePlot(Seurat_obj, features = c("CAF", "Macrophage"))
#'
#' @rdname addTo.Seurat
#' @export
#'
addTo.Seurat  <- function(SpaCET_obj, Seurat_obj)
{
  sliceNum <- length(SpaCET_obj)

  if(sliceNum==1)
  {
    propMat <- SpaCET_obj@results$deconvolution$propMat
    colnames(propMat) <- SpaCET_obj@input$metaData[,"barcode"]

    if(ncol(propMat)!=nrow(Seurat_obj@meta.data)) stop("SpaCET_obj and Seurat_obj have different spot number.")

    Seurat_obj[["propMatFromSpaCET"]] <- Seurat::CreateAssayObject(data=propMat)
    Seurat_obj

  }else{

    for(i in 1:length(SpaCET_obj))
    {
      propMat <- SpaCET_obj[[i]]@results$deconvolution$propMat

      if(i==1)
      {
        propMatComb <- propMat
      }else{
        propMatComb <- cbind(propMatComb,propMat)
      }

    }
    colnames(propMatComb) <- colnames(Seurat_obj@assays$Spatial@counts)
    Seurat_obj[["propMatFromSpaCET"]] <- Seurat::CreateAssayObject(data=propMatComb)
    Seurat_obj
  }
}


rm_duplicates <- function(mat){
  gene_count <- table(rownames(mat))
  gene_dupl <- names(gene_count)[gene_count>1]

  if(length(gene_dupl) > 0){
    gene_unique <- names(gene_count)[gene_count==1]
    gene_unique_index <- which(rownames(mat)%in%gene_unique)

    gene_dupl_index <- c()
    for(gene in gene_dupl)
    {
      gene_dupl_index_gene <- which(rownames(mat)%in%gene)
      mat_dupl_gene <- mat[gene_dupl_index_gene,]
      dupl_sum <- Matrix::rowSums(mat_dupl_gene)
      max_flag <- which(dupl_sum==max(dupl_sum))
      gene_dupl_index <- c(gene_dupl_index,gene_dupl_index_gene[max_flag[1]])
    }

    mat <- mat[sort(c(gene_unique_index,gene_dupl_index)),]
  }

  return(mat)
}
