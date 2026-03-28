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

# Visium hexagonal grid geometry constants
VISIUM_SPOT_SPACING_UM <- 100
VISIUM_HEX_SCALE_X <- 0.5
VISIUM_HEX_SCALE_Y <- 0.5 * sqrt(3)

addVisiumMicrometerCoords <- function(spotCoordinates)
{
  spotCoordinates[["coordinate_x_um"]] <- spotCoordinates[,"array_col"] * VISIUM_HEX_SCALE_X * VISIUM_SPOT_SPACING_UM
  spotCoordinates[["coordinate_y_um"]] <- spotCoordinates[,"array_row"] * VISIUM_HEX_SCALE_Y * VISIUM_SPOT_SPACING_UM
  spotCoordinates[["coordinate_y_um"]] <- max(spotCoordinates[["coordinate_y_um"]]) - spotCoordinates[["coordinate_y_um"]]
  spotCoordinates
}


#' @title Create a SpaCET object from 10X Visium
#' @description Read an ST dataset to create a SpaCET object.
#' @param visiumPath Path to the Space Ranger output folder. See ‘details’ for more information.
#' @param organism Organism of the sample, e.g., human or mouse.
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
create.SpaCET.object.10X <- function(visiumPath, organism="human")
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

  spotCoordinates <- addVisiumMicrometerCoords(spotCoordinates)

  SpaCET_obj <- create.SpaCET.object(
    counts=st.matrix.data,
    spotCoordinates=spotCoordinates,
    metaData=metaData,
    imagePath=imagePath,
    platform=platform,
    organism=organism
  )

  SpaCET_obj
}


#' @title Create a SpaCET object
#' @description Read an ST dataset to create a SpaCET object.
#' @param counts Count matrix with gene name (row) x spot ID (column).
#' @param spotCoordinates Spot coordinate matrix with spot ID (row) x coordinates (column). This matrix should include two columns, i,e., X and Y coordinates, respectively, which represent the position of spots in H&E image.
#' @param imagePath Path to the H&E image file. Can be NA if it is not available.
#' @param platform A character string indicating the platform, i.e., "Visium", "OldST", or "Slide-Seq". "OldST" is the early in situ capturing method from which "Visium" was developed.
#' @param organism Organism of the sample, e.g., human or mouse.
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
create.SpaCET.object <- function(counts, spotCoordinates, metaData=NULL, imagePath=NA, platform, organism="human")
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
      platform=platform,
      organism=organism
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
#' @param organism Organism of the sample, e.g., human or mouse.
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
convert.Seurat <- function(Seurat_obj, platform, visiumPath=NULL, organism="human")
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
        platform=platform,
        organism=organism
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
          platform=platform,
          organism=organism
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


mouse2human_mat <- function(mat) {
  m2h <- read.csv( system.file("extdata",'Mouse2Human_filter.csv',package = 'SpaCET'), row.names=1)
  m2h <- m2h[,c("mouse","human")]

  if (is.null(rownames(mat))) {
    stop("mat must have rownames.")
  }

  # keep only mapped mouse genes
  idx <- match(rownames(mat), m2h$mouse)
  keep <- !is.na(idx)

  mat2 <- mat[keep, , drop = FALSE]
  human <- m2h$human[idx[keep]]

  if (nrow(mat2) == 0) {
    stop("No mouse genes matched the mapping table.")
  }

  # group rows by human gene
  grp <- factor(human)
  lev <- levels(grp)

  # sparse aggregation matrix: one row per human gene
  A <- Matrix::sparseMatrix(
    i = as.integer(grp),
    j = seq_along(grp),
    x = 1,
    dims = c(length(lev), length(grp)),
    dimnames = list(lev, NULL)
  )

  # collapse mouse rows to human rows
  out <- A %*% mat2

  rownames(out) <- lev
  out
}


#' @title Create a SpaCET object from NanoString CosMx
#'
#' @description Read NanoString CosMx SMI output into a SpaCET object.
#'
#' @param cosmxPath Path to the CosMx output folder containing the expression matrix and metadata.
#' @param fov Field of view indices to include. Default NULL loads all FOVs.
#' @param organism Species, either "human" (default) or "mouse".
#'
#' @return A SpaCET object.
#'
#' @examples
#' \dontrun{
#' SpaCET_obj <- create.SpaCET.object.CosMx(cosmxPath = "/path/to/cosmx_output")
#' }
#'
#' @rdname create.SpaCET.object.CosMx
#' @export
#'
create.SpaCET.object.CosMx <- function(cosmxPath, fov=NULL, organism="human")
{
  if(!file.exists(cosmxPath))
  {
    stop("The cosmxPath does not exist. Please input the correct path.")
  }

  # Find expression matrix file
  expr_file <- list.files(cosmxPath, pattern="exprMat_file\\.csv$|tx_file\\.csv$", full.names=TRUE, recursive=TRUE)
  if(length(expr_file) == 0)
  {
    # Try alternative naming patterns
    expr_file <- list.files(cosmxPath, pattern="(expression|counts).*\\.csv$", full.names=TRUE, recursive=TRUE)
  }
  if(length(expr_file) == 0) stop("No expression matrix file found in cosmxPath.")

  message("Reading CosMx expression data...")
  expr_data <- read.csv(expr_file[1], row.names=1, check.names=FALSE)

  # Find metadata file
  meta_file <- list.files(cosmxPath, pattern="metadata_file\\.csv$|metadata\\.csv$", full.names=TRUE, recursive=TRUE)
  if(length(meta_file) > 0)
  {
    message("Reading CosMx metadata...")
    meta_data <- read.csv(meta_file[1], row.names=1, check.names=FALSE)

    # Filter by FOV if requested
    if(!is.null(fov) && "fov" %in% colnames(meta_data))
    {
      keep_cells <- rownames(meta_data)[meta_data$fov %in% fov]
      expr_data <- expr_data[, colnames(expr_data) %in% keep_cells, drop=FALSE]
      meta_data <- meta_data[keep_cells, , drop=FALSE]
    }

    # Extract spatial coordinates
    coord_cols <- intersect(c("CenterX_global_px", "CenterY_global_px",
                               "x_global_px", "y_global_px",
                               "CenterX_local_px", "CenterY_local_px"), colnames(meta_data))
    if(length(coord_cols) >= 2)
    {
      spotCoordinates <- data.frame(
        x = meta_data[, coord_cols[1]],
        y = meta_data[, coord_cols[2]],
        row.names = rownames(meta_data)
      )
    } else {
      stop("Could not find spatial coordinate columns in CosMx metadata.")
    }
  } else {
    stop("No metadata file found in cosmxPath.")
  }

  # Transpose if genes are in rows (CosMx format: cells x genes)
  if(nrow(expr_data) > ncol(expr_data))
  {
    # Likely cells x genes, transpose to genes x cells
    st.matrix.data <- Matrix::Matrix(t(as.matrix(expr_data)), sparse=TRUE)
  } else {
    st.matrix.data <- Matrix::Matrix(as.matrix(expr_data), sparse=TRUE)
  }

  # Ensure cell IDs match between expression and coordinates
  common_cells <- intersect(colnames(st.matrix.data), rownames(spotCoordinates))
  if(length(common_cells) == 0) stop("No matching cell IDs between expression data and coordinates.")

  st.matrix.data <- st.matrix.data[, common_cells]
  spotCoordinates <- spotCoordinates[common_cells, ]

  # Look for composite image
  image_files <- list.files(cosmxPath, pattern="CellComposite.*\\.jpg$|CellComposite.*\\.png$|composite.*\\.png$",
                             full.names=TRUE, recursive=TRUE)
  imagePath <- if(length(image_files) > 0) image_files[1] else NA

  message(paste0("CosMx data loaded: ", nrow(st.matrix.data), " genes x ", ncol(st.matrix.data), " cells"))

  SpaCET_obj <- create.SpaCET.object(
    counts=st.matrix.data,
    spotCoordinates=spotCoordinates,
    metaData=if(exists("meta_data")) meta_data[common_cells, , drop=FALSE] else NULL,
    imagePath=imagePath,
    platform="CosMx",
    organism=organism
  )

  return(SpaCET_obj)
}


#' @title Create a SpaCET object from 10x Xenium
#'
#' @description Read 10x Genomics Xenium output into a SpaCET object.
#'
#' @param xeniumPath Path to the Xenium output folder.
#' @param organism Species, either "human" (default) or "mouse".
#'
#' @return A SpaCET object.
#'
#' @examples
#' \dontrun{
#' SpaCET_obj <- create.SpaCET.object.Xenium(xeniumPath = "/path/to/xenium_output")
#' }
#'
#' @rdname create.SpaCET.object.Xenium
#' @export
#'
create.SpaCET.object.Xenium <- function(xeniumPath, organism="human")
{
  if(!file.exists(xeniumPath))
  {
    stop("The xeniumPath does not exist. Please input the correct path.")
  }

  # Read cell-feature matrix
  matrix_dir <- file.path(xeniumPath, "cell_feature_matrix")
  h5_file <- file.path(xeniumPath, "cell_feature_matrix.h5")

  if(dir.exists(matrix_dir))
  {
    message("Reading Xenium cell-feature matrix from directory...")
    st.matrix.data <- Matrix::readMM(file.path(matrix_dir, "matrix.mtx.gz"))
    genes <- read.csv(file.path(matrix_dir, "features.tsv.gz"), header=FALSE, sep="\t")
    barcodes <- read.csv(file.path(matrix_dir, "barcodes.tsv.gz"), header=FALSE, sep="\t")
    rownames(st.matrix.data) <- genes[, 2]  # Gene symbols
    colnames(st.matrix.data) <- barcodes[, 1]
  } else if(file.exists(h5_file))
  {
    message("Reading Xenium cell-feature matrix from H5...")
    st.matrix.data <- Seurat::Read10X_h5(h5_file)
  } else {
    stop("No cell_feature_matrix found. Expected directory or .h5 file in xeniumPath.")
  }

  # Read cell metadata with spatial coordinates
  cells_file <- file.path(xeniumPath, "cells.csv.gz")
  cells_parquet <- file.path(xeniumPath, "cells.parquet")

  if(file.exists(cells_file))
  {
    message("Reading Xenium cell coordinates...")
    cells <- read.csv(cells_file)
  } else if(file.exists(cells_parquet))
  {
    cells <- as.data.frame(arrow::read_parquet(cells_parquet))
  } else {
    stop("No cells.csv.gz or cells.parquet found in xeniumPath.")
  }

  # Extract coordinates
  rownames(cells) <- cells$cell_id
  spotCoordinates <- data.frame(
    x = cells$x_centroid,
    y = cells$y_centroid,
    row.names = cells$cell_id
  )

  # Align cell IDs
  common_cells <- intersect(colnames(st.matrix.data), rownames(spotCoordinates))
  if(length(common_cells) == 0) stop("No matching cell IDs between expression data and coordinates.")

  st.matrix.data <- st.matrix.data[, common_cells]
  spotCoordinates <- spotCoordinates[common_cells, ]

  message(paste0("Xenium data loaded: ", nrow(st.matrix.data), " genes x ", ncol(st.matrix.data), " cells"))

  SpaCET_obj <- create.SpaCET.object(
    counts=st.matrix.data,
    spotCoordinates=spotCoordinates,
    metaData=cells[common_cells, , drop=FALSE],
    imagePath=NA,
    platform="Xenium",
    organism=organism
  )

  return(SpaCET_obj)
}
