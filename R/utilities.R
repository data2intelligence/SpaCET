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


#' @title Create an SpaCET object from 10X Visium
#' @description Read an ST dataset to create an SpaCET object.
#' @param visiumPath Path to the Space Ranger output folder. See ‘details’ for more information.
#' @return An SpaCET object.
#' @details
#' If users are analyzing an ST data set from 10X Visium platform, they only need to input "visiumPath".
#' Please make sure that "visiumPath" points to the standard output folder of 10X Space Ranger,
#' which have both `filtered_feature_bc_matrix` and `spatial` folders.
#'
#' The "filtered_feature_bc_matrix" folder includes \cr
#' "barcodes.tsv.gz": spot level barcodes; \cr
#' "features.tsv.gz": list of genes; \cr
#' "matrix.mtx.gz": (sparse) matrix of counts.
#'
#' The "spatial" folder includes \cr
#' “tissue_positions_list.csv” : barcodes and spatial information; \cr
#' “tissue_lowres_image.png” : hematoxylin and eosin (H&E) image; \cr
#' “scalefactors_json.json” : scaling factors for adjusting the coordinates.
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

  st.matrix.data <- Matrix::readMM(paste0(visiumPath,"/filtered_feature_bc_matrix/matrix.mtx.gz")) #dgT
  st.matrix.data <- methods::as(st.matrix.data, "dgCMatrix")

  st.matrix.gene <- as.matrix(read.csv(paste0(visiumPath,"/filtered_feature_bc_matrix/features.tsv.gz"),as.is=T,header=F,sep="\t"))
  st.matrix.anno <- as.matrix(read.csv(paste0(visiumPath,"/filtered_feature_bc_matrix/barcodes.tsv.gz"),as.is=T,header=F,sep="\t"))

  rownames(st.matrix.data) <- st.matrix.gene[,2]
  colnames(st.matrix.data) <- st.matrix.anno[,1]


  jsonFile <- jsonlite::fromJSON(paste0(visiumPath,"/spatial/scalefactors_json.json"))

  if(file.exists(paste0(visiumPath,"/spatial/tissue_positions_list.csv")))
  {
    barcode <- read.csv(paste0(visiumPath,"/spatial/tissue_positions_list.csv"),as.is=T,header=FALSE)
  }else{
    barcode <- read.csv(paste0(visiumPath,"/spatial/tissue_positions.csv"),as.is=T,header=TRUE)
  }

  colnames(barcode) <- c("barcode","in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres")
  rownames(barcode) <- paste0(barcode[,"array_row"],"x",barcode[,"array_col"])

  if(file.exists(paste0(visiumPath,"/spatial/tissue_lowres_image.png")))
  {
    imagePath <- paste0(visiumPath,"/spatial/tissue_lowres_image.png")
    barcode[["pxl_row"]] <- round(barcode[,"pxl_row_in_fullres"]*jsonFile$tissue_lowres_scalef,3)
    barcode[["pxl_col"]] <- round(barcode[,"pxl_col_in_fullres"]*jsonFile$tissue_lowres_scalef,3)
  }else{
    imagePath <- paste0(visiumPath,"/spatial/tissue_hires_image.png")
    barcode[["pxl_row"]] <- round(barcode[,"pxl_row_in_fullres"]*jsonFile$tissue_hires_scalef,3)
    barcode[["pxl_col"]] <- round(barcode[,"pxl_col_in_fullres"]*jsonFile$tissue_hires_scalef,3)
  }

  spotCoordinates <- barcode[match(colnames(st.matrix.data),barcode[,"barcode"]),c("pxl_row","pxl_col")]
  colnames(st.matrix.data) <- rownames(spotCoordinates)

  SpaCET_obj <- create.SpaCET.object(
    counts=st.matrix.data,
    spotCoordinates=spotCoordinates,
    imagePath=imagePath,
    platform="Visium"
  )

  SpaCET_obj
}


#' @title Create an SpaCET object
#' @description Read an ST dataset to create an SpaCET object.
#' @param counts Count matrix with gene name (row) x spot ID (column).
#' @param spotCoordinates Spot coordinate matrix with spot ID (row) x coordinates (column). This matrix should include two columns, i,e., X and Y coordinates, respectively, which represent the position of spots in H&E image.
#' @param imagePath Path to the H&E image file. Can be NA if it is not available.
#' @param platform A character string indicating the platform, i.e., "Visium", "OldST", or "Slide-Seq". "OldST" is the early in situ capturing method from which "Visium" was developed.
#' @return An SpaCET object.
#' @details
#' To create an SpaCET object, user need to input four parameters, i.e., "counts", "spotCoordinates", "imagePath", and "platform".
#' However, if analyzing the Visium data, it is more easy to use `create.SpaCET.object.10X` to read ST data.
#'
#' @examples
#' SpaCET_obj <- create.SpaCET.object(counts, spotCoordinates, imagePath, platform)
#'
#' @rdname create.SpaCET.object
#' @export
#'
create.SpaCET.object <- function(counts, spotCoordinates, imagePath, platform)
{
  if(!identical(colnames(counts),rownames(spotCoordinates)))
  {
    stop("The Spot IDs in the counts and spotCoordinates matrices are not identical.")
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

  st.matrix.data <- methods::as(counts, "dgCMatrix")
  st.matrix.data <- rm_duplicates(st.matrix.data)

  SpaCET_obj <- methods::new("SpaCET",
    input=list(
      counts=st.matrix.data,
      spotCoordinates=spotCoordinates,
      image=list(path=imagePath,grob=rg),
      platform=platform
    )
  )

  SpaCET_obj
}


#' @title Filter spatial spots and calculate the QC metrics
#' @description Spots with less than `min.genes` expressed genes would be removed.
#' @param SpaCET An SpaCET object.
#' @param min.genes Minimum number of expressed genes. Default: 1.
#' @return An SpaCET object.
#' @examples
#' SpaCET_obj <- SpaCET.quality.control(SpaCET_obj)
#' @rdname SpaCET.quality.control
#' @export
SpaCET.quality.control  <- function(SpaCET_obj, min.genes=1)
{
  st.matrix.data <- SpaCET_obj@input$counts
  expressed.genes<- Matrix::colSums(st.matrix.data>0)

  remaining.spots <- expressed.genes>=min.genes
  remaining.spots.num <- sum(remaining.spots)

  print(paste0(length(remaining.spots)-remaining.spots.num," spots are removed."))
  print(paste0(remaining.spots.num," spots are kept."))

  st.matrix.data <- st.matrix.data[,remaining.spots]

  SpaCET_obj@input$counts <- st.matrix.data
  SpaCET_obj@input$spotCoordinates <- SpaCET_obj@input$spotCoordinates[colnames(st.matrix.data),]

  SpaCET_obj@results$metrics <- rbind(
    UMI=Matrix::colSums(st.matrix.data),
    Gene=Matrix::colSums(st.matrix.data>0)
  )

  SpaCET_obj
}


rm_duplicates <- function(mat){
  dupl <- duplicated(rownames(mat))
  if (sum(dupl) > 0){
    dupl_genes <- unique(rownames(mat)[dupl])
    mat_dupl <- mat[rownames(mat) %in% dupl_genes,,drop=F]
    mat_dupl_names <- rownames(mat_dupl)
    mat <- mat[!dupl,,drop=F]

    for(gene in dupl_genes){
      mat_dupl_gene <- mat_dupl[mat_dupl_names == gene,]
      dupl_sum <- apply(mat_dupl_gene,1,sum)
      max_flag <- which(dupl_sum==max(dupl_sum))
      mat[gene,] <- mat_dupl_gene[max_flag[1],] # in case two values are max
    }
  }
  return(mat)
}

