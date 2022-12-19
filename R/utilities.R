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


#' @title Create an SpaCET object
#' @description Read an ST dataset to create an SpaCET object.
#' @param counts Count matrix with gene name (row) x spot ID (column)
#' @param spotCoordinates Spot coordinate matrix with spot ID (column) x coordinates (row). This matrix should include two columns, i,e., X and Y coordinates, respectively, which represent the position of spots in H&E image.
#' @param imagePath Path to the H&E image file. Can be NA if it is not available.
#' @param platform A character string indicating the platform, i.e., "Visium", "oldST", or "Slide-Seq".
#' @return An SpaCET object
#' @details
#' To create an SpaCET object, user need to input four parameters, i.e., "counts", "spotCoordinates", "imageFile", and "platform".
#' However, if analyzing Visium data, please use create.SpaCET.object.10X to read data.
#'
#' @examples
#' SpaCET_obj <- create.SpaCET.object(counts, spotCoordinates, imagePath, platform)
#'
#' @rdname create.SpaCET.object
#' @export
#'
create.SpaCET.object <- function(counts,spotCoordinates,imagePath,platform)
{
  library(Matrix)

  st.matrix.data <- as(counts, "dgCMatrix")

  st.matrix.data <- rm_zeroRows(st.matrix.data)
  st.matrix.data <- rm_zeroCols(st.matrix.data)
  st.matrix.data <- rm_duplicates(st.matrix.data)

  edgeMat <- matrix(0,ncol=2,nrow=4)
  rownames(edgeMat) <- c("top","bottom","left","right")
  colnames(edgeMat) <- c("cut_capture_area","cut_tissue_region")

  edgeMat["left","cut_capture_area"]   <- floor( min(spotCoordinates[,1]) )
  edgeMat["right","cut_capture_area"]  <- ceiling( max(spotCoordinates[,1]) )
  edgeMat["bottom","cut_capture_area"] <- floor( min(spotCoordinates[,2]) )
  edgeMat["top","cut_capture_area"]    <- ceiling( max(spotCoordinates[,2]) )


  olp <- intersect(colnames(st.matrix.data),rownames(spotCoordinates))
  spotCoordinates <- spotCoordinates[olp,]
  st.matrix.data <- st.matrix.data[,olp]

  edgeMat["left","cut_tissue_region"]   <- floor( min(spotCoordinates[,1]) )
  edgeMat["right","cut_tissue_region"]  <- ceiling( max(spotCoordinates[,1]) )
  edgeMat["bottom","cut_tissue_region"] <- floor( min(spotCoordinates[,2]) )
  edgeMat["top","cut_tissue_region"]    <- ceiling( max(spotCoordinates[,2]) )


  if(!is.na(imagePath))
  {
    if(platform=="Visium")
    {
      r <- png::readPNG(imagePath)
      rg <- grid::rasterGrob(r, width=grid::unit(1,"npc"), height=grid::unit(1,"npc"))
    }
  }else{
    rg <- NA
  }

  metrics <- rbind(
    UMI=colSums(as.matrix(st.matrix.data)),
    Gene=colSums(as.matrix(st.matrix.data)>0)
  )

  SpaCET_obj <- methods::new("SpaCET",
    input=list(
      counts=st.matrix.data,
      spotCoordinates=spotCoordinates,
      image=list(path=imagePath,grob=rg,edgeMat=edgeMat),
      platform=platform
    ),
    results=list(
      metrics=metrics
    )
  )

  SpaCET_obj
}


#' @title Create an SpaCET object from 10X Visium
#' @description Read an ST dataset to create an SpaCET object.
#' @param visiumPath Path to the SpaCET Ranger output folder. See ‘Details’ for more information.
#' @param resolution A character string indicating the resolution of the H&E image to be used, i.e., "low" or "high".
#' @return An SpaCET object
#' @details
#' If user are analyzing an ST data set from 10X Visium, they only need to input "visiumPath".
#' Please make sure that "visiumPath" points to the standard output folder of 10X SpaCET Ranger,
#' which have both `filtered_feature_bc_matrix` and `spatial` folders.
#'
#' The "filtered_feature_bc_matrix" folder includes \cr
#' "barcodes.tsv.gz": spot level barcodes; \cr
#' "features.tsv.gz": list of genes; \cr
#' "matrix.mtx.gz": (sparse) matrix of counts.
#'
#' The "spatial" folder includes \cr
#' “tissue_positions_list.csv” : barcodes and spatial information; \cr
#' “tissue_hires_image.png” : hematoxylin and eosin (H&E) image; \cr
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
create.SpaCET.object.10X <- function(visiumPath, resolution="low")
{
  platform <- "Visium"

  library(Matrix)

  st.matrix.data <- readMM(paste0(visiumPath,"/filtered_feature_bc_matrix/matrix.mtx.gz")) #dgT
  st.matrix.data <- as(st.matrix.data, "dgCMatrix")

  st.matrix.gene <- as.matrix(read.csv(paste0(visiumPath,"/filtered_feature_bc_matrix/features.tsv.gz"),as.is=T,header=F,sep="\t"))
  st.matrix.anno <- as.matrix(read.csv(paste0(visiumPath,"/filtered_feature_bc_matrix/barcodes.tsv.gz"),as.is=T,header=F,sep="\t"))

  rownames(st.matrix.data) <- st.matrix.gene[,2]
  colnames(st.matrix.data) <- st.matrix.anno[,1]

  library(jsonlite)

  jsonFile <- fromJSON(paste0(visiumPath,"/spatial/scalefactors_json.json"))

  barcode <- read.csv(paste0(visiumPath,"/spatial/tissue_positions_list.csv"),as.is=T,row.names=1,header=F)
  colnames(barcode) <- c("in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres")

  barcode[["pxl_row_in_lowres"]] <- round(barcode[,4]*jsonFile$tissue_lowres_scalef,3)
  barcode[["pxl_col_in_lowres"]] <- round(barcode[,5]*jsonFile$tissue_lowres_scalef,3)
  barcode[["pxl_row_in_hires"]] <- round(barcode[,4]*jsonFile$tissue_hires_scalef,3)
  barcode[["pxl_col_in_hires"]] <- round(barcode[,5]*jsonFile$tissue_hires_scalef,3)

  if(resolution=="low")
  {
    imagePath <- paste0(visiumPath,"/spatial/tissue_lowres_image.png")
    imageRes <- c("pxl_row_in_lowres","pxl_col_in_lowres")
  }else{
    imagePath <- paste0(visiumPath,"/spatial/tissue_hires_image.png")
    imageRes <- c("pxl_row_in_hires","pxl_col_in_hires")
  }

  spotCoordinates <- barcode[,imageRes]

  SpaCET_obj <- create.SpaCET.object(
    counts=st.matrix.data,
    spotCoordinates=spotCoordinates,
    imagePath=imagePath,
    platform=platform
  )

  SpaCET_obj
}


rm_zeroRows <- function(mat){
  mat[rowSums(as.matrix(mat))>0,]
}

rm_zeroCols <- function(mat){
  mat[,colSums(as.matrix(mat))>0]
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
