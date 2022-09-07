#' @title Define the SpaCE object
#' @slot input The input data
#' @slot results The results
#'
setClass("SpaCE",
  slots = c(
    input = "list",
    results = "list"
  )
)


#' @title Create an SpaCE object
#' @description Read an ST dataset to create an SpaCE object.
#' @param counts Count matrix with gene name (row) x spot ID (column)
#' @param spotCoordinates Spot coordinate matrix with spot ID (column) x coordinates (row). This matrix should include two columns, i,e., X and Y coordinates, respectively, which represent the position of spots in H&E image.
#' @param imagePath Path to the H&E image file. Can be NA if it is not available.
#' @param platform A character string indicating the platform, i.e., "Visium", "oldST", or "Slide-Seq".
#' @return An SpaCE object
#' @details
#' To create an SpaCE object, user need to input four parameters, i.e., "counts", "spotCoordinates", "imageFile", and "platform".
#' However, if analyzing Visium data, please use create.SpaCE.object.10X to read data.
#'
#' @examples
#' SpaCE_obj <- create.SpaCE.object(counts,spotCoordinates,imagePath,platform)
#'
#' @rdname create.SpaCE.object
#' @export
#'
create.SpaCE.object <- function(counts,spotCoordinates,imagePath,platform)
{
  library(Matrix)

  st.matrix.data <- as(counts, "dgCMatrix")

  st.matrix.data <- rm_zeroRows(st.matrix.data)
  st.matrix.data <- rm_zeroCols(st.matrix.data)
  st.matrix.data <- rm_duplicates(st.matrix.data)

  olp <- intersect(colnames(st.matrix.data),rownames(spotCoordinates))
  spotCoordinates <- spotCoordinates[olp,]
  st.matrix.data <- st.matrix.data[,olp]

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

  SpaCE_obj <- methods::new("SpaCE",
    input=list(
      counts=st.matrix.data,
      spotCoordinates=spotCoordinates,
      image=list(path=imagePath,grob=rg),
      platform=platform
    ),
    results=list(
      metrics=metrics
    )
  )

  SpaCE_obj
}


#' @title Create an SpaCE object from 10X Visium
#' @description Read an ST dataset to create an SpaCE object.
#' @param visiumPath Path to the Space Ranger output folder. See ‘Details’ for more information.
#' @param resolution A character string indicating the resolution of the H&E image to be used, i.e., "low" or "high".
#' @return An SpaCE object
#' @details
#' If user are analyzing an ST data set from 10X Visium, they only need to input "visiumPath".
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
#' “tissue_hires_image.png” : hematoxylin and eosin (H&E) image; \cr
#' “tissue_lowres_image.png” : hematoxylin and eosin (H&E) image; \cr
#' “scalefactors_json.json” : scaling factors for adjusting the coordinates.
#'
#' @examples
#' visiumPath <- file.path(system.file(package = "SpaCE"), "extdata/Visium_BC")
#' SpaCE_obj <- create.SpaCE.object.10X(visiumPath = visiumPath)
#'
#' @rdname create.SpaCE.object.10X
#' @export
#' @importFrom Matrix readMM
#' @importFrom jsonlite fromJSON
#'
create.SpaCE.object.10X <- function(visiumPath, resolution="low")
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

  SpaCE_obj <- create.SpaCE.object(
    counts=st.matrix.data,
    spotCoordinates=spotCoordinates,
    imagePath=imagePath,
    platform=platform
  )

  SpaCE_obj
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
