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
#' @param spotCoordinates Spot coordinate matrix with coordinates (row) x spot ID (column). This matrix should include two rows, i,e., X and Y coordinates, respectively, which represent the position of spots in H&E image.
#' @param imageFile Path to the H&E image file. Can be "NA" if it is not available.
#' @param platform A character string indicating the platform, i.e., "Visium", "oldST", or "Slide-Seq".
#' @param visiumPath Path to the visium output folder. See ‘Details’ for more information.
#' @return An SpaCE object
#' @details
#' To create an SpaCE object, user need to input four parameters, i.e., "counts", "spotCoordinates", "imageFile", and "platform".
#'
#' However, if user are analyzing Visium data, they only need to input "visiumPath". Please make sure that "visiumPath" point the standard output folders of 10x Visium, which have both `filtered_feature_bc_matrix` and `spatial` folders.
#'
#' The "filtered_feature_bc_matrix" folder includes \cr
#' "barcodes.tsv.gz": spot level barcodes; \cr
#' "features.tsv.gz": list of genes; \cr
#' "matrix.mtx.gz": (sparse) matrix of counts.
#'
#' The "spatial" folder includes \cr
#' “tissue_positions_list.csv” : barcodes and spatial information; \cr
#' “tissue_lowres_image.png” : hematoxylin and eosin (H&E) image; \cr
#' “scalefactors_json.json” : scaling factors for adjusting the coordinates .
#' @examples
#' visiumPath <- system.file("extdata",'Visium_BC',package = 'SpaCE')
#' SpaCE_obj <- create.SpaCE.object(visiumPath = visiumPath)
#' @rdname create.SpaCE.object
#' @export
#' @importFrom Matrix readMM
#' @importFrom jsonlite fromJSON
#'
create.SpaCE.object <- function(counts,spotCoordinates,imageFile,platform=c("Visium","oldST","Slide-Seq"),visiumPath=NULL)
{
  if(!is.null(visiumPath)) platform<-"Visium"

  library(Matrix)
  if(platform=="Visium")
  {
    st.matrix.data <- readMM(paste0(visiumPath,"/filtered_feature_bc_matrix/matrix.mtx.gz"))
    st.matrix.data <- as(st.matrix.data, "dgCMatrix")

    st.matrix.gene <- as.matrix(read.csv(paste0(visiumPath,"/filtered_feature_bc_matrix/features.tsv.gz"),as.is=T,header=F,sep="\t"))
    st.matrix.anno <- as.matrix(read.csv(paste0(visiumPath,"/filtered_feature_bc_matrix/barcodes.tsv.gz"),as.is=T,header=F,sep="\t"))

    rownames(st.matrix.data) <- st.matrix.gene[,2]
    colnames(st.matrix.data) <- st.matrix.anno[,1]

    st.matrix.data <- rm_duplicates(st.matrix.data)
    st.matrix.data <- rm_zeroRows(st.matrix.data)

    library(jsonlite)
    jsonFile <- fromJSON(paste0(visiumPath,"/spatial/scalefactors_json.json"))
    scalef <- jsonFile$tissue_lowres_scalef

    barcode <- read.csv(paste0(visiumPath,"/spatial/tissue_positions_list.csv"),as.is=T,row.names=1,header=F)
    barcode[["comb"]] <- paste0(barcode[,2],"x",barcode[,3])
    barcode[["comb2"]] <- paste0(round(barcode[,4]*scalef,3),"x",round(barcode[,5]*scalef,3))
    barcode[["X"]] <- round(barcode[,4]*scalef,3)
    barcode[["Y"]] <- round(barcode[,5]*scalef,3)


    spotCoordinates <- barcode[colnames(st.matrix.data),c("X","Y")]
    rownames(spotCoordinates) <- barcode[colnames(st.matrix.data),c("comb")]

    colnames(st.matrix.data) <- rownames(spotCoordinates) # barcode to spot id

    imageFile <- paste0(visiumPath,"/spatial/tissue_lowres_image.png")

    deconvolution <- as.matrix(read.csv(paste0(visiumPath,"/propMat_SpaCE.csv"),row.names=1))
    colnames(deconvolution) <- gsub("X","",colnames(deconvolution))

  }else{
    st.matrix.data <- as(counts, "dgCMatrix")

    st.matrix.data <- rm_duplicates(st.matrix.data)
    st.matrix.data <- rm_zeroRows(st.matrix.data)
  }

  UMI <- colSums(as.matrix(st.matrix.data))
  Gene <- colSums(as.matrix(st.matrix.data)>0)
  metrics <- rbind(UMI=UMI,Gene=Gene)

  SpaCE_obj <- methods::new("SpaCE",
    input=list(
      counts=st.matrix.data,
      spotCoordinates=spotCoordinates,
      imageFile=imageFile,
      platform=platform
    ),
    results=list(
      metrics=metrics
    )
  )
  SpaCE_obj
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


rm_zeroRows <- function(mat){
  mat[rowSums(as.matrix(mat))>0,]
}
