#' @title Data Input
#' @description Read ST data set to create an SpaCE object.
#' @param seqPath Path to the folder 'sequencing'.
#' @param imagePath Path to the folder 'image'.
#' @return An SpaCE object
#' @details 
#' Basically, `create.SpaCE` requires two parameters `seqPath` and `imagePath`, which point the standard output folders of 10x Visium data. 
#' 
#' The `seqPath` folder should include \cr
#' "barcodes.tsv.gz": spot level barcodes; \cr
#' "features.tsv.gz": list of genes; \cr
#' "matrix.mtx.gz": (sparse) matrix of counts.
#' 
#' The `imagePath` folder should include \cr
#' “tissue_positions_list.csv” : barcodes and spatial information; \cr
#' “tissue_lowres_image.png” : hematoxylin and eosin (H&E) image; \cr
#' “scalefactors_json.json” : scaling factors for adjusting the coordinates .
#' @examples 
#' 
#' @rdname create.SpaCE
#' @export 
#' @importFrom Matrix readMM
#' @importFrom jsonlite fromJSON
create.SpaCE <- function(seqPath,imagePath)
{
  st.matrix.data <- as.matrix(Matrix::readMM(paste0(seqPath,"/matrix.mtx.gz")))
  st.matrix.gene <- as.matrix(utils::read.csv(paste0(seqPath,"/features.tsv.gz"),as.is=T,header=F,sep="\t"))
  st.matrix.anno <- as.matrix(utils::read.csv(paste0(seqPath,"/barcodes.tsv.gz"),as.is=T,header=F,sep="\t")) 
  
  rownames(st.matrix.data) <- st.matrix.gene[,2]
  colnames(st.matrix.data) <- st.matrix.anno[,1]
  
  jsonFile <- jsonlite::fromJSON(paste0(imagePath,"/scalefactors_json.json"))
  scalef <- jsonFile$tissue_lowres_scalef

  barcode <- utils::read.csv(paste0(imagePath,"/tissue_positions_list.csv"),as.is=T,row.names=1,header=F)
  barcode[["comb"]] <- paste0(barcode[,2],"x",barcode[,3])
  barcode[["comb2"]] <- paste0(barcode[,4]*scalef,"x",barcode[,5]*scalef)
  
  spotID <- barcode[colnames(st.matrix.data),"comb"]
  colnames(st.matrix.data) <- barcode[colnames(st.matrix.data),"comb2"]
  
  methods::new("SpaCE",
      input=list(
        counts=st.matrix.data,
        HEimage=paste0(imagePath,"/tissue_lowres_image.png"),
        spotID=spotID
      ), 
      results=list()
  )
}

