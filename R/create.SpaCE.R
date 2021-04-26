#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param seqPath PARAM_DESCRIPTION
#' @param imagePath PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[Matrix]{externalFormats}}
#'  \code{\link[jsonlite]{toJSON, fromJSON}}
#' @rdname create.SpaCE
#' @export 
#' @importFrom Matrix readMM
#' @importFrom jsonlite fromJSON
create.SpaCE <- function(seqPath,imagePath)
{
  st.matrix.data <- as.matrix(Matrix::readMM(paste0(seqPath,"/matrix.mtx.gz")))
  st.matrix.gene <- as.matrix(read.csv(paste0(seqPath,"/features.tsv.gz"),as.is=T,header=F,sep="\t"))
  st.matrix.anno <- as.matrix(read.csv(paste0(seqPath,"/barcodes.tsv.gz"),as.is=T,header=F,sep="\t")) 
  
  rownames(st.matrix.data) <- st.matrix.gene[,2]
  colnames(st.matrix.data) <- st.matrix.anno[,1]
  
  jsonFile <- jsonlite::fromJSON(paste0(imagePath,"/scalefactors_json.json"))
  scalef <- jsonFile$tissue_lowres_scalef

  barcode <- read.csv(paste0(imagePath,"/tissue_positions_list.csv"),as.is=T,row.names=1,header=F)
  barcode[["comb"]] <- paste0(barcode[,4]*scalef,"x",barcode[,5]*scalef)
  
  colnames(st.matrix.data) <- barcode[colnames(st.matrix.data),"comb"]
  
  new("SpaCE",
      input=list(counts=st.matrix.data,HEimage=paste0(imagePath,"/tissue_lowres_image.png")), 
      results=list()
      )
}

