#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ST PARAM_DESCRIPTION
#' @param cutoffMalignant PARAM_DESCRIPTION
#' @param nClone PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ST.malignant.clone.heatmap
#' @export 
ST.malignant.clone.heatmap <- function(ST, cutoffMalignant, nClone)
{
  deconv.res <- ST@results$fraction
  cnv <- log2(ST@results$CNV)
    
  spot_Malignant <- colnames(deconv.res)[deconv.res["Malignant",]>0.5]
  cnv_malignant <- t(as.matrix(cnv[,spot_Malignant]))
  
  gene_order <- read.csv(system.file("extdata",'gencode_v19_gene_pos.txt', package = 'SpaCE'),as.is=T,sep="\t",header=F)
  gene_order_filter <- gene_order[gene_order[,1]%in%rownames(cnv),]
  gene_order_filter[,2] <- gsub("chr","",gene_order_filter[,2]) #1~22,X
  gene_order_filter[,"chr"] <- gsub("chr","",gene_order_filter[,2]) #1~23
  gene_order_filter[,"chr_label"] <- gsub("chr","",gene_order_filter[,2]) #odd,even
  
  gene_order_filter[gene_order_filter[,2]%in%"X","chr"] <- "23"
  
  chr_midPoint <- c()
  for(i in c(1:23))
  {
    if(i%%2==1)
    {
      gene_order_filter[gene_order_filter[,"chr"]%in%as.character(i),"chr_label"] <- "odd"
    }else{
      gene_order_filter[gene_order_filter[,"chr"]%in%as.character(i),"chr_label"] <- "even"
    }
    
    geneInChr <- gene_order_filter[gene_order_filter[,"chr"]%in%as.character(i),1]
    midGeneInChr <- geneInChr[round(length(geneInChr)/2)] 
    mid <- which(gene_order_filter[,1]==midGeneInChr)
    chr_midPoint <- c(chr_midPoint,mid)
  }

  library(ComplexHeatmap)
  ht_opt$message = FALSE
  
  column_ha = columnAnnotation(
    foo = anno_mark(
      at = chr_midPoint, 
      labels = c(1:22,"X"),
      labels_rot=0,
      link_width = unit(3, "mm")
      ),
    chr = gene_order_filter[,"chr_label"],
    col = list(
      foo = "black",
      chr = c("odd" = "grey", "even" = "black")
      ),
    show_legend = FALSE,
    show_annotation_name=TRUE
    )
  
  ht <- Heatmap(
    cnv_malignant, 
    name = "Log2(\ninferred\nCNV\nvalue)", 
    row_km = nClone, 
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    
    clustering_method_rows = "ward.D",
    
    show_row_names = FALSE,
    show_column_names = FALSE,
    
    top_annotation = column_ha
  )
  
  ht = draw(ht)
  ncluster <- row_order(ht)
  

  visiualVector <- deconv.res["Malignant",spot_Malignant]
  visiualVector[!is.na(visiualVector)] <- NA
  
  for(i in names(ncluster))
  {
    spot <- rownames(cnv_malignant)[ncluster[[i]]]
    visiualVector[spot] <- i
  } 

  ST@results$clone <- visiualVector
  
  ST
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ST PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ST.malignant.clone.scatter
#' @export 
ST.malignant.clone.scatter <- function(ST)
{
  visualSpatial(
    ST@results$clone,
    ST@input$HEimage,
    NULL,
    NULL,
    "Malignant cell clone",
    "Clone",
    is.continuous=FALSE
    )

}



