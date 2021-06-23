#' @title ST Deconvolution step 1
#' @description Malignant cell fraction inference
#' @param ST PARAM_DESCRIPTION
#' @param Ref PARAM_DESCRIPTION
#' @param infercnv_dir PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname inferMalignant
#' @export 
inferMalignant <- function(ST,Ref,infercnv_dir)
{
  MacrophageName <- c("M1","M2","Macrophage other")
  TCD4Name <- c("Treg","T CD4 naive","T helper","T CD4 other")
  
  
  Signature <- Ref$sigGenes[!names(Ref$sigGenes)%in%c(MacrophageName,TCD4Name,"Malignant")]
  
  
  gmt <- list(normalMarkers=unlist(Signature))
  
  gsvaScore <- GSVA::gsva(log2(ST+1),gmt)
  gsvaScore1 <- gsvaScore[1,]
  
  gsvaScore1sort <- sort(gsvaScore1,decreasing=T)
  top5p <- length(gsvaScore1sort)*0.05
  top5p.spots <- names(gsvaScore1sort)[1:top5p]
  
  Content <- names(gsvaScore1)%in%top5p.spots
  names(Content) <- names(gsvaScore1)
  
  
  Content[Content==TRUE] <- "GSVA score the top 0.05"
  Content[Content=="FALSE"] <- "GSVA score the rest 0.95"
 
  
  # create annotation file for infercnv
  annotations <- cbind(names(Content),Content)
  annotations_file <- paste0(infercnv_dir,"annotations.txt")
  write.table(annotations,annotations_file,quote=F,sep="\t",row.names=F,col.names=F)
  
  # create gene order file for infercnv
  gene_order_file <- system.file("extdata",'gencode_v19_gene_pos.txt', package = 'SpaCE')
  
  # create infercnv object
  infercnv_obj = infercnv::CreateInfercnvObject(
    raw_counts_matrix=ST,
    annotations_file=annotations_file,
    gene_order_file=gene_order_file,
    ref_group_names=c("GSVA score the top 0.05"),
    chr_exclude=c("chrY","chrM")
  )
  
  temp_out_dir <- paste0(infercnv_dir,"temp/")
  dir.create(temp_out_dir)
  
  # run infercnv
  infercnv_obj_default = infercnv::run(
    infercnv_obj,
    cutoff=1,  # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir=temp_out_dir,
    cluster_by_groups=TRUE,
    HMM=TRUE,
    BayesMaxPNormal=0,
    denoise=FALSE,
    num_threads=20,
    png_res=200
  )
  
  system(paste0("cp ",temp_out_dir,"infercnv.png ",infercnv_dir))
  system(paste0("cp ",temp_out_dir,"infercnv.references.txt ",infercnv_dir))
  system(paste0("cp ",temp_out_dir,"infercnv.observations.txt ",infercnv_dir))
  system(paste0("cp ",temp_out_dir,"infercnv.observations_dendrogram.txt ",infercnv_dir))
  system(paste0("cp ",temp_out_dir,"infercnv.observation_groupings.txt ",infercnv_dir))
  

  system(paste0("cp ",temp_out_dir,"infercnv.17_HMM_predHMMi6.hmm_mode-samples.png ",infercnv_dir))
  system(paste0("cp ",temp_out_dir,"infercnv.17_HMM_predHMMi6.hmm_mode-samples.references.txt ",infercnv_dir))
  system(paste0("cp ",temp_out_dir,"infercnv.17_HMM_predHMMi6.hmm_mode-samples.observations.txt ",infercnv_dir))

  
  refNorm <- "infercnv.references.txt"
  obsNorm <- "infercnv.observations.txt"
  
  ref <- read.csv(paste0(infercnv_dir,refNorm), as.is=T, sep=" ", row.names=1)
  obs <- read.csv(paste0(infercnv_dir,obsNorm), as.is=T, sep=" ", row.names=1)
  
  colnames(ref) <- gsub("X","",colnames(ref))
  colnames(obs) <- gsub("X","",colnames(obs))
  
  HMMfile <- "infercnv.17_HMM_predHMMi6.hmm_mode-samples.observations.txt"
  HMM <- read.csv(paste0(infercnv_dir,HMMfile), as.is=T, sep=" ", row.names=1)
  
  cnvRegion <- rownames(HMM)[HMM[,1]!=3]
  
  cnvRegionGain <- rownames(HMM)[HMM[,1]>3]
  cnvRegionLoss <- rownames(HMM)[HMM[,1]<3]
  
  
  allSpots <- cbind(ref[cnvRegion,],obs[cnvRegion,])
  allSpots <- allSpots-1
  
  
  cnv <- colSums(abs(allSpots[cnvRegion,]))
  
  temp <- sort(cnv)
  p5 <- mean(temp[1:top5p])
  p95 <- mean(temp[(length(temp)-top5p+1):length(temp)])
  
  cnv[cnv<=p5] <- p5
  cnv[cnv>=p95] <- p95
  
  malProp <- ( cnv-min(cnv) ) / ( max(cnv)-min(cnv) )
  
  malProp <- malProp[colnames(ST)]
  
  malProp[is.nan(malProp)] <- 0
  
  #system(paste0("rm -rf ",temp_out_dir))
  #system(paste0("rm -rf ",infercnv_dir))
  
  malProp
}


#' @title ST Deconvolution step 2
#' @description Hierarchical deconvolution of non-malignant cell fractions
#' @param ST PARAM_DESCRIPTION
#' @param Ref PARAM_DESCRIPTION
#' @param malInput PARAM_DESCRIPTION, Default: TRUE
#' @param malProp PARAM_DESCRIPTION
#' @param Unidentifiable PARAM_DESCRIPTION, Default: TRUE
#' @param MacrophageOther PARAM_DESCRIPTION, Default: TRUE
#' @param TCD4Other PARAM_DESCRIPTION, Default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[stats]{constrOptim}}
#' @rdname SpatialDeconv
#' @export 
#' @importFrom stats constrOptim
SpatialDeconv <- function(
  ST, 
  Ref, 
  malInput=TRUE, 
  malProp,
  Unidentifiable=TRUE, 
  MacrophageOther=TRUE,
  TCD4Other=TRUE
)
{
  MacrophageName <- c("M1","M2","Macrophage other")
  TCD4Name <- c("Treg","T CD4 naive","T helper","T CD4 other")
  
  
  Reference <- Ref$refProfiles
  Signature <- Ref$sigGenes
  
  olpGenes <- intersect(rownames(ST), rownames(Reference))
  
  ST <- ST[olpGenes,]
  Reference <- Reference[olpGenes,]
  
  ST <- t( t(ST)*1e6/colSums(ST) )
  Reference <- t( t(Reference)*1e6/colSums(Reference) )
  
  tempReference <- Reference
  tempSignature <- Signature
  
  
    ###### level 1 deconv ######
    if(malInput==TRUE)
    {
      Level1 <- setdiff(colnames(tempReference),c(MacrophageName,TCD4Name,"Malignant"))
    }else{
      Level1 <- setdiff(colnames(tempReference),c(MacrophageName,TCD4Name))
    }
    
    mixture <- ST[,!is.nan(ST[1,])]
    Reference <- tempReference[,Level1]
    Signature <- tempSignature[c(Level1,"T cell")]
    #Signature <- tempSignature[c(Level1)]
    nSpot <- dim(mixture)[2]
    nCell <- dim(Reference)[2]
    thetaSum <- (1-malProp)-1e-5
    
    Signature <- unique(unlist(Signature))
    Signature <- Signature[Signature%in%olpGenes]
    
    mixture <- mixture[Signature,]
    Reference <- Reference[Signature,]
    
    f <- function(A, x, b){
      sum( (A %*% x - b)^2 )
    }
    
    
    propList <- lapply(
      1:nSpot, 
      FUN=function(i){
        theta <- rep(thetaSum[i]/nCell, nCell)
        if(thetaSum[i]>0.01)
        {
          if(Unidentifiable)
          {
            ppmin <- 0
          }else{
            ppmin <- 1-malProp[i]-2e-5
          }
          
          ppmax <- 1-malProp[i]
          
          ui <- rbind(diag(nCell), rep(1, nCell), rep(-1, nCell))
          ci <- c(rep(0,nCell), ppmin, -ppmax) #ppmin, ppmax
          res <- stats::constrOptim(theta=theta, f=f, grad=NULL, ui=ui, ci=ci, A=Reference, b=mixture[,i])
          prop <- res$par
          names(prop) <- colnames(Reference)
        }else{
          prop <- theta
        }
        return(prop)
      }
    )
    
    propMat <- as.matrix(as.data.frame(propList))
    colnames(propMat) <- colnames(mixture)	
    rownames(propMat) <- colnames(Reference)	
    
    if(malInput==TRUE)
    {
      propMat <- rbind(Malignant=malProp, propMat)
    }
    
    if(Unidentifiable==TRUE)
    {
      propMat <- rbind(propMat, Unidentifiable=1-colSums(propMat))
    }
    
    propMatLevel1 <- propMat
    
    
    ###### level 2 deconv ######
    
    # Macrophage
    cellSpe <- "Macrophage"
    cellsub <- MacrophageName[1:2]
    
    mixture <- ST[,!is.nan(ST[1,])] - tempReference[,setdiff(Level1,cellSpe),drop=F] %*% propMatLevel1[setdiff(Level1,cellSpe),,drop=F]
    Reference <- tempReference[,colnames(tempReference)%in%cellsub,drop=F]
    Signature <- tempSignature[names(tempSignature)%in%cellsub]
    nSpot <- dim(mixture)[2]
    nCell <- dim(Reference)[2]
    thetaSum <- propMatLevel1[cellSpe,]-1e-5
    
    Signature <- unique(unlist(Signature))
    Signature <- Signature[Signature%in%olpGenes]
    
    mixture <- mixture[Signature,]
    Reference <- Reference[Signature,,drop=F]
    
    f <- function(A, x, b){
      sum( (A %*% x - b)^2 )
    }
    
    ui <- rbind(diag(nCell), rep(1, nCell), rep(-1, nCell))
    
    propList <- lapply(
      1:nSpot, 
      FUN=function(i){
        theta <- rep(thetaSum[i]/nCell, nCell)
        if(thetaSum[i]>0.01)
        {
          if(MacrophageOther)
          {
            ppmin <- 0
          }else{
            ppmin <- propMatLevel1[cellSpe,i]-2e-5
          }
          
          ppmax <- propMatLevel1[cellSpe,i]
          
          ci <- c(rep(0,nCell), ppmin, -ppmax)
          res <- stats::constrOptim(theta=theta, f=f, grad=NULL, ui=ui, ci=ci, A=Reference, b=mixture[,i])
          prop <- res$par
        }else{
          prop <- theta
        }
        return(prop)
      }
    )
    
    propMat <- as.matrix(as.data.frame(propList))
    colnames(propMat) <- colnames(mixture)	
    rownames(propMat) <- colnames(Reference)	
    
    if(MacrophageOther)
    {
      propMat <- rbind(propMat, "Macrophage other"=propMatLevel1[cellSpe,]-colSums(propMat))
    }
    
    propMatLevel2M <- propMat
    
    # T CD4
    cellSpe <- "T CD4"
    cellsub <- TCD4Name[1:3]
    
    mixture <- ST[,!is.nan(ST[1,])] - tempReference[,setdiff(Level1,cellSpe),drop=F] %*% propMatLevel1[setdiff(Level1,cellSpe),,drop=F]
    Reference <- tempReference[,colnames(tempReference)%in%cellsub,drop=F]
    Signature <- tempSignature[names(tempSignature)%in%cellSpe] ######
    nSpot <- dim(mixture)[2]
    nCell <- dim(Reference)[2]
    thetaSum <- propMatLevel1[cellSpe,]-1e-5
    
    Signature <- unique(unlist(Signature))
    Signature <- Signature[Signature%in%olpGenes]
    
    mixture <- mixture[Signature,]
    Reference <- Reference[Signature,,drop=F]
    
    f <- function(A, x, b){
      sum( (A %*% x - b)^2 )
    }
    
    ui <- rbind(diag(nCell), rep(1, nCell), rep(-1, nCell))
    
    propList <- lapply(
      1:nSpot, 
      FUN=function(i){
        theta <- rep(thetaSum[i]/nCell, nCell)
        if(thetaSum[i]>0.01)
        {
          if(TCD4Other)
          {
            ppmin <- 0
          }else{
            ppmin <- propMatLevel1[cellSpe,i]-2e-5
          }
          
          ppmax <- propMatLevel1[cellSpe,i]
          
          ci <- c(rep(0,nCell), ppmin, -ppmax)
          res <- stats::constrOptim(theta=theta, f=f, grad=NULL, ui=ui, ci=ci, A=Reference, b=mixture[,i])
          prop <- res$par
        }else{
          prop <- theta
        }
        return(prop)
      }
    )
    
    propMat <- as.matrix(as.data.frame(propList))
    colnames(propMat) <- colnames(mixture)	
    rownames(propMat) <- colnames(Reference)	
    
    if(TCD4Other)
    {
      propMat <- rbind(propMat, "T CD4 other"=propMatLevel1[cellSpe,]-colSums(propMat))
    }
    
    propMatLevel2T <- propMat
    
    
    propMat <- rbind(propMatLevel1,propMatLevel2M,propMatLevel2T)
    propMat[propMat < 0] <- 0
    
  
  
    propMat
  
}

#' @title ST Deconvolution
#' @description Deconvovle tumor ST data hierarchically
#' @param ST An SpaCE object
#' @param tempFilePath Path to a temporary folder for storing infercnv results.
#' @return An SpaCE object
#' @details This function first estimates cancer cell abundance through modeling the segmental copy number variations across ST spots. Then, a constrained regression model calibrates local tissue densities and determines stromal and immune cell lineage hierarchies.
#' @examples 
#' ST <- deconvolve.ST(ST,tempFilePath)
#' @rdname deconvolve.ST
#' @export 
deconvolve.ST <- function(ST,tempFilePath)
{
  infercnv_dir <- paste0(tempFilePath,"/tempFile_SpaCE/")
  #dir.create(infercnv_dir)
  
  st.matrix.data <- ST@input$counts
  st.matrix.data <- t(t(st.matrix.data)*1e6/colSums(st.matrix.data))
  
  load( system.file("extdata",'Ref.rda',package = 'SpaCE') )
  
  malProp <- inferMalignant(
    st.matrix.data,
    Ref,
    infercnv_dir
  )
  
  propMat <- SpatialDeconv(
    st.matrix.data,
    Ref,
    malInput=TRUE,
    malProp=malProp
  )
  
  ST@results$fraction <- propMat
  
  cnv1 <- read.csv(paste0(infercnv_dir,"/infercnv.references.txt"),as.is=T,sep=" ",quote = "\"")
  cnv2 <- read.csv(paste0(infercnv_dir,"/infercnv.observations.txt"),as.is=T,sep=" ",quote = "\"")
  cnv <- as.matrix(cbind(cnv1,cnv2))
  colnames(cnv) <- gsub("X","",colnames(cnv))
  system(paste0("rm -rf ",infercnv_dir))
  
  ST@results$CNV <- cnv[,colnames(propMat)]
  
  ST
}
