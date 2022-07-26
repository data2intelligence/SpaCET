corMat <- function(X,Y,method="pearson")
{
  olp <- intersect(rownames(X),rownames(Y))
  X_olp <- X[olp,,drop=F]
  Y_olp <- Y[olp,,drop=F]

  cc_corr <- psych::corr.test(X_olp,Y_olp,method=method,adjust="none",ci=FALSE)

  cc_corr_r <- round(cc_corr$r[,1],3)
  cc_corr_p <- signif(cc_corr$p[,1],3)

  cc_corr_rp <- data.frame(
    cor_r=cc_corr_r,
    cor_p=cc_corr_p,
    cor_padj=p.adjust(cc_corr_p, method = "BH")
  )

  cc_corr_rp
}

inferMal_cor <- function(st.matrix.data, cancerType)
{
  print("Stage 1 - Step 1. Clustering.")

  # clustering
  set.seed(123)
  library(MUDAN)

  matnorm.info <- normalizeVariance(methods::as(st.matrix.data, "dgCMatrix"),details=TRUE,verbose=FALSE)
  matnorm <- log10(matnorm.info$mat+1)
  pcs <- getPcs(matnorm[matnorm.info$ods,],nGenes=length(matnorm.info$ods),nPcs=30,verbose=FALSE)

  d <- as.dist(1-cor(t(pcs)))
  hc <- hclust(d, method='ward.D')

  cluster_numbers <- 2:9
  clustering <- cutree(hc,k=cluster_numbers)
  clustering <- t(clustering)
  rownames(clustering) <- paste0("c",rownames(clustering))

  # silhouette
  library(factoextra)
  library(NbClust)
  library(cluster)

  v <- c()
  for(i in cluster_numbers)
  {
    clustering0 <- cutree(hc,k=i)
    sil <- silhouette(clustering0, d, Fun=mean)
    v <- c(v, mean(sil[,3]))
  }

  v_diff <- v[1:(length(v)-1)]-v[2:length(v)]
  maxN <- which( v_diff == max(v_diff) ) +1  # find one that has the biggest decrease

  silMat <- cbind(cluster=cluster_numbers,silhouette=v)
  silMat <- cbind(silMat,maxN=cluster_numbers%in%maxN)

  clustering <- clustering[paste0("c",maxN),] # find optimal k


  print("Stage 1 - Step 2. Find tumor clusters.")

  st.matrix.data <- as.matrix(st.matrix.data)
  st.matrix.data.diff <- t(t(st.matrix.data)*1e6/colSums(st.matrix.data))
  st.matrix.data.diff <- log2(st.matrix.data.diff+1)
  st.matrix.data.diff <- st.matrix.data.diff-rowMeans(st.matrix.data.diff)

  malFlag <- TRUE
  load( system.file("extdata",'cancerDictionary.rda',package = 'SpaCE') )

  for(CNA_expr in c("CNA","expr"))
  {
    cancerTypeExists <- grepl(cancerType,names(cancerDictionary[[CNA_expr]]))

    if(sum(cancerTypeExists) > 0 )
    {
      sig <- as.matrix(cancerDictionary[[CNA_expr]][cancerTypeExists][[1]],ncol=1)
    }else if(CNA_expr=="CNA" & sum(cancerTypeExists)==0 ){
      next
    }else{
      sig <- as.matrix(cancerDictionary[[CNA_expr]]["TCGA_PANCAN"],ncol=1)
    }

    cor_sig <- corMat(as.matrix(st.matrix.data.diff),sig)

    stat.df <- data.frame()
    for(i in sort(unique(clustering)))
    {
      cor_sig_clustering <- cor_sig[clustering==i,]

      stat.df[i,"cluster"] <- i
      stat.df[i,"spotNum"] <- nrow(cor_sig_clustering)
      stat.df[i,"mean"] <- mean(cor_sig_clustering[,1])
      stat.df[i,"fraction_spot_padj"] <- sum(cor_sig_clustering[,"cor_r"]>0&cor_sig_clustering[,"cor_padj"]<0.25)/nrow(cor_sig_clustering)
      stat.df[i,"wilcoxTestG0"] <- wilcox.test(cor_sig_clustering[,1],mu=0,alternative="greater")$ p.value
    }

    stat.df[i+1,"cluster"] <- "All"
    stat.df[i+1,"spotNum"] <- nrow(cor_sig)
    stat.df[i+1,"mean"] <- mean(cor_sig[,1])
    stat.df[i+1,"fraction_spot_padj"] <- sum(cor_sig[,"cor_r"]>0&cor_sig[,"cor_padj"]<0.25)/nrow(cor_sig)
    stat.df[i+1,"wilcoxTestG0"] <- wilcox.test(cor_sig[,1],mu=0,alternative="greater")$ p.value

    stat.df[,"mean"] <- round(stat.df[,"mean"],6)
    stat.df[,"fraction_spot_padj"] <- round(stat.df[,"fraction_spot_padj"],6)

    clusterMal <- which(stat.df[1:i,"mean"]>0&stat.df[1:i,"wilcoxTestG0"]<0.05&stat.df[1:i,"fraction_spot_padj"]>=stat.df[i+1,"fraction_spot_padj"])

    if(length(clusterMal)!=0)
    {
      if(sum(cancerTypeExists) > 0)
      {
        print(paste0("                  > Use cancer-type specific ",CNA_expr," signature."))
      }else{
        print(paste0("                  > Use pan-cancer expr signature."))
      }
      break
    }else{
      if(CNA_expr=="expr")
      {
        print(paste0("                  > No malignant cell detected in this ST data set."))
        malFlag <- FALSE
      }
    }

  }


  print("Stage 1 - Step 3. Infer malignant cell.")

  if(malFlag)
  {
    spotMal <- names(clustering)[clustering%in%clusterMal & cor_sig[,"cor_r"]>0]
    malRef <- rowMeans( t( t(st.matrix.data[,spotMal])*1e6/colSums(st.matrix.data[,spotMal]) ) )

    sig <- apply(st.matrix.data.diff[,spotMal,drop=F],1,mean)
    sig <- matrix(sig)
    rownames(sig) <- rownames(st.matrix.data.diff)

    cor_sig <- corMat(as.matrix(st.matrix.data.diff),sig)

    malProp <- cor_sig[,"cor_r"]
    names(malProp) <- rownames(cor_sig)

    malPropSorted <- sort(malProp)
    top5p <- round(length(malPropSorted)*0.05)
    p5 <- malPropSorted[top5p]
    p95 <- malPropSorted[length(malPropSorted)-top5p+1]

    malProp[malProp<=p5] <- p5
    malProp[malProp>=p95] <- p95

    malProp <- ( malProp-min(malProp) ) / ( max(malProp)-min(malProp) )

    list("malRef"=malRef,"malProp"=malProp)
  }else{
    malProp <- rep(0,dim(st.matrix.data.diff)[2])
    names(malProp) <- colnames(st.matrix.data.diff)

    list("malRef"=NULL,"malProp"=malProp)
  }

}


SpatialDeconv <- function(
    ST,
    Ref,
    malProp,
    malRef,
    mode=c("standard","deconvMal","deconvWithSC"),
    Unidentifiable=TRUE,
    MacrophageOther=TRUE,
    coreNo=8
)
{
  Reference <- Ref$refProfiles
  Signature <- Ref$sigGenes
  Tree <- Ref$lineageTree

  olpGenes <- intersect(rownames(ST), rownames(Reference))

  ST <- ST[olpGenes,]
  Reference <- Reference[olpGenes,]

  ST <- t( t(ST)*1e6/colSums(ST) )
  Reference <- t( t(Reference)*1e6/colSums(Reference) )

  ST <- ST[,!is.nan(ST[1,])]

  if(sum(malProp)>0)
  {
    if(is.matrix(malRef))
    {
      olpGenes <- intersect(rownames(ST), rownames(malRef))

      ST <- ST[olpGenes,]
      malRef <- malRef[olpGenes,]

      ST <- t( t(ST)*1e6/colSums(ST) )
      malRef <- t( t(malRef)*1e6/colSums(malRef) )

      ST <- ST[,!is.nan(ST[1,])]

      mixtureMal <- malRef%*%malProp[-1,] # -1 for minus unidentifiable
    }else{
      malRef <- malRef[rownames(ST)]
      malRef <- malRef*1e6/sum(malRef)
      mixtureMal <- matrix(malRef,ncol=1)%*%matrix(malProp,nrow=1)
    }
    mixtureMinusMal <- ST - mixtureMal
  }else{
    mixtureMinusMal <- ST
  }

  tempReference <- Reference
  tempSignature <- Signature

  print("Stage 2 - Level 1. Estimate the major lineage.")

  ###### level 1 deconv ######
  Level1 <- names(Tree)[names(Tree)%in%colnames(tempReference)]

  if(mode!="deconvMal")
  {
    mixture <- mixtureMinusMal
    Reference <- tempReference[,Level1]
    Signature <- tempSignature[c(Level1,"T cell")]

    nSpot <- dim(mixture)[2]
    nCell <- dim(Reference)[2]
    thetaSum <- (1-malProp)-1e-5

    Signature <- unique(unlist(Signature))
    Signature <- Signature[Signature%in%olpGenes]

    mixture <- mixture[Signature,]
    Reference <- Reference[Signature,]

    propList <- parallel::mclapply(
      1:nSpot,
      FUN=function(i){
        theta <- rep(thetaSum[i]/nCell, nCell)

        if(thetaSum[i]>0.01)
        {
          if(Unidentifiable==TRUE)
          {
            ppmin <- 0
          }else{
            ppmin <- 1-malProp[i]-2e-5
          }
          ppmax <- 1-malProp[i]

          ui <- rbind(diag(nCell), rep(1, nCell), rep(-1, nCell))
          ci <- c(rep(0,nCell), ppmin, -ppmax) #ppmin, ppmax


          f0 <- function(A, x, b){
            sum( (A %*% x - b)^2 )
          }
          res <- stats::constrOptim(theta=theta, f=f0, grad=NULL, ui=ui, ci=ci, A=Reference, b=mixture[,i])
          prop <- res$par
          names(prop) <- colnames(Reference)


          bhat <- Reference %*% prop
          f <- function(A, x, b){
            sum( (A %*% x - b)^2 * ( 1 / ( (bhat +1) ) ) )
          }

          res <- stats::constrOptim(theta=theta, f=f, grad=NULL, ui=ui, ci=ci, A=Reference, b=mixture[,i])
          prop <- res$par
          names(prop) <- colnames(Reference)
        }else{
          prop <- theta
        }

        return(prop)
      },
      mc.cores=coreNo
    )

    propMat <- as.matrix(as.data.frame(propList))
    colnames(propMat) <- colnames(mixture)
    rownames(propMat) <- colnames(Reference)

    if(mode=="standard")
    {
      propMat <- rbind(Malignant=malProp, propMat)
    }

    if(Unidentifiable==TRUE)
    {
      if(mode=="standard")
      {
        propMat <- rbind(propMat, Unidentifiable=1-colSums(propMat))
      }else{
        propMat <- t( t(propMat)/colSums(propMat) )
      }
    }

    propMatLevel1 <- propMat

  }else{
    propMatLevel1 <- matrix(1-colSums(malProp),nrow=1)
    rownames(propMatLevel1) <- "Malignant"
    colnames(propMatLevel1) <- colnames(malProp)
  }

  print("Stage 2 - Level 2. Estimate the sub lineage.")

    ###### level 2 deconv ######
    for(cellSpe in names(Tree)[unlist(lapply(Tree,function(x) length(x)>=2))])
    {
      if(!cellSpe%in%rownames(propMatLevel1)) next

      cellsub <- Tree[[cellSpe]]
      cellsub <- setdiff(cellsub,"Macrophage other")

      if( length(setdiff(Level1,cellSpe)) > 0)
      {
        mixture <- mixtureMinusMal - tempReference[,setdiff(Level1,cellSpe),drop=F] %*% propMatLevel1[setdiff(Level1,cellSpe),,drop=F]
      }else{
        mixture <- mixtureMinusMal
      }
      Reference <- tempReference[,colnames(tempReference)%in%cellsub,drop=F]
      Signature <- tempSignature[names(tempSignature)%in%cellsub] ######

      nSpot <- dim(mixture)[2]
      nCell <- dim(Reference)[2]
      thetaSum <- propMatLevel1[cellSpe,]-1e-5

      Signature <- unique(unlist(Signature))
      Signature <- Signature[Signature%in%olpGenes]

      mixture <- mixture[Signature,]
      Reference <- Reference[Signature,,drop=F]

      ui <- rbind(diag(nCell), rep(1, nCell), rep(-1, nCell))

      propList <- parallel::mclapply(
        1:nSpot,
        FUN=function(i){
          theta <- rep(thetaSum[i]/nCell, nCell)

          if(thetaSum[i]>0.01)
          {
            if(cellSpe=="Macrophage")
            {
              if(MacrophageOther)
              {
                ppmin <- 0
              }else{
                ppmin <- propMatLevel1[cellSpe,i]-2e-5
              }
              ppmax <- propMatLevel1[cellSpe,i]
            }else{
              ppmin <- propMatLevel1[cellSpe,i]-2e-5
              ppmax <- propMatLevel1[cellSpe,i]
            }

            ci <- c(rep(0,nCell), ppmin, -ppmax)

            f0 <- function(A, x, b){
              sum( (A %*% x - b)^2 )
            }
            res <- stats::constrOptim(theta=theta, f=f0, grad=NULL, ui=ui, ci=ci, A=Reference, b=mixture[,i])
            prop <- res$par
            names(prop) <- colnames(Reference)

            bhat <- Reference %*% prop
            f <- function(A, x, b){
              sum( (A %*% x - b)^2 * ( 1 / ( (bhat +1) ) ) )
            }

            res <- stats::constrOptim(theta=theta, f=f, grad=NULL, ui=ui, ci=ci, A=Reference, b=mixture[,i])
            prop <- res$par
            names(prop) <- colnames(Reference)

          }else{
            prop <- theta
          }

          return(prop)
        },
        mc.cores=coreNo
      )

      propMat <- as.matrix(as.data.frame(propList))
      colnames(propMat) <- colnames(mixture)
      rownames(propMat) <- colnames(Reference)

      if(cellSpe=="Macrophage")
      {
        if(MacrophageOther&mode=="standard")
        {
          propMat <- rbind(propMat, "Macrophage other"=propMatLevel1[cellSpe,]-colSums(propMat))
        }
      }

      propMatLevel1 <- rbind(propMatLevel1, propMat)
    }

    propMat <- propMatLevel1

    propMat[propMat<0] <- 0
    propMat[propMat>1] <- 1

    propMat
}


#' @title Deconvolve tumor ST data set
#' @description Estimate the fraction of cell lineage and sub lineage.
#' @param SpaCE_obj An SpaCE object.
#' @param cancerType Cancer type of this tumor ST dataset.
#' @param coreNo Core number.
#' @return An SpaCE object
#' @examples
#' SpaCE_obj <- SpaCE.deconvolution(SpaCE_obj, cancerType="BRCA", coreNo=8)
#' @rdname SpaCE.deconvolution
#' @export
#'
SpaCE.deconvolution <- function(SpaCE_obj,cancerType,coreNo=8)
{
  st.matrix.data <- as.matrix(SpaCE_obj@input$counts)

  print("Stage 1. Infer malignant cell fraction.")
  malRes <- inferMal_cor(st.matrix.data,cancerType)

  load( system.file("extdata",'combRef_0.5.rda',package = 'SpaCE') )

  print("Stage 2. Hierarchically deconvolve non-malignant cell fracton.")

  propMat <- SpatialDeconv(
    ST=st.matrix.data,
    Ref=Ref,
    malProp=malRes$malProp,
    malRef=malRes$malRef,
    mode="standard",
    coreNo=coreNo
  )

  SpaCE_obj@results$Ref <- Ref
  SpaCE_obj@results$deconvolution <- propMat
  SpaCE_obj
}
