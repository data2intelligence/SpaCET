#' @title Deconvolve tumor ST data set
#' @description Estimate the cell fraction of cell lineages and sub lineages.
#' @param SpaCET_obj A SpaCET object.
#' @param cancerType Cancer type of the current tumor ST dataset.
#' @param signatureType Indicate the tumor signature type, NULL, CNA or expr. Default: NULL (automatically detect CNA or expr).
#' @param adjacentNormal Indicate whether your sample is normal tissue adjacent to the tumor. If TURE, SpaCET will skip the stage of malignant cell inference. Default: FALSE.
#' @param coreNo Core number in parallel computation.
#' @return A SpaCET object.
#' @examples
#' SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType="BRCA", coreNo=6)
#'
#' @rdname SpaCET.deconvolution
#' @export
#'
SpaCET.deconvolution <- function(SpaCET_obj, cancerType, signatureType=NULL, adjacentNormal=FALSE, coreNo=6)
{
  coreNoDect <- parallel::detectCores(logical = FALSE)
  if(coreNoDect<coreNo)
  {
    message(paste0("Since the number of your physical cores is ",coreNoDect,", coreNo=",coreNoDect," is used automatically."))
    coreNo <- coreNoDect
  }
  if(Sys.info()[['sysname']] == "Windows")
  {
    message("Since Windows does not support > 1 core, coreNo=1 is used automatically.")
    coreNo <- 1
  }

  st.matrix.data <- SpaCET_obj@input$counts
  st.matrix.data <- st.matrix.data[Matrix::rowSums(st.matrix.data)>0,]

  # filter the matrix with ref genes in case the matrix is too big
  load( system.file("extdata",'combRef_0.5.rda',package = 'SpaCET') )

  if(cancerType%in%c("LIHC","CHOL"))
  {
    load( system.file("extdata",'Ref_Normal_LIHC.rda',package = 'SpaCET') )

    olp <- intersect(rownames(Ref$refProfiles), rownames(Ref_Normal$refProfiles))

    Ref$refProfiles <- cbind(Ref$refProfiles[olp,], Ref_Normal$refProfiles[olp,])
    Ref$sigGenes <- append(Ref$sigGenes, Ref_Normal$sigGenes)
    Ref$lineageTree <- append(Ref$lineageTree, Ref_Normal$lineageTree)
  }

  if(ncol(st.matrix.data) > 20000)
  {
    st.matrix.data <- st.matrix.data[rownames(st.matrix.data)%in%rownames(Ref$refProfiles),]
  }

  if(adjacentNormal==TRUE)
  {
    message("Stage 1. Infer malignant cell fraction (skip).")

    malProp <- rep(0,dim(st.matrix.data)[2])
    names(malProp) <- colnames(st.matrix.data)

    malRes <- list("malRef"=NULL,"malProp"=malProp)
  }else{
    message("Stage 1. Infer malignant cell fraction.")
    malRes <- inferMal_cor(
      st.matrix.data,
      cancerType=cancerType,
      signatureType=signatureType
    )
  }

  message(" ")
  message("Stage 2. Hierarchically deconvolve non-malignant cell fraction.")

  if(ncol(st.matrix.data) <= 20000)
  {
    propMat <- SpatialDeconv(
      ST=st.matrix.data,
      Ref=Ref,
      malProp=malRes$malProp,
      malRef=malRes$malRef,
      mode="standard",
      coreNo=coreNo
    )
  }else{
    subNo <- ceiling(ncol(st.matrix.data)/5000)
    for(x in 1:subNo)
    {
      message(paste0("Processing ",x,"/",subNo))

      if(x!=subNo)
      {
        spotSub <- (5000*(x-1)+1):(5000*x)
      }else{
        spotSub <- (5000*(x-1)+1):ncol(st.matrix.data)
      }

      propMatSub <- SpatialDeconv(
        ST=st.matrix.data[,spotSub],
        Ref=Ref,
        malProp=malRes$malProp[spotSub],
        malRef=malRes$malRef,
        mode="standard",
        coreNo=coreNo
      )

      if(x==1)
      {
        propMat <- propMatSub
      }else{
        propMat <- cbind(propMat, propMatSub)
      }
    }
  }

  SpaCET_obj@results$deconvolution$malRes <- malRes
  SpaCET_obj@results$deconvolution$Ref <- Ref
  SpaCET_obj@results$deconvolution$propMat <- propMat
  SpaCET_obj
}


inferMal_cor <- function(st.matrix.data, cancerType, signatureType)
{
  load( system.file("extdata", 'cancerDictionary.rda', package = 'SpaCET') )
  cancerTypes <- unique(c(names(cancerDictionary$CNA),names(cancerDictionary$expr)))
  cancerTypes <- sapply(strsplit(cancerTypes,"_",fixed=T),function(x) return(x[2]))
  if(!cancerType%in%cancerTypes)
  {
    stop("The input cancer type does not match anyone in the build-in dictionary of SpaCET.
         Please make sure you have input the correct cancer type name.
         If yes, it means the dictionary of SpaCET does not include the signature for input cancer type.
         User can set cancerType='PANCAN' to use the pan-cancer expression signature.")
  }

  seq_depth <- Matrix::colSums(st.matrix.data>0)

  #st.matrix.data.diff <- Matrix::t(Matrix::t(st.matrix.data)*1e6/Matrix::colSums(st.matrix.data))
  st.matrix.data.diff <- sweep(st.matrix.data, 2, Matrix::colSums(st.matrix.data), "/") *1e6

  st.matrix.data.diff[is.na(st.matrix.data.diff)] <- 0
  st.matrix.data.diff@x <- log2(st.matrix.data.diff@x+1)
  st.matrix.data.diff <- st.matrix.data.diff-Matrix::rowMeans(st.matrix.data.diff)

  if(ncol(st.matrix.data.diff) < 20000)
  {
    message("Stage 1 - Step 1. Clustering.")

    # clustering
    set.seed(123)
    suppressPackageStartupMessages(
      library(MUDAN)
    )

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
    suppressPackageStartupMessages({
      library(factoextra)
      library(NbClust)
      library(cluster)
    })

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


    message("Stage 1 - Step 2. Find tumor clusters.")

    if(is.null(signatureType))
    {
      if(cancerType=="PANCAN")
      {
        comb_list <- list(c("expr","PANCAN"))
      }else{
        comb_list <- list(c("CNA",cancerType), c("expr",cancerType), c("expr","PANCAN"))
      }

      malFlag <- FALSE
      for(n in 1:length(comb_list))
      {
        CNA_expr <- comb_list[[n]][1]
        cancerType <- comb_list[[n]][2]

        cancerTypeExists <- grepl(cancerType,names(cancerDictionary[[CNA_expr]]))

        if(sum(cancerTypeExists)>0)
        {
          sig <- as.matrix(cancerDictionary[[CNA_expr]][cancerTypeExists][[1]],ncol=1)

          cor_sig <- corMat(as.matrix(st.matrix.data.diff),sig)

          stat.df <- data.frame()
          for(i in sort(unique(clustering)))
          {
            cor_sig_clustering <- cor_sig[clustering==i,]
            seq_depth_clustering <- seq_depth[clustering==i]

            stat.df[i,"cluster"] <- i
            stat.df[i,"spotNum"] <- nrow(cor_sig_clustering)
            stat.df[i,"mean"] <- mean(cor_sig_clustering[,1])
            stat.df[i,"wilcoxTestG0"] <- suppressWarnings(wilcox.test(cor_sig_clustering[,1],mu=0,alternative="greater")$p.value)
            stat.df[i,"fraction_spot_padj"] <- sum(cor_sig_clustering[,"cor_r"]>0&cor_sig_clustering[,"cor_padj"]<0.25)/nrow(cor_sig_clustering)
            stat.df[i,"seq_depth_diff"] <- mean(seq_depth_clustering)-mean(seq_depth)
            stat.df[i,"clusterMal"] <- stat.df[i,"seq_depth_diff"]>0 &
              stat.df[i,"mean"]>0 &
              stat.df[i,"wilcoxTestG0"]<0.05 &
              stat.df[i,"fraction_spot_padj"] >= sum(cor_sig[,"cor_r"]>0&cor_sig[,"cor_padj"]<0.25)/nrow(cor_sig)
          }

          if(sum(stat.df[,"clusterMal"])>0) # find malignant spots.
          {
            message(paste0("                  > Use ",CNA_expr," signature: ",cancerType,"."))
            malFlag <- TRUE
            break
          }
        }

      }

    }else{

      if(signatureType=="CNA")
      {
        for(CNA_expr in c("CNA"))
        {
          cancerTypeExists <- grepl(cancerType,names(cancerDictionary[[CNA_expr]]))

          sig <- as.matrix(cancerDictionary[[CNA_expr]][cancerTypeExists][[1]],ncol=1)

          cor_sig <- corMat(as.matrix(st.matrix.data.diff),sig)

          stat.df <- data.frame()
          for(i in sort(unique(clustering)))
          {
            cor_sig_clustering <- cor_sig[clustering==i,]

            stat.df[i,"cluster"] <- i
            stat.df[i,"spotNum"] <- nrow(cor_sig_clustering)
            stat.df[i,"mean"] <- mean(cor_sig_clustering[,1])
            stat.df[i,"fraction_spot_padj"] <- sum(cor_sig_clustering[,"cor_r"]>0&cor_sig_clustering[,"cor_padj"]<0.25)/nrow(cor_sig_clustering)
            stat.df[i,"wilcoxTestG0"] <- suppressWarnings(wilcox.test(cor_sig_clustering[,1],mu=0,alternative="greater")$ p.value)
            stat.df[i,"clusterMal"] <- stat.df[i,"mean"]>0 &
              stat.df[i,"wilcoxTestG0"]<0.05 &
              stat.df[i,"fraction_spot_padj"] >= sum(cor_sig[,"cor_r"]>0&cor_sig[,"cor_padj"]<0.25)/nrow(cor_sig)
          }

          if(sum(stat.df[,"clusterMal"])>0) # find malignant spots.
          {
            message(paste0("                  > Use ",CNA_expr," signature: ",cancerType,"."))
            malFlag <- TRUE
          }else{
            message(paste0("                  > No malignant cells detected in this tumor ST data set."))
            malFlag <- FALSE
          }
        }
      }

      if(signatureType=="expr")
      {
        for(CNA_expr in c("expr"))
        {
          cancerTypeExists <- grepl(cancerType,names(cancerDictionary[[CNA_expr]]))

          sig <- as.matrix(cancerDictionary[[CNA_expr]][cancerTypeExists][[1]],ncol=1)

          cor_sig <- corMat(as.matrix(st.matrix.data.diff),sig)

          stat.df <- data.frame()
          for(i in sort(unique(clustering)))
          {
            cor_sig_clustering <- cor_sig[clustering==i,]

            stat.df[i,"cluster"] <- i
            stat.df[i,"spotNum"] <- nrow(cor_sig_clustering)
            stat.df[i,"mean"] <- mean(cor_sig_clustering[,1])
            stat.df[i,"fraction_spot_padj"] <- sum(cor_sig_clustering[,"cor_r"]>0&cor_sig_clustering[,"cor_padj"]<0.25)/nrow(cor_sig_clustering)
            stat.df[i,"wilcoxTestG0"] <- suppressWarnings(wilcox.test(cor_sig_clustering[,1],mu=0,alternative="greater")$ p.value)
            stat.df[i,"clusterMal"] <- stat.df[i,"mean"]>0 &
              stat.df[i,"wilcoxTestG0"]<0.05 &
              stat.df[i,"fraction_spot_padj"] >= sum(cor_sig[,"cor_r"]>0&cor_sig[,"cor_padj"]<0.25)/nrow(cor_sig)
          }

          if(sum(stat.df[,"clusterMal"])>0) # find malignant spots.
          {
            message(paste0("                  > Use ",CNA_expr," signature: ",cancerType,"."))
            malFlag <- TRUE
          }else{
            message(paste0("                  > No malignant cells detected in this tumor ST data set."))
            malFlag <- FALSE
          }
        }
      }

    }


    message("Stage 1 - Step 3. Infer malignant cells.")
    top5p <- round(length(seq_depth)*0.05)
    if(malFlag)
    {
      spotMal <- names(clustering)[clustering%in%stat.df[stat.df[,"clusterMal"]==TRUE,"cluster"] & cor_sig[,"cor_r"]>0]
    }else{
      seq_depthSorted <- sort(seq_depth,decreasing = T)
      spotMal <- names(seq_depthSorted)[1:top5p]

      CNA_expr <- "expr";
      cancerType <- "seq_depth"
      stat.df <- NULL
      message(paste0("                  > Use ",CNA_expr," signature: ",cancerType,"."))
    }

    malRef <- Matrix::rowMeans( Matrix::t( Matrix::t(st.matrix.data[,spotMal])*1e6/Matrix::colSums(st.matrix.data[,spotMal]) ) )

    sig <- apply(st.matrix.data.diff[,spotMal,drop=F],1,mean)
    sig <- matrix(sig)
    rownames(sig) <- rownames(st.matrix.data.diff)

    cor_sig <- corMat(as.matrix(st.matrix.data.diff),sig)

    malProp <- cor_sig[,"cor_r"]
    names(malProp) <- rownames(cor_sig)

    malPropSorted <- sort(malProp)
    p5 <- malPropSorted[top5p]
    p95 <- malPropSorted[length(malPropSorted)-top5p+1]

    malProp[malProp<=p5] <- p5
    malProp[malProp>=p95] <- p95

    malProp <- ( malProp-min(malProp) ) / ( max(malProp)-min(malProp) )

    list("sig"=c(CNA_expr, cancerType),"stat.df"=stat.df,"malRef"=malRef,"malProp"=malProp)

  }else{ # spot > 20000
    CNA_expr <- "CNA"

    # first round
    cancerTypeExists <- grepl(cancerType,names(cancerDictionary[[CNA_expr]]))
    sig <- as.matrix(cancerDictionary[[CNA_expr]][cancerTypeExists][[1]],ncol=1)

    subNo <- ceiling(ncol(st.matrix.data.diff)/5000)
    malProp <- c()
    for(x in 1:subNo)
    {
      if(x!=subNo)
      {
        cor_sig <- corMat(as.matrix(st.matrix.data.diff[,(5000*(x-1)+1):(5000*x)]),sig)
      }else{
        cor_sig <- corMat(as.matrix(st.matrix.data.diff[,(5000*(x-1)+1):ncol(st.matrix.data.diff)]),sig)
      }

      malPropSub <- cor_sig[,"cor_r"]
      names(malPropSub) <- rownames(cor_sig)

      malProp <- c(malProp, malPropSub)
    }

    malPropSorted <- sort(malProp)
    top5p <- round(length(malPropSorted)*0.01)
    top5p <- 100
    p5 <- malPropSorted[top5p]
    p95 <- malPropSorted[length(malPropSorted)-top5p+1]

    malProp[malProp<=p5] <- p5
    malProp[malProp>=p95] <- p95

    malProp <- ( malProp-min(malProp) ) / ( max(malProp)-min(malProp) )


    # second round
    spotMal <- colnames(st.matrix.data.diff)[malProp>=1]

    sig <- apply(st.matrix.data.diff[,spotMal,drop=F],1,mean)
    sig <- matrix(sig)
    rownames(sig) <- rownames(st.matrix.data.diff)

    malProp <- c()
    for(x in 1:subNo)
    {
      if(x!=subNo)
      {
        cor_sig <- corMat(as.matrix(st.matrix.data.diff[,(5000*(x-1)+1):(5000*x)]),sig)
      }else{
        cor_sig <- corMat(as.matrix(st.matrix.data.diff[,(5000*(x-1)+1):ncol(st.matrix.data.diff)]),sig)
      }

      malPropSub <- cor_sig[,"cor_r"]
      names(malPropSub) <- rownames(cor_sig)

      malProp <- c(malProp, malPropSub)
    }

    malPropSorted <- sort(malProp)
    top5p <- round(length(malPropSorted)*0.01)
    top5p <- 100
    p5 <- malPropSorted[top5p]
    p95 <- malPropSorted[length(malPropSorted)-top5p+1]

    malProp[malProp<=p5] <- p5
    malProp[malProp>=p95] <- p95

    malProp <- ( malProp-min(malProp) ) / ( max(malProp)-min(malProp) )


    malRef <- Matrix::rowMeans( Matrix::t( Matrix::t(st.matrix.data[,malProp>=1])*1e6/Matrix::colSums(st.matrix.data[,malProp>=1]) ) )

    list("sig"=c(CNA_expr, cancerType),"stat.df"=NULL,"malRef"=malRef,"malProp"=malProp)
  }

}


SpatialDeconv <- function(
    ST,
    Ref,
    malProp,
    malRef,
    mode=c("standard","deconvMal","deconvWithSC","deconvWithSC_alt"),
    Unidentifiable=TRUE,
    MacrophageOther=TRUE,
    coreNo
)
{
  Reference <- Ref$refProfiles
  Signature <- Ref$sigGenes
  Tree <- Ref$lineageTree

  olpGenes <- intersect(rownames(ST), rownames(Reference))

  ST <- ST[olpGenes,]
  Reference <- Reference[olpGenes,]

  ST <- Matrix::t( Matrix::t(ST)*1e6/Matrix::colSums(ST) )
  Reference <- t( t(Reference)*1e6/colSums(Reference) )

  ST <- ST[,!is.nan(ST[1,])]

  if(sum(malProp)>0)
  {
    if(is.matrix(malRef)|is.data.frame(malRef))
    {
      olpGenes <- intersect(rownames(ST), rownames(malRef))

      ST <- ST[olpGenes,]
      malRef <- malRef[olpGenes,]

      ST <- Matrix::t( Matrix::t(ST)*1e6/Matrix::colSums(ST) )
      malRef <- t( t(malRef)*1e6/colSums(malRef) )

      ST <- ST[,!is.nan(ST[1,])]

      mixtureMal <- malRef%*%malProp[colnames(malRef),] # -1 for minus unidentifiable
    }else{
      malRef <- malRef[rownames(ST)]
      malRef <- malRef*1e6/sum(malRef)
      mixtureMal <- matrix(malRef,ncol=1)%*%matrix(malProp,nrow=1)
      colnames(mixtureMal) <- names(malProp)
    }

    olpSpots <- intersect(colnames(ST), colnames(mixtureMal))
    ST <- ST[,olpSpots]
    mixtureMal <- mixtureMal[,olpSpots]

    mixtureMinusMal <- ST - mixtureMal
  }else{
    mixtureMinusMal <- ST
  }

  tempReference <- Reference
  tempSignature <- Signature


  ###### level 1 deconv ######
  Level1 <- names(Tree)[names(Tree)%in%colnames(tempReference)]

  if(mode!="deconvMal")
  {
    message("Stage 2 - Level 1. Estimate the major lineage.")

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

    propList <- pbmcapply::pbmclapply(
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

    if(mode%in%c("standard","deconvWithSC_alt"))
    {
      propMat <- rbind(Malignant=malProp[colnames(propMat)], propMat)
    }

    if(Unidentifiable==TRUE)
    {
      if(mode=="standard")
      {
        propMat <- rbind(propMat, Unidentifiable=1-colSums(propMat))
      }else if(mode=="deconvWithSC_alt"){
        propMat_nonMal <- propMat[-1,]
        propMat_nonMal <- t( t(propMat_nonMal)/colSums(propMat_nonMal) )
        propMat_nonMal <- t( t(propMat_nonMal)*(1-propMat[1,]) )
        propMat[-1,] <- propMat_nonMal
      }else{ # deconvWithSC
        propMat <- t( t(propMat)/colSums(propMat) )
      }
    }

    propMatLevel1 <- propMat

  }else{
    propMatLevel1 <- matrix(1-colSums(malProp),nrow=1)
    rownames(propMatLevel1) <- Level1
    colnames(propMatLevel1) <- colnames(malProp)
  }

  if(mode!="deconvMal")
  {
    message("Stage 2 - Level 2. Estimate the sub lineage.")
  }

  ###### level 2 deconv ######
  for(cellSpe in names(Tree)[unlist(lapply(Tree,function(x) length(x)>=2))])
  {
      if(!cellSpe%in%rownames(propMatLevel1)) next

      message(paste0("                  > ",cellSpe,":"))

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

      propList <- pbmcapply::pbmclapply(
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

      if(mode=="standard"&MacrophageOther&cellSpe=="Macrophage")
      {
        propMat <- rbind(propMat, "Macrophage other"=propMatLevel1[cellSpe,]-colSums(propMat))
      }

      propMatLevel1 <- rbind(propMatLevel1, propMat)
  }

  propMat <- propMatLevel1

  propMat[propMat<0] <- 0
  propMat[propMat>1] <- 1

  propMat
}


corMat <- function(X,Y,method="pearson")
{
  olp <- intersect(rownames(X),rownames(Y))

  cc_corr <- psych::corr.test(
    X[olp,,drop=F],
    Y[olp,,drop=F],
    method=method,adjust="none",ci=FALSE)

  cc_corr_r <- round(cc_corr$r[,1],3)
  cc_corr_p <- signif(cc_corr$p[,1],3)

  cc_corr_rp <- data.frame(
    cor_r=cc_corr_r,
    cor_p=cc_corr_p,
    cor_padj=p.adjust(cc_corr_p, method = "BH")
  )

  cc_corr_rp
}
