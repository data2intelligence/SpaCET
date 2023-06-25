#' @title Deconvolve malignant cell fraction
#' @description Explore different malignant cell states in tumor ST dataset.
#' @param SpaCET_obj An SpaCET object.
#' @param Malignant Indicates the name of malignant cell type in the major lineage layer from the deconvolution results. Default: "Malignant".
#' @param malignantCutoff Fraction cutoff for defining spots with high abundant malignant cells. Default: 0.7.
#' @param coreNo Core number in parallel.
#' @return An SpaCET object
#' @examples
#' SpaCET_obj <- SpaCET.deconvolution.malignant(SpaCET_obj)
#'
#' @rdname SpaCET.deconvolution.malignant
#' @export
SpaCET.deconvolution.malignant <- function(SpaCET_obj, Malignant="Malignant", malignantCutoff=0.7, coreNo=8)
{
  if(length(SpaCET_obj@results$deconvolution$Ref$lineageTree[[Malignant]])>1)
  {
    stop("Your deconvolution results have included multiple malignant cell states. We do not recommend deconvolve malignant cell fraction further.")
  }

  if(is.null(SpaCET_obj@results$deconvolution$propMat))
  {
    stop("Please do the complete deconvolution first by using SpaCET.deconvolution.")
  }else{
    res_deconv <- SpaCET_obj@results$deconvolution$propMat
  }

  if(length(Malignant)>1 | length(Malignant)==0)
  {
    stop("Please input the only one major malignant cell type.")
  }else{
    if(!Malignant%in%rownames(res_deconv))
    {
      stop("The input malignant cell type does not exist in the deconvolution results. Please check whether you input correct the name of malignant cell type. Of note, R language is case sensitive generally.")
    }
  }

  if(malignantCutoff>1 | malignantCutoff<0)
  {
    stop("Please input a value within 0~1 for the cutoff of malignant spots.")
  }

  coreNoDect <- parallel::detectCores()
  if(coreNoDect<coreNo) coreNo <- coreNoDect


  st.matrix.data <- as.matrix(SpaCET_obj@input$counts)
  st.matrix.data <- st.matrix.data[rowSums(st.matrix.data)>0,]

  st.matrix.data.mal <- st.matrix.data[,res_deconv[Malignant,]>=malignantCutoff]
  st.matrix.data.mal.CPM <- t( t(st.matrix.data.mal)*1e5/colSums(st.matrix.data.mal) )
  st.matrix.data.mal.log <- log2(st.matrix.data.mal.CPM+1)

  # clustering
  set.seed(123)
  suppressPackageStartupMessages(
    library(MUDAN)
  )

  matnorm.info <- normalizeVariance(methods::as(st.matrix.data.mal, "dgCMatrix"),details=TRUE,verbose=FALSE)
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

  maxN <- which( v == max(v) )+1

  silMat <- cbind(cluster=cluster_numbers,silhouette=v)
  silMat <- cbind(silMat,maxN=cluster_numbers%in%maxN)

  clustering <- apply(clustering,1:2,function(x) LETTERS[x])


  Content <- as.character(clustering[maxN-1,])
  names(Content) <- colnames(clustering)

  states <- sort(unique(Content))

  print(paste0("Identify ",length(states)," malignant cell states"))

  refProfiles <- data.frame()
  sigGenes <- list()
  lineageTree <- list()

  refProfiles[rownames(st.matrix.data.mal.CPM),"Malignant"] <- rowMeans(st.matrix.data.mal.CPM)

  for(i in states)
  {
    refProfiles[rownames(st.matrix.data.mal.CPM),paste0("Malignant cell state ",i)] <- rowMeans(st.matrix.data.mal.CPM[,Content==i])

    tempMarkers <- c()
    for(j in setdiff(states,i))
    {
      library(limma)
      TT <- as.numeric(Content%in%c(i))
      WT <- as.numeric(Content%in%c(j))
      design <- cbind(TT,WT)
      fit <- lmFit(st.matrix.data.mal.log,design)
      cont.matrix <- makeContrasts(TTvsWT=TT-WT,levels=design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2)
      res <- topTable(fit2,coef=1,number=nrow(st.matrix.data.mal.log))

      res <- res[order(res[,"t"],decreasing=T),]

      tempMarkers <- c(tempMarkers, rownames(res)[1:500])
    }
    tempMarkers <- table(tempMarkers)

    sigGenes[[paste0("Malignant cell state ",i)]] <- names(tempMarkers)[tempMarkers==1]
  }

  lineageTree[[Malignant]] <- paste0("Malignant cell state ",states)
  Refnew <- list(refProfiles=refProfiles, sigGenes=sigGenes, lineageTree=lineageTree)

  load( system.file("extdata",'combRef_0.5.rda',package = 'SpaCET') )


  knownCellTypes <- names(SpaCET_obj@results$deconvolution$Ref$lineageTree)
  knownCellTypes <- setdiff(knownCellTypes,Malignant)

  if("Unidentifiable"%in%rownames(SpaCET_obj@results$deconvolution$propMat))
  {
    knownCellFractions <- c(knownCellTypes,"Unidentifiable")
  }else{
    knownCellFractions <- knownCellTypes
  }

  propMat <- SpatialDeconv(
    ST=st.matrix.data,
    Ref=Refnew,
    malProp=res_deconv[knownCellFractions,],
    malRef=Ref$refProfiles[,knownCellTypes],
    mode="deconvMal",
    coreNo=coreNo
  )

  propMat<- rbind(res_deconv,propMat[!rownames(propMat)%in%"Malignant",])
  SpaCET_obj@results$deconvolution$propMat <- propMat

  SpaCET_obj
}


#' @title Deconvolve ST data set with matched scRNAseq data
#' @description Estimate the fraction of cell lineage and sub lineage.
#' @param SpaCET_obj An SpaCET object.
#' @param sc_includeMalignant Logical. Indicate whether the single cell data includes malignant cells. If no, please input a cancer type and then SpaCET will predict the malignant cell fraction based on its build-in reference.
#' @param cancerType Cancer type of the current tumor ST dataset.
#' @param sc_counts Single cell count matrix with gene name (row) x cell ID (column).
#' @param sc_annotation Single cell annotation matrix. This matrix should include two columns, i,e., cellID and cellType. Each row represents a single cell.
#' @param sc_lineageTree Cell lineage tree. This should be organized by using a list, and the name of each element are major lineages while the value of elements are the corresponding sublineages. If a major lineage does not have any sublineages, the value of this major lineage should be itself.
#' @param sc_nCellEachLineage Cell count each lineage. Default: 100. If a cell type is comprised of >100 cells, only 100 cells per cell identity are randomly selected to generate cell type reference.
#' @param coreNo Core number.
#' @return An SpaCET object
#' @examples
#' SpaCET_obj <- SpaCET.deconvolution.matched.scRNAseq(SpaCET_obj, sc_counts, sc_annotation, sc_lineageTree)
#'
#' @rdname SpaCET.deconvolution.matched.scRNAseq
#' @export
SpaCET.deconvolution.matched.scRNAseq <- function(SpaCET_obj, sc_includeMalignant=TRUE, cancerType, sc_counts, sc_annotation, sc_lineageTree, sc_nCellEachLineage=100, coreNo=8)
{
  coreNoDect <- parallel::detectCores()
  if(coreNoDect<coreNo) coreNo <- coreNoDect

  if(!identical(rownames(sc_annotation),as.character(sc_annotation[,"cellID"])))
  {
    rownames(sc_annotation) <- as.character(sc_annotation[,"cellID"])
  }

  if(ncol(sc_counts)!=nrow(sc_annotation))
  {
    stop("The cell number in the count and annotation matrix is not identical.")
  }

  if(!identical(sort(colnames(sc_counts)), sort(rownames(sc_annotation))))
  {
    stop("The cell IDs in the count and annotation matrix are not matched.")
  }

  if(length(sc_lineageTree)==0)
  {
    stop("Your lineage tree is empty. Please build a lineage tree as a list.")
  }

  allCellTypes <- unlist(sc_lineageTree)
  allCellTypes_flag <- allCellTypes%in%unique(sc_annotation[,2])
  if(sum(!allCellTypes_flag)>0)
  {
    stop(paste0(
      "These cell types (i.e., ",
      paste0(allCellTypes[!allCellTypes_flag], collapse = ", "),
      ") in the lineage tree do not match with the cell typpes in sc_annotation. Please double-check your lineage tree."
    ))
  }

  set.seed(123)
  idx <- split(sc_annotation[,1], sc_annotation[,2])
  c_keep <- lapply(idx, function(i) {
    n <- length(i)
    if (n > sc_nCellEachLineage)
      n <- sc_nCellEachLineage
    sample(i, n)
  })

  sc_counts <- sc_counts[,unlist(c_keep)]
  sc_annotation <- sc_annotation[unlist(c_keep),]

  sc_counts <- sc_counts[rowSums(sc_counts)>0,]

  print("1. Generate the reference from the matched scRNAseq data.")
  Ref <- generateRef(
    sc.matrix.data = sc_counts,
    sc.matrix.anno = sc_annotation,
    sc.matrix.tree = sc_lineageTree,
    coreNo=coreNo
  )

  print("2. Hierarchically deconvolve the Spatial Transcriptomics dataset.")

  st.matrix.data <- SpaCET_obj@input$counts
  st.matrix.data <- st.matrix.data[Matrix::rowSums(st.matrix.data)>0,]

  if(sc_includeMalignant)
  {
    propMat <- SpatialDeconv(
      ST=st.matrix.data,
      Ref=Ref,
      malProp=rep(0,ncol(st.matrix.data)),
      malRef=NULL,
      mode="deconvWithSC",
      Unidentifiable=TRUE,
      MacrophageOther=FALSE,
      coreNo=coreNo
    )
  }else{
    print("Stage 1. Infer malignant cell fraction.")
    malRes <- inferMal_cor(st.matrix.data,cancerType)

    print("Stage 2. Deconvolve non-malignant cell fracton.")
    propMat <- SpatialDeconv(
      ST=st.matrix.data,
      Ref=Ref,
      malProp=malRes$malProp,
      malRef=malRes$malRef,
      mode="deconvWithSC_alt",
      Unidentifiable=TRUE,
      MacrophageOther=FALSE,
      coreNo=coreNo
    )

  }

  SpaCET_obj@results$deconvolution$Ref <- Ref
  SpaCET_obj@results$deconvolution$propMat <- propMat
  SpaCET_obj
}


generateRef <- function(
    sc.matrix.data = sc.matrix.data,
    sc.matrix.anno = sc.matrix.anno,
    sc.matrix.tree = sc.matrix.tree,
    coreNo = coreNo
)
{
  sc.matrix.data.norm <- t(t(sc.matrix.data)*1e5/colSums(sc.matrix.data))
  sc.matrix.data.log2 <- log2(sc.matrix.data.norm+1)

  cellTypes_level_1 <- names(sc.matrix.tree)
  cellTypes_level_1_toBeSplit <- names(sc.matrix.tree)[sapply(sc.matrix.tree,function(x) length(x)!=1)]

  refProfiles <- data.frame()
  sigGenes <- list()

  for(cellType in cellTypes_level_1)
  {
    refProfiles[rownames(sc.matrix.data.norm),cellType] <- apply(
      sc.matrix.data.norm[,sc.matrix.anno[,"bio_celltype"]%in%sc.matrix.tree[[cellType]],drop=F],
      1,mean)

    run_limma <- function(cellTypeOther)
    {
      library(limma)
      TT <- as.numeric(sc.matrix.anno[,"bio_celltype"]%in%sc.matrix.tree[[cellType]])
      WT <- as.numeric(sc.matrix.anno[,"bio_celltype"]%in%sc.matrix.tree[[cellTypeOther]])
      design <- cbind(TT,WT)
      fit <- lmFit(sc.matrix.data.log2,design)
      cont.matrix <- makeContrasts(TTvsWT=TT-WT,levels=design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2)
      res <- topTable(fit2,coef=1,number=nrow(sc.matrix.data.log2))
      res <- res[order(res[,"t"],decreasing=T),]
      res <- res[1:500,]

      rownames(res)[res[,"logFC"]>0.25&res[,"adj.P.Val"]<0.01]
    }

    markers <- parallel::mclapply(setdiff(cellTypes_level_1,cellType), run_limma, mc.cores=coreNo)
    temp <- table(unlist(markers))

    sigGenes[[cellType]] <- names(temp)[temp>=length(cellTypes_level_1)-1]


    if(cellType%in%cellTypes_level_1_toBeSplit)
    {
      cellTypeSubs <- sc.matrix.tree[[cellType]]

      for(cellTypeSub in cellTypeSubs)
      {
        refProfiles[rownames(sc.matrix.data.norm),cellTypeSub] <- apply(
          sc.matrix.data.norm[,sc.matrix.anno[,"bio_celltype"]%in%cellTypeSub,drop=F],
          1,mean)

        library(limma)
        TT <- as.numeric(sc.matrix.anno[,"bio_celltype"]%in%cellTypeSub)
        WT <- as.numeric(sc.matrix.anno[,"bio_celltype"]%in%setdiff(sc.matrix.tree[[cellType]],cellTypeSub))
        design <- cbind(TT,WT)
        fit <- lmFit(sc.matrix.data.log2,design)
        cont.matrix <- makeContrasts(TTvsWT=TT-WT,levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2)
        res <- topTable(fit2,coef=1,number=nrow(sc.matrix.data.log2))
        res <- res[order(res[,"t"],decreasing=T),]
        res <- res[1:500,]

        markers <- rownames(res)[res[,"logFC"]>0.25&res[,"adj.P.Val"]<0.01]
        sigGenes[[cellTypeSub]] <- markers
      }
    }
  }

  Ref <- list(refProfiles=refProfiles, sigGenes=sigGenes, lineageTree=sc.matrix.tree)
  Ref
}
