#' @title Deconvolve malignant cell fraction
#' @description Explore different cancer cell state in tumor ST dataset.
#' @param SpaCE_obj An SpaCE object.
#' @param malignantCutoff Fraction cutoff for defining spots with high abundant malignant cells.
#' @param coreNo Core number in parallel.
#' @return An SpaCE object
#' @examples
#' SpaCE_obj <- SpaCE.deconvolution.malignant(SpaCE_obj)
#' @rdname SpaCE.deconvolution.malignant
#' @export
SpaCE.deconvolution.malignant <- function(SpaCE_obj, malignantCutoff=0.7, coreNo=8)
{
  st.matrix.data <- as.matrix(SpaCE_obj@input$counts)
  res_deconv <- SpaCE_obj@results$deconvolution

  st.matrix.data.mal <- st.matrix.data[,res_deconv["Malignant",]>=malignantCutoff]
  st.matrix.data.mal.CPM <- t( t(st.matrix.data.mal)*1e6/colSums(st.matrix.data.mal) )
  st.matrix.data.mal.log <- log2(st.matrix.data.mal.CPM+1)

  # clustering
  set.seed(123)
  library(MUDAN)

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
  suppressPackageStartupMessages(library(factoextra))
  library(NbClust)
  library(cluster)

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

  print(paste0("Identify ",length(states)," cancer cell states"))

  refProfiles <- data.frame()
  sigGenes <- list()

  refProfiles[rownames(st.matrix.data.mal.CPM),"Malignant"] <- rowMeans(st.matrix.data.mal.CPM)

  for(i in states)
  {
    refProfiles[rownames(st.matrix.data.mal.CPM),paste0("Cancer cell state ",i)] <- rowMeans(st.matrix.data.mal.CPM[,Content==i])

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

    sigGenes[[paste0("Cancer cell state ",i)]] <- names(tempMarkers)[tempMarkers==1]
  }

  lineageTree <- list("Malignant"=paste0("Cancer cell state ",states))
  Refnew <- list(refProfiles=refProfiles, sigGenes=sigGenes, lineageTree=lineageTree)

  load( system.file("extdata",'combRef_0.5.rda',package = 'SpaCE') )

  propMat <- SpatialDeconv(
    ST=st.matrix.data,
    Ref=Refnew,
    malProp=res_deconv[2:14,],
    malRef=Ref$refProfiles[,1:12],
    mode="deconvMal",
    coreNo=coreNo
  )

  propMat<- rbind(res_deconv,propMat[!rownames(propMat)%in%"Malignant",])
  SpaCE_obj@results$deconvolution <- propMat

  SpaCE_obj
}


#' @title Deconvolve ST data set with matched scRNAseq
#' @description Estimate the fraction of cell lineage and sub lineage.
#' @param SpaCE_obj An SpaCE object.
#' @param sc_counts Single cell count matrix with gene name (row) x cell ID (column).
#' @param sc_annotation Single cell annotation matrix. This matrix should include two columns, i,e., cellID and cellType. Each row represents a single cell.
#' @param sc_lineageTree Cell lineage tree. This should be organized by using a list, and the name of each element are major lineages while the value of elements are the corresponding sublineages. If a major lineage does not have any sublineages, the value of this major lineage should be itself.
#' @param coreNo Core number.
#' @return An SpaCE object
#' @examples
#' SpaCE_obj <- SpaCE.deconvolution.matched.scRNAseq(SpaCE_obj, sc_counts, sc_annotation, sc_lineageTree)
#' @rdname SpaCE.deconvolution.matched.scRNAseq
#' @export
SpaCE.deconvolution.matched.scRNAseq <- function(SpaCE_obj, sc_counts, sc_annotation, sc_lineageTree, coreNo=8)
{
  st.matrix.data <- as.matrix(SpaCE_obj@input$counts)

  print("1. Generate the reference from the matched scRNAseq data.")
  Ref <- generateRef(
    sc.matrix.data = sc_counts,
    sc.matrix.anno = sc_annotation,
    sc.matrix.tree = sc_lineageTree
  )

  print("2. Hierarchically deconvolve the Spatial Transcriptomics dataset.")

  propMat <- SpatialDeconv(
    ST=st.matrix.data,
    Ref=Ref,
    malProp=rep(0,ncol(st.matrix.data)),
    malRef=NULL,
    mode="deconvWithSC",
    Unidentifiable=TRUE,
    coreNo=coreNo
  )

  SpaCE_obj@results$Ref <- Ref
  SpaCE_obj@results$deconvolution <- propMat
  SpaCE_obj

}

generateRef <- function(
    sc.matrix.data = sc.matrix.data,
    sc.matrix.anno = sc.matrix.anno,
    sc.matrix.tree = sc.matrix.tree
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

    markers <- c()
    for(cellTypeOther in setdiff(cellTypes_level_1,cellType))
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

      temp <- rownames(res)[res[,"logFC"]>0.25&res[,"adj.P.Val"]<0.01]
      markers <- c(markers,temp)
    }
    temp <- table(markers)

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
