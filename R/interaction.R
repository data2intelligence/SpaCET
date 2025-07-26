#' @title Cell-cell colocalization analysis
#' @description Calculate the cell-type pair colocalization by using Spearman correlation analysis.
#' @param SpaCET_obj A SpaCET object.
#' @return A SpaCET object.
#' @examples
#' SpaCET_obj <- SpaCET.CCI.colocalization(SpaCET_obj)
#'
#' @rdname SpaCET.CCI.colocalization
#' @export
#'
SpaCET.CCI.colocalization <- function(SpaCET_obj)
{
  res_deconv <- SpaCET_obj@results$deconvolution$propMat
  res_deconv <- res_deconv[!rownames(res_deconv)%in%c("Unidentifiable","Macrophage other"),]
  res_deconv <- round(res_deconv,2)
  overallFraction <- rowMeans(res_deconv)

  cc_corr <- psych::corr.test(t(res_deconv),t(res_deconv),method="spearman",adjust="none",ci=FALSE)

  cc_corr_r <- cc_corr$r
  cc_corr_p <- cc_corr$p

  cc_corr_r.m <- reshape2::melt(cc_corr_r)
  cc_corr_p.m <- reshape2::melt(cc_corr_p)

  summary_df <- data.frame(
    cell_type_1 = as.character(cc_corr_r.m[,1]),
    cell_type_2 = as.character(cc_corr_r.m[,2]),
    fraction_product = overallFraction[as.character(cc_corr_r.m[,1])]*overallFraction[as.character(cc_corr_r.m[,2])],
    fraction_rho = round(cc_corr_r.m[,3],3),
    fraction_pv = cc_corr_p.m[,3]
    )
  rownames(summary_df) <- paste0(summary_df[,1],"_",summary_df[,2])
  summary_df[is.na(summary_df[,"fraction_rho"] ),"fraction_rho"] <- 0
  summary_df[is.na(summary_df[,"fraction_pv"] ),"fraction_pv"] <- 1

  Ref <- SpaCET_obj@results$deconvolution$Ref
  reff <- Ref$refProfiles[unique(unlist(Ref$sigGenes[names(Ref$sigGenes)%in%c(names(Ref$lineageTree),"T cell")])),]
  reff <- reff-rowMeans(reff)

  cc_corr <- psych::corr.test(reff,reff,method="spearman",adjust="none",ci=FALSE)

  cc_corr_r <- cc_corr$r
  cc_corr_p <- cc_corr$p

  cc_corr_r.m <- reshape2::melt(cc_corr_r)
  cc_corr_p.m <- reshape2::melt(cc_corr_p)

  summary_df2 <- data.frame(
    cell_type_1 = as.character(cc_corr_r.m[,1]),
    cell_type_2 = as.character(cc_corr_r.m[,2]),
    reference_rho = round(cc_corr_r.m[,3],3),
    reference_pv = cc_corr_p.m[,3]
  )
  rownames(summary_df2) <- paste0(summary_df2[,1],"_",summary_df2[,2])

  summary_df[rownames(summary_df2),"reference_rho"] <- summary_df2[,"reference_rho"]
  summary_df[rownames(summary_df2),"reference_pv"] <- summary_df2[,"reference_pv"]

  summary_df <- summary_df[ summary_df[,1] != summary_df[,2], ] #remove same cell type

  SpaCET_obj@results$CCI$colocalization <- summary_df

  SpaCET_obj
}


#' @title Cell-cell colocalization visualization
#' @description Visualize the cell-type pair colocalization.
#' @param SpaCET_obj A SpaCET object.
#' @return A ggplot object.
#' @examples
#' SpaCET.visualize.colocalization(SpaCET_obj)
#'
#' @rdname SpaCET.visualize.colocalization
#' @export
#'
SpaCET.visualize.colocalization <- function(SpaCET_obj)
{
  summary_df <- SpaCET_obj@results$CCI$colocalization

  summary_df[summary_df[,"fraction_product"]>0.02,"fraction_product"] <- 0.02

  if("Malignant cell state A"%in%rownames(SpaCET_obj@results$deconvolution$propMat) & "Unidentifiable"%in%rownames(SpaCET_obj@results$deconvolution$propMat))
  {
    states <- rownames(SpaCET_obj@results$deconvolution$propMat)[grepl("Malignant cell state",rownames(SpaCET_obj@results$deconvolution$propMat))]
    ctOrder <- c(states,unlist(SpaCET_obj@results$deconvolution$Ref$lineageTree))
  }else{
    if("Unidentifiable"%in%rownames(SpaCET_obj@results$deconvolution$propMat))
    {
      ctOrder <- c("Malignant",unlist(SpaCET_obj@results$deconvolution$Ref$lineageTree))
    }else{
      ctOrder <- unlist(SpaCET_obj@results$deconvolution$Ref$lineageTree)
    }
  }
  ctOrder <- unique(ctOrder)

  summary_df <- summary_df[summary_df[,1]%in%ctOrder,]
  summary_df <- summary_df[summary_df[,2]%in%ctOrder,]

  summary_df[,1] <- factor(summary_df[,1],levels=rev(ctOrder))
  summary_df[,2] <- factor(summary_df[,2],levels=ctOrder)

  p1 <- ggplot(summary_df, aes( cell_type_2, cell_type_1)) +
    geom_point(aes(colour = fraction_rho, size=fraction_product), na.rm = TRUE) +
    scale_colour_gradient2(low = "blue", high = "red", mid = "white", na.value = NA,
                           midpoint = 0, limit = c(-0.6,0.6), space = "Lab",
                           name="Rho",oob = scales::squish)+
    scale_size(range = c(0, 6))+
    ggtitle("Cell-cell colocalization")+
    xlab("Cell lineages")+
    ylab("Cell lineages")+
    labs(size = "Fraction\nproduct")+
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white",colour = "white", size = 0.5, linetype = "solid"),
      strip.text = element_blank(),
      axis.text.x = element_text(size = 11, colour = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 11, colour = "black")
    )


  summary_df <- summary_df[as.character(summary_df[,1])<as.character(summary_df[,2]),]

  summary_df <- summary_df[as.character(summary_df[,1])%in%unlist(SpaCET_obj@results$deconvolution$Ref$lineageTree),]
  summary_df <- summary_df[as.character(summary_df[,2])%in%unlist(SpaCET_obj@results$deconvolution$Ref$lineageTree),]

  summary_df <- cbind(summary_df,label=paste0(summary_df[,1],"_",summary_df[,2]))
  summary_df[summary_df[,"fraction_product"]<0.0005 ,"label"] <- ""
  summary_df[summary_df[,"fraction_rho"]<0.1 ,"label"] <- ""

  summary_df[is.na(summary_df[,"fraction_rho"]),c("fraction_rho")] <- 0

  res <- lm(summary_df$fraction_rho ~ summary_df$reference_r)
  summary_df <- cbind(summary_df,Residual=abs(residuals(res)))

  res_weight <- lm(summary_df$fraction_rho ~ summary_df$reference_r, weight=summary_df$fraction)

  p2 <- ggplot(summary_df, aes( reference_rho, fraction_rho, label=label )) +
    geom_point(aes(size=fraction_product),color="#856aad")+
    geom_smooth(method="lm",formula="y~x",color="darkgrey",mapping = aes(weight = fraction_product))+
    scale_size(range = c(0, 6))+
    ggrepel::geom_text_repel(max.overlaps=Inf)+
    ggtitle("Correlation of cell fractions and cell reference profiles")+
    xlab("Cor (Reference profiles)")+
    ylab(" \nCor (Cell fractions)")+
    labs(size = "Fraction\nProduct")+
    ylim(min(summary_df[,"fraction_rho"])-0.02,max(summary_df[,"fraction_rho"])+0.02)+
    xlim(min(summary_df[,"reference_rho"])-0.02,max(summary_df[,"reference_rho"])+0.02)+
    theme_bw()+
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.background = element_blank(),
      panel.grid = element_blank(),
      strip.text = element_blank(),
      axis.text = element_text(size = 11, colour = "black"),
      axis.title = element_text(size = 12, colour = "black")
    )

  library(patchwork)
  p1+p2
}


#' @title Ligand-Receptor interaction enrichment analysis
#' @description Calculate the overall intensity of L-R interactions within ST spots.
#' @param SpaCET_obj A SpaCET object.
#' @param coreNo Core number in parallel computation.
#' @return A SpaCET object.
#' @examples
#' SpaCET_obj <- SpaCET.CCI.LRNetworkScore(SpaCET_obj, coreNo=6)
#' SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "LRNetworkScore", spatialFeatures=c("Network_Score","Network_Score_pv"))
#'
#' @rdname SpaCET.CCI.LRNetworkScore
#' @export
#'
SpaCET.CCI.LRNetworkScore <- function(SpaCET_obj, coreNo=6)
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

  st.matrix.data <- Matrix::t(Matrix::t(st.matrix.data)*1e6/Matrix::colSums(st.matrix.data))
  st.matrix.data@x <- log2(st.matrix.data@x+1)

  spotNum <- ncol(st.matrix.data)

  ###### preprocess L-R network ######
  LRdb <- read.csv(system.file("extdata",'Ramilowski2015.txt',package = 'SpaCET'),as.is=T,sep="\t")
  LRdb <- data.frame(L=LRdb[,2],R=LRdb[,4],stringsAsFactors=FALSE)

  Ls <- unique(LRdb[,"L"])
  Rs <- unique(LRdb[,"R"])

  # intersect(Ls,Rs)
  # [1] "DLK2"
  LRdb <- LRdb[LRdb[,"R"]!="DLK2",]

  ###### filter
  LRdb <- LRdb[LRdb[,1]%in%rownames(st.matrix.data)&LRdb[,2]%in%rownames(st.matrix.data),]
  rownames(LRdb) <- paste0(LRdb[,1],"_",LRdb[,2])

  ###### transform to matrix
  Ls <- unique(LRdb[,"L"])
  Rs <- unique(LRdb[,"R"])

  LRdb_mat <- matrix(0,nrow=length(Ls),ncol=length(Rs))
  rownames(LRdb_mat) <- Ls
  colnames(LRdb_mat) <- Rs

  for(i in 1:dim(LRdb)[1])
  {
    LRdb_mat[LRdb[i,1],LRdb[i,2]] <- 1
  }

  ###### permute L-R network
  message("Step 1. Permute Ligand-Receptor network.")

  LRdb_bg <- function(i)
  {
    LRdb_rand <- BiRewire::birewire.rewire.bipartite(LRdb_mat,verbose=FALSE)

    LRdb_rand.m <- as.matrix(reshape2::melt(LRdb_rand))
    LRdb_rand.m <- LRdb_rand.m[LRdb_rand.m[,3]==1,1:2]
    colnames(LRdb_rand.m) <- c("L","R")

    c(t(LRdb_rand.m))
  }

  set.seed(123456)
  LRdb_bg_List <- pbmcapply::pbmclapply(1:1000,LRdb_bg,mc.cores=coreNo)
  LRdb_rand_comb <- matrix(unlist(LRdb_bg_List),ncol=2,byrow=TRUE)


  ###### Calculate L-R NS
  message("Step 2. Calculate L-R network score.")
  LRcoexpr <- function(LRdb,spot)
  {
    Ls <- LRdb[,1]
    Rs <- LRdb[,2]
    mean( spot[Ls] * spot[Rs] )
  }

  LRcoexpr1000 <- function(LRdb,spot)
  {
    Ls <- LRdb[,1]
    Rs <- LRdb[,2]
    colMeans(matrix(spot[Ls] * spot[Rs], ncol=1000))
  }

  LRNetworkScore <- function(i)
  {
    spot <- st.matrix.data[,i]

    LR_raw <- LRcoexpr(LRdb,spot=spot)

    LR1000 <- LRcoexpr1000(LRdb_rand_comb,spot=spot)

    score <- LR_raw/mean(LR1000)
    pv <- (sum(LR1000>=LR_raw)+1)/1001

    list(LR_raw,score,pv)
  }

  LRNetworkScoreList <- pbmcapply::pbmclapply(1:spotNum,LRNetworkScore,mc.cores=coreNo)

  LRNetworkScoreMat <- matrix(unlist(LRNetworkScoreList),ncol=ncol(st.matrix.data))
  colnames(LRNetworkScoreMat) <- colnames(st.matrix.data)
  rownames(LRNetworkScoreMat) <- c("Raw_expr","Network_Score","Network_Score_pv")
  ###### Calculate L-R network score ######

  SpaCET_obj@results$CCI$LRNetworkScore <- LRNetworkScoreMat
  SpaCET_obj
}


#' @title Ligand-Receptor analysis for a co-localized cell-type pair
#' @description Test the co-expression of ligands and receptors within the same ST spot for the co-localized cell-type pairs.
#' @param SpaCET_obj A SpaCET object.
#' @param cellTypePair A pair of cell-types.
#' @return A SpaCET object.
#' @examples
#' SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("CAF","Macrophage M2"))
#'
#' @rdname SpaCET.CCI.cellTypePair
#' @export
#'
SpaCET.CCI.cellTypePair <- function(SpaCET_obj, cellTypePair)
{
  if(length(cellTypePair)!=2)
  {
    stop("Please input a pair of cell-types.")
  }

  if(sum(cellTypePair%in%rownames(SpaCET_obj@results$deconvolution$propMat))!=2)
  {
    stop("Please input the correct cell-type name. Of note, R language is case sensitive generally.")
  }

  cellTypePair <- sort(cellTypePair)

  res_deconv <- SpaCET_obj@results$deconvolution$propMat
  LRNetworkScoreMat <- SpaCET_obj@results$CCI$LRNetworkScore

  if(is.null(SpaCET_obj@results$CCI$interaction$groupMat))
  {
    groupMat <- data.frame()
    testRes <- data.frame()
  }else{
    groupMat <- SpaCET_obj@results$CCI$interaction$groupMat
    testRes <- SpaCET_obj@results$CCI$interaction$testRes
  }

  summary_df <- SpaCET_obj@results$CCI$colocalization
  rho <- summary_df[summary_df[,1]==cellTypePair[1]&summary_df[,2]==cellTypePair[2],"fraction_rho"]
  pv1 <-  summary_df[summary_df[,1]==cellTypePair[1]&summary_df[,2]==cellTypePair[2],"fraction_pv"]

  testRes[paste0(cellTypePair[1],"_",cellTypePair[2]),"colocalization_rho"] <- rho
  testRes[paste0(cellTypePair[1],"_",cellTypePair[2]),"colocalization_pv"] <- pv1


  cutoff1 <- summary(res_deconv[cellTypePair[1],])[5]
  cutoff2 <- summary(res_deconv[cellTypePair[2],])[5]

  cutoff11 <- quantile(res_deconv[cellTypePair[1],],0.85)
  cutoff22 <- quantile(res_deconv[cellTypePair[2],],0.85)

  Content <- res_deconv[1,]
  for(i in 1:length(Content))
  {
    if(res_deconv[cellTypePair[1],i]>cutoff11 & res_deconv[cellTypePair[2],i]>cutoff22)
    {
      Content[i] <- "Both"
    }else if(res_deconv[cellTypePair[1],i]>cutoff11 & res_deconv[cellTypePair[2],i]<cutoff2){
      Content[i] <- cellTypePair[1]
    }else if(res_deconv[cellTypePair[2],i]>cutoff22 & res_deconv[cellTypePair[1],i]<cutoff1){
      Content[i] <- cellTypePair[2]
    }else{
      Content[i] <- NA
    }
  }

  groupMat[paste0(cellTypePair[1],"_",cellTypePair[2]),colnames(res_deconv)] <- Content

  if(sum(Content%in%"Both")>5)
  {
    fg.df <- data.frame(group=Content,value=LRNetworkScoreMat[2,],stringsAsFactors=FALSE)
    fg.df <- fg.df[!fg.df[,"group"]%in%NA,]
    fg.df[fg.df[,1]%in%cellTypePair,1] <- "Single"

    cohend_res <- psych::cohen.d(fg.df, group="group", alpha=.05, std=TRUE)
    cd1 <- signif(cohend_res$cohen.d["value","effect"],2)

    n1 <- sum(fg.df[,1]%in%"Both")
    n2 <- sum(!fg.df[,1]%in%"Both")
    pv2 <- signif(wilcox.test(fg.df[fg.df[,1]%in%"Both",2],fg.df[!fg.df[,1]%in%"Both",2])$p.value,2)

    testRes[paste0(cellTypePair[1],"_",cellTypePair[2]),"groupCompare_cohen.d"] <- cd1
    testRes[paste0(cellTypePair[1],"_",cellTypePair[2]),"groupCompare_pv"] <- pv2

    if(rho > 0 & pv1 < 0.05 & cd1 < 0 & pv2 < 0.05)
    {
      message(paste0("Based on colocalization analysis and L-R enrichment analysis, ",cellTypePair[1]," and ",cellTypePair[2], " have potential intercellular interaction in the current tissue."))
      testRes[paste0(cellTypePair[1],"_",cellTypePair[2]),"Interaction"] <- TRUE
     }else{
      message("Based on colocalization analysis and L-R enrichment analysis, the intercellular interaction is not significant for the current cell-type pair. Please check other cell-type pairs.")
      testRes[paste0(cellTypePair[1],"_",cellTypePair[2]),"Interaction"] <- FALSE
     }
  }else{
    message("The colocalization analysis is not significant for the current cell-type pair. Please check other cell-type pairs.")
    testRes[paste0(cellTypePair[1],"_",cellTypePair[2]),"Interaction"] <- FALSE
  }

  SpaCET_obj@results$CCI$interaction$groupMat <- groupMat
  SpaCET_obj@results$CCI$interaction$testRes <- testRes

  SpaCET_obj
}


#' @title Cell-type pair visualization
#' @description Visualize the interaction analysis of a co-localized cell-type pair.
#' @param SpaCET_obj A SpaCET object.
#' @param cellTypePair A pair of cell types.
#' @return A ggplot object.
#' @examples
#' SpaCET.visualize.cellTypePair(SpaCET_obj, cellTypePair=c("CAF","Macrophage M2"))
#'
#' @rdname SpaCET.visualize.cellTypePair
#' @export
#'
SpaCET.visualize.cellTypePair <- function(SpaCET_obj, cellTypePair)
{
  if(length(cellTypePair)!=2)
  {
    stop("Please input a pair of cell-types.")
  }

  if(sum(cellTypePair%in%rownames(SpaCET_obj@results$deconvolution$propMat))!=2)
  {
    stop("Please input the correct cell-type name. Of note, R language is case sensitive generally.")
  }
  cellTypePair <- sort(cellTypePair)


  if(!paste0(cellTypePair[1],"_",cellTypePair[2])%in%rownames(SpaCET_obj@results$CCI$interaction$testRes))
  {
    stop("Please run SpaCET.CCI.cellTypePair first.")
  }else{
    res_deconv <- SpaCET_obj@results$deconvolution$propMat
    groupMat <- SpaCET_obj@results$CCI$interaction$groupMat
    testRes <- SpaCET_obj@results$CCI$interaction$testRes

    rho <- testRes[paste0(cellTypePair[1],"_",cellTypePair[2]),"colocalization_rho"]
    pv1 <- testRes[paste0(cellTypePair[1],"_",cellTypePair[2]),"colocalization_pv"]
    pv1 <- signif(as.numeric(pv1),3)
    pv1 <-  ifelse(pv1==0," < 2.2e-16",paste0(" = ",pv1))

    cd1 <- testRes[paste0(cellTypePair[1],"_",cellTypePair[2]),"groupCompare_cohen.d"]
    pv2 <- testRes[paste0(cellTypePair[1],"_",cellTypePair[2]),"groupCompare_pv"]

    Content <- unlist(groupMat[paste0(cellTypePair[1],"_",cellTypePair[2]),colnames(res_deconv)])
    visiualVector <- Content
    spotID <- names(visiualVector)
    names(visiualVector) <- paste0(SpaCET_obj@input$spotCoordinates[,1],"x",SpaCET_obj@input$spotCoordinates[,2])

    if(which(sort(unique(visiualVector))=="Both")==1)
    {
      icolors <- c("green","red","blue")
    }else if(which(sort(unique(visiualVector))=="Both")==2){
      icolors <- c("red","green","blue")
    }else{
      icolors <- c("red","blue","green")
    }

    p1 <- visualSpatial(
      visiualVector,
      SpaCET_obj@input$image,
      SpaCET_obj@input$platform,
      scaleType="color-discrete",
      colors=icolors,
      limits=NULL,
      pointSize=1,
      pointAlpha=1,
      titleName="Spatial distribution of two cell-types",
      legendName="Spot",
      legend.position="none",
      imageBg=TRUE,
      imageSize = "CaptureArea",
      spotID=spotID
    )

    # scatter plot

    Content[is.na(Content)] <- "Other"

    fg.df <- as.data.frame(cbind(
      x=res_deconv[cellTypePair[1],],
      y=res_deconv[cellTypePair[2],]
    ))
    fg.df <- cbind(fg.df,group=Content)
    fg.df[,3] <- factor(fg.df[,3],levels=c("Both",cellTypePair,"Other"))

    p2 <- ggplot(fg.df,aes(x=x, y=y)) +
      geom_point(aes(colour=group),size=0.5)+
      scale_colour_manual(values=c("green","red","blue","grey"))+
      ggtitle(paste0("Spearman Rho = ",rho,", P ",pv1))+
      xlab(paste0("Cell fraction (",cellTypePair[1],")"))+
      ylab(paste0("Cell fraction (",cellTypePair[2],")"))+
      theme_bw()+
      theme(
        plot.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size=14,hjust = 0.5),
        axis.title = element_text(size=14,colour = "black"),
        axis.text = element_text(size=13,colour = "black"),
        legend.title=element_blank(),
        legend.position=c(.75,.85)
      )+
      geom_smooth(method = "lm", formula = "y ~ x", color="orange")


    # boxplot

    LRNetworkScoreMat <- SpaCET_obj@results$CCI$LRNetworkScore
    LRNetworkScoreMat[2,LRNetworkScoreMat[2,]> 2] <- 2

    fg.df <- data.frame(group=Content,value=LRNetworkScoreMat[2,],stringsAsFactors=FALSE)
    fg.df <- fg.df[!fg.df[,"group"]%in%"Other",]
    fg.df[fg.df[,1]%in%cellTypePair,1] <- "Single"

    xlab <- paste0(cellTypePair[1]," - ",cellTypePair[2])
    ylab <- "LR network score"

    p3 <- ggplot(fg.df,aes(x=group, y=value)) +
      geom_jitter(aes(colour = group),size=0.3)+
      geom_boxplot(aes(colour = group), outlier.shape = NA, alpha = 0.8,size=0.8)+
      scale_colour_manual(values=c("green","purple"))+
      xlab(xlab)+
      ylab(ylab)+
      ggtitle(paste0("Cohen's d=",cd1,", P=",pv2))+
      theme_bw()+
      theme(
        plot.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size=14,hjust = 0.5),
        axis.title = element_text(size=14,colour = "black"),
        axis.text = element_text(size=13,colour = "black"),
        legend.position="none"
      )

    library(patchwork)
    p1+p2+p3
  }
}


#' @title Identify tumor-stroma interface
#' @description Identify the spots at the tumor-stroma interface.
#' @param SpaCET_obj A SpaCET object.
#' @param Malignant Indicates the name of malignant cell type in the major lineage layer from the deconvolution results. Default: "Malignant".
#' @param MalignantCutoff Malignant cell fraction cutoff for tumor boundary. Default: 0.5.
#' @return A SpaCET object.
#' @examples
#' SpaCET_obj <- SpaCET.identify.interface(SpaCET_obj)
#' SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "Interface")
#'
#' @rdname SpaCET.identify.interface
#' @export
#'
SpaCET.identify.interface <- function(SpaCET_obj, Malignant="Malignant", MalignantCutoff=0.5)
{
  if(!grepl("visium", tolower(SpaCET_obj@input$platform)))
  {
    stop("This function is only applicable to 10X Visium data.")
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

  if(MalignantCutoff>1 | MalignantCutoff<0)
  {
    stop("Please input a value within 0~1 for the cutoff of malignant spots.")
  }

  Content <- res_deconv[Malignant,]
  names(Content) <- colnames(res_deconv)
  Content[Content>=MalignantCutoff] <- "Tumor"
  Content[Content!="Tumor"] <- "Stroma"

  spot_Non_Malignant <- names(Content)[Content=="Stroma"]
  Content_new <- Content

  for(spot in spot_Non_Malignant)
  {
    sxy <- unlist(strsplit(spot,"x"))
    sx <- as.numeric(sxy[1])
    sy <- as.numeric(sxy[2])

    spot_neighbor <- c(
      paste0(sx,"x",sy-2),
      paste0(sx,"x",sy+2),
      paste0(sx-1,"x",sy-1),
      paste0(sx+1,"x",sy+1),
      paste0(sx-1,"x",sy+1),
      paste0(sx+1,"x",sy-1)
    )

    spot_neighbor_cellType <- Content[names(Content)%in%spot_neighbor]

    if(length(unique(spot_neighbor_cellType))==1)
    {
      if(unique(spot_neighbor_cellType)=="Stroma")
      {
        Content_new[spot] <- "Stroma"
      }else{
        Content_new[spot] <- "Interface"
      }
    }else{
      Content_new[spot] <- "Interface"
    }

  }

  SpaCET_obj@results$CCI$interface <- matrix(
    Content_new, nrow=1, byrow=TRUE,
    dimnames = list("Interface",colnames(res_deconv))
    )

  SpaCET_obj
}



#' @title Combine interaction spots to interface
#' @description Demonstrate the spatial position of interaction spots at the tumor microenvironment.
#' @param SpaCET_obj A SpaCET object.
#' @param cellTypePair A pair of cell types.
#' @return A SpaCET object.
#' @examples
#' SpaCET_obj <- SpaCET.combine.interface(SpaCET_obj,cellTypePair=c("CAF", "Macrophage M2"))
#' SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "Interface", spatialFeature = "Interaction")
#'
#' @rdname SpaCET.combine.interface
#' @export
#'
SpaCET.combine.interface <- function(SpaCET_obj, cellTypePair)
{
  if(!grepl("visium", tolower(SpaCET_obj@input$platform)))
  {
    stop("This function is only applicable to 10X Visium data.")
  }

  if(length(cellTypePair)!=2)
  {
    stop("Please input a pair of cell-types.")
  }

  if(sum(cellTypePair%in%rownames(SpaCET_obj@results$deconvolution$propMat))!=2)
  {
    stop("Please input the correct cell-type name. Of note, R language is case sensitive generally.")
  }
  cellTypePair <- sort(cellTypePair)

  if(!"Interface"%in%rownames(SpaCET_obj@results$CCI$interface))
  {
    stop("Please run SpaCET.identify.interface first.")
  }

  if(!paste0(cellTypePair[1],"_",cellTypePair[2])%in%rownames(SpaCET_obj@results$CCI$interaction$testRes))
  {
    stop("Please run SpaCET.CCI.cellTypePair first.")
  }else{
    M2CAF <- SpaCET_obj@results$CCI$interaction$groupMat[paste0(cellTypePair[1],"_",cellTypePair[2]),]
    M2CAF <- names(M2CAF)[M2CAF%in%c("Both")]

    stromal <- colnames(SpaCET_obj@results$CCI$interface)[SpaCET_obj@results$CCI$interface["Interface",]=="Stroma"]
    M2CAF <- M2CAF[M2CAF%in%stromal]

    Content_new_new <- SpaCET_obj@results$CCI$interface["Interface",]
    Content_new_new[names(Content_new_new)%in%M2CAF] <- "Interaction"

    gname <- paste0("Interface&",cellTypePair[1],"_",cellTypePair[2])
    if(gname%in%rownames(SpaCET_obj@results$CCI$interface))
    {
      SpaCET_obj@results$CCI$interface[gname,] <- Content_new_new
    }else{
      SpaCET_obj@results$CCI$interface <- rbind(SpaCET_obj@results$CCI$interface,Content_new_new)
      rownames(SpaCET_obj@results$CCI$interface)[nrow(SpaCET_obj@results$CCI$interface)] <- gname
    }

    SpaCET_obj
  }

}



#' @title Distance to tumor border
#' @description Calculate the distance of cell-cell interactions to tumor-immune interface.
#' @param SpaCET_obj A SpaCET object.
#' @param cellTypePair A pair of cell types.
#' @param nPermutation Permutation number.
#' @return A ggplot object.
#' @examples
#' SpaCET.distance.to.interface(SpaCET_obj, cellTypePair=c("CAF", "Macrophage M2"))
#'
#' @rdname SpaCET.distance.to.interface
#' @export
#'
SpaCET.distance.to.interface <- function(SpaCET_obj, cellTypePair=c("CAF","Macrophage M2"), nPermutation = 1000)
{
  if(!grepl("visium", tolower(SpaCET_obj@input$platform)))
  {
    stop("This function is only applicable to 10X Visium data.")
  }

  if(length(cellTypePair)!=2)
  {
    stop("Please input a pair of cell-types.")
  }

  if(sum(cellTypePair%in%rownames(SpaCET_obj@results$deconvolution$propMat))!=2)
  {
    stop("Please input the correct cell-type name. Of note, R language is case sensitive generally.")
  }
  cellTypePair <- sort(cellTypePair)
  testRes <- SpaCET_obj@results$CCI$interaction$testRes

  if(testRes[paste0(cellTypePair[1],"_",cellTypePair[2]),"Interaction"] == FALSE)
  {
    if(is.na(testRes[paste0(cellTypePair[1],"_",cellTypePair[2]),"groupCompare_pv"]))
    {
      message("The colocalization analysis is not significant for the current cell-type pair. Please check other cell-type pairs.")
    }else{
      message("Based on colocalization analysis and L-R enrichment analysis, the intercellular interaction is not significant for the current cell-type pair. Please check other cell-type pairs.")
    }
  }else{

    interface <- SpaCET_obj@results$CCI$interface
    spot_Malignant_border <- colnames(interface)[interface[1,]%in%"Interface"]
    spot_Non_Malignant_core <- colnames(interface)[interface[1,]%in%"Stroma"]

    groupCellTypes <- SpaCET_obj@results$CCI$interaction$groupMat [paste0(cellTypePair[1],"_",cellTypePair[2]), ,drop=F]

    M2CAF <- colnames(groupCellTypes)[groupCellTypes[1,]%in%c("Both")]
    M2CAF <- M2CAF[M2CAF%in%spot_Non_Malignant_core] #filter by malignant

    M2_CAF <- colnames(groupCellTypes)[groupCellTypes[1,]%in%cellTypePair]
    M2_CAF <- M2_CAF[M2_CAF%in%spot_Non_Malignant_core] #filter by malignant

    calSpotSpotDistance <- function(spot1, spot2)
    {
      spot1_xy <- as.numeric(unlist(strsplit(spot1,"x")))
      spot2_xy <- as.numeric(unlist(strsplit(spot2,"x")))

      sqrt((spot1_xy[1]-spot2_xy[1])^2+(spot1_xy[2]-spot2_xy[2])^2)
    }

    calSpotVecDistance <- function(spot1, vec2)
    {
      min(sapply(vec2,calSpotSpotDistance,spot1=spot1))
    }

    calVecVecDistance <- function(vec1,vec2)  # vec1 M2CAF, vec2 border
    {
      sapply(vec1,calSpotVecDistance,vec2=vec2)
    }

    d0 <- mean( calVecVecDistance(M2CAF,spot_Malignant_border) )
    db <- calVecVecDistance(M2_CAF,spot_Malignant_border)

    dd <- function(i)
    {
      set.seed(i)
      mean( db[sample(M2_CAF,length(M2CAF))] )
    }

    d_list <- sapply(1:nPermutation, dd)
    d_vec <- unlist(d_list)

    pv <- (sum(d_vec<=d0)+1)/(nPermutation+1)

    fg.df <- data.frame(value=d_vec)

    ggplot(fg.df,aes(x=value)) +
      geom_histogram(aes(y=..density..), colour="grey1", fill="gainsboro", alpha=0.5)+
      geom_density(color="grey3",size=0.6)+
      ggtitle(paste0("Permutation\n P = ",signif(pv,3)))+
      ylab("Density")+
      xlab("Distance to Tumor-Stroma Interface")+
      theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size=14,colour = "black"),
        axis.text = element_text(size=13,colour = "black"),
        legend.position="none",
        #panel.border = element_blank(),
        axis.line.y.left = element_line(color = 'black'),
        axis.line.x.bottom = element_line(color = 'black')
      )+
      geom_vline(xintercept=d0, color = "green", linetype="dashed",size=1.3)
  }
}
