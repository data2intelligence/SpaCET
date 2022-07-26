#' @title Cell-cell colocalization analysis
#' @description Calculate the cell-type pair colocalization by using Spearman correlation analysis.
#' @param SpaCE_obj An SpaCE object.
#' @return An SpaCE object
#' @examples
#' SpaCE_obj <- SpaCE.CCI.colocalization(SpaCE_obj)
#' @rdname SpaCE.CCI.colocalization
#' @export
#'
SpaCE.CCI.colocalization <- function(SpaCE_obj)
{
  res_deconv <- SpaCE_obj@results$deconvolution
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


  load( system.file("extdata",'combRef_0.5.rda',package = 'SpaCE') )
  reff <- Ref$refProfiles[unique(unlist(Ref$sigGenes[names(Ref$sigGenes)%in%c(names(Ref$lineageTree),"T cell")])),]
  reff <- reff-rowMeans(reff)

  cc_corr <- psych::corr.test(reff,reff,method="pearson",adjust="none",ci=FALSE)

  cc_corr_r <- cc_corr$r
  cc_corr_p <- cc_corr$p

  cc_corr_r.m <- reshape2::melt(cc_corr_r)
  cc_corr_p.m <- reshape2::melt(cc_corr_p)

  summary_df2 <- data.frame(
    cell_type_1 = as.character(cc_corr_r.m[,1]),
    cell_type_2 = as.character(cc_corr_r.m[,2]),
    reference_r = round(cc_corr_r.m[,3],3),
    reference_pv = cc_corr_p.m[,3]
  )
  rownames(summary_df2) <- paste0(summary_df2[,1],"_",summary_df2[,2])

  summary_df[rownames(summary_df2),"reference_r"] <- summary_df2[,"reference_r"]
  summary_df[rownames(summary_df2),"reference_pv"] <- summary_df2[,"reference_pv"]

  summary_df <- summary_df[ summary_df[,1] != summary_df[,2], ] #remove same cell type

  SpaCE_obj@results$CCI.colocalization <- summary_df

  SpaCE_obj
}


#' @title Cell-cell colocalization visualization
#' @description Visualize the cell-type pair colocalization.
#' @param SpaCE_obj An SpaCE object.
#' @return An ggplot object
#' @examples
#' SpaCE.visualize.colocalization(SpaCE_obj)
#' @rdname SpaCE.visualize.colocalization
#' @export
#'
SpaCE.visualize.colocalization <- function(SpaCE_obj)
{
  summary_df <- SpaCE_obj@results$CCI.colocalization

  summary_df[summary_df[,"fraction_product"]>0.02,"fraction_product"] <- 0.02

  ctOrder <- c("Malignant",unlist(SpaCE_obj@results$Ref$lineageTree))

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
    ggtitle(" ")+
    xlab(" ")+
    ylab(" ")+
    labs(size = "Fraction\nproduct")+
    theme_bw() +
    theme(
      panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"),
      strip.text = element_blank(),
      axis.text.x = element_text(size = 11, colour = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 11, colour = "black")
    )


  summary_df <- summary_df[as.character(summary_df[,1])<as.character(summary_df[,2]),]

  summary_df <- summary_df[!as.character(summary_df[,1])%in%c("Malignant"),]
  summary_df <- summary_df[!as.character(summary_df[,2])%in%c("Malignant"),]

  summary_df <- cbind(summary_df,label=paste0(summary_df[,1],"_",summary_df[,2]))
  summary_df[summary_df[,"fraction_product"]<0.0005 ,"label"] <- NA

  summary_df[is.na(summary_df[,"fraction_rho"]),c("fraction_rho")] <- 0

  res <- lm(summary_df$fraction_rho ~ summary_df$reference_r)
  summary_df <- cbind(summary_df,Residual=abs(residuals(res)))

  res_weight <- lm(summary_df$fraction_rho ~ summary_df$reference_r, weight=summary_df$fraction)

  p2 <- ggplot(summary_df, aes( reference_r, fraction_rho, label=label )) +
    geom_point(aes(size=fraction_product),color="#856aad")+
    geom_smooth(method="lm",color="darkgrey",mapping = aes(weight = fraction_product))+
    scale_size(range = c(0, 6))+
    ggrepel::geom_text_repel()+
    ggtitle(" ")+
    xlab("Cor (Reference profiles)")+
    ylab("Cor (Cell fractions)")+
    labs(size = "Fraction\nProduct")+
    ylim(min(summary_df[,"fraction_rho"])-0.1,max(summary_df[,"fraction_rho"])+0.1)+
    xlim(-0.6,1)+
    theme_bw()+
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      strip.text = element_blank(),
      axis.text = element_text(size = 16, colour = "black"),
      axis.title = element_text(size = 19, colour = "black"),
      panel.border = element_blank(),
      axis.line.y.left = element_line(color = 'black'),
      axis.line.x.bottom = element_line(color = 'black')
    )

  library(patchwork)
  p1+p2
}


#' @title Ligand-Receptor interaction enrichment analysis
#' @description Calculate the overall intensity of L-R interactions at each ST spot.
#' @param SpaCE_obj An SpaCE object.
#' @return An SpaCE object
#' @examples
#' SpaCE_obj <- SpaCE.CCI.LRNetworkScore(SpaCE_obj)
#' @rdname SpaCE.CCI.LRNetworkScore
#' @export
#'
SpaCE.CCI.LRNetworkScore <- function(SpaCE_obj)
{
  st.matrix.data <- as.matrix(SpaCE_obj@input$counts)

  st.matrix.data <- t(t(st.matrix.data)*1e5/colSums(st.matrix.data))
  st.matrix.data <- log2(st.matrix.data+1)
  spotNum <- ncol(st.matrix.data)

  ###### permute L-R network ######
  LRdb <- read.csv(system.file("extdata",'LR.txt',package = 'SpaCE'),as.is=T,sep="\t")
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

  ###### permute
  LRdb_string <- c()

  set.seed(123456)
  for(i in 1:1000)
  {
    LRdb_rand <- BiRewire::birewire.rewire.bipartite(LRdb_mat,verbose=FALSE)

    LRdb_rand.m <- as.matrix(reshape2::melt(LRdb_rand))
    LRdb_rand.m <- LRdb_rand.m[LRdb_rand.m[,3]==1,1:2]
    colnames(LRdb_rand.m) <- c("L","R")

    LRdb_string <- c(LRdb_string,c(t(LRdb_rand.m)))
  }

  LRdb_rand_comb <- matrix(LRdb_string,ncol=2,byrow=TRUE)
  ###### permute L-R network ######


  LRcoexpr <- function(LRdb,spot)
  {
    Ls <- LRdb[,1]
    Rs <- LRdb[,2]
    mean( spot[Ls] * spot[Rs] )
  }

  LRNetworkScore <- function(i)
  {
    spot <- st.matrix.data[,i]

    LR_raw <- LRcoexpr(LRdb,spot=spot)

    LR1000 <- LRcoexpr(LRdb_rand_comb,spot=spot)

    score <- LR_raw/LR1000

    list(LR_raw,score)
  }

  LRNetworkScoreList <- parallel::mclapply(1:spotNum,LRNetworkScore,mc.cores=8)

  LRNetworkScoreMat <- matrix(unlist(LRNetworkScoreList),ncol=ncol(st.matrix.data))
  colnames(LRNetworkScoreMat) <- colnames(st.matrix.data)
  rownames(LRNetworkScoreMat) <- c("Raw_expr","Network_Score")


  SpaCE_obj@results$LRNetworkScore <- LRNetworkScoreMat

  SpaCE_obj
}


#' @title Ligand-Receptor analysis for a co-localized cell-type pair
#' @description Test the co-expression of ligands and receptors within the same ST spot for the co-localized cell-type pairs.
#' @param SpaCE_obj An SpaCE object.
#' @param cellTypePair Cancer type of this tumor ST dataset.
#' @return An SpaCE object
#' @examples
#' SpaCE.CCI.cellTypePair(SpaCE_obj,cellTypePair=c("CAF","Macrophage M2"))
#' @rdname SpaCE.CCI.cellTypePair
#' @export
#'
SpaCE.CCI.cellTypePair <- function(SpaCE_obj,cellTypePair)
{
  res_deconv <- SpaCE_obj@results$deconvolution

  summary_df <- SpaCE_obj@results$CCI.colocalization
  rho <- summary_df[summary_df[,1]==cellTypePair[1]&summary_df[,2]==cellTypePair[2],"fraction_rho"]
  pv <-  summary_df[summary_df[,1]==cellTypePair[1]&summary_df[,2]==cellTypePair[2],"fraction_pv"]
  pv <- signif(as.numeric(pv),3)
  pv <-  ifelse(pv==0," < 2.2e-16",paste0(" = ",pv))


  cutoff1 <- summary(res_deconv[cellTypePair[1],])[5]
  cutoff2 <- summary(res_deconv[cellTypePair[2],])[5]

  cutoff11 <- quantile(res_deconv[cellTypePair[1],],0.85)
  cutoff22 <- quantile(res_deconv[cellTypePair[2],],0.85)

  Content <- res_deconv["Unidentifiable",]
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

  visiualVector <- Content
  names(visiualVector) <- paste0(SpaCE_obj@input$spotCoordinates[1,],"x",SpaCE_obj@input$spotCoordinates[2,])
  p2 <- visualSpatial(
    visiualVector,
    SpaCE_obj@input$image,
    SpaCE_obj@input$platform,
    scaleType="color-discrete",
    c("green","red","blue"),
    pointSize=1,
    pointAlpha=1,
    limits,
    "",
    "Spot",
    imageBg=TRUE)


  Content[is.na(Content)] <- "None"

  fg.df <- as.data.frame(cbind(
    x=res_deconv[cellTypePair[1],],
    y=res_deconv[cellTypePair[2],]
  ))
  fg.df <- cbind(fg.df,group=Content)

  fg.df <- fg.df[order(fg.df[,3],decreasing=T),]

  p1 <- ggplot(fg.df,aes(x=x, y=y)) +
    geom_point(aes(colour=group),size=0.5)+
    scale_colour_manual(values=c("green","red","blue","grey"))+
    ggtitle(paste0("Spearman Rho = ",rho,", P ",pv))+
    xlab(paste0("Cell fraction (",cellTypePair[1],")"))+
    ylab(paste0("Cell fraction (",cellTypePair[2],")"))+
    guides(colour=guide_legend(title=""))+
    theme_bw()+
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size=14,hjust = 0.5),
      axis.title = element_text(size=14,colour = "black"),
      axis.text = element_text(size=13,colour = "black"),
      legend.position = "none",
      panel.border = element_blank(),
      axis.line.y.left = element_line(color = 'black'),
      axis.line.x.bottom = element_line(color = 'black')
    )+
    geom_smooth(method = "lm", formula = "y ~ x", color="orange")


  # boxplot

  LRNetworkScoreMat <- SpaCE_obj@results$LRNetworkScore

  #LRNetworkScoreMat[3,] <- -log10(LRNetworkScoreMat[3,])

  fg.df <- data.frame(group=Content,value=LRNetworkScoreMat[2,],stringsAsFactors=FALSE)
  fg.df <- fg.df[!fg.df[,"group"]%in%"None",]

  n1 <- sum(fg.df[,1]%in%"Both")
  n2 <- sum(!fg.df[,1]%in%"Both")

  pv <- signif(wilcox.test(fg.df[fg.df[,1]%in%"Both",2],fg.df[fg.df[,1]%in%cellTypePair,2])$p.value,2)

  fg.df[fg.df[,1]%in%cellTypePair,1] <- "Single"

  cohend_res <- psych::cohen.d(fg.df, group="group", alpha=.05, std=TRUE)
  cd1 <- signif(cohend_res$cohen.d["value","effect"],2)

  ylab <- "-Log10 (pv)"

  p3 <- ggplot(fg.df,aes(x=group, y=value)) +
    geom_jitter(aes(colour = group),size=0.3)+
    geom_boxplot(aes(colour = group), outlier.shape = NA, alpha = 0.8,size=0.8)+
    scale_colour_manual(values=c("green","purple"))+
    ylab(ylab)+
    ggtitle(paste0("Cohen's d=",cd1,", P=",pv))+
    theme_bw()+
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size=14,hjust = 0.5),
      axis.text = element_text(size=13,colour = "black"),
      axis.title = element_text(size=14,colour = "black"),
      axis.title.x = element_blank(),
      panel.border = element_blank(),
      axis.line.y.left = element_line(color = 'black'),
      legend.position="none"
    )

  library(patchwork)
  p1+p2+p3
}

