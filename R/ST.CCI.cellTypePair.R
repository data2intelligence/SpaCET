#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ST PARAM_DESCRIPTION
#' @param cellTypePair PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[cowplot]{plot_grid}}
#' @rdname ST.CCI.cellTypePair.scatter
#' @export 
#' @importFrom cowplot plot_grid
ST.CCI.cellTypePair.scatter <- function(ST,cellTypePair)
{
  res_deconv <- ST@results$fraction

  cellTypePairName<- paste0(cellTypePair,collapse="-")
  Content <- colSums(res_deconv[cellTypePair,])
  
  for(i in 1:length(Content))
  {
    if(res_deconv[cellTypePair[2],i]>summary(res_deconv[cellTypePair[2],])[5]&res_deconv[cellTypePair[1],i]>summary(res_deconv[cellTypePair[1],])[5])
    {
      Content[i] <- cellTypePairName
    }else{
      if(res_deconv[cellTypePair[2],i]>summary(res_deconv[cellTypePair[2],])[5])
      {
        Content[i] <- cellTypePair[2]
      }else if(res_deconv[cellTypePair[1],i]>summary(res_deconv[cellTypePair[1],])[5]){
        Content[i] <- cellTypePair[1]
      }else{
        Content[i] <- "rest"
      }
      
    }
    
  }
  
  
  library(ppcor)
  ppcor_res <- pcor.test(res_deconv[cellTypePair[1],],res_deconv[cellTypePair[2],],res_deconv["Malignant",], method="spearman")
  ppcor_rho <- round(ppcor_res[1,"estimate"],3)
  ppcor_pv <-  signif(ppcor_res[1,"p.value"],2)
  
  fg.df <- data.frame(
    x=res_deconv[cellTypePair[1],],
    y=res_deconv[cellTypePair[2],],
    group=Content
  )
  
  fg.df[,3] <- factor(fg.df[,3], levels=c(cellTypePairName,cellTypePair,"rest")) 
  
  p1 <- ggplot(fg.df,aes(x=x, y=y)) + 
    geom_point(aes(colour=group),size=0.5)+
    scale_colour_manual(values=c("green","red","blue","grey"))+
    ggtitle(paste0("Partial Correlation: rho = ",ppcor_rho,", p = ",ppcor_pv))+
    xlab(paste0("Cell fraction (",cellTypePair[1],")"))+
    ylab(paste0("Cell fraction (",cellTypePair[2],")"))+
    guides(colour=guide_legend(title=""))+
    theme_bw()+ 
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size=14,colour = "black"),
      axis.text = element_text(size=13,colour = "black"),
      legend.position = "right"
    )+
    geom_smooth(method = "lm",color="orange")
  
  colVec <- c("green","red","blue")
  names(colVec) <- c(cellTypePairName,cellTypePair)
  
  p2 <- visualSpatial(
    Content[Content!="rest"],
    ST@input$HEimage,
    cols=colVec[sort(names(colVec))],
    colsLimits=NULL,
    "Spatial distribution of cell co-localization",
    "",
    legendPosition="none",
    is.continuous=FALSE
  )
  
  cowplot::plot_grid(p1,p2, ncol = 2, align = "h",rel_widths=c(1.2, 1))
  
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ST PARAM_DESCRIPTION
#' @param cellTypePair PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[cowplot]{plot_grid}}
#' @rdname ST.CCI.cellTypePair.boxplot
#' @export 
#' @importFrom cowplot plot_grid
ST.CCI.cellTypePair.boxplot <- function(ST,cellTypePair)
{
  res_deconv <- ST@results$fraction
  
  cellTypePairName<- paste0(cellTypePair,collapse="-")
  Content <- colSums(res_deconv[cellTypePair,])
  
  for(i in 1:length(Content))
  {
    if(res_deconv[cellTypePair[2],i]>summary(res_deconv[cellTypePair[2],])[5]&res_deconv[cellTypePair[1],i]>summary(res_deconv[cellTypePair[1],])[5])
    {
      Content[i] <- cellTypePairName
    }else{
      if(res_deconv[cellTypePair[2],i]>summary(res_deconv[cellTypePair[2],])[5])
      {
        Content[i] <- cellTypePair[2]
      }else if(res_deconv[cellTypePair[1],i]>summary(res_deconv[cellTypePair[1],])[5]){
        Content[i] <- cellTypePair[1]
      }else{
        Content[i] <- "rest"
      }
      
    }
    
  }
  
  
  bb <- data.frame(
    x=Content,
    y=ST@results$LRInteraction[1,],
    z=-log10(ST@results$LRInteraction[2,])
  )
  
  bb <- bb[bb[,"x"]%in%c(cellTypePairName,cellTypePair),]
  
  pv1 <- signif(wilcox.test(bb[bb[,1]%in%cellTypePairName,2],bb[bb[,1]%in%c(cellTypePair),2])$p.value,2)
  pv2 <- signif(wilcox.test(bb[bb[,1]%in%cellTypePairName,3],bb[bb[,1]%in%c(cellTypePair),3])$p.value,2)
  
  bb[bb[,1]%in%cellTypePair,1] <- paste0(paste0(cellTypePair,collapse="/"),"\ndominated")   
  bb[bb[,1]%in%cellTypePairName,1] <- paste0(cellTypePairName,"\nco-localization")
  
  
  p3 <- ggplot(bb,aes(x=x, y=y)) + 
    geom_boxplot(aes(colour = x), alpha = 0.6,size=0.8)+
    scale_colour_manual(values=c("green","purple"))+
    ggtitle(paste0("Wilcoxon test: P = ",pv1))+
    ylab("Network Score")+
    theme_bw()+ 
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size=12,colour = "black"),
      axis.title = element_text(size=15,colour = "black"),
      axis.title.x = element_blank(),
      legend.position="none"
    )
  
  p4 <- ggplot(bb,aes(x=x, y=z)) + 
    geom_boxplot(aes(colour = x), alpha = 0.6,size=0.8)+
    scale_colour_manual(values=c("green","purple"))+
    ggtitle(paste0("Wilcoxon test: P = ",pv2))+
    ylab("-Log10pv")+
    theme_bw()+ 
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size=12,colour = "black"),
      axis.title = element_text(size=15,colour = "black"),
      axis.title.x = element_blank(),
      legend.position="none"
    )

  cowplot::plot_grid(p3,p4, ncol = 2)
  
}





