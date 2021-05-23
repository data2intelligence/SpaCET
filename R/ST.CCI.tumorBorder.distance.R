#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ST An SpaCE object
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
#' @rdname ST.CCI.tumorBorder.distance
#' @export 
#' @importFrom cowplot plot_grid
ST.CCI.tumorBorder.distance <- function(ST,cellTypePair)
{
    res_deconv <- ST@results$fraction
    res_deconv_a <- res_deconv[1:11,]
    
    Content <- sapply(1:dim(res_deconv_a)[2],function(x) names(sort(res_deconv_a[,x],decreasing=T))[1])
    names(Content) <- ST@input$spotID
    
    Content[Content!="Malignant"] <- "Non-Malignant"
    
    spot_Malignant <- names(Content)[Content=="Malignant"]
    spot_Non_Malignant <- names(Content)[Content=="Non-Malignant"]
    
    Content_new <- Content

    for(spot in spot_Malignant)
    {
      sxy <- unlist(strsplit(spot,"x"))
      sx <- as.numeric(sxy[1])
      sy <- as.numeric(sxy[2])
      
      spot_neighbor <- c(
        #paste0(sx,"x",sy),
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
        if(unique(spot_neighbor_cellType)=="Malignant")
        {
          Content_new[spot] <- "Malignant_core"
        }else{
          Content_new[spot] <- "Malignant_border"
        }
      }else{
        Content_new[spot] <- "Malignant_border"
      }
    }
    
    names(Content_new) <- colnames(res_deconv_a)
    
    p1 <- visualSpatial(
      Content_new[Content_new=="Malignant_border"],
      ST@input$HEimage,
      cols=c("yellow"),
      colsLimits=NULL,
      "Tumor Border",
      "",
      legendPosition="none",
      is.continuous=FALSE
    )
    
    
    
    
    
    
    
    names(Content_new) <- ST@input$spotID
    
    spot_Malignant_border <- names(Content_new)[Content_new=="Malignant_border"]
    spot_Malignant_core <- names(Content_new)[Content_new=="Malignant_core"]
    
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
    names(Content) <- ST@input$spotID
    
    M2CAF <- names(Content)[Content%in%cellTypePairName]
    
    M2CAF <- setdiff(M2CAF,spot_Malignant) #filter by malignant
    
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
    
    calVecVecDistance <- function(vec1,vec2) # vec1 M2CAF, vec2 border
    {
      mean(sapply(vec1,calSpotVecDistance,vec2=vec2))
    }
    
    d0 <- calVecVecDistance(M2CAF,spot_Malignant_border)
    
    d_vec <- c()
    for(i in 1:1000)
    {
      temp <- sample(spot_Non_Malignant,length(M2CAF))
      
      d_vec <- c(d_vec, calVecVecDistance(temp,spot_Malignant_border) )
    }
    
    no <- sum(d_vec<d0)
    if(no==0)
    {
      pv <- "< 0.001"
    }else{
      pv<- paste0("= ",no/1000)
    }

    fg.df <- data.frame(value=d_vec)
    library(ggplot2)
    p2 <- ggplot(fg.df,aes(x=value)) + 
      geom_histogram(aes(y=..density..), colour="black", fill="grey")+
      ylab("Density")+
      ggtitle(paste0("P ",pv))+
      theme_bw()+ 
      theme(
        plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size=17,colour = "black"),
        axis.text = element_text(size=15,colour = "black"),
        legend.position="none"
      )+
      geom_vline(xintercept=d0, color = "green", linetype="dashed",size=2)
    
    cowplot::plot_grid(p1,p2, ncol = 2, align = "h",rel_widths=c(1, 1))
    
}  
