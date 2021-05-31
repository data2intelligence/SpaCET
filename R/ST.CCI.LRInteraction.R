#' @title Ligand-Receptor interaction
#' @description Calculate the Ligand-Receptor interaction with network score and p value
#' @param ST An SpaCE object
#' @return An SpaCE object
#' @details SpaCE computes an L-R interactions network score for each spot as the sum of expression products between L-R pairs, normalized by the average L-R network score from 1,000 random networks with the same degrees. The p-value is calculated as the fraction of random network scores that exceeded the original score.
#' @examples 
#' ST <- ST.CCI.LRInteraction(ST)
#' @rdname ST.CCI.LRInteraction
#' @export 
#' @importFrom BiRewire birewire.rewire.bipartite
#' @importFrom reshape2 melt
ST.CCI.LRInteraction <- function(ST)
{
  st.matrix.data <- ST@input$counts
  
  st.matrix.data <- t(t(st.matrix.data)*1e5/colSums(st.matrix.data))
  st.matrix.data <- log2(st.matrix.data+1)
  st.matrix.data <- st.matrix.data[rowSums(st.matrix.data)>0,]
  spotNum <- ncol(st.matrix.data)
  
  
  ###### filter and permute L-R network ######
  
  LR2015 <- read.csv(system.file("extdata",'LR.txt', package = 'SpaCE'),as.is=T,sep="\t")
  LR2015 <- data.frame(L=LR2015[,2],R=LR2015[,4],stringsAsFactors=FALSE)
  
  Ls <- unique(LR2015[,"L"])
  Rs <- unique(LR2015[,"R"])
  
  LR2015 <- LR2015[LR2015[,"R"]!="DLK2",]
  
  LR2015 <- LR2015[LR2015[,1]%in%rownames(st.matrix.data)&LR2015[,2]%in%rownames(st.matrix.data),]
  rownames(LR2015) <- paste0(LR2015[,1],"_",LR2015[,2])
  
  LR_list <- list("0"=LR2015)
  

  Ls <- unique(LR2015[,"L"])
  Rs <- unique(LR2015[,"R"])
  
  LR2015_mat <- matrix(0,nrow=length(Ls),ncol=length(Rs))
  rownames(LR2015_mat) <- Ls
  colnames(LR2015_mat) <- Rs
  
  for(i in 1:dim(LR2015)[1])
  {
    LR2015_mat[LR2015[i,1],LR2015[i,2]] <- 1
  }
  
  set.seed(123456)
  for(i in 1:1000)
  {
    LR2015_rand <- BiRewire::birewire.rewire.bipartite(LR2015_mat,verbose=FALSE)
    
    LR2015_rand.m <- reshape2::melt(LR2015_rand)
    LR2015_rand.m <- LR2015_rand.m[LR2015_rand.m[,3]==1,1:2]
    colnames(LR2015_rand.m) <- c("L","R")
    
    LR_list[[as.character(i)]] <- as.matrix(LR2015_rand.m)
  }
  
  ###### filter and permute L-R network ######
  

  ###### LR_coexpr_raw ######
  
  LRcal <- function(j,spot)
  {
    sum( spot[ LR_list[[j]][,1] ] * spot[ LR_list[[j]][,2] ] )
  }
  
  statSig <- function(i)
  {
    spot <- st.matrix.data[,i]
    
    LR0 <- LRcal(as.character("0"),spot)
    
    LR1000 <- sapply(as.character(1:1000),LRcal,spot=spot)
    
    score <- LR0/mean(LR1000)
    pv <- (sum(LR1000>LR0)+1)/1001
    
    list(score,pv)
  }
  
  
  library(parallel)
  statSigList <- mclapply(1:spotNum,statSig,mc.cores=5)

  statSigMat <- matrix(unlist(statSigList),byrow=TRUE,ncol=2)
  rownames(statSigMat) <- colnames(st.matrix.data)

  ST@results$LRInteraction <- t(statSigMat)
  rownames(ST@results$LRInteraction ) <- c("Network Score", "P value")
  
  ST
}


#' @title Ligand-Receptor interaction plot
#' @description Plot the Ligand-Receptor interaction with network score and p value
#' @param ST An SpaCE object
#' @return A ggplot2 object
#' @details This function show the L-R network scores and p values for ST spots.
#' @examples 
#' ST.CCI.LRInteraction.plot(ST)
#' @rdname ST.CCI.LRInteraction.plot
#' @export 
#' @importFrom cowplot plot_grid
ST.CCI.LRInteraction.plot <- function(ST)
{
  Content <- ST@results$LRInteraction[1,]
  Content[Content==0] <- sort(setdiff(Content,0))[1] 
  Content <- log2(Content)
  
  Content[Content>0.5] <- 0.5
  Content[Content< -0.5] <- -0.5
  
  p1 <- visualSpatial(
    Content,
    ST@input$HEimage,
    cols=c("#081d58","#253494","blue","cyan","yellow"),
    colsLimits=NULL,
    "",
    "Log2\n(Network\nScore)"
  )
  
  Content <- ST@results$LRInteraction[2,]
  Content <- -log10(Content)
  
  p2 <- visualSpatial(
    Content,
    ST@input$HEimage,
    cols=c("blue","cyan","yellow"),
    colsLimits=NULL,
    "",
    "-Log10pv"
  )
  
  cowplot::plot_grid(p1,p2, nrow=1, rel_widths=c(1, 1))
  
}
