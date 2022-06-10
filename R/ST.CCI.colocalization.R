#' @title Cell colocalization
#' @description Calculate the cell-cell colocalization
#' @param ST An SpaCE object
#' @return An SpaCE object
#' @details This function calculates the Spearman correlation of cell types as colocalization.
#' @examples 
#' 
#' @rdname ST.CCI.colocalization
#' @export 
#' @importFrom reshape2 melt
#' @importFrom psych corr.test
ST.CCI.colocalization <- function(ST)
{
  MacrophageName <- c("M1","M2","Macrophage other")
  TCD4Name <- c("Treg","T CD4 naive","T helper","T CD4 other")
  
  ST.deconv.res <- ST@results$fraction
  
  ST.deconv.res <- ST.deconv.res[!rownames(ST.deconv.res)%in%c("Unidentifiable","Macrophage","T CD4",MacrophageName[3],TCD4Name[4]),]
  ST.deconv.res <- round(ST.deconv.res,2)
  
  
  cc_corr <- psych::corr.test(
    t(ST.deconv.res),
    t(ST.deconv.res),
    method="spearman",adjust="none",ci=FALSE)
  
  cc_corr_r <- round(cc_corr$r,3)
  cc_corr_p <- round(-log10(cc_corr$p),3)
  
  summary_df <- reshape2::melt(cc_corr_r)
  summary_df2 <- reshape2::melt(cc_corr_p)
  summary_df <- cbind(summary_df,log10pv=summary_df2[,3])
  
  summary_df[is.na(summary_df[,"log10pv"]),c("value","log10pv")] <- 0
  if(sum(summary_df[,"log10pv"]>5)>0)
  {	
    summary_df[summary_df[,"log10pv"]>5,"log10pv"] <- 5
  }
  
  summary_df <- summary_df[!(summary_df[,1]==summary_df[,2]),]
  
  
  for(i in 1:dim(summary_df)[1])
  {
    if(summary_df[i,1]!="Malignant"&summary_df[i,2]!="Malignant"&summary_df[i,1]!=summary_df[i,2])
    {
      ppcor_res <- ppcor::pcor.test(
        ST.deconv.res[summary_df[i,1],],
        ST.deconv.res[summary_df[i,2],],
        ST.deconv.res["Malignant",],
        method="spearman")
      summary_df[i,3] <- ppcor_res[1,"estimate"]
      summary_df[i,4] <- -log10(ppcor_res[1,"p.value"])
    }
  }
  
  ST@results$colocalization <- summary_df
  
  ST
}


#' @title Cell colocalization plot
#' @description visualize the cell-cell colocalization
#' @param ST An SpaCE object
#' @return A gggplot2 object
#' @details This function visualizes the cell colocalization in point and network plots.
#' @examples 
#' 
#' @rdname ST.CCI.colocalization.plot
#' @export 
#' @importFrom cowplot plot_grid
ST.CCI.colocalization.plot <- function(ST)
{
  ST@results$colocalization -> summary_df
 
  if(sum(summary_df[,"log10pv"]>5)>0)
  {	
    summary_df[summary_df[,"log10pv"]>5,"log10pv"] <- 5
  }
  
  summary_df[,1] <- factor(summary_df[,1],levels=rev(levels(summary_df[,1])))
  
  p1 <- ggplot(summary_df, aes( Var2, Var1)) + 
    geom_point(aes(colour = value,size=log10pv),na.rm = TRUE) + 
    scale_colour_gradient2(low = "blue", high = "red", mid = "white", na.value = NA,
                           midpoint = 0, limit = c(-0.6,0.6), space = "Lab", 
                           name="Rho")+
    ggtitle("")+
    labs(size = "-Log10pv")+
    theme_bw() + 
    theme(
      strip.text = element_blank(),
      axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(colour = "black"),
      axis.title = element_blank()
    )
  

  summary_df <- summary_df[order(abs(summary_df[,3])),]
  graph <- as_tbl_graph(summary_df[as.character(summary_df[,1])>as.character(summary_df[,2]),]) 
  
  # plot using ggraph
  p2 <- ggraph(graph, layout = 'circle') +
    geom_edge_fan2(aes(colour = value, width = log10pv), show.legend = TRUE)+
    scale_edge_colour_gradient2(
      name = "Rho",
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      space = "Lab",
      guide = "edge_colourbar"
    )+
    scale_edge_width(name = "-Log10pv")+
    geom_node_point(size = 5, color="grey") +
    geom_node_label(aes(label = name), size = 5, color="black") +
    ggtitle(" ")+
    xlab("")+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title.y = element_blank()
    )
  
  
  cowplot::plot_grid(p1, p2, labels = "AUTO", ncol = 2)
  
}

