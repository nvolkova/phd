plot_decomposition <- function(decomposition, mm, intnames, col) {
  for (i in 1:nrow(decomposition))
    decomposition[i,] = decomposition[i,] / sum(decomposition[i,])
  new.cont.mat <- t(decomposition[intersect(intnames,rownames(mm)),])
  for (y in colnames(new.cont.mat)) {
    new.cont.mat[,y] = new.cont.mat[,y] * rowSums(mm)[y]
  }
  m.new.cont.mat <- melt(new.cont.mat)
  colnames(m.new.cont.mat) = c("Signature","Sample","Contribution")
  order(rowSums(mm),decreasing = F) -> sampleorder
  names(sampleorder) <- row.names(mm)[sampleorder]
  plot = ggplot(m.new.cont.mat, aes(x = factor(Sample,levels=names(sampleorder)), 
                                    y = Contribution, fill = factor(Signature,levels=colnames(decomposition)), order = Sample)) + 
    geom_bar(stat = "identity", colour = "black") + 
    labs(x = "", y = "Absolute contribution \n (no. mutations)") + 
    theme_bw() + 
    scale_fill_manual(name="",values=col) + 
    theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) + 
    theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,size=6))
  plot
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}