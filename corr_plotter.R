#function for corr plotting
corr_plotter <- function(df,x,y, save = 'F', transform = 'T', savedirectory = './'){
  #change colnames
  colnames(df) <- gsub('-','_',colnames(df))
  #edit inputted colnames
  y <- gsub('-','_',y)
  
  #slice for what we need
  df <- df[,c(x,y)]
  
  if (transform == 'T'){
    df <- log2(df + 1) #log2 transform
    
    #z-score transform input
    df[,x] <- scale(df[,x])
    df[,y] <- scale(df[,y])
  } else {
    #z-score transform input
    df[,x] <- scale(df[,x])
    df[,y] <- scale(df[,y])
  }
  
  
  # Plotting
  gg <- ggplot(df, aes(x=df[,x], y=df[,y])) + geom_jitter(aes(), width = 0.01) + geom_smooth(method=lm, se=FALSE) +
    xlab(x) + ylab(y)
  gg <- gg + theme(plot.title = element_blank()) + stat_cor(method = 'spearman')

  # Save options
  if (save == 'F'){
    gg
  } else if (save == 'T'){
    savedirectory = './' #where you want to save, defaults to current folder
    fname = paste(x,"_",y,".png",sep="")
    print(paste0("Saving...",fname))
    ggsave(path = savedirectory, 
           plot = gg, width = 4, height = 4, dpi = 300, filename = fname)
  }
  
}
