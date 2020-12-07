violin_plotter <- function(df, dds, desired_parameter, save = 'F', savedirectory = './'){
  # Processing
  # change colnames
  colnames(df) <- gsub('-','.',colnames(df))
  # edit inputted colnames
  desired_parameter <- gsub('-','.',desired_parameter)
  
  if (desired_parameter %in% c('viral_load','age')){
    df <- df[which(df$immune_influx != 'negative'),] #omit SARS-CoV-2 negative patients
    df <- df[which(!is.na(df[,desired_parameter])),] #drop NA records
    y <- desired_parameter
    pval <- round(wilcox.test(df[,desired_parameter]~df[,"immune_influx"])$p.value, digits = 5)
    fold_change <- round(mean(df[which(df$immune_influx == 'high'),desired_parameter]) / mean(df[which(df$immune_influx == 'low'),desired_parameter]), digits = 5)
    color_choice <- c('#FFD700','#FF0000')
    if(desired_parameter == 'age'){
      y_label <- 'Age (years)'
    } else if (desired_parameter == 'viral_load'){
      y_label <- 'N1 CT value'
    }
  } else { #for when desired_parameter is a gene
    df[,desired_parameter] <- log2(df[,desired_parameter] + 1) #log2 transform
    df_ref <- df[which(df$infection_status == 'negative'),] #set reference
    df <- df[which(df$immune_influx != 'high'),] #omit SARS-CoV-2 positive high immune influx patients
    res <- results(dds,contrast=c("immune_influx","low","negative"))
    rownames(res) <- gsub('-','.',rownames(res)) #change res rownames to match
    pval <- formatC(res[rownames(res) == desired_parameter, 'padj'], format = 'e', digits = 2) #format pval
    fold_change <- as.character(round(2^(res[rownames(res) == desired_parameter, 'log2FoldChange']),4))
    y_label <- paste("Gene expression", "(Z-score normalized count)", sep = '\n')
    color_choice <- c('#16A085','#C0392B')
    
    #z-score transformation
    ref_mean <- mean(df_ref[,desired_parameter])
    ref_sd <- sd(df_ref[,desired_parameter])
    
    y <- (df[,desired_parameter] - ref_mean) / (ref_sd)
  }
  

  # Plotting
  dp <- ggplot(df, aes_string(x="immune_influx", y=y, fill="immune_influx")) + 
    geom_violin(trim=FALSE, alpha = 0.8)+
    geom_boxplot(width=0.2, fill="white")+
    labs(title=desired_parameter,x= paste('', '\n','p = ',pval,'\n', 'FC = ',fold_change, sep=""), y = y_label)
  dp <- dp + scale_fill_manual(values=color_choice)
  dp <- dp + theme_classic() + theme(legend.position = 'none', axis.line = element_line(colour = "black", size=1),
                                     axis.ticks = element_line(size=1))
  
  # Save options
  if (save == 'F'){
    dp
  } else if (save == 'T'){
    savedirectory = './' #where you want to save, defaults to current folder
    fname = paste(desired_parameter,".png",sep="")
    print(paste0("Saving...",fname))
    ggsave(path = savedirectory,
           plot = dp, width = 1.5, height = 3, dpi = 300, filename = fname)
  }
  
  
}