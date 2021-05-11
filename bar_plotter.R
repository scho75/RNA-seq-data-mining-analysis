#plotter
bar_plotter <- function(df, dds, desired_gene, save = 'F', savedirectory = './'){
  
  #edit condition to drop early IFNbeta timepoints
  df$condition <- as.character(df$condition)
  df <- df[which(!df$condition %in% c('IFNB_4h','IFNB_6h')),]
  df$condition <- factor(df$condition, levels = c('Mock_24h','SARS.CoV.2','Mock_12h','IAV','IAVdNS1','IFNB_12h'))
  
  #edit input to match res
  desired_gene <- gsub('\\.','-',desired_gene)
  #get fold change and pval info
  res <- results(dds,contrast=c("condition","SARS.CoV.2","Mock_24h"))
  sars_pval <- formatC(res[rownames(res) == desired_gene, 'padj'], format = 'e', digits = 2)
  sars_fold_change <- as.character(round(2^(res[rownames(res) == desired_gene, 'log2FoldChange']),4))
  res <- results(dds,contrast=c("condition","IAV","Mock_12h"))
  iav_pval <- formatC(res[rownames(res) == desired_gene, 'padj'], format = 'e', digits = 2)
  iav_fold_change <- as.character(round(2^(res[rownames(res) == desired_gene, 'log2FoldChange']),4))
  res <- results(dds,contrast=c("condition","IAVdNS1","Mock_12h"))
  iavdns1_pval <- formatC(res[rownames(res) == desired_gene, 'padj'], format = 'e', digits = 2)
  iavdns1_fold_change <- as.character(round(2^(res[rownames(res) == desired_gene, 'log2FoldChange']),4))
  
  #assign x axis label
  x_label <- paste('SARS-CoV-2 vs Mock_24h: FC = ', sars_fold_change, '\n', '   p = ', sars_pval, '\n',
                   'IAV vs Mock_12h: FC = ', iav_fold_change, '\n', '   p = ', iav_pval, '\n',
                   'IAVdNS1 vs Mock_12h: FC = ', iavdns1_fold_change, '\n', '   p = ', iavdns1_pval, sep="")

  
  #revert input edit
  desired_gene <- gsub('\\-','.',desired_gene)
  
  #zscore transformation to reference i.e. Mock_24h for SARS-CoV-2 samples and Mock_12h for IAV, IAVdNS1 and IFNb
  ref_mean_24h <- mean(df[which(df$condition == 'Mock_24h'),desired_gene])
  ref_sd_24h <- sd(df[which(df$condition == 'Mock_24h'),desired_gene])
  ref_mean_12h <- mean(df[which(df$condition == 'Mock_12h'),desired_gene])
  ref_sd_12h <- sd(df[which(df$condition == 'Mock_12h'),desired_gene])

  sars_exp <- (df[1:6,desired_gene] - ref_mean_24h) / (ref_sd_24h)
  iav_exp <- (df[7:20,desired_gene] - ref_mean_12h) / (ref_sd_12h)

  #y <- c(sars_exp,iav_exp)
  y <- df[,desired_gene]
  
  #plotting
  gg <- ggplot(df, aes(x = condition, y = y, color = condition, fill = condition)) +
    geom_bar(stat='summary', fun = 'mean', alpha = 0.4, color = 'black') +
    geom_beeswarm(size=4,cex=2) +
    labs(title= desired_gene, x= x_label, y= paste('Gene expression', '(normalized count)', sep='\n'))
  gg <- gg + scale_color_manual(values=c('#16A085','#C0392B', '#62C01D', '#FF8C00','#AE19F1','#0000FF')) +
    scale_fill_manual(values=c('#16A085','#C0392B','#62C01D', '#FF8C00','#AE19F1','#0000FF'))
  gg <- gg + theme_bw() + theme(legend.position = 'none', plot.title = element_text(hjust=0.5),
                                panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size=1),
                                axis.ticks = element_line(size=1), axis.text.x = element_text(angle = 45, hjust=1))
  
  #save options
  if (save == 'F'){
    gg
  } else if (save == 'T'){
    savedirectory = savedirectory #where you want to save, defaults to current folder
    fname = paste(desired_gene,".png",sep="")
    print(paste0("Saving...",fname))
    ggsave(path = savedirectory, 
           plot = gg, width = 4, height = 5, dpi = 300, filename = fname)
  }
  
}
