
#-------------------------------------- 
# Libraries
#--------------------------------------
library(stringr)
library(ggplot2)
library(dplyr)
library(plyr)
library(DESeq2)
library(ggbeeswarm)

#-------------------------------------- 
# Fig 1D
# Imbalanced Host Response to SARS-CoV-2 Drives Development of COVID-19 (2020)
# GSE147507, NHBE data
#--------------------------------------
#-------------------------------------- 
# Data processing
#-------------------------------------- 
#read in count matrix
df <- as.matrix(read.table("./GSE147507_RawReadCounts_Human.tsv.gz", header=TRUE))

#make metadata df
#samp types = A549, A549.ACE2, NHBE, Calu3, COVID19Lung, HealthyLungBiopsy
#conditions = Mock, SARS.CoV.2, RSV, IAV, HPIV3, IAVdNS1, IFNB_4h, IFNB_6h, IFNB_12h, SARS.CoV.2_Rux
metadata_df <- data.frame(matrix(ncol=3,nrow=78, dimnames=list(NULL, c("samp_id", "samp_type", "condition"))))
metadata_df$samp_id <- colnames(df)
metadata_df$samp_type <- str_extract(metadata_df$samp_id,"A549[.]*[A-Z]*2?|NHBE|Calu3|Lung")
metadata_df$condition <- str_extract(metadata_df$samp_id,"Mock|SARS.CoV.2([_]*[Rux])*|RSV|IAV(dNS1)*|HPIV3|IAVdNS1|IFNB_4h|IFNB_6h|IFNB_12h|Healthy|COVID19")
metadata_df$series <- str_extract(metadata_df$samp_id,"Series[0-9]+")
metadata_df <- transform(metadata_df, row.names= samp_id, samp_id =NULL)

#slice for NHBE data
metadata_df <- metadata_df[which(metadata_df$samp_type %in% c('NHBE')),]
#mock and sars-cov-2 conditions were done as series 1 experiment and mock, IAV, IAVdNS1, ifnb done as series9 experiment
#as this creates two distinct mock conditions one at 24h and another at 12h for two independent experiments, we will take this into account
metadata_df$condition <- paste(metadata_df$condition, metadata_df$series, sep = '_')

#assign factor levels
metadata_df$condition <- factor(metadata_df$condition, levels = c('Mock_Series1','SARS.CoV.2_Series1',
                                                                  'Mock_Series9','IAV_Series9','IAVdNS1_Series9',
                                                                  'IFNB_4h_Series9','IFNB_6h_Series9','IFNB_12h_Series9'))
#rename for ease
metadata_df$condition <- mapvalues(metadata_df$condition, from = c('Mock_Series1','SARS.CoV.2_Series1',
                                          'Mock_Series9','IAV_Series9','IAVdNS1_Series9',
                                          'IFNB_4h_Series9','IFNB_6h_Series9','IFNB_12h_Series9'),
                                      to = c('Mock_24h','SARS.CoV.2','Mock_12h','IAV','IAVdNS1','IFNB_4h','IFNB_6h','IFNB_12h'))


#subset our count matrix to only have samples with metadata
df <- df[,which(colnames(df) %in% rownames(metadata_df))]

#-------------------------------------- 
# DESeq2
#-------------------------------------- 
# Deseq
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = metadata_df,
                              design = ~ condition)

# keep genes that have more than 10 total counts when all samples are summed
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#run Deseq
dds <- DESeq(dds)

#extract norm counts
ncount_df <- counts(dds, normalized=TRUE)

#merge normalized data with clinical data
df_master <- merge(data.frame(t(ncount_df)), metadata_df, by='row.names')

#-------------------------------------- 
# Fig. 1D - Bar plots for MHC class I genes
#-------------------------------------- 
source('bar_plotter.R')

interested <- c('NLRC5','HLA.A','HLA.B','HLA.C','TAP1', 'PSMB9')

for (i in interested){
  bar_plotter(df_master, dds = dds, desired_gene = i, save = 'T', savedirectory = './NHBE_interested_normcounts')
}
