#Processing and figures for Yoo et al., Cell Host & Microbe 2021 

#-------------------------------------- 
# Libraries
#-------------------------------------- 
library(ggplot2)
library(dplyr)
library(DESeq2)
library(ComplexHeatmap)
library(ggpubr)

#-------------------------------------- 
# Fig 1A-C, Fig S1A
# In vivo antiviral host transcriptional response to SARS-CoV-2 by viral load, sex and age (2020)
# GSE152075, NP swab data
#--------------------------------------
#-------------------------------------- 
# Data processing
#-------------------------------------- 

#read in count matrix and metadata csv derived from "GSE152075_series_matrix.txt" file
df <- as.matrix(read.table("./GSE152075_raw_counts_GEO.txt.gz", header=TRUE))
metadata <- read.table("./GSE152075_metadata.csv", header=TRUE, sep = ',')

#make metadata for SARS-CoV-2 infection status
df_clinical <- data.frame(matrix(ncol=2,nrow=484, dimnames=list(NULL, c("samp_id", "infection_status"))))
df_clinical$samp_id <- colnames(df)
df_clinical$infection_status <- ifelse(grepl("POS",df_clinical$samp_id) == T, 'positive','negative')
df_clinical <- merge(df_clinical, metadata, by = 'samp_id')
df_clinical <- transform(df_clinical, row.names= samp_id, samp_id =NULL)
df_clinical <- df_clinical[which(df_clinical$viral_load != 'Unknown'),] #drop records with unknown viral load
#convert to numeric
df_clinical$viral_load <- as.numeric(as.character(df_clinical$viral_load))
#add viral load groupings
df_clinical$viral_load_groupings <- ifelse(is.na(df_clinical$viral_load), 'negative', 
                                           ifelse( df_clinical$viral_load <= 19, 'high',
                                                   ifelse( df_clinical$viral_load <= 24, 'mid', 'low')))
df_clinical$viral_load_groupings <- factor(df_clinical$viral_load_groupings, levels = c('negative','low','mid','high'))

#add age groups
df_clinical$age <- as.numeric(as.character(df_clinical$age))
df_clinical$age_groups <- ifelse(is.na(df_clinical$age), 'Unknown',
                                 ifelse(df_clinical$age <= 29, '0-29',
                                        ifelse(df_clinical$age <= 59, '30-59','60-100')))

# Make sure rownames of metadata is in the same order as count data matrix columns
df <- df[,rownames(df_clinical)]

#-------------------------------------- 
# DESeq2
#-------------------------------------- 
# seq batch was not included in design because all SARS-CoV-2 negative patients i.e. seq batch C, R, S, T, U are exclusively only in negative patients
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = df_clinical,
                              design = ~ infection_status)

# keep genes that have more than 10 total counts when all samples are summed
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# run
dds <- DESeq(dds)

#extract norm counts
ncount_df <- counts(dds, normalized=TRUE)
#merge normalized data with clinical data
df_master <- merge(data.frame(t(ncount_df)),df_clinical, by='row.names')
# define immune influx parameter based on PTPRC aka CD45 median expression
df_master$immune_influx <- ifelse(df_master$infection_status == 'negative', 'negative',
                                  ifelse(df_master$PTPRC < median(df_master$PTPRC), 'low','high'))
df_master$immune_influx <- factor(df_master$immune_influx, levels = c('negative','low','high'))

#-------------------------------------- 
# Add in immune influx parameter to dds
#-------------------------------------- 
all(rownames(colData(dds)) == df_master$Row.names) #check order of rows are the same
colData(dds) <- cbind(colData(dds)[1:8], df_master$immune_influx)
names(colData(dds))[9] <- 'immune_influx'

design(dds) <- formula(~ immune_influx)
dds <- DESeq(dds)

#extract norm counts
ncount_df <- counts(dds, normalized=TRUE)
#merge normalized data with metadata
df_master <- merge(data.frame(t(ncount_df)), data.frame(colData(dds)), by='row.names')

#-------------------------------------- 
# Fig. 1A - Heatmapping with ComplexHeatmap
#-------------------------------------- 
source('heatmap_from_matrix_count.R')
heatmap_mhc <- c('NLRC5','B2M','HLA-A','HLA-B','HLA-C','HLA-E','HLA-F','TAP1', 'PSMB9','PTPRC')

heatmap_from_matrix_count(matrix_count = ncount_df, 
                          df = df_master, 
                          gene_list = heatmap_mhc)

#-------------------------------------- 
# Fig. 1B - Violin plots for age and viral load
#-------------------------------------- 
source('violin_plotter.R')
violin_plotter(df_master, dds = dds, desired_parameter = 'age', save = 'T', savedirectory = './')
violin_plotter(df_master, dds = dds, desired_parameter = 'viral_load',  save = 'T', savedirectory = './')

#-------------------------------------- 
# Fig. 1C - Violin plotting for MHC class I genes
#-------------------------------------- 
fig1c_mhc <- c('NLRC5','HLA.A','HLA.B','HLA.C','TAP1', 'PSMB9')

for (i in fig1c_mhc){
  violin_plotter(df_master, dds = dds, desired_parameter = i, save = 'T', savedirectory = './')
}

#-------------------------------------- 
# Fig. S1A - Correlation plots
#--------------------------------------
source('corr_plotter.R')

figs1a_mhc <- c('NLRC5','B2M','HLA.A','HLA.B','HLA.C','TAP1', 'PSMB9','STAT1','IRF1')
for (i in figs1a_mhc){
  corr_plotter(df_master, x = 'PTPRC', i, save = 'T', savedirectory = './')
}


