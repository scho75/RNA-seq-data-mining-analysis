#Processing and figures for Yoo et al., Nat Commun 2021 

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

#-------------------------------------- 
# Fig. S1A - Correlation plots
#--------------------------------------
source('corr_plotter.R')

figs1a_mhc <- c('NLRC5','B2M','HLA.A','HLA.B','HLA.C','TAP1', 'PSMB9','STAT1','IRF1')
for (i in figs1a_mhc){
  corr_plotter(df_master, x = 'PTPRC', i, save = 'T', savedirectory = './figs1a')
}

#-------------------------------------- 
# Use infection status negative patients' CD45 expression distribution as ref and exclude patients in infected group with high CD45
#-------------------------------------- 
q <- quantile(df_master[which(df_master$infection_status == 'negative'),'PTPRC'], probs=c(.25,.75), na.rm=FALSE)
iqr <- IQR(df_master[which(df_master$infection_status == 'negative'),'PTPRC'])
upper_bound <- q[2] + 1.5*iqr #129
lower_bound <- q[1] - 1.5*iqr #-65.2

pat_ids <- df_master[!(df_master$PTPRC > upper_bound), 'Row.names'] #retrieve patient ids that fit criteria
#edit original df
df <- df[,pat_ids]
#edit clinical df
df_clinical <- df_clinical[which(rownames(df_clinical) %in% pat_ids),]
df <- df[,rownames(df_clinical)]

#dds
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

#-------------------------------------- 
# Fig. 1A - Heatmapping with ComplexHeatmap
#-------------------------------------- 
source('heatmap_from_matrix_count.R')
heatmap_mhc <- c('NLRC5','B2M','HLA-A','HLA-B','HLA-C','HLA-E','HLA-F','TAP1', 'PSMB9')

heatmap_from_matrix_count(matrix_count = ncount_df, 
                          df = df_master, 
                          gene_list = heatmap_mhc)

#-------------------------------------- 
# Fig. 1B - Violin plots for age and CD45
#-------------------------------------- 
source('violin_plotter.R')
violin_plotter(df_master, dds = dds, desired_parameter = 'age', save = 'T', savedirectory = './fig1b')
violin_plotter(df_master, dds = dds, desired_parameter = 'PTPRC', save = 'T', savedirectory = './fig1b')
#gender %s
table(df_master[which(df_master$infection_status == 'negative'),'gender'])
table(df_master[which(df_master$infection_status == 'positive' & df_master$gender != 'not collected'),'gender'])
#median age values
median(df_master[which(df_master$infection_status == 'negative'),'age'])
median(df_master[which(df_master$infection_status == 'positive' & df_master$age != 'NA'),'age'])

#-------------------------------------- 
# Fig. 1C - Violin plotting for MHC class I genes
#-------------------------------------- 
fig1c_mhc <- c('NLRC5','B2M','HLA.A','HLA.B','HLA.C','TAP1', 'PSMB9')

for (i in fig1c_mhc){
  violin_plotter(df_master, dds = dds, desired_parameter = i, save = 'T', savedirectory = './fig1c')
}
