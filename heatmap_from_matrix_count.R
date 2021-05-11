heatmap_from_matrix_count <- function(matrix_count, df, gene_list){
  #takes matrix count input as patients (cols) x genes (rows)
  #IMPORTANT: matrix count col patients are in same order as row patients of df
  
  #slice for genes we need
  matrix_count <- matrix_count[rownames(matrix_count) %in% c(gene_list),]
  #transpose
  matrix_count <- as.matrix(t(matrix_count))
  #transform to z-score of log2 norm count
  matrix_count <- log2(matrix_count + 1)
  matrix_count <- scale(matrix_count)
  #transpose back to use in heatmap
  matrix_count <- as.matrix(t(matrix_count))
  
  anno_df <- df[,c('infection_status','age_groups','gender','viral_load_groupings'), drop=FALSE]
  mat_colors <- list(infection_status = c('negative' = '#32CD32', 'positive' = '#FF0000'),
                     age_groups = c('Unknown' = '#000000','0-29' = '#D8BFD8','30-59' = '#BA55D3','60-100' = '#800080'),
                     gender = c('M' = '#1E90FF','F' = '#F08080','not collected' = '#000000'),
                     viral_load_groupings = c('negative' = '#32CD32','low' = '#FFD700','mid' = '#FF8C00','high' = '#FF0000'))
  
  #plot
  ht <- Heatmap(matrix_count, name = 'z-score', na_col='black', 
                column_split = df$infection_status,
                column_names_gp = gpar(fontsize = 0.5),
                row_names_gp = gpar(fontsize = 10),
                top_annotation = HeatmapAnnotation(df = anno_df,
                                                   col = mat_colors
                ),
                cluster_rows = FALSE,
                row_order = gene_list,
                cluster_column_slices = FALSE,
                column_gap = unit(2,"mm"),
                border = FALSE,
                column_title = "Patient IDs",
                row_title = "Genes")
  draw(ht)
  
}
