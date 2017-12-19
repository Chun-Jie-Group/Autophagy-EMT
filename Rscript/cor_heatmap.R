#autophagy and EMT
gene_list_path <- "/project/liucj/projects/6.autophagy/01_autophagy_gene_list/"
expr_path <- "/data/TCGA/TCGA_data/pancan33_expr.rds.gz"
EMT_score_path <- "/home/shimw/TCGA/pancan33_EMT_score.rds.gz"

autophagy_gene <- readr::read_rds(file.path(gene_list_path,"nature_gene_list.rds.gz"))
autophagy_gene <-autophagy_gene[-which(autophagy_gene$gene_symbol=="ATP6V1G3"),]
expr <- readr::read_rds(expr_path)
emt_score <-  readr::read_rds(EMT_score_path)

#filter autophagy expr
filter_autophagy_gene=function(.x){
  .x %>%
    dplyr::filter(symbol%in%autophagy_gene$gene_symbol)
}

purrr::map(.x=expr$expr,filter_autophagy_gene) %>%
tibble::tibble("cancer_types"=expr$cancer_types,"expr"=.) -> filter_autophagy_expr
readr::write_rds(filter_autophagy_expr,path ="/home/shimw/TCGA/pancan33_filter_autophagy_expr.rds.gz",compress = "gz")


#dataframe of autophagy gene expr and emt score

purrr::map2(emt_score$EMT_score,filter_autophagy_expr$expr,function(x,y){

 # t_y <- tibble::as.tibble(t(y[,-c(1,2)]))
  #names(t_y) <- t(y)[1,]
  
#  t_y %>%
  dplyr::select(y, -entrez_id) %>%
    as.data.frame() %>%
    tibble::column_to_rownames(.,var = "symbol")%>%
    t() %>% tibble::as.tibble()%>%
    
  
    purrr::map(.,function(m){
      m = ifelse(m==0,0.0001,m)
      m = log2(m)
      cor.test(m,x$EMT_score,method = "pearson")%>%
        broom::tidy() %>% 
        dplyr::select(coef = estimate, pval = p.value) %>% 
        tibble::as_tibble()
        
    }) %>%
    dplyr::bind_rows()%>%
    tibble::add_column(.,"symbol"=y$symbol,.before = 1)
}) %>%
tibble::tibble("cancer_types"=expr$cancer_types,"autophagy_emt"=.) -> autophagy_emt_cor

#heatmap
##as.matrix
library(pheatmap)
purrr::map2(autophagy_emt_cor$cancer_types,autophagy_emt_cor$autophagy_emt,function(x,y){
  y%>%
    tibble::add_column(.,"tumor_type"=x)
})%>%
  dplyr::bind_rows() %$%
  .[,-3]%>%
  tidyr::spread(tumor_type,coef) -> cor_matrix
  heatmap_matrix<- as.matrix(cor_matrix[,-1])
rownames(heatmap_matrix) <- cor_matrix$symbol



#red <- factoextra::hcut(t(heatmap_matrix), k = 2, hc_func = 'hclust', hc_method = 'ward.D', hc_metric = 'pearson', stand = T)


#res <- factoextra::hcut(heatmap_matrix, k = 2, hc_func = 'hclust', hc_method = 'ward.D', hc_metric = 'pearson', stand = T)
heatmap_matrix[res$order,red$order]

 Heatmap(
   heatmap_matrix[res$order,red$order],
   col = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red"), space = "RGB"),
   name="Correlation",
   #top_annotation  = ha,
   show_row_names = T, 
   show_column_names = T,
   cluster_columns = FALSE,
   cluster_rows = FALSE,
  # clustering_distance_columns = "pearson",
  # clustering_method_columns = "ward.D",
  # clustering_distance_rows = "pearson",
  # clustering_method_rows = "ward.D",
  # row_names_gp = gpar(fontsize = 5)
 )


#pheatmap(heatmap_matrix,col = bluered(100),fontsize_row =4,main = "Correlation of the expression of key genes in autophagy with EMT score by tumor type",filename = "/home/shimw/TCGA/correlation_heatmap.pdf",width = 8.5,height = 11)


