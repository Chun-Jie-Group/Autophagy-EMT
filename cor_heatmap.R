#autophagy and EMT
gene_list_path <- "/project/liucj/projects/6.autophagy/01_autophagy_gene_list/"
expr_path <- "/data/TCGA/TCGA_data/pancan33_expr.rds.gz"
EMT_score_path <- "/home/shimw/TCGA/pancan33_EMT_score.rds.gz"

autophagy_gene <- readr::read_rds(file.path(gene_list_path,"nature_gene_list.rds.gz"))
autophagy_gene <-autophagy_gene[-which(autophagy_gene$gene_symbol=="ATP6V1G3"),]
expr <- readr::read_rds(expr_path)
emt_score <-  readr::read_rds(EMT_score_path)

#filter autophage expr
filter_autophagy_gene=function(.x){
  .x %>%
    dplyr::filter(symbol%in%autophagy_gene$gene_symbol)
}

purrr::map(.x=expr$expr,filter_autophagy_gene) %>%
tibble::tibble("cancer_types"=expr$cancer_types,"expr"=.) -> filter_autophage_expr

#dataframe of autophagy gene expr and emt score

purrr::map2(emt_score$EMT_score,filter_autophage_expr$expr,function(x,y){

  t_y <- tibble::as.tibble(t(y[,-c(1,2)]))
  names(t_y) <- t(y)[1,]
  
  t_y %>%
  
    purrr::map(.,function(m){
      m = ifelse(m==0,0.0001,m)
      m = log2(m)
      cor.test(m,x$EMT_score,method = "pearson")%>%
        broom::tidy() %>% 
        dplyr::select(coef = estimate, pval = p.value) %>% 
        tibble::as_tibble()
        
    }) %>%
    dplyr::bind_rows()%>%
    tibble::add_column(.,"symbol"=names(t_y),.before = 1)
}) %>%
tibble::tibble("cancer_types"=expr$cancer_types,"autophage_emt"=.) -> autophage_emt_cor

#heatmap
##as.matrix
library(pheatmap)
purrr::map2(autophage_emt_cor$cancer_types,autophage_emt_cor$autophage_emt,function(x,y){
  y%>%
    tibble::add_column(.,"tumor_type"=x)
})%>%
  dplyr::bind_rows() %$%
  .[,-3]%>%
  tidyr::spread(tumor_type,coef) -> cor_matrix
  heatmap_matrix<- as.matrix(cor_matrix[,-1])
rownames(heatmap_matrix) <- cor_matrix$symbol

pheatmap(heatmap_matrix,col = bluered(100),fontsize_row =4,main = "Correlation of the expression of key genes in autophage with EMT score by tumor type",filename = "/home/shimw/TCGA/correlation_heatmap.pdf",width = 8.5,height = 11)
