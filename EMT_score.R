#EMT_score
expr <- readr::read_rds("/data/TCGA/TCGA_data/pancan33_expr.rds.gz")
signature_list <- readr::read_csv("/home/shimw/signature.csv",col_names =T)

filter_makers_gene=function(.x){
    .x %>%
    dplyr::filter(symbol%in%signature_list$gene)
}

purrr::map(.x=expr$expr,filter_makers_gene) %>%
  tibble::tibble("cancer_types"=expr$cancer_types,"expr"=.) -> filter_gene_expr



M_list <- signature_list$gene[signature_list$E_VS_M=="M"]
E_list <- signature_list$gene[signature_list$E_VS_M=="E"]

M_score <- function(m){
  m=ifelse(m==0,0.0001,m)
  sum(log2(m))/55
}

E_score <- function(m){
  m=ifelse(m==0,0.0001,m)
  sum(log2(m))/22
}


EMT_score <- function(x){
  x %>%
  dplyr::filter(symbol%in%M_list) %>%
  dplyr::select(3:length(names(x))) %>%
  sapply(M_score) -> M
  
  x %>%
  dplyr::filter(symbol%in%E_list) %>%
  dplyr::select(3:length(names(x))) %>%
  sapply(E_score) -> E  
  
  tibble::tibble("barcode"=names(M-E),"EMT_score"=M-E)
}
 lapply(filter_gene_expr$expr,EMT_score) %>%
 tibble::tibble("cancer_type"=filter_gene_expr$cancer_types,"EMT_score"=.)->result
 readr::write_rds(result,path ="/home/shimw/TCGA/pancan33_EMT_score.rds.gz",compress = "gz")


