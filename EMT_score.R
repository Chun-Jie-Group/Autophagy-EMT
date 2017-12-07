expr <- readr::read_rds("/data/TCGA/TCGA_data/pancan33_expr.rds.gz")
signature_list <- readr::read_csv("/home/shimw/signature.csv",col_names =T)

filter_makers_gene=function(.x){
    .x %>%
    dplyr::filter(symbol%in%signature_list$X1)
}

purrr::map(.x=expr$expr,filter_makers_gene) %>%
  tibble::tibble("cancer_types"=expr$cancer_types,"expr"=.) -> filter_gene_expr



M_list <- signature_list$X1[signature_list$`E vs M`=="M"]
E_list <- signature_list$X1[signature_list$`E vs M`=="E"]

M_score <- function(m){
  m=ifelse(m==0,0.0001,m)
  sum(log2(m))/55
}

E_score <- function(m){
  m=ifelse(m==0,0.0001,m)
  sum(log2(m))/22
}

x=filter_gene_expr$expr[[1]]

EMT_score <- function(x){
  x %>%
  dplyr::filter(symbol%in%M_list) %>%
  dplyr::select(3:length(names(x))) %>%
  sapply(M_score) -> M
  
  x %>%
  dplyr::filter(symbol%in%E_list) %>%
  dplyr::select(3:length(names(x))) %>%
  sapply(E_score) -> E  
  
  tibble::tibble("barcode"=names(M-E),"EMT score"=M-E)
}
 lapply(filter_gene_expr$expr,EMT_score) %>%
 tibble::tibble("cancer type"=filter_gene_expr$cancer_types,"EMT score"=.)->result
 readr::write_rds(result,path ="/home/shimw/TCGA/pancan33_EMT_score.rds.gz",compress = "gz")


