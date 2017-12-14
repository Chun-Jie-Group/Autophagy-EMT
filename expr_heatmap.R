path<-"/home/shimw/TCGA/"
pdf_path<- "/home/shimw/TCGA/heatmap/"
filter_autophage_expr <- readr::read_rds(file.path(path,"pancan33_filter_autophage_expr.rds.gz"))
EMT_sco <- readr::read_rds(file.path(path,"pancan33_EMT_score.rds.gz"))

purrr::map2(filter_autophage_expr$expr,EMT_sco$EMT_score,function(x,y){

  dplyr::select(x, -symbol, -entrez_id) %>%
 purrr::map(function(m){
   ifelse(m==0,0.0001,m)
 }) %>%
    dplyr::bind_rows()%>%
    log2() %>%
  t() %>%
    scale(center = T,scale = T) %>% t() %>%
    tibble::as.tibble() %>%
    tibble::add_column(.,"symbol"=x$symbol,.before = 1) -> temp
  names(temp)[-1] = y$EMT_score
    return(temp)
}) %>%
  tibble::tibble("cancer_type"=filter_autophage_expr$cancer_types,"expr"=.) -> filter_autophage_scale

library(ComplexHeatmap)
library(circlize)
dev.off()
pdf(file.path(pdf_path, "all_expr_heatmap.pdf"), width=8.5, height = 11)
purrr::walk2(filter_autophage_cpm$expr,filter_autophage_cpm$cancer_type,function(x,y){
  ha=HeatmapAnnotation(
    df = data.frame(EMT_score=sort(as.numeric(colnames(x[,-1])))),
  
    col = list(EMT_score=colorRamp2(c(-4, 0, 4), c("blue", "white", "red"), space = "RGB"))
  )
  
  as.data.frame(x) %>%
  tibble::column_to_rownames(var = "symbol") %$% 
    .[,sort(names(x[,-1]))]%>%
    data.matrix() %>%
    Heatmap(
      .,
      col=colorRamp2(c(-2, 0, 2), c("blue", "white", "yellow"), space = "RGB"),
      name="Expression",
      top_annotation = ha,
      column_title = y,
      column_title_side = "bottom",
      show_row_names = T, 
      show_column_names = FALSE, 
      show_row_dend=F,
      show_column_dend=F,
      cluster_columns = FALSE,
      cluster_rows = FALSE,
      row_names_gp = gpar(fontsize = 5),
      width = 8
      
    ) ->ht
  draw(ht)

})
 
dev.off()


