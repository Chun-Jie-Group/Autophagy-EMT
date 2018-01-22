# Library -----------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)

options(digits=12)
# Path --------------------------------------------------------------------
path<-"/home/shimw/TCGA/"
pdf_path<- "/home/shimw/github/EMT/pdf/"
color_path <- "/data/TCGA/TCGA_data/PanCanAtlasTumors_color_coded_by_organ_system_20170302.tsv"


# Load data ---------------------------------------------------------------
pancan_color <- readr::read_tsv(color_path)
filter_autophagy_expr <- readr::read_rds(file.path(path,"pancan33_filter_autophagy_expr.rds.gz"))
EMT_sco <- readr::read_rds(file.path(path,"pancan33_EMT_score.rds.gz"))
pcc <- pancan_color %>% dplyr::pull(color)
names(pcc) <- pancan_color %>% dplyr::pull(cancer_types)

# calculate 
purrr::map2(filter_autophagy_expr$expr,EMT_sco$EMT_score,function(x,y){

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
  tibble::tibble("cancer_type"=filter_autophagy_expr$cancer_types,"expr"=.) -> filter_autophagy_scale


# all cancer in a table ----------------------------------------------
purrr::map2(filter_autophagy_scale$cancer_type,filter_autophagy_scale$expr,function(x,y){
  tibble::tibble("cancer_type"=x,"EMT_score"=as.double(colnames(y[,-1])))
})%>%
  dplyr::bind_rows()->cancer_type_sort
cancer_type_sort=dplyr::arrange(cancer_type_sort,EMT_score)


purrr::map(filter_autophagy_scale$expr,function(x){
  x[,-1]
}) %>%dplyr::bind_cols()%>%
  tibble::add_column(.,"symbol"=filter_autophagy_scale$expr[[1]]$symbol,.before = 1)->all_autophagy_scale

all_autophagy_scale%>%
  as.data.frame()%>%
  tibble::column_to_rownames(var = "symbol") %$%
  .[,order(as.numeric(names(all_autophagy_scale[,-1])))]%>% 
  data.matrix() ->all_autophagy_scale_matrix


pdf(file.path(pdf_path, "all_expr_heatmap.pdf"), width=11, height = 10)

  ha=HeatmapAnnotation(
    df = data.frame(EMT_score=cancer_type_sort$EMT_score,Cancer=cancer_type_sort$cancer_type),
    annotation_legend_param = list(EMT_score=list(title = "EMT score")),
    col = list(EMT_score=colorRamp2(c(-4, 0, 4), c("blue",  "white", "yellow"), space = "RGB"),Cancer=pcc)
  )
  all_autophagy_scale_matrix%>%
    Heatmap(
      .,
      col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red"), space = "RGB"),
      name="Expression",
      top_annotation = ha,
      show_row_names = F, 
      show_column_names = FALSE, 
      show_row_dend=F,
      show_column_dend=F,
      clustering_distance_rows = "pearson",
      clustering_method_rows = "ward.D",
      cluster_columns = FALSE,
      cluster_rows = TRUE,
      row_names_gp = gpar(fontsize = 5),
      width = 8,
      row_title = "mRNA Expression",
      row_title_rot = 90,
      show_heatmap_legend = F
    ) ->ht
  draw(ht)


dev.off()

