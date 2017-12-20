library(magrittr)
gene_list_path <- "/project/liucj/projects/6.autophagy/01_autophagy_gene_list/"
expr_path <- "/data/TCGA/TCGA_data/pancan33_expr.rds.gz"
EMT_score_path <- "/home/shimw/TCGA/pancan33_EMT_score.rds.gz"



autophagy_gene <- readr::read_rds(file.path(gene_list_path,"nature_gene_list.rds.gz"))
#one gene's expr is 0,so I get rid of it
autophagy_gene <-autophagy_gene[-which(autophagy_gene$gene_symbol=="ATP6V1G3"),]
expr <- readr::read_rds(expr_path)
emt_score <-  readr::read_rds(EMT_score_path)

#filter autophagy expr
filter_autophagy_gene=function(.x){
  .x %>%
    dplyr::filter(symbol%in%autophagy_gene$gene_symbol)%>%
    dplyr::select(symbol,entrez_id,matches("TCGA-.{2}-.{4}-0[1-9].-.{3}-.{4}-.{2}"))
  #get the tumor expr
  
}

purrr::map(.x=expr$expr,filter_autophagy_gene) %>%
  tibble::tibble("cancer_types"=expr$cancer_types,"expr"=.) -> filter_autophagy_expr
readr::write_rds(filter_autophagy_expr,path ="/home/shimw/TCGA/pancan33_filter_autophagy_expr.rds.gz",compress = "gz")

#dataframe of autophagy gene expr and emt score


purrr::map2(emt_score$EMT_score,filter_autophagy_expr$expr,function(x,y){

  dplyr::select(y, -entrez_id) %>%
    as.data.frame() %>%
    tibble::column_to_rownames(.,var = "symbol")%>%
    t() %>% tibble::as.tibble()%>%
    
    
    purrr::map(.,function(m){
      m = ifelse(m==0,0.0001,m)
      m = log2(m)
      cor.test(m,x$EMT_score,method = "pearson")%>%
        broom::tidy() %>% 
        dplyr::select(coef = estimate, p.value = p.value) %>% 
        tibble::as_tibble()%>%
        dplyr::mutate(fdr = p.adjust(p.value, method = "fdr"))
      
    }) %>%
    dplyr::bind_rows()%>%
    tibble::add_column(.,"symbol"=y$symbol,.before = 1)
    
  
})%>%
  #create a dataframe for ggplot
  purrr::map2(.,filter_autophagy_expr$cancer_types,function(x,y){
    x %>%
      tibble::add_column(.,"tumor_type"=y,.before = 1)
    
  })%>%
   dplyr::bind_rows() %>%
  #filter p value and fdr
   dplyr::mutate(p.value = -log10(p.value)) %>% 
   dplyr::mutate(p.value = ifelse(p.value > 15, 15, p.value)) %>% 
   dplyr::filter(fdr <= 0.05) %>%
   dplyr::mutate(fdr = -log10(fdr)) %>% 
   
   dplyr::mutate(fdr = ifelse(fdr > 15, 15, fdr)) -> autophagy_emt_cor

#order the gene by median
#  tidyr::spread(autophagy_emt_cor[,-4],tumor_type,coef) %>%
#    as.data.frame() %>% 
#    tibble::remove_rownames()%>%
#    tibble::column_to_rownames(.,var = "symbol") %>%
#    as.matrix() %>%
#    apply(1, median) %>%
#    sort()%>%names() -> gene_rank

#order the gene by positive number
autophagy_emt_cor%>%
  dplyr::select(tumor_type,symbol,coef)%>%
  tidyr::spread(tumor_type,coef)%>%
  dplyr::rowwise() %>%
  dplyr::do(
    symbol = .$symbol,
    rank =  unlist(.[-1], use.names = F) %>% sum(na.rm = TRUE),
  ) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::arrange(rank) -> gene_rank
  
#order the tunor type 
autophagy_emt_cor %>%
  dplyr::select(tumor_type,symbol,coef)%>%
  group_by(tumor_type)%>%
  dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum(abs(.)))) %>%
  dplyr::arrange(dplyr::desc(coef)) -> tumor_rank
  




  
  CPCOLS <-c ("#FF0000", "#FFFFFF", "#00EE00")
  
  autophagy_emt_cor%>%
    ggplot(aes(x = tumor_type, y = symbol)) + 
    geom_point(aes(size = fdr, col = coef)) +
    scale_color_gradient2(
      low = CPCOLS[3],
      mid = CPCOLS[2],
      high = CPCOLS[1],
      midpoint = 0,
      na.value = "white",
      breaks=seq(-1, 1,length.out = 9),
      limit = c(-1,1),

      name = "Correlation"
    ) +
    scale_size_continuous(
      limit = c(-log10(0.05),15),
      #range = c(0.5,5),
      #breaks =  c(-log10(0.05), 2, 3, 4, 5),
      range = c(1, 4),
      breaks = c(-log10(0.05), 5, 10, 15),
      labels = c("0.05", latex2exp::TeX("$10^{-5}$"), latex2exp::TeX("$10^{-10}$"), latex2exp::TeX("$< 10^{-15}$")),
     #labels = c("0.05" ,latex2exp::TeX("$10^{-2}$"), latex2exp::TeX("$10^{-3}$"), latex2exp::TeX("$ 10^{-4}$"),latex2exp::TeX("$ < 10^{-5}$")),
      name = "FDR"
    )+
    scale_y_discrete(limit = gene_rank$symbol)+
    scale_x_discrete(limit = tumor_rank$tumor_type) +
    theme(
      panel.background = element_rect(colour = "black", fill = "white"),
      panel.grid = element_line(colour = "grey", linetype = "dashed"),
      panel.grid.major = element_line(
        colour = "grey",
        linetype = "dashed",
        size = 0.2
      ),
      axis.title = element_blank(),
      axis.ticks = element_line(color = "black"),
      
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 6),
      legend.text = element_text(size = 8),
      legend.text.align = 0.5,
      legend.title = element_text(size = 10),
      legend.key = element_rect(fill = "white", colour = "black")
    ) ->p
    
  ggsave(
    filename = "point_plot.pdf",
    plot = p,
    device = "pdf",
    width = 8.5,
    height = 11,
    path = "/home/shimw/github/EMT/pdf"
  )
 
  
  


  
  
  