# Library -----------------------------------------------------------------
library(magrittr)
library(GSVA)
library(ggplot2)
library(igraph)
library(ggraph)

# Path --------------------------------------------------------------------
gene_list_path <- "/project/liucj/projects/6.autophagy/01_autophagy_gene_list/"
expr_path <- "/data/TCGA/TCGA_data/pancan33_expr.rds.gz"
EMT_score_path <- "/home/shimw/TCGA/pancan33_EMT_score.rds.gz"


# Load data ---------------------------------------------------------------
expr <- readr::read_rds(expr_path)
emt_score <-  readr::read_rds(EMT_score_path)
autophagy_gene <- readr::read_rds(file.path(gene_list_path,"nature_gene_list.rds.gz"))


# Calculate gsva for every sample -----------------------------------------
gene_sets <- list(atg_lys = autophagy_gene$gene_symbol)
fn_gsva <- function(.y, gene_sets = gene_sets){
  .y %>% 
    
    tidyr::drop_na()%>%
    dplyr::select( -entrez_id) %>%
    dplyr::select(symbol,matches("TCGA-.{2}-.{4}-0[1-9].-.{3}-.{4}-.{2}"))-> .d
  

  
  .d_mat <- as.matrix(.d[,-1])
  rownames(.d_mat) <- .d$symbol
  
   # as.data.frame() %>%
   # tibble::column_to_rownames(var = "symbol")%>%
   # as.matrix() -> .d_mat

    
  .es_dif <- gsva(.d_mat, gset.idx.list = gene_sets, method = "gsva", mx.diff = TRUE, verbose = FALSE, parallel.sz = 1)
  
  .es_dif %>% 
    as.matrix() %>% t()%>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "barcode")%>%
    tibble::as_tibble() ->.d_es
     
}

emt_score%>% 
  dplyr::mutate(
    gsva = purrr::map(
      .x = expr$expr,
      .f = function(.x) {
        fn_gsva(.x, gene_sets = gene_sets)
      }
    )
  ) -> emt_gsva

dplyr::mutate(emt_gsva[,1],
emt_with_gsva=purrr::map2(emt_gsva$EMT_score,emt_gsva$gsva,function(.x,.y){
    tibble::tibble(EMT_score=.x$EMT_score,gsva=.y$atg_lys)
  }) 
)->emt_with_gsva


#calculate the correlate------------------------------------------
correlation=purrr::map(emt_with_gsva$emt_with_gsva,function(.x){
  cor.test(.x$EMT_score,.x$gsva)->cor
  return(cor$estimate[[1]])
})%>%unlist()

Pvalue=purrr::map(emt_with_gsva$emt_with_gsva,function(.x){
  cor.test(.x$EMT_score,.x$gsva)->cor
  cor$p.value[[1]]
})%>%unlist()
emt_with_gsva%>%
  dplyr::mutate(
    correlation=correlation,
    fdr=p.adjust(Pvalue, method = "fdr")
  ) -> cor_gsva


#save the table--------------------------------------------------------------------------
readr::write_rds(cor_gsva,path ="/home/shimw/TCGA/pancan33_EMT_with_gsva_cor.rds.gz",compress = "gz")

#draw the picture ---------------------------------------------------------

#point----------
cor_gsva%>%
  dplyr::filter(abs(correlation)>=0.25)%>%
  dplyr::select(cancer_type,emt_with_gsva)->emt_with_gsva_filter

  purrr::map2(emt_with_gsva_filter$cancer_type,emt_with_gsva_filter$emt_with_gsva,function(.x,.y){
      .y%>%
      tibble::add_column(tumor_type=.x)
    })%>%
    dplyr::bind_rows()%>%
    ggplot(aes(x=EMT_score,y=gsva)) +geom_point(aes(colour=tumor_type),alpha = .7,size = 3)+
    scale_color_manual(values=c("#FAD2D9", "#9EDDF9", "#B2509E","#97D1A9","#ED1C24","#F8AFB3","#A084BD","#D97D25","#6E7BA2","#DAF1FC","#00A99D"))+
    theme(
      panel.background = element_rect(colour = "black", fill = "white"),
      panel.grid = element_line(colour = "grey", linetype = "dashed"),
      panel.grid.major = element_line(
        colour = "grey",
        linetype = "dashed",
        size = 0.2
      ),
      axis.ticks = element_line(color = "black")
    )+
    labs(
      x = "EMT Score",
      y = "GSVA",
      colour="Cancer"
    )->p
  ggsave(
    filename = "gsva_point.pdf",
    plot = p,
    device = "pdf",
    width = 10,
    height = 8,
    path = "/home/shimw/github/EMT/pdf"
  )


#bar------------
#purrr::walk2(emt_with_gsva$cancer_type,emt_with_gsva$emt_with_gsva,function(.x,.y){
  #.y%>%
 # ggplot(aes(x=EMT_score,y=gsva)) +geom_point()
#})



cor_gsva%>%
    dplyr::arrange(desc(correlation))%>%
    dplyr::mutate(fdr = -log10(fdr)) %>% 
    dplyr::mutate(fdr = ifelse(fdr > 15, 15, fdr))%>%
    dplyr::filter(fdr>-log10(0.05))->cor_gsva_fdr
  
cor_gsva_fdr%>%
  ggplot(aes(cancer_type,correlation))+geom_point(aes(size=fdr))+
  scale_x_discrete(limits=cor_gsva_fdr$cancer_type)+
  theme(
    axis.text.x = element_text(angle = 90),
    panel.background = element_rect(colour = "black", fill = "white"),
   #  panel.grid = element_line(colour = "grey", linetype = "dashed"),
   # panel.grid.major = element_line(
   #    colour = "grey",
   #    linetype = "dashed",
   #    size = 0.2
   #  ),
    panel.grid =element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.text.align = 0.5
  )+
  scale_size_continuous(
    limit = c(-log10(0.05),15),

    range = c(1, 4),
    breaks = c(-log10(0.05), 5, 10, 15),
    labels = c("0.05", latex2exp::TeX("$10^{-5}$"), latex2exp::TeX("$10^{-10}$"), latex2exp::TeX("$< 10^{-15}$")),

    name = "FDR"
  )+
  labs(
    x = "Cancer Type",
    y = "Correlation"
  )-> p
  ggsave(
    filename = "gsva_cor_bub.pdf",
    plot = p,
    device = "pdf",
    width = 5,
    height = 4,
    path = "/home/shimw/github/EMT/pdf"
  )


  
  cor_gsva_fdr%>%dplyr::select(cancer_type,correlation,fdr)%>%
    tibble::add_column("EMT"="EMT",.before = 1)%>%
    graph_from_data_frame(directed = FALSE)->graph_cors
  ggraph(graph_cors) +
    geom_edge_link(aes(edge_alpha = abs(correlation), edge_width = fdr, edge_colour = correlation))+
    guides(edge_alpha = "none")+
    scale_edge_colour_gradientn(limits = c(-0.3, 0.5), colors = c("firebrick2", "dodgerblue2"))+
    geom_node_point(color = "#FFEFDB", size = 5)+
    geom_node_text(aes(label = name), repel = TRUE,size=3) +
    labs(
      edge_colour = "Correlation",
      edge_width = "FDR"
    )+
    scale_edge_width_continuous(breaks = c(-log10(0.05), 5, 10, 15),
                                labels = c("0.05", latex2exp::TeX("$10^{-5}$"), latex2exp::TeX("$10^{-10}$"), latex2exp::TeX("$< 10^{-15}$")))+
    theme_graph(base_family = "sans")->p
  ggsave(
    filename = "gsva_cor_net.pdf",
    plot = p,
    device = "pdf",
    width = 8,
    height = 7,
    path = "/home/shimw/github/EMT/pdf"
  )
