# Library -----------------------------------------------------------------
library(magrittr)
library(GSVA)
library(ggplot2)

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

#draw the picture ---------------------------------------------------------

#point----------
cor_gsva%>%
  dplyr::filter(abs(correlation)>=0.3)%>%
  dplyr::select(cancer_type,emt_with_gsva)->emt_with_gsva_filter

  purrr::map2(emt_with_gsva_filter$cancer_type,emt_with_gsva_filter$emt_with_gsva,function(.x,.y){
      .y%>%
      tibble::add_column(tumor_type=.x)
    })%>%
    dplyr::bind_rows()%>%
    ggplot(aes(x=EMT_score,y=gsva)) +geom_point(aes(colour=tumor_type))+
    scale_color_manual(values=c("#FAD2D9", "#97D1A9", "#A084BD","#6E7BA2","#00A99D"))+
    theme(
      panel.background = element_rect(colour = "black", fill = "white"),
      panel.grid = element_line(colour = "grey", linetype = "dashed"),
      panel.grid.major = element_line(
        colour = "grey",
        linetype = "dashed",
        size = 0.2
      ),
      axis.ticks = element_line(color = "black")
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
    dplyr::mutate(fdr = ifelse(fdr > 15, 15, fdr))->cor_gsva
cor_gsva%>%
  ggplot(aes(cancer_type,correlation))+geom_point(aes(size=fdr))+
  scale_x_discrete(limits=cor_gsva$cancer_type)+
  theme(
    axis.text.x = element_text(angle = 90),
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    axis.ticks = element_line(color = "black"),
    legend.text.align = 0.5
  )+
  scale_size_continuous(
    limit = c(-log10(0.05),15),

    range = c(1, 4),
    breaks = c(-log10(0.05), 5, 10, 15),
    labels = c("0.05", latex2exp::TeX("$10^{-5}$"), latex2exp::TeX("$10^{-10}$"), latex2exp::TeX("$< 10^{-15}$")),
    
    name = "FDR"
  )-> p
  ggsave(
    filename = "gsva_cor_bar.pdf",
    plot = p,
    device = "pdf",
    width = 8,
    height = 8,
    path = "/home/shimw/github/EMT/pdf"
  )

  









