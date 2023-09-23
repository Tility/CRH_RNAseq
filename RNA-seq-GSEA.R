## GSEA enrichment and Figure ##
library(ggplot2) 
library(xlsx)
library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)

### read data =======
setwd("/Users/liangtingting/hp/lab_file/WMM/RNAseq_CRH/Data/GSEA/")

geneList1 = read.table("/Users/liangtingting/hp/lab_file/WMM/RNAseq_CRH/Data/GSEA/H-vs-C-all.gene.txt", header=TRUE)
geneList = geneList1[,c('gene_id','log2FoldChange')]

### check duplicated
table(duplicated(geneList$gene_id))
##FALSE 
##14834 
rownames(geneList) <- geneList$gene_id
geneList_1 = arrange(geneList, desc(log2FoldChange))
geneList1 <- geneList_1
head(geneList1)

id <- geneList1$log2FoldChange
names(id) <- rownames(geneList1)
head(id)
##THBS2    POSTN     RND3    CRLF1    LAMA4   COL1A2 
##6.993519 6.193616 6.018307 5.524683 5.470309 5.252303 

# get all human gene sets
#all_sets = msigdbr(species = "Homo sapiens" ) 
load(file = 'human_geneset.Rdata')

egmt2 <- GSEA(id, TERM2GENE= all_sets[,c('gs_name','gene_symbol')] , 
              minGSSize = 1,
              eps = 0, 
              pvalueCutoff = 0.99,
              verbose=FALSE)

gsea_results3 <- egmt2@result 
write.xlsx(gsea_results3,file = 'gsea_all_humanpathyway.xlsx')

### |NES|>1，NOM pvalue<0.05，FDR（padj）<0.25
g2 <- gsea_results3[gsea_results3$pvalue<0.05 & gsea_results3$p.adjust<0.25 & abs(gsea_results3$NES)>1,]
g2 <- g2[order(g2$NES,decreasing = T),]
write.xlsx(g2,file = 'gsea_humanpathyway_p_0.05_2.xlsx')

###### save Rdata result 
save(egmt2,file = 'gsea_result.Rdata')

load(file = 'gsea_result.Rdata')

### Pathway Visualization
g3 <- c('REACTOME_RND3_GTPASE_CYCLE',
        'REACTOME_RHO_GTPASES_ACTIVATE_PKNS',
        'REACTOME_RHOJ_GTPASE_CYCLE',
        'REACTOME_RHO_GTPASE_CYCLE',
        'GOBP_POSITIVE_REGULATION_OF_SMALL_GTPASE_MEDIATED_SIGNAL_TRANSDUCTION',
        'KEGG_ECM_RECEPTOR_INTERACTION',
        'KEGG_FOCAL_ADHESION'
)


source('/Users/liangtingting/hp/lab_file/WMM/RNAseq_CRH/gsea_plot2_modify.R')

for (i in 1:length(g3)) {
  
  pdf(paste0('/Users/liangtingting/hp/lab_file/WMM/RNAseq_CRH/Figure/',g3[i],'.pdf'),width = 6, height = 6)
  p = gseaplot2_mod_tt(egmt2,geneSetID = g3[i],
                       title = NULL,
                       color = "green",
                       color_hmp = c("#39518F","#7D2125"),
                       base_size = 10,
                       rel_heights = c(1.5, 0.5, 1),
                       subplots = 1:3,
                       pvalue_table = F,
                       fontface = "plain", 
                       pCol = "black", 
                       pvalSize = 6,
                       pvalX = 0.9, pvalY = 0.8,
                       ES_geom = "line"
  )
  print(p)
  dev.off()
  
}
sessionInfo()
##R version 4.1.2 (2021-11-01)
##Platform: x86_64-apple-darwin17.0 (64-bit)
##Running under: macOS 13.5.1

##Matrix products: default
##LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

##locale:
##[1] zh_CN.UTF-8/zh_CN.UTF-8/zh_CN.UTF-8/C/zh_CN.UTF-8/zh_CN.UTF-8

##attached base packages:
##[1] stats     graphics  grDevices utils     datasets  methods   base     

##other attached packages:
##[1] enrichplot_1.14.2     msigdbr_7.5.1         clusterProfiler_4.2.2 forcats_0.5.2         stringr_1.4.1         dplyr_1.0.10          purrr_0.3.5          
##[8] readr_2.1.3           tidyr_1.2.1           tibble_3.1.8          tidyverse_1.3.2       xlsx_0.6.5            ggplot2_3.4.0        

##loaded via a namespace (and not attached):
##[1] shadowtext_0.1.2       readxl_1.4.1           backports_1.4.1        fastmatch_1.1-3        plyr_1.8.8             igraph_1.3.5          
##[7] lazyeval_0.2.2         splines_4.1.2          BiocParallel_1.28.3    usethis_2.1.6          GenomeInfoDb_1.30.1    digest_0.6.30         
##[13] yulab.utils_0.0.5      htmltools_0.5.3        GOSemSim_2.20.0        viridis_0.6.2          GO.db_3.14.0           fansi_1.0.3           
##[19] magrittr_2.0.3         memoise_2.0.1          googlesheets4_1.0.1    tzdb_0.3.0             remotes_2.4.2          Biostrings_2.62.0     
##[25] graphlayouts_0.8.3     modelr_0.1.10          timechange_0.1.1       prettyunits_1.1.1      colorspace_2.0-3       blob_1.2.3            
##[31] rvest_1.0.3            ggrepel_0.9.2          haven_2.5.1            callr_3.7.3            crayon_1.5.2           RCurl_1.98-1.10       
##[37] jsonlite_1.8.4         scatterpie_0.1.8       ape_5.6-2              glue_1.6.2             polyclip_1.10-4        gtable_0.3.1          
##[43] gargle_1.2.1           zlibbioc_1.40.0        XVector_0.34.0         pkgbuild_1.3.1         BiocGenerics_0.40.0    scales_1.2.1          
##[49] DOSE_3.20.1            DBI_1.1.3              miniUI_0.1.1.1         Rcpp_1.0.10            viridisLite_0.4.1      xtable_1.8-4          
##[55] tidytree_0.4.1         gridGraphics_0.5-1     bit_4.0.5              stats4_4.1.2           profvis_0.3.7          htmlwidgets_1.5.4     
##[61] httr_1.4.4             fgsea_1.20.0           RColorBrewer_1.1-3     ellipsis_0.3.2         urlchecker_1.0.1       pkgconfig_2.0.3       
##[67] rJava_1.0-6            farver_2.1.1           dbplyr_2.2.1           utf8_1.2.2             ggplotify_0.1.0        tidyselect_1.2.0      
##[73] rlang_1.0.6            reshape2_1.4.4         later_1.3.0            AnnotationDbi_1.56.2   munsell_0.5.0          cellranger_1.1.0      
##[79] tools_4.1.2            cachem_1.0.6           downloader_0.4         cli_3.6.0              generics_0.1.3         RSQLite_2.2.20        
##[85] devtools_2.4.5         broom_1.0.1            fastmap_1.1.0          ggtree_3.2.1           babelgene_22.9         processx_3.8.0        
##[91] bit64_4.0.5            fs_1.5.2               tidygraph_1.2.2        KEGGREST_1.34.0        ggraph_2.1.0           nlme_3.1-160          
##[97] mime_0.12              aplot_0.1.8            DO.db_2.9              xml2_1.3.3             compiler_4.1.2         rstudioapi_0.14       
##[103] png_0.1-8              treeio_1.18.1          reprex_2.0.2           tweenr_2.0.2           stringi_1.7.8          ps_1.7.2              
##[109] lattice_0.20-45        Matrix_1.5-1           vctrs_0.5.2            pillar_1.8.1           lifecycle_1.0.3        data.table_1.14.4     
##[115] bitops_1.0-7           patchwork_1.1.2        httpuv_1.6.6           qvalue_2.26.0          R6_2.5.1               promises_1.2.0.1      
##[121] gridExtra_2.3          IRanges_2.28.0         sessioninfo_1.2.2      MASS_7.3-58.1          assertthat_0.2.1       pkgload_1.3.1         
##[127] xlsxjars_0.6.1         withr_2.5.0            S4Vectors_0.32.4       GenomeInfoDbData_1.2.7 parallel_4.1.2         hms_1.1.2             
##[133] grid_4.1.2             ggfun_0.0.8            googledrive_2.0.0      ggforce_0.4.1          Biobase_2.54.0         shiny_1.7.3           
##[139] lubridate_1.9.0 