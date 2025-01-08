library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggnetwork)
library(network)
library(GGally)
library(GO.db)
library(org.Hs.eg.db)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(wesanderson)
library(topGO)
library(EnsDb.Hsapiens.v79)
library(statmod)
library(patchwork)
library(multienrichjam)
library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(GOSemSim)
library(ggupset)
library(circlize)
library(httr)
library(jsonlite)
library(rmarkdown)
library(CeTF)

tr <- tr.CCvsctrl

genes <- tr$genes %>% 
  tibble::rownames_to_column("ENSEMBL")

tr.genes <- tr$table %>% 
  tibble::rownames_to_column("ENSEMBL") %>% 
  left_join(., genes) %>% 
  relocate(Symbol, .after = "ENSEMBL")

data("TFs")

TF <- data.frame(ENSEMBL = TFs, ID = "TF")

TF.gene <- tr.genes %>% 
  left_join(., TF) %>% 
  mutate(ID = replace_na(ID, "target_gene"))
