# install.packages('tidyverse')
# install.packages('Seurat')
# install.packages('ggsci')
# install.packages('future')
# install.packages('progressr')
# install.packages('RcppAnnoy')
# install.packages('cowplot')
# install.packages('munsell')
# install.packages('irlba')
# install.packages('patchwork')
# install.packages('plotly')
# install.packages('sctransform')
# install.packages('tensor')
# install.packages('goftest')

library(tidyverse)
library(Seurat)
library(ggsci)
library(readr)
library(dplyr)

# AML329_dem_tbl <- read_tsv("https://ftp.ncbi.nlm.nGSM3587940_AML329-D0.dem.txt.gz", show_col_types = FALSE)
AML329_dem_tbl <- read_tsv("GSM3587940_AML329-D0.dem.txt.gz", show_col_types = FALSE)
dim(AML329_dem_tbl)
AML329_dem_tbl[1:8, 1:8]

# AML329_annot_tbl <- read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3587nnn/GSM3587941/suppl/GSM3587941_AML329-D0.anno.txt.gz", show_col_types = FALSE)
AML329_annot_tbl <- read_tsv("GSM3587941_AML329-D0.anno.txt.gz", show_col_types = FALSE)
dim(AML329_annot_tbl)
head(AML329_annot_tbl)
class(AML329_dem_tbl)

AML329_dem_mat <- AML329_dem_tbl %>% column_to_rownames("Gene") %>% as.matrix()
AML329_dem_mat[1:8, 1:3]

AML329_dem_mat %>% colSums() %>% summary()
