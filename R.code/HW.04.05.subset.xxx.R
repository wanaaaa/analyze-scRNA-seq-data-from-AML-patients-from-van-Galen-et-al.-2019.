library(readr)
library(tidyverse)
library(Seurat)
library(ggsci)
library(readr)
library(dplyr)

AML419_dem_tbl <- read_tsv("../Data/hwData/GSM3587950_AML419A-D0.dem.txt.gz", show_col_types = FALSE)
AML419_annot_tbl <- read_tsv("../Data/hwData/GSM3587951_AML419A-D0.anno.txt.gz", show_col_types = FALSE)

AML419_annot_tbl$CellType[400:700]
unique(AML419_annot_tbl$CellType)

AML419_annot_tbl$PredictionRefined
  