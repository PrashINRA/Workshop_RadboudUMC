# Workshop_RadboudUMC
workshop_scripts

# This interactive workshop is designed to get scRNA-seq data (GSE57872) and process it for biological meannigful insights



**Step1:download and unzip the GSE57872 data**

```console
# make a workspace- open terminal and do
mkdir -p GSE57872 && cd GSE57872

# download and unzip the data: also in terminal
wget -O - "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE57872&format=file&file=GSE57872%5FGBM%5Fdata%5Fmatrix%2Etxt%2Egz" \
  | gunzip > GSE57872_GBM_data_matrix.txt
```
**Step 2: Score per cell**

```{r}
library(UCell)
# keep only genes present
DDR_early_use <- intersect(DDR_early, rownames(Seu))
RS_use        <- intersect(RS,        rownames(Seu))
Seu <- AddModuleScore_UCell(Seu, features = list(DDR_early = DDR_early_use, RS = RS_use))
#Seu is a seurat object containing single-cell data from SSc+SjS patients and all the metadata
# Columns created would be: 'UCell_DDR_early', 'UCell_RS'

#Adjust for replication stress
library(dplyr)
library(broom)
meta <- Seu@meta.data %>%
  mutate(CTs = factor(CTs), Organ = factor(Organ))

adj <- meta %>%
  group_by(CTs, Organ) %>%
  do({
    m <- lm(UCell_DDR_early ~ UCell_RS + S.Score + G2M.Score + nCount_RNA + percent.mt, data = .)
    tibble(adj_DDR = resid(m))
  }) %>% ungroup()

meta$DDR_adj <- adj$adj_DDR

```
  
**Map adjusted score**
```{r}
library(scales)
meta <- meta %>%
  group_by(CTs, Organ) %>%
  mutate(DDR_adj_01 = rescale(DDR_adj, to = c(0,1))) %>%
  ungroup()

# Binary call (tune thresholds as you like)
meta$DDR_call <- ifelse(meta$DDR_adj_01 >= 0.8, 1L, 0L)  # top 20% = high damage

```
This will create an adjuste DDR scores (removing confounding from replication stress). 

