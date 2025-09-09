# Workshop_RadboudUMC
workshop_scripts

# This interactive workshop is designed to get scRNA-seq data (GSE57872) and process it for biological meannigful insights



**Step1:download the GSE57872 data**

```{r}
DDR_early <- c(
  "CDKN1A","GADD45A","GADD45B","RRM2B","DDB2","XPC","MDM2","PMAIP1","BBC3",
  "SESN1","TP53INP1","POLH",        # early p53 targets / damage responders
  "PARP1","H2AFX","ATM","ATR","CHEK1","CHEK2","TP53BP1","MDC1"  # sensors/adaptors
)
#For replication stress
RS <- c(
  "MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","PCNA","RRM2",
  "RPA1","RPA2","RPA3","CLSPN","TIMELESS","TIPIN","ETAA1","TOPBP1",
  "CDC45","TYMS","TK1","POLA1","POLE","RAD51","BRCA1","FANCD2"  # ATR–CHK1 axis / fork protection
)


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

