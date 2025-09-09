# Workshop_RadboudUMC
workshop_scripts

# This interactive workshop is designed to get scRNA-seq data (GSE57872) and process it for biological meannigful insights



**Step1:download and unzip the GSE57872 data**

```console
# make a workspace- open terminal and do
mkdir -p ~/GSE120221 && cd ~/GSE120221 && wget -O GSE120221_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE120221&format=file" && tar -xvf GSE120221_RAW.tar

```
**Step 2: make a Seurat object**

```{r}

# load packages
pkgs <- c( 'ggplot2', 'Seurat','dplyr','rstatix')
sapply(pkgs, library, character.only = T)
theme_set(theme_bw())
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

