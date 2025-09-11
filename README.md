# Workshop_RadboudUMC
workshop_scripts

# This interactive workshop is designed to get scRNA-seq data (GSE120221) and process it to obtain biologically meannigful insights

**Step1:download and unzip the GSE120221 data**

```console
# make a workspace- open terminal and type
mkdir -p ~/GSE120221 && cd ~/GSE120221 && wget -O GSE120221_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE120221&format=file" && tar -xvf GSE120221_RAW.tar

```
**Step 2: make a Seurat object**

```{r}

# load packages

pkgs <- c( 'ggplot2', 'Seurat','dplyr','stringr')
sapply(pkgs, library, character.only = T)
theme_set(theme_bw())   ##Install the packages if you haven't already

setwd("/Users/prashant/GSE120221")  # adjust with you path

# List triplets
genes   <- sort(Sys.glob("GSM*_genes_*.tsv.gz"))
bars    <- sort(Sys.glob("GSM*_barcodes_*.tsv.gz"))
mtx     <- sort(Sys.glob("GSM*_matrix_*.mtx.gz"))

stopifnot(length(genes)==length(bars), length(bars)==length(mtx), length(mtx) > 0)

# Derive sample keys from the trailing token (e.g., A, B, C1, …)
key_from <- function(x) sub(".*_(.+)\\.(mtx|tsv)\\.gz$", "\\1", x)

samples <- unique(key_from(mtx))
objs <- list()

for (k in samples) {
  g <- genes[ key_from(genes) == k ]
  b <- bars [ key_from(bars)  == k ]
  m <- mtx  [ key_from(mtx)   == k ]
  stopifnot(length(g)==1, length(b)==1, length(m)==1)
  
  
  mat <- tryCatch(
    ReadMtx(mtx = m, features = g, cells = b, feature.column = 2),
    error = function(e) ReadMtx(mtx = m, features = g, cells = b, feature.column = 1)
  )
  
  objs[[k]] <- CreateSeuratObject(counts = mat, min.cells=30, project = paste0("BM_", k))
}

# Merge into one Seurat object; prefix cell barcodes with sample keys
seu <- if (length(objs) == 1) {
  objs[[1]]
} else {
  merge(x = objs[[1]], y = objs[-1], add.cell.ids = names(objs), 
        project = "test")
}

##Make a sample metadta
seu$Sample <- factor(seu$orig.ident)
levels(factor(seu$Sample))

rm(list=setdiff(ls(), c('seu')))
gc()

##Next we want to see sequencing depth per sample before QC or normalization

ggplot(seu@meta.data, aes(x = Sample, y = nCount_RNA, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5) +scale_y_log10()+
  labs(y = "Sequencing depth (UMIs per cell)", x = "Sample") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = 'none')
```
  
**Step 3 QC:Cell level filtering**
```{r}
seu <- JoinLayers(seu, assay = 'RNA')
counts <- GetAssayData(seu, assay = 'RNA', layer = 'counts')
counts[1:10,1:3]
genes_per_cell <- Matrix::colSums(counts>0) # count a gene only if it has non-zero reads mapped.
counts_per_cell <- Matrix::colSums(counts)
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')

MIN_GENES_PER_CELL <- #What should be the value here?
MAX_GENES_PER_CELL <- #What should be the value here?  

# now replot with the thresholds being shown:
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
abline(h=MIN_GENES_PER_CELL, col='magenta')  # lower threshold
abline(h=MAX_GENES_PER_CELL, col='blue') # upper threshold

```
This will create a library complexity plot (~cell vs sorted genes-per-cell). 

**Step 3 QC:filter out cells with high % of MT-genes**

```{r}
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
mito_genes <- grep("^mt-", rownames(counts) , ignore.case=T, value=T)
mito_gene_read_counts = Matrix::colSums(counts[mito_genes,])
pct_mito = mito_gene_read_counts / counts_per_cell * 100
plot(sort(pct_mito), xlab = "cells sorted by percentage mitochondrial counts", ylab = 
       "percentage mitochondrial counts")

MAX_PCT_MITO <- #what should be cutoff for %MT

plot(sort(pct_mito))
abline(h=MAX_PCT_MITO, col='red')

rm(list=setdiff(ls(), c('seu')))
gc()

#Subset the Seurat object with required cutoffs
dim(seu)
seu <- subset(seu, subset = nFeature_RNA >  & nFeature_RNA <  & percent.mt< )
dim(seu)  ##Check the data dimension before and after QC
```
**Normalize data and check normaize seq depths**

```{r}
seu <- NormalizeData(seu, scale.factor = median(seu$nCount_RNA))

# Extract normalized data
norm_mat <- GetAssayData(seu, slot = "data")

# Sum across genes per cell
seu$norm_counts_per_cell <- Matrix::colSums(norm_mat)

# Plot, same as before
ggplot(seu@meta.data, aes(x = Sample, y = norm_counts_per_cell, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(y = "Normalized sequencing depth (log-normalized counts per cell)",
       x = "Sample") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = 'none')

##remove the sample - 'BM_G", why?

seu <- subset(seu, subset = Sample != "BM_G")

seu$Sample <- droplevels(seu$Sample)

table(seu$Sample)

```

**Now we are downsampling  the seurat object to make it smaller and processable in local machines with ~16GB RAM**

```{r}

set.seed(100) 

# sample 1000 cells per Sample group (or all cells if fewer available)
n_cells <- 1000
cells_keep <- unlist(lapply(split(colnames(seu), seu$Sample), function(cells) {
  if (length(cells) > n_cells) {
    sample(cells, n_cells)
  } else {
    cells
  }
}))

# subset Seurat object
seu <- subset(seu, cells = cells_keep)

# check result
table(seu$Sample)

rm(list=setdiff(ls(), c('seu')))
gc()
```

**Process the downsampled data**

```{r}
seu <- seu |> NormalizeData(scale.factor = median(seu$nCount_RNA)) |>
      FindVariableFeatures(nfeatures = 3000) |>
      ScaleData(features = rownames(seu))

##Visualiza variable features

library(ggrepel)

# Get HVF info and mark hypervariable
variance.data <- as_tibble(HVFInfo(seu), rownames = "Gene") %>%
  mutate(hypervariable = Gene %in% VariableFeatures(seu))

# Rank HVGs by variance (or standardized variance)
top20 <- variance.data %>%
  filter(hypervariable) %>%
  arrange(desc(variance.standardized)) %>%
  slice(1:20)

# Plot with labels
variance.data %>%
  ggplot(aes(x = log(mean), y = log(variance), color = hypervariable)) +
  geom_point() +
  geom_text_repel(
    data = top20,
    aes(label = Gene),
    color = "firebrick",
    size = 3,
    max.overlaps = 50
  ) +
  scale_color_manual(values = c("black", "chartreuse2")) +
  labs(x = "log(mean expression)", y = "log(variance)", color = "Hypervariable")

##What does this plot suggest?

```

**Further processing for Dimension Reduction, clustering and it's optimization**


```{r}
seu <- seu  |> 
          RunPCA(features = VariableFeatures(seu), ndims.print = 1:3) |>
           FindNeighbors(reduction = 'pca', dims = 1:30, annoy.metric='cosine')

##Malke Elbow plot to asses PCs

ElbowPlot(seu, ndims=50)


#Optimize resolution for clustering

library('clustree')
BMH <- seu  ##Make a copy of seurat object to test resolutions
BMH <- FindClusters(BMH, resolution =  seq(from = 0.1, to = 1, by = 0.1))
clust_seurat <- BMH@meta.data %>% dplyr::select(dplyr::contains("RNA_snn_res."))
clustree::clustree(clust_seurat, prefix="RNA_snn_res.")

rm(list=setdiff(ls(), c('seu')))
gc()

#choose less noisy and biologically meaningful resolution
seu <- FindClusters(seu, resolution = 0.3) |>
        RunUMAP( reduction = 'pca', dims = 1:30, n.components = 3L, metric ='cosine', umap.method = 'uwot')

```

**Visualizations**

```{r}
library(RColorBrewer)
n <-18
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(cols, n))
clust_cols <- cols[1:n]

my_levels <- paste0('C', 1:n)
levels(seu@active.ident)<- my_levels
levels(seu$seurat_clusters)<- my_levels

seu$Celltypes <- seu$seurat_clusters

##Plot UMAP
DimPlot(seu, group.by = 'Celltypes', cols = clust_cols,label = T, reduction = 'umap')+NoLegend()

#or
devtools::install_github("pwwang/scplotter")
scplotter::CellDimPlot(seu, reduction = 'umap', group_by='Celltypes',label=T,label_insitu = TRUE,cols=clust_cols, 
                       theme = "theme_blank")

```

**How to annotate cells**

```{r}
##Find DEGs per cluser

cl.markers <-FindAllMarkers(seu, min.pct = 0.2, min.diff.pct = 0.2,only.pos = T, 
                            verbose = T,slot = 'data')

##Filter genes on the basis of logFC and adjP-values

markers_cl <- cl.markers %>%
  dplyr::filter(avg_log2FC > 1, p_val_adj < 1e-3) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE)

#Get top genes for plotting 
top5_cl <- markers_cl %>%
  dplyr::filter(avg_log2FC > 0.5, p_val_adj < 1e-3) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 5)

top5_cl <- unique(top5_cl$gene)

##Get avg expression of top genes in each cluster 
avg_exp <- AverageExpression(seu, assays = 'RNA', features = top5_cl, group.by = 'Celltypes')$RNA

# Create the heatmap
pheatmap::pheatmap(avg_exp,
                   scale = "row",
                   cluster_rows = T, cluster_cols = T,
                   show_rownames = TRUE,angle_col = 45, fontsize_row = 7,
                   show_colnames = TRUE,treeheight_row = 5,treeheight_col = 5,
                   color = colorRampPalette(c("lightblue", "white", 'firebrick'))(100),
                   main = "Top DEGs per Cell Type")

##This will create a heatmap for DEGs per cluster and will help you to annotate cell-types

#Use scillus to find functionality of each cluster
devtools::install_github("xmc811/Scillus", ref = "development")

Scillus::plot_cluster_go(markers_cl, ont='BP',org = 'human', cluster_name = 'C1')

levels(seu$Celltypes) <- c("Naive_TCs", "CD8Ts", "Erythroids", "GMP", "MEP", "BCs", 
                           "Erythroids", "CD4T", "HSPCs", "T_progenitors", "NKs", "pDCs", "Pre_BCs",
                           "Monocytes", "Monocytes", "PlasmaCs", "MEP", "Pro_BCs")


scplotter::CellDimPlot(seu, reduction = 'umap', group_by='Celltypes',label=T,label_insitu = TRUE,cols=clust_cols, 
                       theme = "theme_blank")


rm(list=setdiff(ls(), c('seu', 'cl.markers', 'clust_cols')))
gc()
```

**How to check contribution of each donor on celltypes?**

