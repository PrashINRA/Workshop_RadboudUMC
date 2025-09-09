# Workshop_RadboudUMC
workshop_scripts

# This interactive workshop is designed to get scRNA-seq data (GSE120221) and process it for biological meannigful insights

**Step1:download and unzip the GSE120221 data**

```console
# make a workspace- open terminal and do
mkdir -p ~/GSE120221 && cd ~/GSE120221 && wget -O GSE120221_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE120221&format=file" && tar -xvf GSE120221_RAW.tar

```
**Step 2: make a Seurat object**

```{r}

# load packages
pkgs <- c( 'ggplot2', 'Seurat','dplyr','stringr')
sapply(pkgs, library, character.only = T)
theme_set(theme_bw())

setwd("/home/prashant/GSE120221")  # adjust with you path

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
  
  objs[[k]] <- CreateSeuratObject(counts = mat, project = paste0("GSE120221_", k))
}

# Merge into one Seurat object; prefix cell barcodes with sample keys
seu <- if (length(objs) == 1) {
  objs[[1]]
} else {
  merge(x = objs[[1]], y = objs[-1], add.cell.ids = names(objs), 
        project = "GSE120221")
}

##Make a sample metadta
seu$Sample <- factor(seu$orig.ident)
seu$Sample <- sub("^GSE120221_", "BM_", seu$Sample)
levels(factor(seu$Sample))

rm(list=setdiff(ls(), c('seu')))
gc()
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
```

