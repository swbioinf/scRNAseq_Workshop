# Catchup script

# Rebuild the seurat_object used in analysis in one go, without all the plotting and exploration.
# Run up to the appropriate section heading. 
# Handy if you have accidentally overwritten/brokern the object at some point.

library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)


## Load
pbmc.data <- Read10X_h5("data/filtered_feature_bc_matrix.h5")
metadata <- read.table("data/metadata.txt")
seurat_object <- CreateSeuratObject(counts = pbmc.data ,
                                    assay = "RNA", project = 'pbmc')
seurat_object  <- AddMetaData(object = seurat_object, metadata = metadata)




## QC

seurat_object$percent.mt <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & percent.mt < 5)



## Normalise
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 1e4)


## PCA and UMAP
seurat_object <- FindVariableFeatures(seurat_object, selection.method = 'vst', nfeatures = 2000)
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)

seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
seurat_object <- RunUMAP(seurat_object, dims = 1:10)

# save intermediate file
#qs2::qs_save(seurat_object, "data/seurat_object_preprocessed.qs2")

# ---------------
# Part 2
# ---------------

library(qs2)
library(Seurat)
library(harmony)
library(tidyverse)
library(SingleCellExperiment)
library(SingleR)
library(celldex)
library(patchwork)
library(limma)
library(edgeR)

# Harmony

# Load file
seurat_object <- qs2::qs_read("data/seurat_object_preprocessed.qs2")

seurat_object<- FindNeighbors(seurat_object, reduction="pca", dims=1:10)
seurat_object <- FindClusters(seurat_object, resolution=0.6)
seurat_object$pca_clusters <- seurat_object$seurat_clusters

seurat_object <- RunHarmony(seurat_object, c("stim", "ind"), reduction="pca",reduction.save="harmony")
seurat_object <- RunUMAP(seurat_object, reduction="harmony", dims=1:10, reduction.name="umap_harmony")
seurat_object <- FindNeighbors(seurat_object, reduction="harmony", dims=1:10)
seurat_object <- FindClusters(seurat_object, resolution=0.6)
seurat_object$harmony_clusters <- seurat_object$seurat_clusters

## Clustering
seurat_object <- qs2::qs_read("data/seurat_object_preprocessed.harmony.qs2")
min_res <- 0.4
max_res <- 2
interval <- 0.2
seurat_object <- FindClusters(seurat_object, resolution = seq(min_res, max_res, interval), verbose = F)

Idents(seurat_object) <- seurat_object$RNA_snn_res.0.6

## Cluster markers
seurat_object <- qs2::qs_read("data/seurat_object_preprocessed.harmony.clustered.qs2")

new.cluster.ids <- c(
  "0" = "CD4 T cells",
  "1" = "CD14+ Monocytes",
  "2" = "NK cells",
  "3" = "B cells",
  "4" = "FCGR3A+ Monocytes",
  "5" = "Megakaryocytes",
  "6" = "Ambiguous",
  "7" = "Dendritic cells"
)

seurat_object@meta.data$markers <- new.cluster.ids[as.character(seurat_object$RNA_snn_res.0.6)]

## Single R

seurat_object <- qs2::qs_read("data/seurat_object_preprocessed.harmony.clustered.markers.qs2")
Idents(seurat_object) <- seurat_object$markers

sce <- as.SingleCellExperiment(seurat_object)
ref.set <- celldex::HumanPrimaryCellAtlasData()

pred.cnts <- SingleR::SingleR(test = sce, ref = ref.set, labels = ref.set$label.main)
lbls.keep <- table(pred.cnts$labels)>10
seurat_object$SingleR.labels <- ifelse(lbls.keep[pred.cnts$labels], pred.cnts$labels, 'Other')


## DE
seurat_object <- qs2::qs_read("data/seurat_object_preprocessed.harmony.clustered.markers.singler.qs2")
seurat_object$sample_name          <- paste(seurat_object$stim, seurat_object$ind, sep=".")
seurat_object$sample_celltype_name <- paste(gsub(" ","",seurat_object$cell), seurat_object$ind, seurat_object$stim, sep=".")

total_per_gene <- rowSums(GetAssayData(seurat_object, assay='RNA', layer='counts'))
seurat_object  <- seurat_object[total_per_gene >= 50, ]

Idents(seurat_object) <- seurat_object$cell
seurat_object_celltype <- seurat_object[, seurat_object$cell == "CD14+ Monocytes"]

pseudobulk_annotation_table <- FetchData(seurat_object_celltype,
                                         vars = c('sample_name', 'ind', 'stim', 'cell', 'sample_celltype_name')) |>
  as_tibble() |>
  group_by(across(everything())) |>
  summarise(n_cells = n(), .groups = "drop")

Idents(seurat_object_celltype) <- seurat_object_celltype$sample_name
pseudobulk_matrix_list <- AggregateExpression(seurat_object_celltype, slot = 'counts', assays = 'RNA')
pseudobulk_matrix      <- pseudobulk_matrix_list[['RNA']]
pseudobulk_matrix      <- pseudobulk_matrix[, pseudobulk_annotation_table$sample_name]

dge      <- DGEList(pseudobulk_matrix)
dge      <- calcNormFactors(dge)
ind      <- as.character(pseudobulk_annotation_table$ind)
stim     <- pseudobulk_annotation_table$stim
design   <- model.matrix(~0 + stim + ind)
vm       <- voom(dge, design = design, plot = FALSE)
fit      <- lmFit(vm, design = design)
contrasts <- makeContrasts(stimstim - stimctrl, levels = coef(fit))
fit      <- contrasts.fit(fit, contrasts)
fit      <- eBayes(fit)
de_result <- arrange(topTable(fit, n = Inf, adjust.method = "BH"), adj.P.Val)

pseudobulk_annotation_table.allcelltypes <- FetchData(seurat_object,
                                                       vars = c('sample_name', 'ind', 'stim', 'cell', 'sample_celltype_name')) %>%
  as_tibble() %>%
  group_by(across(everything())) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 20)

Idents(seurat_object) <- seurat_object$sample_celltype_name
pseudobulk_matrix_list.allcelltypes <- AggregateExpression(seurat_object, slot = 'counts', assays = 'RNA')
pseudobulk_matrix.allcelltypes      <- pseudobulk_matrix_list.allcelltypes[['RNA']]
pseudobulk_matrix.allcelltypes      <- pseudobulk_matrix.allcelltypes[, pseudobulk_annotation_table.allcelltypes$sample_celltype_name]

dge           <- DGEList(pseudobulk_matrix.allcelltypes)
dge           <- calcNormFactors(dge)
norm_matrix   <- edgeR::cpm(dge)
lognorm_matrix <- log2(norm_matrix + 1) 

 

