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


## Harmony

seurat_object<- FindNeighbors(seurat_object, reduction="pca", dims=1:10)
seurat_object <- FindClusters(seurat_object, resolution=0.5)
seurat_object$pca_clusters <- seurat_object$seurat_clusters

library(harmony)
seurat_object <- RunHarmony(seurat_object, c("stim", "ind"), reduction="pca",reduction.save="harmony")
seurat_object <- RunUMAP(seurat_object, reduction="harmony", dims=1:10, reduction.name="umap_harmony")
seurat_object <- FindNeighbors(seurat_object, reduction="harmony", dims=1:10)
seurat_object <- FindClusters(seurat_object, resolution=0.5)
seurat_object$harmony_clusters <- seurat_object$seurat_clusters


## Clustering

resolution = 2
seurat_object <- FindClusters(seurat_object, resolution = seq(0.1, resolution, 0.1)) # Slow, can be skipped if already past the cluster resolution

Idents(seurat_object) <- seurat_object$RNA_snn_res.0.5



## Cluster makrers

genes_markers <- list(Naive_CD4_T = c("IL7R", "CCR7"))
seurat_object <- AddModuleScore(object = seurat_object, features = genes_markers, ctrl = 5, name = "Naive_CD4_T",  search = TRUE)
seurat_object$cell_label = NA
seurat_object$cell_label[seurat_object$Naive_CD4_T1 > 1] = "Naive_CD4_T"
Idents(seurat_object) = seurat_object$cell_label



Idents(seurat_object) <- seurat_object$RNA_snn_res.0.5


## Single R

library(SingleCellExperiment)
library(SingleR)
library(celldex)

sce <- as.SingleCellExperiment(seurat_object)
ref.set <- celldex::HumanPrimaryCellAtlasData()

pred.cnts <- SingleR::SingleR(test = sce, ref = ref.set, labels = ref.set$label.main)
lbls.keep <- table(pred.cnts$labels)>10
seurat_object$SingleR.labels <- ifelse(lbls.keep[pred.cnts$labels], pred.cnts$labels, 'Other')


## DE
# Only run when up to this point, as it will subset the main object.
# rest of code should be run from this section, it builds new objects, and doesn't change seurat_object
seurat_object$sample_name          <- paste(seurat_object$stim, seurat_object$ind, sep=".")
seurat_object$sample_celltype_name <- paste(gsub(" ","",seurat_object$cell), seurat_object$ind, seurat_object$stim, sep=".")


total_per_gene <- rowSums(GetAssayData(seurat_object, assay='RNA', layer='counts'))
seurat_object<- seurat_object[total_per_gene >= 50, ] 

 

