
# Useful plot for haemonglobin exercise.

seurat_object$pc.haemoglobin <- PercentageFeatureSet(seurat_object, features = c("HBA1","HBA2","HBB"))
VlnPlot(seurat_object, 'pc.haemoglobin')

FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "pc.haemoglobin")
ggplot(seurat_object@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, col=pc.haemoglobin)) +
  geom_point() +
  theme_bw()


