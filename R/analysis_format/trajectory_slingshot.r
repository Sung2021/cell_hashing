library(slingshot)
library(SingleCellExperiment)

sce <- slingshot(obj.srt@reductions$umap@cell.embeddings, 
                 clusterLabels = obj.srt@meta.data$method_time, reducedDim = 'umap')
sce
sce %>% str()

plot(reducedDims(sce))

obj.srt@meta.data[1:3,]

sce <- as.SingleCellExperiment(obj.srt, assay = "RNA")
sce <- slingshot(sce, reducedDim = 'UMAP', 
                 clusterLabels = obj.srt@meta.data$method_time, start.clus = '4',
                 approx_points = 150)
plot(reducedDims(sce)$UMAP, col=obj.srt@meta.data$method_time)
lines(SlingshotDataSet(sce), lwd=2)

