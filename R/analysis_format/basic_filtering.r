## 4. initial analysis related to QC
obj.srt@meta.data %>% ggplot(aes(time, nCount_RNA)) + geom_boxplot() + 
  geom_jitter(size=0.1, alpha=0.5) +
  geom_hline(yintercept = 1000, color='red')

obj.srt[["percent.mt"]] <- PercentageFeatureSet(obj.srt, pattern = "^mt-")
obj.srt@meta.data %>% ggplot(aes(time, percent.mt)) + geom_boxplot() + 
  geom_jitter(size=0.1, alpha=0.5)

obj.srt@meta.data %>% filter(nCount_RNA > 1000) %>% dim()
## 4.1. filter out bad quality cells 
obj.srt <- subset(obj.srt, subset = nCount_RNA > 1000  & percent.mt < 20)
obj.srt@meta.data[1:3,]
obj.srt@meta.data$time <- factor(obj.srt@meta.data$time, levels = c('naive','48hr_','72hr_','day5_','day30'))
obj.srt@meta.data$time %>% table()

## 4.2. initial seurat obj flow
obj.srt <- FindVariableFeatures(obj.srt, selection.method = 'vst',nfeatures = 3000) 
all.genes <- rownames(obj.srt)
obj.srt <- ScaleData(obj.srt, features = all.genes)
obj.srt <- NormalizeData(obj.srt)
obj.srt <- RunPCA(obj.srt, features = VariableFeatures(object = obj.srt), npcs = 10)
DimPlot(obj.srt, reduction = "pca", group.by = 'time')
obj.srt@meta.data[,c(12:19)] <- 'NA'
obj.srt <- FindNeighbors(obj.srt, dims = 1:10)
obj.srt <- FindClusters(obj.srt, resolution = c(0.1,0.2,0.3,0.4,0.5,1))
obj.srt <- RunTSNE(obj.srt, dims = 1:10)
obj.srt <- RunUMAP(obj.srt, dims = 1:10)
Idents(obj.srt) <- 'time'
DimPlot(obj.srt, group.by = 'time', label = T)
DimPlot(obj.srt, group.by = 'time', label = T, split.by = 'time')

DimPlot(obj.srt, group.by = 'RNA_snn_res.0.1', label = T, label.box = F, pt.size = 0.1)
DimPlot(obj.srt, group.by = 'RNA_snn_res.0.5', label = T, label.box = F, pt.size = 0.1)
DimPlot(obj.srt, group.by = 'RNA_snn_res.1', label = T, label.box = F, pt.size = 0.1)
