library(monocle3)
library(SingleCellExperiment)
library(SeuratWrappers)

cds <- as.cell_data_set(obj.srt)
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

## Step 6: Order cells
cds <- order_cells(cds)

ps_tim <- pseudotime(cds, reduction_method = 'UMAP')
plot_cells(cds, color_cells_by = 'sample')
plot_cells(cds, color_cells_by = 'pseudotime',
           graph_label_size = 0, show_trajectory_graph = T)
obj.srt[['ps_tim']] <- ps_tim
FeaturePlot(obj.srt, features = 'ps_tim')
VlnPlot(obj.srt, features = 'ps_tim', group.by = 'RNA_snn_res.0.5', sort = T)
VlnPlot(obj.srt, features = 'ps_tim', group.by = 'sample', sort = T)

cds %>% saveRDS('2022.cell_hashing/cutoff.2000.5.activated.only.22.06.02.cds.22.06.03.rds')
