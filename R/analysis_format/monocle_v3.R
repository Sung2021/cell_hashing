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

####################################################
####################################################
### monocle trajectory analysis 

library(monocle3)
library(SingleCellExperiment)
library(SeuratWrappers)

cds <- as.cell_data_set(obj.srt)
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

## Step 6: Order cells
cds <- order_cells(cds)

ps_tim <- pseudotime(cds, reduction_method = 'UMAP')
plot_cells(cds, color_cells_by = 'Identity', 
           graph_label_size = 0, show_trajectory_graph = T,
           cell_size = 0.8)+ 
  scale_color_manual(values=c('#999999','#A93226', '#56B4E9', '#D68910','#BB8FCE'))
plot_cells(cds, color_cells_by = 'pseudotime',
           graph_label_size = 0, show_trajectory_graph = T,
           cell_size = 0.8)


cds = readRDS('2022.Qiang.Tle3.scRNA/no.trimming.wt_only.1000.2.5.2022.06.15_cds.rds')

plot_cells(cds, color_cells_by = 'pseudotime',
           graph_label_size = 0, show_trajectory_graph = T,
           cell_size = 0.8, trajectory_graph_color = 'light green', 
           trajectory_graph_segment_size = 2, label_branch_points = F,
           label_roots = F, label_leaves = F)

cds %>% saveRDS('2022.Qiang.Tle3.scRNA/no.trimming.wt_only.1000.2.5.2022.06.15_cds.rds')



