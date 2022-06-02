genes <- all.markers %>% group_by(cluster) %>% top_n(1, avg_log2FC)
VlnPlot(obj.srt, genes$gene, stack = TRUE, sort = F, flip = TRUE, 
        group.by = 'RNA_snn_res.0.5') +
  theme(legend.position = "none")
