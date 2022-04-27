
## 5 sets including naive
my_colors <- c("#D2B4DE","#A9CCE3", "#2E86C1", "#F1C40F", "#D35400") # Create vector of colors
names(my_colors) <- levels(obj.srt@meta.data$time) 


meta <- obj.srt@meta.data
meta$hhx_cluster_name %>% table()
meta$cluster.0425 <- factor(meta$hhx_cluster_name)
levels(meta$cluster.0425)
new.levels <-  c('naive','activated','trans_I','trans_II','memory_CD4',
                 'Tfh','non_Tfh_I','non_Tfh_II','non_Tfh_III',
                 'MP_I','MP_II')
levels(meta$cluster.0425) <- new.levels

reduction <- obj.srt@reductions$umap@cell.embeddings %>% data.frame(check.names = F)
reduction %>% colnames()

#####################################################################################
my_colors <- c("#81C784",
               '#FBC02D',
               '#F57C00','#FF9800',
               '#F06292',
               '#0277BD',
               '#B71C1C','#D32F2F','#EF5350',
               '#039BE5','#81D4FA') # Create vector of colors
names(my_colors) <- levels(obj.srt@meta.data$cluster.0425) # Extract all levels of both data

reduction %>% ggplot(aes(UMAP_1, UMAP_2, color=meta$cluster.0425)) + 
  geom_point(size=0.5, alpha=0.6) + theme_classic() +
  scale_color_manual(values = my_colors)
