### color palette for samples
my_colors <- c("#D2B4DE","#A9CCE3", "#2E86C1", "#F1C40F", "#D35400") # Create vector of colors
names(my_colors) <- levels(obj.srt@meta.data$time) # Extract all levels of both data

obj.srt@meta.data %>% ggplot(aes(RNA_snn_res.1, fill=time)) + 
  geom_bar() + scale_fill_manual(values = my_colors) + theme_classic()

obj.srt@meta.data %>% ggplot(aes(time, fill=time)) + 
  geom_bar() + scale_fill_manual(values = my_colors) + theme_classic() +RotatedAxis()

reduction <- obj.srt@reductions$umap@cell.embeddings %>% data.frame(check.names = F)
reduction %>% colnames()
meta <- obj.srt@meta.data
reduction %>% ggplot(aes(UMAP_1, UMAP_2, color=meta$time)) + 
  geom_point(size=0.1) + theme_classic()+ scale_color_manual(values = my_colors)

reduction %>% ggplot(aes(UMAP_1, UMAP_2, color=meta$time)) + 
  geom_point(size=0.1) + theme_classic()+ scale_color_manual(values = my_colors)



## cells from each sample distribution-pheatmap
df <- table(meta$time, meta$RNA_snn_res.1)
df <- as.data.frame.matrix(df)
colsum <- colSums(df)
colsum <- replicate(5, colsum) %>% t()
df2 <- (df/colsum)*100
# my.colors <- c(colorRampPalette(colors = c("white","orange"))(20), colorRampPalette(colors = c("orange","red"))(30))
df2 %>% pheatmap::pheatmap(cluster_rows = F, cluster_cols = F, 
                           color = colorRampPalette(colors = c("#EBF5FB","#2874A6"))(100))


pheatmap::pheatmap(df, cluster_rows = F, cluster_cols = F,
                   color=colorRampPalette(c('#EAF2F8',
                                            '#99CCFF',
                                            '#000099'))(1000))

## reorder the cluster
meta <- obj.srt@meta.data
naive <- c(3,6)
hr48 <- c(7)
hr72 <- c(4)
d5 <- c(1,9)
d30 <- c(0,15)
t_late <- c(13)
t_early <- c(11,16,5,10,14,8,12,2)

new.order <- c(naive,hr48,hr72,d5,d30,t_late,t_early)
meta$time_cluster <- factor(meta$RNA_snn_res.1, 
                            levels = new.order)
obj.srt[['time_cluster']] <- meta$time_cluster
meta %>% ggplot(aes(time_cluster, fill=time)) + 
  geom_bar() + scale_fill_manual(values = my_colors) + theme_classic()

meta %>% ggplot(aes(time_cluster, fill=time)) + 
  geom_bar(position = 'fill') + scale_fill_manual(values = my_colors) + 
  theme_classic() +coord_flip() +scale_x_discrete(limits = rev(levels(meta$time_cluster)))


df <- table(meta$time, meta$time_cluster)
df <- as.data.frame.matrix(df)
colsum <- colSums(df)
colsum <- replicate(5, colsum) %>% t()
df2 <- (df/colsum)*100
# my.colors <- c(colorRampPalette(colors = c("white","orange"))(20), colorRampPalette(colors = c("orange","red"))(30))
df2 %>% pheatmap::pheatmap(cluster_rows = F, cluster_cols = F, 
                           color = colorRampPalette(colors = c("#EBF5FB","#2874A6"))(100))

## cells from each sample distribution-pheatmap
meta <- obj.srt@meta.data
df <- table(meta$tag, meta$RNA_snn_res.0.5)
df <- as.data.frame.matrix(df)
colsum <- colSums(df)
colsum <- replicate(nrow(df), colsum) %>% t()
df2 <- (df/colsum)*100
# my.colors <- c(colorRampPalette(colors = c("white","orange"))(20), colorRampPalette(colors = c("orange","red"))(30))
df2 %>% pheatmap::pheatmap(cluster_rows = F, cluster_cols = F, 
                           color = colorRampPalette(colors = c("#EBF5FB","#2874A6"))(100))
df2 %>% pheatmap::pheatmap(cluster_rows = F, cluster_cols = F, 
                           color = colorRampPalette(colors = c("#EBF5FB","#2874A6"))(100),
                           display_numbers = T)

