## cells from each sample distribution-pheatmap
meta <- obj.srt@meta.data
df <- table(meta$sample, meta$RNA_snn_res.0.5)
df <- as.data.frame.matrix(df)
colsum <- colSums(df)
colsum <- replicate(nrow(df), colsum) %>% t()
df2 <- (df/colsum)*100
df2 %>% pheatmap::pheatmap(cluster_rows = F, cluster_cols = F, 
                           display_numbers = T, number_format = "%.0f", fontsize = 15,
                           color = colorRampPalette(colors = c("#EBF5FB","#2874A6"))(100))
meta %>% ggplot(aes(RNA_snn_res.0.5, fill=RNA_snn_res.0.5)) + geom_bar() +
  theme_classic() +theme(legend.position = 'none')

pheatmap::pheatmap(df, cluster_rows = F, cluster_cols = F,
                   color=colorRampPalette(c('#EAF2F8',
                                            '#99CCFF',
                                            '#000099'))(1000))
