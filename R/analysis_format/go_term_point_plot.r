## Howard request
df <- read.csv('Tle3_integration_2022/Great_analysis/Great_analysis.22.06.15.csv')
df[1:3,]
cluster <- 'Cluster3'
df2 <- df %>% filter(group ==cluster) 
df2$order <- rownames(df2)
dev.off()
df2 %>% data.frame() %>% gridExtra::grid.table()
range(-log10(df2$Hyper.FDR.Q.Val))
df2$p_val <- (-log10(df2$Hyper.FDR.Q.Val) - min(-log10(df2$Hyper.FDR.Q.Val)))
df2
df2$rank <- df2$Hyper.Observation.Gene.Hits %>% rank()
col.val <- c('#66FF00','#FF3300')
df2 %>% ggplot(aes(Hyper.Observation.Gene.Hits, rank,
                   color= p_val)) + geom_point(size=5) +
                      ggtitle(cluster) +
  scale_color_gradient(low = col.val[1],high = col.val[2])
ggsave('test.png', width = 5, height = 4, units = 'in')

########################################################
