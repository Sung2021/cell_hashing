library(phateR)


df.phate <- phate(t(obj.srt@assays$RNA@data),
                  t=5) 
p1 <- ggplot(df.phate) +
  geom_point(aes(PHATE1, PHATE2, color=obj.srt@meta.data$time))
df.phate <- phate(t(obj.srt@assays$RNA@data),
                  t=7) 
p2 <- ggplot(df.phate) +
  geom_point(aes(PHATE1, PHATE2, color=obj.srt@meta.data$time))

df.phate <- phate(t(obj.srt@assays$RNA@data),
                  t=3) 
p3 <- ggplot(df.phate) +
  geom_point(aes(PHATE1, PHATE2, color=obj.srt@meta.data$time))
cowplot::plot_grid(p3,p1,p2)

df.phate <- phate(t(obj.srt@assays$RNA@data),
                  t=10) 
p4 <- ggplot(df.phate) +
  geom_point(aes(PHATE1, PHATE2, color=obj.srt@meta.data$time))

df.phate <- phate(t(obj.srt@assays$RNA@data),
                  t=20) 
p5 <- ggplot(df.phate) +
  geom_point(aes(PHATE1, PHATE2, color=obj.srt@meta.data$time))
cowplot::plot_grid(p3,p1,p2,p4,p5,ncol = 3)


df.phate.t20 <- df.phate
df.phate.t20 %>% saveRDS('2022.cell_hashing/imm_timecourse/phate.t20.2022.04.12.rds')




ggplot(df.phate) +
  geom_point(aes(PHATE1, PHATE2, color=obj.srt@meta.data$hhx_cluster_name), size=0.2) + theme_classic() +
  facet_wrap(.~obj.srt@meta.data$hhx_cluster_name)

meta <- obj.srt@meta.data
gene <- 'Runx1'
meta[,gene] <- obj.srt@assays$RNA@data[gene,] 
ggplot(df.phate) +
  geom_point(aes(PHATE1, PHATE2, color=df.magic$result[,gene]), size=0.2) + theme_classic() +
  scale_color_gradient(low = '#FEF9E7',
                       high='#A93226')


df.magic <- Rmagic::magic(t(obj.srt@assays$RNA@data)) 

df.magic$result %>% dim()
df.magic$result[1:3,1:3]

df.magic$result %>% ggplot(aes(Cd226,Tigit)) + geom_point()
