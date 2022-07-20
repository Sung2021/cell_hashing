df.magic <- Rmagic::magic(t(obj.srt@assays$RNA@data)) 

df.magic$result %>% dim()
df.magic$result %>% ggplot(aes(Klf2,Cxcr5, color=obj.srt@meta.data$RNA_snn_res.0.5)) + 
geom_point(size=1, alpha=0.5) + facet_grid(obj.srt@meta.data$dpi~obj.srt@meta.data$mut) + theme_bw()


df.magic <- readRDS('2022.cell_hashing/icos_mt_set/icos_mt.magic.22.07.20.rds')

## scatter plot by gene1 and gene2

g1 ='Cxcr5'
g2='Klf2'
## basic
df.magic %>% ggplot(aes(get(g1),get(g2))) + geom_point() + xlab(g1) + ylab(g2)

