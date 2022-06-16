
library(Seurat)
library(dplyr)

##### Marker gene analysis #########
sample <- 'imm_timecourse'
res <- 'time_cluster'
obj.tmp <- obj.srt
Idents(obj.tmp) <- res
all.markers <- FindAllMarkers(obj.tmp, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = log2(1.5))
all.markers[1:3,]
all.markers %>% dim()
all.markers %>% select(cluster) %>% table()

all.markers %>% write.csv(paste0('2022.cell_hashing/',sample, '_',res,'.csv'))

all.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC) %>% 
  write.csv(paste0('2022.cell_hashing/',sample, '_',res,'.top20.csv'))


### find markers between two sets
markers <- FindMarkers(obj.srt, ident.1 = c('a05','a30'), ident.2 = c('i05','i30')
                       ,min.pct = 0.25, only.pos = FALSE)
markers %>% dim()
markers %>% write.csv('2022.cell_hashing/arm_imm.day30.comparsion.csv')
markers %>% ggplot(aes(avg_log2FC, -log10(p_val_adj))) + 
  geom_point(size=1, alpha=0.5) + geom_hline(yintercept = -log10(0.05), color='red') +
  theme_classic()

markers %>% ggplot(aes(avg_log2FC, -log10(p_val_adj), label=rownames(markers))) + 
  geom_point(size=1, alpha=0.5) + geom_hline(yintercept = -log10(0.05), color='red') +
  theme_classic() + geom_text()

## label genes in the plot
genes <- c('Nkg7','Gzmb','Lgals3','Cxcr3','Batf',
           'Mki67','Ifng','Tigit','Pdcd1', 'Icos',
           'Il12rb2','Cdk6','Tbx21','Tnfsf8','Irf7',
           'Bcl2','Nfkbiz','Il21','Runx2','Stat4',
           'Tox2','Cxcr5','Ezh2','Ikzf2','Btla',
           'Il16','Cd4','Lef1','Jarid2','Junb',
           'Il10ra','Tmtc2','Aff3','Bach2','Il6ra',
           'Il7r','Foxp1','Id3','Btg1','Ccr7',
           'Sell','Id2','Cd226')

markers$gene <- ''
markers[genes,]$gene <- genes
markers %>% ggplot(aes(avg_log2FC, -log10(p_val_adj), label=gene)) + 
  geom_point(size=1, alpha=0.1) + geom_hline(yintercept = -log10(0.05), color='red') +
  theme_classic() + geom_text(vjust = 0.5, nudge_y = -0.5) + geom_vline(xintercept = c(-log2(1.5),log2(1.5)), color='red')



### find markers between two sets
Idents(obj.srt) <- 'cell_type'
markers <- FindMarkers(obj.srt, ident.1 = 'Tem', ident.2 = 'Tcm' ,min.pct = 0.25, only.pos = FALSE)
markers %>% dim()
markers %>% filter(avg_log2FC > 0) %>% head()

markers %>% ggplot(aes(avg_log2FC, -log10(p_val_adj))) + 
  geom_point(size=1, alpha=0.5) + geom_hline(yintercept = -log10(0.05), color='red') +
  geom_vline(xintercept = c(-log2(1.5),log2(1.5)), color='red') +
  xlim(c(-2,2)) +
  theme_classic()

genes <- c((markers %>% filter(avg_log2FC > 0) %>% top_n(20, avg_log2FC) %>% rownames()),
           (markers %>% filter(avg_log2FC < 0) %>% top_n(-40, avg_log2FC) %>% rownames()))
DoHeatmap(obj.srt, features = genes, group.by = 'cell_type') + viridis::scale_fill_viridis() + 
  theme(axis.text.y = element_text(size = 7))



## PPT figure update 22.06.16

genes <- c('Il7r','Cxcr3','Ccr7','Sell','Tcf7','Gzma','Klrg1','Cx3cr1','Bhlhe40', 'Gzmb','Ccl5', 'Zeb2')
FeaturePlot(obj.srt, features = genes, pt.size = 0.1, ncol =3 )
FeaturePlot(obj.srt, features = 'pseudotime')

tmp.genes <- read.csv('~/Desktop/tmp.csv', row.names = 1)
blue <- tmp.genes$Tcm %>% as.vector()
red <- tmp.genes$Tem %>% as.vector()

markers <- read.csv('2022.Qiang.Tle3.scRNA/no.trimming.wt_only.1000.2.5.Tcm.Tem.mk.csv',
                    row.names = 1)
markers[1:5,]
markers$color <- 'no significant'
genes <- markers %>% filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05) %>% rownames()
markers[genes,]$color <- 'sig'
genes <- markers %>% filter(avg_log2FC < -log2(1.5) & p_val_adj < 0.05) %>% rownames()
markers[genes,]$color <- 'sig'
markers[blue,]$color <- 'blue'
markers[red,]$color <- 'red'
markers %>% ggplot(aes(avg_log2FC, -log10(p_val_adj), color=color)) + 
  geom_point(size=1, alpha=0.5) + geom_hline(yintercept = -log10(0.05), color='red') +
  geom_vline(xintercept = c(-log2(1.5),log2(1.5)), color='red') +
  xlim(c(-2,2)) +
  theme_classic() + scale_color_manual(values = c('blue','grey','red','black'))

markers$gene <- ''
genes <- c(blue,red)
markers[genes,]$gene <- genes
markers %>% ggplot(aes(avg_log2FC, -log10(p_val_adj), label=gene)) + 
  geom_point(size=1, alpha=0.1, aes(color=color)) + geom_hline(yintercept = -log10(0.05), color='red') +
  theme_classic() + geom_text(vjust = 0.5, nudge_y = -0.5) + 
  geom_vline(xintercept = c(-log2(1.5),log2(1.5)), color='red') + 
  scale_color_manual(values = c('blue','grey','red','black')) +xlim(c(-2,2))

genes <- c(tmp.genes[,1],tmp.genes[,2],tmp.genes[,3])
DoHeatmap(obj.srt, features = genes, group.by = 'Identity') + 
  theme(axis.text.y = element_text(size = 10)) 

