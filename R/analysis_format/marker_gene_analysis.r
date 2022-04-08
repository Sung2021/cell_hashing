
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

