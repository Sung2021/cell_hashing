
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
