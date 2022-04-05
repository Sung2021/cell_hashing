## 1. prepare the gene expression dataset

obj.srt <- obj.early
gene.exp <- obj.srt@assays$RNA@scale.data[VariableFeatures(obj.srt),]
gene.exp %>% dim()
## 2. WGCNA correlation run
test.cor <- WGCNA::cor(gene.exp %>% t(), method = 'pearson')
test.cor <- WGCNA::cor(gene.exp %>% t())
test.cor <- test.cor %>% data.frame()

## 
test.cor <- read.csv('rds/2022.cell_hashing/correlation.early.all.csv',
                     row.names = 1)
test.cor[1:3,1:3]

test.cor %>% ggplot(aes(Cxcr5)) + geom_density()


gene1 <- 'Cxcr5'
gene2 <- 'Tox2'
gene1 <- 'Tigit'
gene2 <- 'Cd226'
test.cor %>% ggplot(aes(get(gene1), get(gene2), label=rownames(test.cor))) + geom_point(size=1, alpha=0.2) + 
  xlab(gene1) +ylab(gene2) + 
  xlim(c(-0.25,0.5)) +ylim(c(-0.25,0.5)) +geom_text(size=4)


naive_like <- c('Id3','Ccr7','Il7r','Lef1','Tcf7','Slamf6','Bcl2')
Th1 <- c('Cxcr6','Tbx21','Id2','Gzmb','Ifng','Ly6c2') 
Tfh <- c('Cxcr5','Bcl6','Ascl2','Pdcd1','Icos','Il21','Il4','Tigit')
Treg <- c('Foxp3','Il2ra','Il10','Cd81','Cd74','Klrg1')

genes <- Tfh[Tfh %in% rownames(test.cor)]
test.cor[genes.row, genes] %>% write.csv('~/Desktop/correlation.tmp.csv')
test.cor %>% select(Tox) %>% top_n(100)
test.cor %>% select(Tox) %>% ggplot(aes(Tox)) + geom_density()

