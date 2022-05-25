## sample statistics
obj.srt@meta.data$sample %>% table()
meta <- obj.srt@meta.data

## umi
meta %>% ggplot(aes(sample,nCount_RNA)) + geom_boxplot() + 
  geom_jitter(size=0.1, alpha=0.2) + geom_hline(yintercept = c(200,1000,5000,10000),color='red') +
  RotatedAxis()
## number of detected genes
meta %>% ggplot(aes(sample,nFeature_RNA)) + geom_boxplot() + 
  geom_jitter(size=0.1, alpha=0.2) + geom_hline(yintercept = c(200,1000,7000),color='red') +
  RotatedAxis()
## number of mitochondrial population
meta %>% ggplot(aes(sample,percent.mt)) + geom_boxplot() + 
  geom_jitter(size=0.1, alpha=0.2) + geom_hline(yintercept = c(2.5,5,10),color='red') +
  RotatedAxis()

## quality comparisons
meta %>% ggplot(aes(nCount_RNA,percent.mt,
                    color=sample)) + geom_point(size=1,alpha=0.5) + 
  facet_wrap(.~sample, ncol = 5) + 
  geom_hline(yintercept = c(2.5,5,10),color='red') +
  geom_vline(xintercept = c(1000,5000,10000),color='red') + 
  RotatedAxis()

meta %>% ggplot(aes(nCount_RNA,percent.mt,
                    color=sample)) + geom_point(size=1,alpha=0.5) + 
  facet_wrap(.~sample, ncol = 5) + 
  geom_hline(yintercept = c(2.5,5,10),color='red') +
  geom_vline(xintercept = c(1000,5000,10000),color='red') + 
  ylim(c(-2.5,20)) + xlim(c(0,60000)) +theme_bw() +
  RotatedAxis()

## check by combinations  
mt.filter <- 5
read.filter <- 2000
meta %>% filter(!(sample %in% c('ARM_d5','ARM_d30'))) %>% 
  filter(percent.mt <= mt.filter) %>% filter(nCount_RNA >= read.filter) %>% 
  ggplot(aes(nCount_RNA,percent.mt,
             color=sample)) + geom_point(size=1,alpha=0.5) + 
  facet_wrap(.~sample, nrow = 1) + 
  geom_hline(yintercept = c(mt.filter),color='red') +
  geom_vline(xintercept = c(read.filter),color='red') + 
  ylim(c(-1,10)) + xlim(c(0,60000)) +theme_bw() +
  RotatedAxis()
## number of cells in each sample 
meta %>% filter(!(sample %in% c('ARM_d5','ARM_d30'))) %>% 
  filter(percent.mt <=mt.filter) %>% filter(nCount_RNA >= read.filter) %>% 
  select(sample) %>% table()
