### functions

## color gene manually in tSNE 
tsnegene <- function(gene, obj.srt=obj.srt){
  meta <- obj.srt@meta.data
  meta[,gene] <- obj.srt@assays$RNA@data[gene,] 
  reduction <- obj.srt@reductions$tsne@cell.embeddings %>% data.frame(check.names = F)
  p <- reduction %>% ggplot(aes(tSNE_1, tSNE_2, color=meta[,gene])) + 
    geom_point(size=0.1) + theme_classic()+ 
    scale_color_gradient(low = '#FEF9E7',
                         high= '#A93226') +ggtitle(gene)
  print(p)
}

tsnegene('Runx1')

umapgene <- function(gene){
  meta <- obj.srt@meta.data
  meta[,gene] <- obj.srt@assays$RNA@data[gene,] 
  reduction <- obj.srt@reductions$umap@cell.embeddings %>% data.frame(check.names = F)
  p <- reduction %>% ggplot(aes(UMAP_1, UMAP_2, color=meta[,gene])) + 
    geom_point(size=0.1) + theme_classic()+ 
    scale_color_gradient(low = '#FEF9E7',
                         high='#A93226') +ggtitle(gene)
  print(p)
}

p1 <- umapgene('Aff3')
p2 <- umapgene('Il21')
p3 <- umapgene('Il4')
cowplot::plot_grid(p1,p2,p3,ncol=1)


## color gene manually in diffusion map
# dm : diffusion map object 
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2])
                  
dmgene <- function(gene){
  meta[,gene] <- obj.srt@assays$RNA@data[gene,] 
  p <- ggplot(tmp, aes(x = DC1, y = DC2, colour = meta[,gene])) +
    geom_point(size= 1, alpha=0.5) +
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    scale_color_gradient(low = '#FEF9E7',
                         high='#A93226') +ggtitle(gene)
  print(p)
}

dmgene('Tigit')


dm_pseudotime_gene <- function(gene, color='time'){
  meta <- obj.srt@meta.data
  meta[,gene] <- obj.srt@assays$RNA@data[gene,] 
  tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                    DC2 = eigenvectors(dm)[, 2])
  meta <- cbind(meta,tmp)
  meta %>% ggplot(aes(DC1, get(gene), color=get(color))) + 
    geom_point(size=1, alpha=0.5) + theme_classic() +ylab(gene)
}
dm_pseudotime_gene('Cxcr5', color='time')



naive <- c('3','6')
early_activated <- c('7','10','5','14')
early_transitional <- c('8')
later_transtional <- c('12','2','4')
mature_Tfh <- c('0')
non_Tfh_eff_memory_transitional <- c('13')
non_Tfh_eff_memory_CD4 <- c('1')
non_Tfh_eff_super_memory <- c('15')
Tfh_mp2 <- c('16')
Tfh_mp3 <- c('11')

umap_groups <- function(groups){
  meta <- obj.srt@meta.data
  reduction <- obj.srt@reductions$umap@cell.embeddings %>% data.frame(check.names = F)
  meta$select <- 'NA'
  meta[meta$hhx_cluster %in% c(groups),]$select <- 'selected'
  meta$select <- factor(meta$select)
  p <- reduction %>% ggplot(aes(UMAP_1, UMAP_2, color=meta$select)) + 
    geom_point(size=0.1) + theme_classic() +scale_color_manual(values = c('#EBEDEF',
                                                                          '#EC7063')) +
    theme(legend.position = "None")
  print(p)
}

umap_groups(mature_Tfh)
umap_groups(Tfh_mp2)


umap_factor <- function(group){
  meta <- obj.srt@meta.data
  reduction <- obj.srt@reductions$umap@cell.embeddings %>% data.frame(check.names = F)
  meta$select <- 'NA'
  meta[meta$hhx_cluster %in% c(3,6),]$select <- 'selected'
  p <- reduction %>% ggplot(aes(UMAP_1, UMAP_2, color=meta[,]$select)) + 
    geom_point(size=0.1) + theme_classic()+ 
   ggtitle(group)
  print(p)
}
umap_factor('time')


## gene expression line plot by group
gene_hhx_cluster_time <- function(gene){
  meta <- obj.srt@meta.data
  meta[,gene] <- obj.srt@assays$RNA@data[gene,] 
  df <- meta %>%
    group_by(hhx_cluster_name, time) %>%
    summarize(mean_size = mean(get(gene), na.rm = TRUE),
              sd= sd(get(gene)))
  df %>% ggplot(aes(x=time, y=mean_size, group=hhx_cluster_name, 
                    color=hhx_cluster_name)) + 
    geom_line() + geom_point() + RotatedAxis() +ylab(gene) +
    geom_errorbar(aes(ymin=mean_size-sd, ymax=mean_size+sd), width=.1) + facet_wrap(.~hhx_cluster_name)
}

gene_hhx_cluster_time('Jmjd1c') # Kdm6a
gene_hhx_cluster_time('Kdm6a')


genes <- c('Tox','Aff3','Runx1','Cd226','Cxcr5','Tigit','Bcl6','Nr4a1','Egr1','Icos','Pdcd1','Ctla4','Bach2','Ikzf2')
genes <- c('Tox','Aff3','Runx1','Bach2')
genes <- c('Jarid2','Tmtc2','Hells','Cdk1')
genes <- c('Ccne2','Cdk2')
genes <- c('Sesn3','Depdc5')
genes <- c('Jmjd1c','Kdm6a')
genes <- c('Nsd2','Hdac7','Ezh2')
genes <- c('Nr3c1','Nr4a1','Nr4a2','Nr4a3')
genes <- c('Fkbp5','Tsc22d2','Tsc22d3')
cowplot::plot_grid(plotlist=lapply(genes,gene_hhx_cluster_time),ncol = 2)

