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

