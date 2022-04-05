####### diffusion map #########

## load libraries
library(destiny)
library(scater)
library(ggbeeswarm)

## prepare data
obj.srt <- input.data
deng <- log10(obj.srt@assays$RNA@counts+1)
deng <- deng[VariableFeatures(obj.srt),]
dm <- DiffusionMap(data.frame(t(deng), check.names = F), n_pcs = 50)
dm@eigenvectors[1:3,]

## generate data frame with diffusion map
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2])

ggplot(tmp, aes(x = DC1, y = DC2, colour = obj.srt@meta.data$time_cluster)) +
  geom_point(size= 1, alpha=0.5) +
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +facet_wrap(.~obj.srt@meta.data$time_cluster)

## save data
dm %>% saveRDS('2022.cell_hashing/imm.time_course.2022.04.05.diffusionmap.rds')

## color them
gene <- c('Cxcr5')
meta[,gene] <- obj.srt@assays$RNA@data[gene,] 
ggplot(tmp, aes(x = DC1, y = DC2, colour = meta[,gene])) +
  geom_point(size= 1, alpha=0.5) +
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  scale_color_gradient(low = '#FEF9E7',
                        high='#A93226') +ggtitle(gene) +theme_classic()
                        
## color gene in diffusion map
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
