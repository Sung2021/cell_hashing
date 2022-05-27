library(clusterProfiler)
library(org.Mm.eg.db)
geneset <- all.markers[all.markers$cluster == '1',] %>% rownames()
genes_to_convert <- bitr(geneset, fromType = "SYMBOL", 
                         toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
gene.ego <- enrichGO(gene     = genes_to_convert$ENTREZID, 
                     ont      = "BP",
                     OrgDb = "org.Mm.eg.db",
                     pvalueCutoff = 0.01, qvalueCutoff = 0.01, readable = T)
barplot(gene.ego)

geneset <- all.markers[all.markers$cluster == '1',] %>% rownames()
genes_to_convert <- bitr(geneset, fromType = "SYMBOL", 
                         toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
gene.kk <- enrichKEGG(genes_to_convert$ENTREZID,
                      organism = "mmu", 
                      keyType = "kegg", pvalueCutoff = 0.05, 
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.2, 
                      use_internal_data = FALSE)


gene.kk@result[1:5,]
barplot(gene.kk) + theme(axis.text.y = element_text(size = 8)) + ggtitle('KEGG')


gene.ego@result[5:6,c("Description",'geneID')]


