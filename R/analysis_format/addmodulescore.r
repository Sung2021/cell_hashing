module <-  c('Il7r','Ly6e','Ltb','Cxcr3','Ly6a',
             'Ccr7','Sell','Rack1','Myc','Tcf7')
obj.srt <- AddModuleScore(obj.srt, features = list(module), name = 'Tcm' )
module <- c('Bhlhe40','Gzma','Gzmb','Klrc1','Klrg1',
            'S1pr5','Ccl5','Cx3cr1','Gzmk','Klre1')
obj.srt <- AddModuleScore(obj.srt, features = list(module),name = 'Tem' )

FeaturePlot(obj.srt, features = c('Tcm','Tem'), pt.size = 0.2)
VlnPlot(obj.srt, features = c('Tcm','Tem'),
        stack = T, flip = T, group.by = 'RNA_snn_res.0.5') +theme(legend.position = 'none')
