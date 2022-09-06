ggvenn::ggvenn(list(c11=c11,
                    c16=c16,
                    c9=c9))
Tcell.list = list(Tcm.enr=Tcm.enr,
                  Tem.enr=Tem.enr,
                  Tem.ko.up=Tem.ko.up,
                  Tem.ko.dn=Tem.ko.dn,
                  Tcm.ko.up=Tcm.ko.up,
                  Tcm.ko.dn=Tcm.ko.dn)
ggvenn::ggvenn(Tcell.list)



c11_over_c16 <- setdiff(c11,c16)
all.markers %>% filter(gene %in% c11_over_c16) %>% write.csv('~/Desktop/c11_over-c16.csv')

c16_over_c11 <- setdiff(c16,c11)

all.markers %>% filter(gene %in% c16_over_c11) %>% write.csv('~/Desktop/c16_over-c11.csv')


venn.output <- gplots::venn(list(c11=c11,
                  c16=c16,
                  c9=c9), intersections=TRUE)

venn.output %>% str()
venn.output <- venn.output %>% list()
lengths(attributes(venn.output)$intersections)
attributes(venn.output)$intersections$`c11:c9`
