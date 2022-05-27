obj.srt@meta.data %>% group_by(RNA_snn_res.0.5, sample) %>% count() %>% 
  ggplot(aes(RNA_snn_res.0.5,n)) + geom_bar(aes(fill=sample),
                                                              stat = 'identity',
                                                              position = position_dodge()) +
  theme_bw()
