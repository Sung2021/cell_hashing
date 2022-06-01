# Create the knee plot 
# Convert to dgCMatrix, which is a compressed, sparse matrix format
mat<- as(t(obj.srt@assays$RNA@counts), "dgCMatrix")
tot_counts <- rowSums(mat)  ### UMI per cell
Rank <- rank(-tot_counts)
df <- cbind(tot_counts, Rank)
df <- data.frame(df)
df[1:3,]
ggplot(df, aes(tot_counts,Rank)) + geom_point() +
  scale_x_log10() + scale_y_log10() + annotation_logticks() +
  labs(y = "Barcode rank", x = "Total UMI count")

ggplot(df, aes(tot_counts,Rank)) + geom_point() +
  scale_x_log10() +  annotation_logticks() +
  labs(y = "Barcode rank", x = "Total UMI count") + 
  geom_vline(xintercept = c(300,400,1000), color='red')
