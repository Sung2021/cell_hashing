## https://plotly.com/r/3d-scatter-plots/
library(plotly)

mtcars$am[which(mtcars$am == 0)] <- 'Automatic'
mtcars$am[which(mtcars$am == 1)] <- 'Manual'
mtcars$am <- as.factor(mtcars$am)

fig <- plot_ly(mtcars, x = ~wt, y = ~hp, z = ~qsec, color = ~am, colors = c('#BF382A', '#0C4B8E'))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Weight'),
                                   yaxis = list(title = 'Gross horsepower'),
                                   zaxis = list(title = '1/4 mile time')))

fig

###############################################
df.tmp <- obj.monocle@reducedDimS %>% t() %>% data.frame()
df.tmp[1:3,]
obj.meta <- obj.all@meta.data
obj.meta[1:3,]
ggplot(df.tmp, aes(X1,X2, color=obj.all@meta.data$RNA_snn_res.0.5)) + 
  geom_point(size=0.2)
plot_ly(df.tmp, 
               x=~X1, y= ~X2, z= ~X3,
               color = obj.meta$group.v2, size = 1)
plot_ly(df.tmp, 
        x=~X1, y= ~X2, z= ~X3,
        color = obj.meta$tag, alpha = 0.8, size = 1)

plot_ly(df.tmp, 
        x=~X3, y= ~X4, z= ~X1,
        color = obj.meta$group.v2, alpha = 0.9, size = 1)
