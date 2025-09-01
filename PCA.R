install.packages("ggplot2") 
library(ggplot2)
pca <- read.csv('C:\\Users\\wutianqin\\Desktop\\ppt\\plink.eigenvec.csv',header = T)
gxz <- read.csv('C:\\Users\\wutianqin\\Desktop\\ppt\\plink.eigenval.csv',header = T)
xlab <- paste0("PC1(", round(gxz[1,2] * 100, 2), "%)")
ylab <- paste0("PC1(", round(gxz[2,2] * 100, 2), "%)")
ggplot(pca, aes(x= PC1, y= PC2, color=group))+    
geom_point(size=2)+    
theme_bw()+    
theme(panel.grid = element_blank())+    
labs(x=xlab,y=ylab)+    
stat_ellipse(level=0.95, linetype = 2)+
scale_color_manual(values = c("#FF9900","#68cc99","#EF566B","#3E91D2"))
