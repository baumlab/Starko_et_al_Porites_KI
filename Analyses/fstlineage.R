
setwd("C:/Users/james/My Drive/Documents/BOSTON/Davies/Baum/Code/FST")
install.packages("rgr")
library(rgr)
data=read.csv(file="fst.lineage3pop.csv")
mat=as.matrix(read.csv(file="fst.lineage3pop.csv")[, -1])
apply(mat, 2, as.numeric)
remove.na(mat, iftell = TRUE)
rownames(mat) <- data$ï..
colnames(mat)=colnames(data[-1])
colnames(data[-1])
colnames(mat)
mat
heatmap(mat)

library(ggplot2)

# Dummy data
library(reshape2)

data=setNames(melt(mat), c('X', 'Y', 'Z'))

df = data[complete.cases(data),]


# Heatmap 
sessionInfo()
library(ggplot2)
ggplot(df, aes(X, Y, fill= Z)) + 
  geom_tile((aes(fill = Z)), colour="black")+
  geom_text(aes(label = Z),size=8 ) +
  theme_classic(base_size = 22)+
  scale_fill_gradient2(low = "white", high="orange", midpoint = .28) +
  guides(fill=guide_legend(title="Fst"))+
  xlab(label = "")+ ylab(label="")

?scale_fill_gradient
