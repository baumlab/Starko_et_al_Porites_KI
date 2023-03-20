setwd("C:/Users/james/My Drive/Documents/BOSTON/Davies/Baum/Code/Demo")
#mij - migration rate from population j to population i.
library(reshape2)
library(dplyr)
library(gridExtra)
library(ggplot2)

#12
pop12=read.table(file="boot.pop1pop2.out.params")
colnames(pop12) <- paste0('pop1_pop2', colnames(pop12))
meltpop12=melt(pop12)
p <- ggplot(meltpop12, aes(factor(variable), value)) 
p + geom_boxplot()+  facet_wrap(~variable, scale="free")
mig12=pop12[,1:4]
meltmig12=melt(mig12)
p <- ggplot(meltmig12, aes(factor(variable), value)) 
p + geom_boxplot() + ylim(0,1)+scale_y_continuous(trans='log2') +facet_wrap(~variable)
#more going from 2 to 1 
#13
pop13=read.table(file="boot.pop1pop3.out.params")
colnames(pop13) <- paste0('pop1_pop3', colnames(pop13))
meltpop13=melt(pop13)
p <- ggplot(meltpop13, aes(factor(variable), value)) 
p + geom_boxplot()+  facet_wrap(~variable, scale="free")
mig13=pop13[,1:4]
meltmig13=melt(mig13)
p <- ggplot(meltmig13, aes(factor(variable), value)) 
p + geom_boxplot() + ylim(0,1)+scale_y_continuous(trans='log2') +facet_wrap(~variable)
#more 3 to 1
#23
pop23=read.table(file="boot.pop2pop3.out.params")
colnames(pop23) <- paste0('pop2_pop3', colnames(pop23))
meltpop23=melt(pop23)
p <- ggplot(meltpop23, aes(factor(variable), value)) 
p + geom_boxplot()+  facet_wrap(~variable, scale="free")
mig23=pop23[,1:4]
meltmig23=melt(mig23)
p <- ggplot(meltmig23, aes(factor(variable), value)) 
p + geom_boxplot() + ylim(0,1)+scale_y_continuous(trans='log2') +facet_wrap(~variable)
#more 3 going into 2

allmig=rbind(meltmig12,meltmig13,meltmig23)
p <- ggplot(allmig, aes(factor(variable), value)) 
p + geom_boxplot() + ylim(0,1)+scale_y_continuous(trans='log2')+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#larger proportion of genome shows reduced gene flow 

#only asymmetrical 

names=c("","","","TemperateCore-SubtropD","","","","TemperateEdge-SubtropD","","","","TemperateCore-SubtropL","","","","TemperateEdge-SubtropD","","","","SubtropL-SubtropD")


# Adding column based on other column:
allmig$Status=
ifelse(grepl("pop1_pop2",allmig$variable),"lineage1_lineage2", ifelse(grepl("pop1_pop3",allmig$variable),"lineage1_lineage3","lineage2_lineage3"))

allmig$Status_f=factor(allmig$Status, levels=c("lineage1_lineage2","lineage1_lineage3","lineage2_lineage3"))
allmig$GI=
  ifelse(grepl("i",allmig$variable),"GI", "BG") # had to switch
p <- ggplot(allmig, aes(factor(variable), value, fill=Status)) 
p + geom_boxplot() + ylim(0,1)+scale_y_continuous(trans='log2',labels = scales::scientific_format(
  digits = 3,
  scale = 1,
  prefix = "",
  suffix = "",
  decimal.mark = ".",
  trim = TRUE))+
  scale_x_discrete(labels=names)+
 theme_classic(base_size = 22)+
  scale_fill_manual(values=c("#DFBE99","#729EA1","#DB5375"))+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+ylab("migration(log)")+xlab(NULL)+facet_grid(~Status_f*GI, scales = "free_x")

P12=cbind(mean(pop12[,ncol(pop12)]), 1-mean(pop12[,ncol(pop12)]))
colnames(P12)=c("GI","BG")
P13=cbind(mean(pop13[,ncol(pop13)]), 1-mean(pop13[,ncol(pop13)]))
colnames(P13)=c("GI","BG")
P23=cbind(mean(pop23[,ncol(pop23)]), 1-mean(pop23[,ncol(pop23)]))
colnames(P23)=c("GI","BG")

PlotP12=ggplot(melt(P12), aes(x="", y=value, fill=Var2)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
# geom_text(aes(label = Var2), color = "white", size=6)+
  scale_fill_manual(values=c("#b3b3b3ff","white"))

PlotP13=ggplot(melt(P13), aes(x="", y=value, fill=Var2)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  # geom_text(aes(label = Var2), color = "white", size=6)+
  scale_fill_manual(values=c("#b3b3b3ff","white"))

PlotP23=ggplot(melt(P23), aes(x="", y=value, fill=Var2)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  # geom_text(aes(label = Var2), color = "white", size=6)+
  scale_fill_manual(values=c("#b3b3b3ff","white"))

grid.arrange(PlotP12, PlotP13, PlotP23, nrow = 1)


