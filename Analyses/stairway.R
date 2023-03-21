setwd("C:/Users/james/My Drive/Documents/BOSTON/Davies/Baum/Code/Demo")
#popN=read.table(file="popN.final.summary",header=T)
#popN$population="popN"
pop2=read.table(file="pop2two-epoch.final.summary",header=T)
pop2$population="pop2"
pop1=read.table(file="pop1two-epoch.final.summary",header=T)
pop1$population="pop1"

pop3=read.table(file="pop3two-epoch.final.summary", header=T)
pop3$population="pop3"




all=rbind(pop1,pop3,pop2)
library(ggplot2)
ggplot(data=all, aes(x=year, y=Ne_median, color=population)) +geom_line(size=1.5)+
geom_ribbon(aes(ymin=all$Ne_12.5.,ymax=all$Ne_87.5., fill = factor(population)), alpha = 0.3,colour = NA,show_guide = FALSE)+
 theme_bw(base_size = 22)+ scale_y_continuous(name="Ne", labels = scales::comma)+#scale_x_continuous(name="Years", labels=c("0","100k","200k",
  scale_x_continuous(name="Years",limits =  c(0,1.25e+06), labels=c("0","500k","1000k","1500k"))+                                                                                                                         # "300k"),limits =  c(0,3e+05))+ 
  scale_fill_manual(values=c("#dfbe99ff","#729EA1","#DB5375"))+
  scale_color_manual(values=c("#dfbe99ff","#729EA1","#DB5375"), labels=c("PKir-1","PKir-2","PKir-3"))+theme(legend.title=element_blank())


pop2=read.table(file="pop2_35two-epoch.final.summary",header=T)
pop2$population="pop2"
pop1=read.table(file="pop1_35two-epoch.final.summary",header=T)
pop1$population="pop1"

pop3=read.table(file="pop3_35two-epoch.final.summary", header=T)
pop3$population="pop3"  






