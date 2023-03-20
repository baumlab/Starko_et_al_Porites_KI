setwd("C:/Users/james/My Drive/Documents/BOSTON/Davies/Baum/Code/Hetero")

#glf <- read.table(file = 'third.het.beagle.gz', header=TRUE)[,-c(1:3)]
#glf <- read.table(file = 'lineage_ref.beagle.gz', header=TRUE)[,-c(1:3)]
glf <- read.table(file = 'Hetero_1.beagle.gz', header=TRUE)[,-c(1:3)]

glf <- array(c(t(as.matrix(glf))), c(3, ncol(glf)/3, nrow(glf)))
EMstep <- function (sfs, GL) rowMeans(prop.table(sfs * GL, 2))

SFS <- matrix(1/3,3,dim(glf)[2])

maxiter <- 200
tol <- 1e-8

for(sample in 1:dim(glf)[2])
{
  for (iter in 1:maxiter)
  {
    upd <- EMstep (SFS[,sample], glf[,sample,])
    if (sqrt(sum((upd - SFS[,sample])^2)) < tol)
      break;
    SFS[,sample] <- upd
  }
  if (iter == maxiter) warning("increase maximum number of iterations")
}


poplist.names=c("pop2","pop1","pop3","pop3","pop1","pop3","pop1","pop2","pop2","pop2","pop1","pop3","pop3","pop2","pop2","pop3","pop1","pop1","pop2","pop2","pop2","pop1","pop3","pop2","pop2","pop3","pop1","pop3","pop2","pop2","pop2","pop3","pop1","pop1","pop1","pop3","pop1","pop2","pop1","pop3","pop1","pop3","pop1","pop1","pop1","pop2","pop1","pop1","pop1","pop3","pop1","pop3","pop1","pop2","pop2","pop3","pop3","pop1","pop2","pop1","pop1","pop1","pop2","pop1","pop1","pop2","pop1")
pop1pos=which(poplist.names %in% "pop1")
pop2pos=which(poplist.names %in% "pop2")
pop3pos=which(poplist.names %in% "pop3")


# 3pop lineages
print(c('pop1',round(summary(SFS[2,pop1pos]),4)),quote=F)
print(c('pop2',round(summary(SFS[2,pop2pos]),4)),quote=F)
print(c('pop3',round(summary(SFS[2,pop3pos]),4)),quote=F)

#
sfs=SFS[2,pop1pos]
Site=replicate(29, "pop1")
pop1=data.frame(sfs, Site)

sfs=SFS[2,pop2pos]
Site=replicate(21, "pop2")
pop2=data.frame(sfs, Site)

sfs=SFS[2,pop3pos]
Site=replicate(17, "pop3")
pop3=data.frame(sfs, Site)



All=rbind(pop1,pop2,pop3)

library(ggplot2)
# Basic violin plot
colors=c("#ff0088","#caff70", "#690582ff")
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
p <- ggplot(All, aes(x=Site, y=sfs,fill=Site)) + 
  geom_violin()+
  scale_fill_manual(values=colors)+
  theme_bw(base_size = 22)+
  ylab(label = "Expected Het")+xlab(label="")+
  theme(axis.text.x = element_text(angle = 45, vjust = , hjust=1))
p+stat_summary(fun.y=mean,size=2, 
               geom="point", color="black")

t.test(pop2$sfs,pop3$sfs)
#t = 2.0786, df = 32.383, p-value = 0.04566
t.test(pop1$sfs,pop2$sfs)

