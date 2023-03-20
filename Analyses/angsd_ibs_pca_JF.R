getwd()

setwd("C:/Users/james/My Drive/Documents/BOSTON/Davies/Baum/Code/PCA")
# bams=read.table("bams.csv")[,1] # list of bam files
# goods=c(1:length(bams))
# 
# # reading table of pairs of replicates (tab-delimited) - skip if there are no clones
# clonepairs=read.table("clonepairs.tab",sep="\t")
# repsa= clonepairs[,1]
# repsb= clonepairs[,2]
# # removing "b" replicates
# goods=which(!(bams %in% repsb))

#--------------------
# loading individual to population correspondences
bams=data.frame(read.csv("C:/Users/james/My Drive/Documents/BOSTON/Davies/Baum/Code/bamscl", header = FALSE)) # list of bam files
colnames(bams)<- "Adapter"
samples=data.frame(read.csv("C:/Users/james/My Drive/Documents/BOSTON/Davies/Baum/Code/Adapter2Sample_MetaData_Combined.csv", header=T))
colnames(samples)[1]<-c("Adapter")
# samples <- data.frame(lapply(samples, function(x) {
#   gsub("tr0", "nosymbio.fastq.bt2.bam", x)
# }))

new<-merge(bams,samples, by ="Adapter")
new<-new[ order(match(new$Adapter, bams$Adapter)), ]

#i2p=read.table("inds2pops",sep="\t") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
#row.names(i2p)=i2p[,1]
#i2p=i2p[goods,]
#Nnew=new[!(new$Site=="Oura Bay" | new$Site== "Sekisei"),]

#Nsites=Nnew$Site
site=new$Clade



# settign up colors for plotting
# palette(rainbow(length(unique(Nsites))))
# colors=as.numeric(as.factor(Nsites))
# colpops=as.numeric(as.factor(sort(unique(Nsites))))



# palette(rainbow(length(unique(site))))
# colors=as.numeric(as.factor(site))
# colpops=as.numeric(as.factor(sort(unique(site))))
#--------------------
#run above code 'loading invidividual.." first
# covariance / PCA 

library(vegan)
# choose either of the following two covarince matrices:
co = as.matrix(read.table("myresult.covMat")) # covariance based on single-read resampling
#co = as.matrix(read.table("ok.covar")) # covariance by ngsCovar
dimnames(co)=list(new$RAD_sampleID)
# PCoA and CAP (constranied analysis of proximities)  
conds=data.frame(site)
pp0=capscale(as.dist(1-cov2cor(co))~1) # PCoA
pp=capscale(as.dist(1-cov2cor(co))~site,conds) # CAP

# significance of by-site divergence, based on 1-correlation as distance
adonis(as.dist(1-cov2cor(co))~site,conds)

# eigenvectors: how many are interesting?
plot(pp0$CA$eig) 

# k-means clustering (instead of admixture-based)
kmeans.axes <- sum(diff(summary(eigenvals(pp0))[2,])*-1 > 0.01) #based on eigen plot
# This function finds the total number of eigenvectors in which the difference between last eigenvalue and the one preceding is greater than 0.01
# Standardized method of selecting eigenvectors that contribute the most explanatory power in clustering
# 2 axes
#ks <- kmeans(pp0$CA$u[,1:kmeans.axes], 3, nstart = 100, iter.max = 1000)



axes2plot=c(1,2)  
#quartz()
cc=pp0
plot(cc,choices=axes2plot,type="n") # choices - axes to display
points(cc,choices=axes2plot,pch=19,col=colors)
#ordihull(cmd,choices= axes2plot,groups= conds$grp,draw="polygon",col=1+as.numeric(unique(as.factor(conds$grp))),label=T)
ordispider(cc,choices= axes2plot,groups=conds$site,col="grey80")
ordiellipse(cc,choices= axes2plot,groups= conds$site,draw="polygon",col=colpops)

# unscaled, to identify outliers
n2identify=2
cmd=pp0
cmd$CA$u
# unscaled, to identify outliers
plot(cmd$CA$u[,axes2plot],pch=19,col=colors)
ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=colpops,label=T)
identify(cmd$CA$u[,axes2plot],labels=rownames(co),n=n2identify,cex=0.7)
#
library(dplyr)
library(ggplot2)

pca_s <- as.data.frame(cmd$CA$u[,axes2plot])
#colors=c("#006400FF", "#caff70ff","#6af2ffff","#00bfffff","#5087c1ff","#325fa2ff","#7fabd3ff")
colors=c("#DFBE99","#729EA1","#DB5375")
#square, circle, triangle 0,1,2
eigenvals=pp0$CA$eig
#calc variance explained
varexplained=eigenvals/sum(eigenvals)

ggplot(pca_s, aes(MDS1, MDS2)) +
  geom_point(size=3, aes(color=as.factor(new$Clade), shape=as.factor(new$Clade)), show.legend = NULL)  +
  theme_classic(base_size = 20) +
 # geom_polygon(data = hull_cyl, alpha = 0.5, color='black')+
  stat_ellipse(geom="polygon",level=.99999, alpha = 1/3, aes(fill = as.factor(new$Clade)),show.legend = NULL)+
  scale_fill_manual(values=c("#DFBE99","#729EA1","#DB5375"))+
  xlab(paste0("MDS1 (",formatC((varexplained[1]*100), digits = 2, format = "f"),"%)")) +
  ylab(paste0("MDS2 (",formatC((varexplained[2]*100), digits = 2, format = "f"),"%)")) +
  scale_shape_manual(values=c(15, 16, 17))+
  scale_colour_manual(values=colors, name="", labels = c("lineage1", "lineage2", "lineage3", guide=FALSE)) #
  #scale_fill_manual(values=colors,name = "site" )
  #scale_shape_discrete(name = "Condition", labels = c("High Flow + Heat", "High Flow + No Heat", "Low Flow + Heat","Low Flow + No Heat"))+
  #ggtitle("PCoA Host",expression(paste(italic("Ex Situ"))))

axes2plot=c(2,3)  
pca_s <- as.data.frame(cmd$CA$u[,axes2plot])


#plotting 2 and 3
ggplot(pca_s, aes(MDS2, MDS3)) +
  #geom_point(size=3, aes(color=as.factor(new$Clade), shape=as.factor(new$mortality_cat)), show.legend = NULL)  +
  geom_point(size=3, aes(color=as.factor(new$Clade), shape=as.factor(new$Region))) +
   theme_classic(base_size = 20) +
 # stat_ellipse(geom="polygon",level=.99999, alpha = 1/3, aes(fill = as.factor(new$Clade)),show.legend = NULL)+
  scale_fill_manual(values=c("#DFBE99","#729EA1","#DB5375"))+
  xlab(paste0("MDS3 (",formatC((varexplained[3]*100), digits = 2, format = "f"),"%)")) +
  ylab(paste0("MDS4 (",formatC((varexplained[4]*100), digits = 2, format = "f"),"%)")) +
  scale_shape_manual(values=1:nlevels(as.factor(new$site_name))) +
  scale_colour_manual(values=colors, name="", labels = c("lineage1", "lineage2", "lineage3", guide=FALSE)) #


### Add site 

ggplot(pca_s, aes(MDS1, MDS2)) +
  geom_point(size=3, aes(color=site))  +
  theme_classic(base_size = 20) +
  # geom_polygon(data = hull_cyl, alpha = 0.5, color='black')+
  #stat_ellipse(geom="polygon",level=0.65, alpha = 1/3, aes(fill = heat),show.legend = NA)+
  stat_ellipse(geom="polygon",level=.75, alpha = 1/4, aes(fill = site),show.legend = NA)+
  xlab(paste0("MDS1")) +
  ylab(paste0("MDS2")) +
  coord_cartesian(ylim = c(-.15, .25), xlim=c(-.1,.2))+
  
  scale_fill_manual(values=colors,guide=FALSE)+
  scale_colour_manual(values=colors, name="", labels = c("Amakusa", "Tatsukushi", "Kushima","Kushimoto","Oura Bay","Sekisei Lagoon", "Shirahama")) #+
#scale_fill_manual(values=colors,name = "site" )
#scale_shape_discrete(name = "Condition", labels = c("High Flow + Heat", "High Flow + No Heat", "Low Flow + Heat","Low Flow + No Heat"))+
#ggtitle("PCoA Host",expression(paste(italic("Ex Situ"))))

axes2plot=c(2,3)  
pca_s <- as.data.frame(cmd$CA$u[,axes2plot])


#plotting 2 and 3
ggplot(pca_s, aes(MDS2, MDS3)) +
  geom_point(size=3, aes(color=site))  +
  theme_classic(base_size = 20) +
  # geom_polygon(data = hull_cyl, alpha = 0.5, color='black')+
  #stat_ellipse(geom="polygon",level=0.65, alpha = 1/3, aes(fill = heat),show.legend = NA)+
  stat_ellipse(geom="polygon",level=.75, alpha = 1/4, aes(fill = site),show.legend = NA)+
  xlab(paste0("MDS2")) +
  ylab(paste0("MDS3")) +
  #coord_cartesian(ylim = c(-.15, .15), xlim=c(-.15,.25))+
  
  scale_fill_manual(values=colors,guide=FALSE)+
  scale_colour_manual(values=colors, name="", labels = c("Amakusa", "Tatsukushi", "Kushima","Kushimoto","Oura Bay","Sekisei Lagoon", "Shirahama")) #+
#scale_fill_manual(values=colors,name = "site" )
#scale_shape_discrete(name = "Condition", labels = c("High Flow + Heat", "High Flow + No Heat", "Low Flow + Heat","Low Flow + No Heat"))+
#ggtitle("PCoA Host",expression(paste(italic("Ex Situ"))))





#-------------
# t-SNE:  machine learning to identify groups of samples 
# based on genotypes' correlations

library(Rtsne)
library(vegan)
library(adegenet)
quartz()

# perplexity:  expected number fo neighbors. Set to 0.5x N(samples per pop)
perp=15
rt = Rtsne(as.dist(1-cov2cor(co)), perplexity=perp,max_iter=2,is_distance=T)
for (i in 1:250){
	rt = Rtsne(as.dist(1-cov2cor(co)), perplexity=perp,max_iter=10,Y_init=rt$Y,is_distance=T)
	plot(rt$Y,col=colors,pch=16,cex=0.8,main=i*10)
}
ordispider(rt$Y,groups=site,col="grey80",alpha=0.01)
ordiellipse(rt$Y,groups= site,draw="polygon",col=colpops,label=T)

#-------------
# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)

#First look for clones
ma = as.matrix(read.table("../myresult.withclones.ibsMat"))
bams=data.frame(read.csv("../bams", header = FALSE)) # list of bam files
colnames(bams)<- "Adapter"
samples=data.frame(read.csv("../Adapter2Sample.csv"))
new<-merge(bams,samples, by ="Adapter")
new<-new[ order(match(new$Adapter, bams$Adapter)), ]

#bamz=read.table("bams.csv")
dimnames(ma)=list(new$Sample)
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.7)

#check w/o clones
ma = as.matrix(read.table("myresult.ibsMat"))
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.5)  # this shows how similar clones are

#without clones reran to generate myresult.ibsMat with new bamscl file instead
ma = as.matrix(read.table("myresult.ibsMat"))
bams=data.frame(read.csv("bamscl.csv", header = FALSE)) # list of bam files
colnames(bams)<- "Adapter"
samples=data.frame(read.csv("Adapter2Sample.csv"))
new<-merge(bams,samples, by ="Adapter")
new<-new[ order(match(new$Adapter, bams$Adapter)), ]

#bamz=read.table("bams.csv")
dimnames(ma)=list(new$Sample)
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.7)






#### N sites only
#--------------------
# covariance / PCA 

library(vegan)
# choose either of the following two covarince matrices:
co = as.matrix(read.table("myresult.covMat")) # covariance based on single-read resampling
#co = as.matrix(read.table("ok.covar")) # covariance by ngsCovar

dimnames(co)=list(new$Sample)

# PCoA and CAP (constranied analysis of proximities)  
conds=data.frame(site)
pp0=capscale(as.dist(1-cov2cor(co))~1) # PCoA
pp=capscale(as.dist(1-cov2cor(co))~site,conds) # CAP

# significance of by-site divergence, based on 1-correlation as distance
adonis(as.dist(1-cov2cor(co))~site,conds)

# eigenvectors: how many are interesting?
plot(pp0$CA$eig) 

axes2plot=c(1,2)  
quartz()
cc=pp0
plot(cc,choices=axes2plot,type="n") # choices - axes to display
points(cc,choices=axes2plot,pch=19,col=colors)
#ordihull(cmd,choices= axes2plot,groups= conds$grp,draw="polygon",col=1+as.numeric(unique(as.factor(conds$grp))),label=T)
ordispider(cc,choices= axes2plot,groups=conds$site,col="grey80")
ordiellipse(cc,choices= axes2plot,groups= conds$site,draw="polygon",col=colpops)

########################################


