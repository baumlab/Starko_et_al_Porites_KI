
setwd("C:/Users/james/My Drive/Documents/BOSTON/Davies/Baum/Code/Admixture")

# assembling the input table
dir=("C:/Users/james/My Drive/Documents/BOSTON/Davies/Baum/Code/Admixture/") # path to input files
#ngsAdmix
#inName="mydata_k4_r6.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
#inName="myresult2.6.Q"
inName="myresult.3.Q"


samples=data.frame(read.csv("../Adapter2Sample_MetaData_Combined_v2.csv", header=T))
colnames(samples)[1]="Adapter"

row.names(samples)=samples$Adapter

bams=data.frame(read.csv("../bamscl", header = FALSE)) # list of bam files
colnames(bams)<- "Adapter"


#------------

#npops=as.numeric(sub(".+(\\d+)\\..+","\\1",inName))
npops=3
tbl=read.table(paste(dir,inName,sep=""),header=F)
row.names(tbl)= bams$Adapter

#new<-merge(tbl,i2p,by="row.names")
new1<-transform(merge(tbl,samples,by=0), row.names=Row.names, Row.names=NULL)
#write.csv(new1, file="bamscl_adaptors.csv")
#tbl=cbind(new1[c(1:3,13,8)])
tbl=cbind(new1[c(1:3,14,8)])
head(tbl)
colnames(tbl)=c("V1","V2","V3","pop","ind")
head(tbl)

tbl$pop=as.factor(tbl$pop)
# putting populaitons in desired order (edit pop names as needed or skip to plot them alphabetically)
# tbl$pop=factor(tbl$pop,levels=c("O","K"))
#red=1, blue=3, tan=2
colors=c("#DFBE99","#729EA1","#DB5375")
#new colors: blush lineage 3, tan lineage 1, cadet lineage 2
tbl$pop=factor(tbl$pop, levels=c("VL1","VL2","VL3","L1","L2","M1","M2","M3","M5","VH1","VH2","VH3"))
#Download below from https://github.com/z0on/2bRAD_denovo/blob/master/2bRAD_README.txt
source("C:/Users/james/My Drive/Documents/BOSTON/Davies/Range expansion/Code/plot_admixture_v6_function.R")
ords=plotAdmixture(data=tbl,npops=npops,grouping.method="distance", hshift=-1.5, vshift=-.07,colors = colors)



# recording cluster affiliations
cluster.admix=apply(tbl[,1:npops],1,function(x) {return(paste(which(x>0.85),collapse=".")) })
cluster.assign=as.data.frame(cluster.admix)
write.csv(cluster.assign, file="cluster.assign.csv")
#write.csv(cluster.assign, file="cluster.assign.HALF.ngsADMIX.csv")
save(cluster.admix,file=paste(inName,"_clusters.RData",sep=""))
