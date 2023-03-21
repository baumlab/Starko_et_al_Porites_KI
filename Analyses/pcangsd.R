setwd("C:/Users/james/My Drive/Documents/BOSTON/Davies/Baum/Code/PCAngsd")

library(RcppCNPy)
#install.packages("bigutilsr")
library(bigutilsr)
filename="PCANGSD_1.pcadapt.zscores.npy"

zscores <- npyLoad(filename)
K <- ncol(zscores)

# COuld mess around with below 
#if you wanted to look at only one pc you could do
#d2 <- (zscores - median(zscores))^2 for the column of interest
if (K == 1) {
  d2 <- (zscores - median(zscores))^2
} else {
  d2 <- dist_ogk(zscores)
}

outputname=sub("\\.pcadapt.*", "", filename)
write.table(d2, file=paste0(outputname, ".pcadapt.test.txt"), quote=F, row.names=F, col.names=F)
write.table(pchisq(d2, df=K, lower.tail=F), file=paste0(outputname, ".pcadapt.pval.txt"), quote=F, row.names=F, col.names=F)


# load pvals
pval=read.table(file ="PCANGSD_1.pcadapt.pval.txt")

## read positions (hg38)
p<-read.table("Sites4PCAngsd.sites",colC=c("factor","integer"))

names(p)<-c("chr","pos")

## make Manhattan plot
plot(-log10(pval$V1),col=p$chr,xlab="Chromosomes",main="Manhattan plot")


df=cbind(p,pval$V1)
colnames(df)=c("chr","pos","pval")
df$pos=as.numeric(df$pos)

## Make a prettier Manhattan plot
library(dplyr)
don <- df %>% 
# Compute chromosome size
group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate(BPcum=pos+tot)
axisdf = don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

library(ggplot2)
ggplot(don, aes(x=BPcum, y=-log10(pval))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  #scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
#################################################
#Now splitting the analysis into pairwise comps

library(RcppCNPy)
#install.packages("bigutilsr")
library(bigutilsr)

test=function(filename){
zscores <- npyLoad(filename)
K <- ncol(zscores)
# COuld mess around with below 
#if you wanted to look at only one pc you could do
#d2 <- (zscores - median(zscores))^2 for the column of interest
if (K == 1) {
  d2 <- (zscores - median(zscores))^2
} else {
  d2 <- dist_ogk(zscores)
}
outputname=sub("\\.pcadapt.*", "", filename)
write.table(d2, file=paste0(outputname, ".pcadapt.test.txt"), quote=F, row.names=F, col.names=F)
write.table(pchisq(d2, df=K, lower.tail=F), file=paste0(outputname, ".pcadapt.pval.txt"), quote=F, row.names=F, col.names=F)}

filename="PCANGSD_pop1pop2.pcadapt.zscores.npy"
test(filename)
filename="PCANGSD_pop1pop3.pcadapt.zscores.npy"
test(filename)
filename="PCANGSD_pop2pop3.pcadapt.zscores.npy"
test(filename)

## read positions
p<-read.table("Sites4PCAngsd.sites",colC=c("factor","integer"))
names(p)<-c("chr","pos")

# load pvals 12
pval=read.table(file ="PCANGSD_pop1pop2.pcadapt.pval.txt")
df=cbind(p,pval$V1)
colnames(df)=c("chr","pos","pval")
df$pos=as.numeric(df$pos)
pop12_output=df
#load pvals 13
pval=read.table(file ="PCANGSD_pop1pop3.pcadapt.pval.txt")
df=cbind(p,pval$V1)
colnames(df)=c("chr","pos","pval")
df$pos=as.numeric(df$pos)
pop13_output=df
#load pvals 23
pval=read.table(file ="PCANGSD_pop2pop3.pcadapt.pval.txt")
df=cbind(p,pval$V1)
colnames(df)=c("chr","pos","pval")
df$pos=as.numeric(df$pos)
pop23_output=df

#Convert to q values first to account for multiple loci
#BiocManager::install("qvalue")
library(qvalue)
pop12_output$qval <- qvalue(pop12_output$pval)$qvalues
alpha <- 0.05
pop12_outliers=subset(pop12_output, qval < alpha)
nrow(pop12_outliers)

pop13_output$qval <- qvalue(pop13_output$pval)$qvalues
alpha <- 0.05
pop13_outliers=subset(pop13_output, qval < alpha)
nrow(pop13_outliers)

pop23_output$qval <- qvalue(pop23_output$pval)$qvalues
alpha <- 0.05
pop23_outliers=subset(pop23_output, qval < alpha)
nrow(pop23_outliers)


#Multiple test correction (because we've split it up into pairwise comps)
#12
library(stats)
 fdr <- matrix(ncol=1, nrow=length(pop12_outliers$qval))
 for(i in 1:length(pop12_outliers$qval)){
   fdr[i,]=p.adjust(pop12_outliers$qval[i], method = "BH", n=3)
 }
 pop12_outliers$qval=fdr
 pop12_outliers=subset(pop12_outliers, qval < 0.05)

#13
 fdr <- matrix(ncol=1, nrow=length(pop13_outliers$qval))
 for(i in 1:length(pop13_outliers$qval)){
   fdr[i,]=p.adjust(pop13_outliers$qval[i], method = "BH", n=3)
 }
 pop13_outliers$qval=fdr
 pop13_outliers=subset(pop13_outliers, qval < 0.05)
 #23
 fdr <- matrix(ncol=1, nrow=length(pop23_outliers$qval))
 for(i in 1:length(pop23_outliers$qval)){
   fdr[i,]=p.adjust(pop23_outliers$qval[i], method = "BH", n=3)
 }
 pop23_outliers$qval=fdr
 pop23_outliers=subset(pop23_outliers, qval < 0.05)
 


###Annotation###

##Annotating Bayescan (#function is below )
#Just need to do below once then read in genes file below
#  genes=read.table(file="C:/Users/james/My Drive/Documents/BOSTON/Davies/Baum/Code/PCAngsd/genes2.txt")
#  names(genes)=c("contig","start","end","gene")
#  test=sub(".*ID= *(.*?) *;Name.*", "\\1", genes$gene)
#  genes$gene=test
# # # adding gene annotations, in the same order as genes.txt (skip if you don't have such file):
#  gnames=read.table(file="C:/Users/james/My Drive/Documents/BOSTON/Davies/Baum/Code/PCAngsd/Plutea_iso2gene.tab", header=F, sep='\t')
#  library(dplyr)
#   names(gnames)=c("gene","protein")
#  genes=merge(genes,gnames,by="gene",all.x=T)
#  genes$protein=as.character(genes$protein)
#  genes$protein[is.na(genes$protein)]="unknown"
#  nrow(genes[genes$protein!="unknown",])
#  write.csv(genes, file="P.lutea.genes.annotated")
# #add +/- 1k
#  library(tidyverse)
#  genes=genes %>%
#    mutate(start = start - 1000)%>%
#    mutate(end = end + 1000)
#  write.csv(genes, file="P.lutea.genes2kb.annotated")
 
genes=read.csv(file="P.lutea.genes2kb.annotated")


#install.packages("tidyr")
library(tidyr)
#install.packages("BiocManager")
#BiocManager::install("GenomicRanges")
library(GenomicRanges)

test=function(a){
  filename1 <- deparse(substitute(a));
  genes=read.csv(file="P.lutea.genes2kb.annotated",  row.names = 1);
  #a=  a[ -c(2) ]; only add these two lines if there aren't already just two columns
  #a=separate(data = a, col = V1, into = c("contig", "pos"), sep = "_");
  colnames(a)=c("contig", "pos","qval")
  #or for fishers (no qval)
  #colnames(a)=c("contig", "pos")
  query <- GRanges(seqnames=genes$contig, ranges=IRanges(genes$start, genes$end));
  ref = GRanges(seqnames=a$contig, ranges=IRanges(as.numeric(a$pos), width=1));
  toprint=table(!is.na(findOverlaps(query, ref)));
  print(toprint);
  olaps=(findOverlaps(query, ref));
  Results=cbind(genes[queryHits(olaps),], a[subjectHits(olaps),]);
  write.csv(Results, file = paste0(filename1, ".overlap.csv"),  row.names = FALSE);
}


test(pop12_outliers)
test(pop13_outliers)
test(pop23_outliers)

genes12=read.csv(file="pop12_outliers.overlap.csv")
genes13=read.csv(file="pop13_outliers.overlap.csv")
genes23=read.csv(file="pop23_outliers.overlap.csv")

library(dplyr)                                        
uniq.12=genes12 %>% distinct(gene, .keep_all = TRUE)
uniq.13=genes13 %>% distinct(gene, .keep_all = TRUE)
uniq.23=genes23 %>% distinct(gene, .keep_all = TRUE)

uniq.int.12.13=data.frame(intersect(uniq.12$gene,uniq.13$gene))
#the interesting intersection
uniq.int.23.13=data.frame(intersect(uniq.23$gene,uniq.13$gene))

#reformat for output

colnames(uniq.12)=c("gene.ref","chromosome.ref","start.ref","end.ref","protein","chromosome","pos","pval","qval")
write.csv(uniq.12,file="uniq.12.final.csv")
colnames(uniq.13)=c("gene.ref","chromosome.ref","start.ref","end.ref","protein","chromosome","pos","pval","qval")
write.csv(uniq.13,file="uniq.13.final.csv")
colnames(uniq.23)=c("gene.ref","chromosome.ref","start.ref","end.ref","protein","chromosome","pos","pval","qval")
write.csv(uniq.23,file="uniq.23.final.csv")



