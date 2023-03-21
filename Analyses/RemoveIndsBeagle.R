args[1] beagle file
#args[2] list of column names that you want removed
#need to download caroline package, if you dont have it do
#wget https://cran.r-project.org/src/contrib/caroline_0.7.6.tar.gz
#R CMD INSTALL -l ~/R/ caroline_0.7.6.tar.gz  (or wherever your r packages are downloaded)


args = commandArgs(trailingOnly=TRUE)


beagle=read.table(file=args[1],header=T,check.names = FALSE)
#beagle=read.table(file="Sites4PCAngsd.beagle",header=T,check.names = FALSE)
library(dplyr)
fileB=args[2]
#fileB="pop1"
library(readr)
#make sure you add the lib path that contains the caroline package
.libPaths(c("~/R/", .libPaths()))
library(caroline)
nameremoves=(scan(fileB, character(), quote = ''))

test=paste("^",nameremoves,"$", collapse="|", sep = "")
test

#df_new <- beagle %>% select(-matches(test))
output=beagle[, -grep(test, colnames(beagle))]
output1=sub("[[:punct:]].*", "", colnames(output))
colnames(output)=output1
write.delim(output, file=paste0("no_",args[2], ".beagle"), quote = FALSE, row.names = FALSE, sep = "\t")
#write.delim(output, file=paste0("no_","pop1", ".beagle"), quote = FALSE, row.names = FALSE, sep = "\t")

