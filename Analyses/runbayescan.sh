#run demoangsd.sh

#create vcfs for each comparison


#create .spids for each vcf

# create a file called vcf2bayescan.spid containing this text:
echo "############
# VCF Parser questions
PARSER_FORMAT=VCF
# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION=true
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=DIPLOID
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION=true
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=false
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION=./bspops
# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION=false
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=
# GESTE / BayeScan Writer questions
WRITER_FORMAT=GESTE_BAYE_SCAN
# Specify which data type should be included in the GESTE / BayeScan file  (GESTE / BayeScan can only analyze one data type per file):
GESTE_BAYE_SCAN_WRITER_DATA_TYPE_QUESTION=SNP
############" >vcf2bayescan_n2.spid


echo -e '#!/bin/bash \n#$ -V # inherit the submission environment \n#$ -cwd # start job in submission directory
#$ -N bayescan  # job name, anything you want
#$ -P davieslab
#$ -pe omp 28
#$ -l h_rt=12:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be \n
module load java/1.8.0_181
java -jar /projectnb/davieslab/jfifer/Japan_rad/Fst_Outlier/PGDSpider_2.0.7.1/PGDSpider2-cli.jar -inputfile AllSites.vcf -outputfile Best.bayescan -spid vcf2bayescan.spid
module load bayescan
BayeScan2.1_linux64bits Best.bayescan -threads=20' >bayescan.sh


#run bayescan.sh
#For some reason below code only works by adding one blank line to the top of the _fst.txt file
/projectnb/davieslab/jfifer/Japan_rad/Bayescan/removeBayescanOutliers.pl bayescan=Best.bayesca_fst.txt vcf=./AllSites.vcf FDR=0.5 mode=delete >nooutSites.vcf

while read p; do
if [[ $p =~ "#" ]]; then
echo "$p"
fi; done <AllSites.vcf

awk 'NR==FNR{a[$1]=1;next}{if(!a[FNR])print$0}' deletelistfile largetextfile > resultfile

