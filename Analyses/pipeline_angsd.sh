#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ADMIX_FILT  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be



# angsd settings:
# -minMapQ 20 : only highly unique mappings (prob of erroneous mapping =< 1%)
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option)
# -maxDepth : highest total depth (sum over all samples) to assess; set to 10x number of samples
# -minInd : the minimal number of individuals the site must be genotyped in. Reset to 50% of total N at this stage.

#FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 780 -minInd 39"
# T O   D O :
#TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

# in the following line, -r argument is one chromosome or contig to work with (no need to do this for whole genome as long as the chosen chromosome or contig is long enough, ~1 Mb. When mapping to a real genome, consider chr1:1-1000000 )
# (look up lengths of your contigs in the header of *.sam files)
#/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out dd

# summarizing results (using modified script by Matteo Fumagalli)
#module load R
#Rscript /projectnb/davieslab/jfifer/Japan_rad/plotQC.R prefix=dd
# proportion of sites covered at >5x:
#cat quality.txt

#minInd 80% .8 x 78= 62
#FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 62 -snp_pval 1e-5 -minMaf 0.05"
#TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

# Starting angsd with -P the number of parallel processes. Funny but in many cases angsd runs faster on -P 1
#/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out myresult_withclones

# how many SNPs?
#NSITES=`zcat myresult.mafs.gz | wc -l`
#echo $NSITES
#9178


#without clones 67 samples
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 53 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

# Starting angsd with -P the number of parallel processes. Funny but in many cases angsd runs faster on -P 1
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -b bamscl -GL 1 $FILTERS $TODO -P 1 -out myresult


#module load python2
#python ~/bin/pcangsd/pcangsd.py -beagle myresult.beagle.gz -admix -o pcangsd -inbreed 2 -kinship -selection -threads 12
module load bcftools
for i in *.bcf; do bcftools view $i> ${i/%.bcf}.vcf; done
gunzip myresult.vcf.gz
module load admixture/1.3.0
module load plink/1.90b6.4

plink --vcf myresult.vcf --make-bed --allow-extra-chr 0 --out myresult --const-fid 0
for K in `seq 1 5`; \
do admixture --cv myresult.bed $K | tee myresult_${K}.out; done

# which K is least CV error?
grep -h CV myresult_*.out

#CV error (K=1): 0.59427
#CV error (K=2): 0.46569
#CV error (K=3): 0.38799
#CV error (K=4): 0.42006
#CV error (K=5): 0.46003


###Changing the filter
for minIndDep in 1 2 3 4 5; do echo -e "FILTERS=\"-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 62 -snp_pval 1e-5 -minMaf 0.05 -setMinDepthInd $minIndDep \" \n \
TODO=\"-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2\" \n \
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -b bamscl -GL 1 \$FILTERS \$TODO -P 1 -out myresult_dep$minIndDep \n \
module load bcftools \n \
for i in *.bcf; do bcftools view \$i> \${i/%.bcf}.vcf; done \n \
gunzip myresult_dep$minIndDep.vcf.gz \n \
module load admixture/1.3.0 \n \
module load plink/1.90b6.4 \n \
plink --vcf myresult_dep$minIndDep.vcf --make-bed --allow-extra-chr 0 --out myresult_dep$minIndDep --const-fid 0 \n \
for K in \`seq 1 5\`; \n \
do admixture --cv myresult_dep$minIndDep.bed \$K | tee myresult_dep$minIndDep._\${K}.out; done" >>filter_tests.sh; done

#grep -h CV myresult_*.out


### Get reads and mapping rates
#mapping rates
grep overall maps.e >mappingrates
cut -f1 -d ' ' mappingrates >mappingrates.tmp; paste mappingrates.tmp bams > mappingrates

#reads
for i in *.tr0; do awk '{s++}END{print s/4}' $i >>libsizes.tmp; paste libsizes.tmp bams >libsizes
