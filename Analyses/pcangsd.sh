#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N setup  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
##$ -l h_rt=200:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be

#the following filters are what I used for bayescan minus the -snp_pval because no need to hard call
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 35 -dosnpstat 1 -doHWE 1 -sb_pval 1e-3 -maxHetFreq 0.5 -hetbias_pval 1e-3 -skipTriallelic 1 -minInd 53  -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -b ../bamscl -GL 1 $FILTERS $TODO -P 1 -out Sites4PCAngsd


PCANGSD=/usr4/bi594/jfifer/bin/pcangsd/pcangsd.py
module load python3
python $PCANGSD -beagle Sites4PCAngsd.beagle.gz -pcadapt -o PCANGSD_1

zcat Sites4PCAngsd.mafs.gz| tail -n +2 | cut -f 1,2 > Sites4PCAngsd.sites

#scp PCANGSD_1.pcadapt.zscores.npy and sites file and use R script pcangsd.R

#create wanted pairwise comps

#Rscript RemoveIndsBeagle.R fileA fileB
module load R

Rscript RemoveIndsBeagle.R Sites4PCAngsd.beagle pop1
Rscript RemoveIndsBeagle.R Sites4PCAngsd.beagle pop2
Rscript RemoveIndsBeagle.R Sites4PCAngsd.beagle pop3

gzip no_*beagle

#this is fast no need to run as a job
#need to add -minMaf 0 otherwise site numbers will change
PCANGSD=/usr4/bi594/jfifer/bin/pcangsd/pcangsd.py
module load python3
python $PCANGSD -beagle no_pop1.beagle.gz -pcadapt -minMaf 0 -o PCANGSD_pop2pop3
python $PCANGSD -beagle no_pop2.beagle.gz -pcadapt -minMaf 0 -o PCANGSD_pop1pop3
python $PCANGSD -beagle no_pop3.beagle.gz -pcadapt -minMaf 0 -o PCANGSD_pop1pop2
