#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N setup  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
##$ -l h_rt=200:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be

#doesnt make a difference if you do a reference or not
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 35 -dosnpstat 1 -doHWE 1 -sb_pval 1e-3 -maxHetFreq 0.5 -hetbias_pval 1e-3 -skipTriallelic 1 -minInd 53"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doBcf 1 -doPost 1 -doGlf 2"
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -b ../bamscl -GL 1 $FILTERS $TODO -P 1 -out Hetero_1



