#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N demosetup  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
##$ -l h_rt=200:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be
##$ -pe omp 24


#FILTERS="-minMapQ 30 -minQ 35 -minInd 53 -doHWE 1 -sb_pval 1e-3 -hetbias_pval 1e-3 -skipTriallelic 1 -maxHetFreq 0.5"
#TODO="-doMajorMinor 4 -ref $GENOME_FASTA -doMaf 1 -dosnpstat 1 -doGeno 11 -doPost 2 -doBcf 1"
#/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -b ../bamscl -GL 1 $FILTERS $TODO -P 1 -out AllSitesDemo

#zcat AllSitesDemo.mafs.gz | tail -n +2 | cut -f 1,2 > AllSitesDemo.sites

#cut -d" " -f2- ../Bayescan/Best.baye_fst.txt | tail -n +2> aaa
#find line number in vcf file that contains header
#grep -n "#CHROM" ../Bayescan/AllSites.vcf
#3018 #add 1
#tail --lines=+3019 ../Bayescan/AllSites.vcf | cut -f 1,2 | paste --delimiters "\t" - aaa > baye_fst_pos.txt

#awk '$5<0.5 {print $1"\t"$2}' ./baye_fst_pos.txt > bayeOuts
#grep -vf bayeOuts AllSitesDemo.sites > sites2do_filtBaye
#/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd sites index sites2do_filtBaye

export GENOME_FASTA=/projectnb/davieslab/jfifer/Baum/Genome_ref/Host/plut_final_2.1.fasta
TODO="-doSaf 1 -anc $GENOME_FASTA -ref $GENOME_FASTA"
# In the following lines, set minInd to 75-90% of each pop's sample size
#Trying with 80%
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -sites sites2do_filtBaye -b ../bamscl_pop1 -GL 1 -P 1 -minInd 23 $TODO -out pop1.out # 29
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -sites sites2do_filtBaye -b ../bamscl_pop2 -GL 1 -P 1 -minInd 16 $TODO -out pop2.out # 21
/projectnb/davieslab/jfifer/Japan_rad/angsd/angsd -sites sites2do_filtBaye -b ../bamscl_pop3 -GL 1 -P 1 -minInd 13 $TODO -out pop3.out #17


/projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS pop2.out.saf.idx >pop2.sfs
/projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS pop1.out.saf.idx >pop1.sfs
/projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS pop3.out.saf.idx >pop3.sfs


export GENOME_FASTA=/projectnb/davieslab/jfifer/Baum/Genome_ref/Host/plut_final_2.1.fasta

f=(*.saf.idx)
for B in `seq 1 100`; do :; for ((i = 0; i < ${#f[@]}; i++)); do       for ((j = i + 1; j < ${#f[@]}; j++)); do           echo "/projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS ${f[i]} ${f[j]} \
-ref "\$GENOME_FASTA" -anc "\$GENOME_FASTA"  -bootstrap 6 -P 1 -resample_chr 1 -seed  "\$RAN_SEED" > boot.${f[i]/%.out.saf.idx/}${f[j]/%.out.saf.idx/}/${f[i]/%.out.saf.idx/}${f[j]/%.out.saf.idx/}_$B " >>b100.txt;       done;   done; done
#3pop


###Post boots
#Pop1 29 (*2+1) 59  (*2*.8) 46
#Pop2 21 43 33
#Pop3 17 35 27
#two pop
SFSIZE="59 35" # 2N+1 for each population. In this case we assume that we have sampled 10 diploid individuals from each `p1` and `p2`.
for B in `seq 1 100`; do
echo $SFSIZE > pop1pop3_${B}.sfs;
tail -5 pop1pop3_${B} | awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.3f", a[i]/NR; printf "\t"};printf "\n"}' >> pop1pop3_${B}.sfs;
done

SFSIZE="59 43" # 2N+1 for each population. In this case we assume that we have sampled 10 diploid individuals from each `p1` and `p2`.
for B in `seq 1 100`; do
echo $SFSIZE > pop1pop2_${B}.sfs;
tail -5 pop1pop2_${B} | awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.3f", a[i]/NR; printf "\t"};printf "\n"}' >> pop1pop2_${B}.sfs;
done

SFSIZE="43 35" # 2N+1 for each population. In this case we assume that we have sampled 10 diploid individuals from each `p1` and `p2`.
for B in `seq 1 100`; do
echo $SFSIZE > pop2pop3_${B}.sfs;
tail -5 pop2pop3_${B} | awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.3f", a[i]/NR; printf "\t"};printf "\n"}' >> pop2pop3_${B}.sfs;
done

cp /projectnb/davieslab/jfifer/Japan_rad/AFS-analysis-with-moments/multimodel_inference/allmodels_unfolded allmodels_unfolded

NREPS=6 # number of random restarts per model per bootstrap rep
>mods
for i in `seq 1 $NREPS`;do
cat allmodels_unfolded >>mods;
done

sed -i -e 's/^/\/projectnb\/davieslab\/jfifer\/Japan_rad\/AFS-analysis-with-moments\/multimodel_inference\/py2\//' mods

#Calculating mutation rate per generation per  2bRAD-sequenced genome fraction
#using the 0.138% per Ma substitution rate (1.38e-9) in Prada et al 2014
#using the equation for growth rate using the equation r= 2.13+0.248X  in Klein & Loya 1991 https://www.researchgate.net/profile/Yossi-Loya/publication/250215112_Skeletal_growth_and_density_patterns_of_two_Pontes_corals_from_the_Gulf_of_Eilat_Red_Sea/links/0a85e536b78f5a1cce000000/Skeletal-growth-and-density-patterns-of-two-Pontes-corals-from-the-Gulf-of-Eilat-Red-Sea.pdf
# and reproductive age of lutea (8cm diameter) in Harriot et al 1983 https://link.springer.com/content/pdf/10.1007/BF00304727.pdf
#16=2.13+0.248X
#x=56 years
#OK thats crazy instead lets use the growth rate of 1.3 + 0.3 cm/year calculated for lobata by Patzold 1984 https://link.springer.com/content/pdf/10.1007/BF00263758.pdf
#6 years gen time
#mutation rate per gen (1.38e-9*6)= 8.28e-9
#size of ref genome
#grep "length=" AllSitesDemo.vcf | sed 's/.*th=//' | sed 's/>//' | awk '{ sum += $1 } END { print sum }'
#552020673
#number of new mutations per genome per generation
#8.28e-9 x 552020673= 4.57
#fraction of genome sequenced
#wc -l sites2do_filtBaye
#360608
#360608/552020673 =0.00065
#u=4.57*0.00065=0.0030


#CHANGE MUTATION AND GEN TIME TO REFLECT PORITES
CONTRAST=pop1pop3  # name of population comparison, should match the leading part of the bootstapped SFS names
ARGS="pop1 pop3 46 27 0.003 0.006" # pop1, pop2, projection for pop1, projection for pop2, mutation rate (per genotyped portion of the genome per generation), generation time$
CONTRAST=pop1pop2
ARGS="pop1 pop2 46 33 0.003 0.006"
CONTRAST=pop2pop3
ARGS="pop2 pop3 33 27 0.003 0.006"

rm -f modsel args

>modsel
for B in `seq 1 10`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat mods | wc -l`
>args
>${CONTRAST}.stdout
for i in `seq 1 $NMODELS`; do echo "$INFILE $ARGS >>${CONTRAST}.stdout & PID=\$! && sleep 360m && kill -9 \$PID " >>args; done; paste mods args -d " " >>modsel; done


CONTRAST=pop1pop3
#CONTRAST=pop1pop2
#CONTRAST=pop2pop3
#have to add -a to each grep because for some reason it reading as binary, I think because of weird ^@ characters?
#need to remove weird characters will fuck up lines further down if you don't
sed -i 's/\x00//g' ${CONTRAST}.stdout

grep -a RESULT ${CONTRAST}.stdout -A 4 | grep -aE "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep -a RESULT  >${CONTRAST}.res
#removing any lines that have only partial output
for i in *.res; do awk -F' ' '$6 != ""' $i > tmp && mv tmp $i; done
#go through .res in nano as well and make sure there aren't any lines that were cut off (delete those lines)

rm winner*

#Use /projectnb/davieslab/jfifer/Japan_rad/Demographic_analysis_II/North_Only/Unlinked_Analyses/modSel_summary_JF.R to create .winmod
### To bootstrap winning model

NREPS=6 # number of random restarts per model per bootstrap rep
>winner.mods
for i in `seq 1 $NREPS`;do
awk '{print $2}' *.winmod | sed 's/$/.py/'  >>winner.mods;
done


sed -i -e 's/^/\/projectnb\/davieslab\/jfifer\/Japan_rad\/AFS-analysis-with-moments\/multimodel_inference\/py2\//' winner.mods

CONTRAST=pop1pop3  # name of population comparison, should match the leading part of the bootstapped SFS names
ARGS="pop1 pop3 46 27 0.003 0.006" # pop1, pop2, projection for pop1, projection for pop2, mutation rate (per genotyped portion of the genome per generation), generation time$
#CONTRAST=pop1pop2
#ARGS="pop1 pop2 46 33 0.003 0.006"
#CONTRAST=pop2pop3
#ARGS="pop2 pop3 33 27 0.003 0.006"

>winner.modsel
for B in `seq 1 100`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat winner.mods | wc -l`
>winner.args
>${CONTRAST}.winboots
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS >>${CONTRAST}.winboots & PID=\$! && sleep 360m && kill -9 \$PID " >>winner.args;
done;
paste winner.mods winner.args -d " " >>winner.modsel;
done

#########after booting



CONTRAST=pop1pop3
#CONTRAST=pop1pop2
#CONTRAST=pop2pop3

#have to add -a to each grep because for some reason it reading as binary, I think because of weird ^@ characters?
#need to remove any lines that dont have all elements (just do in nano)
sed -i 's/\x00//g' ${CONTRAST}.winboots
grep -a RESULT ${CONTRAST}.winboots -A 4 | grep -aE "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep -a RESULT >${CONTRAST}.winboots.res
for i in *winboots.res; do awk -F' ' '$6 != ""' $i > tmp && mv tmp $i; done


