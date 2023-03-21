export GENOME_FASTA=/projectnb/davieslab/jfifer/Baum/Genome_ref/Host/plut_final_2.1.fasta
/projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS dadi ../pop1.out.saf.idx ../pop2.out.saf.idx -sfs ../pop1.sfs -sfs ../pop2.sfs -ref $GENOME_FASTA -anc $GENOME_FASTA >12dadiout
/projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS dadi ../pop1.out.saf.idx ../pop3.out.saf.idx -sfs ../pop1.sfs -sfs ../pop3.sfs -ref $GENOME_FASTA -anc $GENOME_FASTA >13dadiout
/projectnb/davieslab/jfifer/Japan_rad/angsd/misc/realSFS dadi ../pop2.out.saf.idx ../pop3.out.saf.idx -sfs ../pop2.sfs -sfs ../pop3.sfs -ref $GENOME_FASTA -anc $GENOME_FASTA >23dadiout

#num individuals

/projectnb/davieslab/jfifer/Japan_rad/admixed_pops_analysis/4pop_analysis/Demographic_Analysis_II/Unlinked_analyses/Stairwell_analyses/realsfs2dadi.pl 12dadiout 29 21 >pop12_dadi.data
/projectnb/davieslab/jfifer/Japan_rad/admixed_pops_analysis/4pop_analysis/Demographic_Analysis_II/Unlinked_analyses/Stairwell_analyses/realsfs2dadi.pl 13dadiout 29 17 >pop13_dadi.data
/projectnb/davieslab/jfifer/Japan_rad/admixed_pops_analysis/4pop_analysis/Demographic_Analysis_II/Unlinked_analyses/Stairwell_analyses/realsfs2dadi.pl 23dadiout 42 17 >pop23_dadi.data

module load python2
export PYTHONPATH=$PYTHONPATH:$HOME/moments

PROJECTION=58 # 2 x (number of individuals in population 0)
/projectnb/davieslab/jfifer/Japan_rad/admixed_pops_analysis/4pop_analysis/Demographic_Analysis_II/Unlinked_analyses/Stairwell_analyses/1dAFS.py pop12_dadi.data pop0 $PROJECTION
mv 1dsfss pop1.swout

PROJECTION=42 # 2 x (number of individuals in population 1)
/projectnb/davieslab/jfifer/Japan_rad/admixed_pops_analysis/4pop_analysis/Demographic_Analysis_II/Unlinked_analyses/Stairwell_analyses/1dAFS.py pop12_dadi.data pop1 $PROJECTION
mv 1dsfss pop2.swout

PROJECTION=34 # 2 x (number of individuals in population 1)
/projectnb/davieslab/jfifer/Japan_rad/admixed_pops_analysis/4pop_analysis/Demographic_Analysis_II/Unlinked_analyses/Stairwell_analyses/1dAFS.py pop13_dadi.data pop1 $PROJECTION
mv 1dsfss pop3.swout

#create pop1.blueprint
#nseq is the same as projection
#L is wc -l sites2do_filtBaye
#mu should be same mut rate used in moments analysis (i.e. incorporates fraction of genome sequenced). This is just my intuition since stairway assums mutation rate to be for the whole genome, number of mutations should be smaller if only looking at a fraction
#^this is wrong, mutation rate should be per gen per base. I ran with the 0.003 mutation rate and it was scaled to like 5 years.
#SFS is taken from the .swout file

module load java
java -cp /usr4/bi594/jfifer/bin/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder pop1.blueprint18


#add this to the top of pop1.blueprint.sh
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sw_test1  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
##$ -l h_rt=200:00:00 #maximum run time
#$ -M james.e.fifer@gmail.com #your email
#$ -m be
##$ -pe omp 24
