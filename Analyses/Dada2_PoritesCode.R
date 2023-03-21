##This is the code used to assign ASVs for Porites samples from KI
##This code was adapted from https://benjjneb.github.io/dada2/ITS_workflow.html by Sam Starko and Danielle Claar 

# Load necessary libraries and programs
library(dada2)
library(phyloseq)
library(Biostrings)
# Make sure you have access to cutadapt on your machine
# See https://cutadapt.readthedocs.io/en/stable/installation.html for installation information
cutadapt <- "./dada2/cutadapt.exe" # Where is cutadapt located on your machine
system2(cutadapt, args = "--version") # Check that cutadapt works

# https://benjjneb.github.io/dada2/ITS_workflow.html

# Create file lists
path <- "./dada2/PoritesOnly" # Path to where sequences are stored
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)


# Identify primers used for this project
FWD <- "GTGAATTGCAGAACTCCGTG" # Note that these are only the ITS2 primers (do not include the Illumina adapters, which are already removed)
REV <- "CCTCCGCTTACTTATATGCTT" # Note that these are only the ITS2 primers (do not include the Illumina adapters, which are already removed)

# Verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), 
               Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# Pre-filter the sequences just to remove those with Ns, and trim reverse reads to deal with low quality reverse reads from run4
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN,  maxN = 0, # dada does not allow Ns, so must be zero
              multithread = TRUE, trimLeft = c(0,53), #Trim reverse reads at position 54 due to issues with the quality of reverse reads for run4
              compress = FALSE) # Turn off compression (default is TRUE), cutadapt doesn't work on compressed files

#Good to this point

# Count the number of times the primers appear in the forward and reverse read, while considering all possible primer orientations.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, 
                         ShortRead::sread(ShortRead::readFastq(fn)), 
                         fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, 
                                fn = fnFs.filtN[[2]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, 
                                fn = fnRs.filtN[[2]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, 
                                fn = fnFs.filtN[[2]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, 
                                fn = fnRs.filtN[[2]]))
# Check to make sure that something shows up for both forward and reverse - this can be a good "spell check", the first time I did it I didn't get any in forward, due to a typo.


# Use cutadapt to remove primers
path.cut <- file.path(path, "cutadapt") # Append the filepath to the filenames - cutadapt
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs)) # Append the filepath to the filenames - forward reads for cutadapt
fnRs.cut <- file.path(path.cut, basename(fnRs)) # Append the filepath to the filenames - reverse reads for cutadapt

FWD.RC <- dada2:::rc(FWD) # Create object that is reverse complement of forward reads
REV.RC <- dada2:::rc(REV) # Create object that is reverse complement of reverse reads
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 


# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# Sanity check to make sure that cutadapt worked and all primers are gone
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))




# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_R")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Look at quality profiles of trimmed reads
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])


# Assigning the filenames for the output of the filtered reads to be stored as fastq.gz files.
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

# ref_seqs <- readDNAStringSet(filepath = "ITS2db_fromSymPortal2.fasta")
ref_seqs <- readDNAStringSet(filepath = "./Dada2/PoritesOnly/ITS2db_trimmed_derep_dada.fasta")
min(width(ref_seqs))
max(width(ref_seqs))


out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, 
                     maxN = 100, # Standard filtering parameter
                     maxEE = c(2, 2),
                     minLen = 100, # Enforce a minLen here, to get rid of spurious short sequences; was 50 before
                     rm.phix = TRUE, # Remove any phiX sequences
                     compress = TRUE, 
                     multithread = TRUE, # on windows, set multithread = FALSE
                     verbose = TRUE, # Sets the maximum number of “expected errors” allowed in a read
                     trimRight = 80 # The number of nucleotides to remove from the end of each read - chosen based on declining quality scors in plot
                     )

head(out)
out2<-out[out[,2]==0,]
head(out2)
##Remove file names that were filtered out during above filtering step
filtFs2<-filtFs[out[,2]!=0]
filtRs2<-filtRs[out[,2]!=0]

# Learn error rates
# Please ignore all the “Not all sequences were the same length.” messages in the next couple sections. We know they aren’t, and it’s OK!
errF <- learnErrors(filtFs2, multithread = TRUE)
errR <- learnErrors(filtRs2, multithread = TRUE)
# Visualize the estimated error rates as a sanity check
plotErrors(errF, nominalQ = TRUE)

# Dereplicate identical reads
derepFs <- derepFastq(filtFs2, verbose = TRUE)
derepRs <- derepFastq(filtRs2, verbose = TRUE)
# Name the derep-class objects by the sample names #SS - I did not include this, since sample names did not parse properly
names(derepFs) <- sample.names[out[,2]!=0]
names(derepRs) <- sample.names[out[,2]!=0]

# Sample Inference - apply the core sample inference algorithm to the dereplicated data.
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                      trimOverhang=TRUE, verbose=TRUE)

# Construct amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
# Inspect distribution of sequence lengths:
hist(nchar(getSequences(seqtab.nochim)))

# Track reads through the pipeline - inspect the # of reads that made it through each step in the pipeline to verify everything worked as expected.
getN <- function(x) sum(getUniques(x))
track <- cbind(out[out[,2]!=0], sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

# Assign taxonomy (Porites matches only)
sym.ref <- "./Dada2/PoritesOnly/ITS2db_trimmed_derep_dada.fasta"

# Try with in-house reference database down to Genus
taxa <- assignTaxonomy(seqtab.nochim, sym.ref, multithread = TRUE, 
                       tryRC = TRUE, minBoot = 80, verbose = TRUE)
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
tp <- data.frame(taxa.print)
tp$Genus
unique(tp$Genus)

# Import to phyloseq
samdf <- read.csv("./dada2/PoritesOnly/Porites_mapping_Jan17.csv",
                    header = TRUE) # Read in sample data
rownames(samdf) <- samdf$SampleID
head(samdf)

# We now construct a phyloseq object directly from the dada2 outputs.

samdf.ps<-sample_data(samdf)
sample_names(samdf.ps)<-samdf$fastq_file_name
#sample_names(samdf.ps)<-substr(sample_names(samdf.ps),1,nchar(sample_names(samdf.ps)))
sample_names(samdf.ps)<-gsub("_R1_001.fastq", "", x = sample_names(samdf.ps))
sample_names(samdf.ps) %in% (otu_table(seqtab.nochim, taxa_are_rows=FALSE) %>% sample_names())
ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),  samdf.ps,
               tax_table(taxa))

# It is more convenient to use short names for our ASVs 
# (e.g. ASV21) rather than the full DNA sequence 
# when working with some of the tables and visualizations 
# from phyloseq, but we want to keep the full DNA 
# sequences for other purposes like merging with other 
# datasets or indexing into reference databases like the 
# Earth Microbiome Project. For that reason we’ll store 
# the DNA sequences of our ASVs in the refseq slot of the 
# phyloseq object, and then rename our taxa to a short string. 
# That way, the short new taxa names will appear 
# in tables and plots, and we can still recover the DNA sequences 
# corresponding to each ASV as needed with refseq(ps).
dna <- Biostrings::DNAStringSet(taxa_names(ps1))
taxa_names(ps1) <- paste0("ASV", seq(ntaxa(ps1)))
names(dna) <- taxa_names(ps1)
psHost <- merge_phyloseq(ps1, dna)
psHost

#This phyloseq object ("psHost") contains all ASVs (host and symbiont)
save(psHost, file="./Dada2/PoritesOnly/Porites_ASV_Jan17_2022.RData")

#Isolate and export ASVs that were annotated as Porites
Por.host<-which(tax_table(psHost)[,"Genus"]=="g__Porites")
Por.Clad<-which(tax_table(psHost)[,"Genus"]=="g__Cladocopium")
write.csv(dna[Por.host], "Porites_ASVs_NoSubset_Jan17_2022.csv")
write.csv(dna[Por.Clad], "Cladocopium_ASVs_NoSubset_Jan17_2022.csv")

#Convert csv to fasta in texteditor and open in geneious (or similar) to confirm that all sequences are in fact Porites; if need be, add additional reference sequences and reassign taxonomy

#Create phyloseq object of JUST Porites sequences
ps_otu <- otu_table(psHost) %>% as.data.frame()
ps_otu2 <- ps_otu[,Por.host]
ps_otu2 %>% dim()

Porites_ASV <- phyloseq(otu_table(ps_otu2, taxa_are_rows = FALSE), sample_data(psHost), tax_table(psHost))
#save
save(Porites_ASV, file="./Dada2/PoritesOnly/Porites_Host_ASV_Jan22_2022.RData")

#Create phyloseq object of JUST Cladocopium 
ps_otu3 <- otu_table(psHost) %>% as.data.frame()
ps_otu4 <- ps_otu[,Por.Clad]
ps_otu4 %>% dim()

Clad_ASV <- phyloseq(otu_table(ps_otu4, taxa_are_rows = FALSE), sample_data(psHost), tax_table(psHost))
#save
save(Clad_ASV, file="./Dada2/PoritesOnly/Cladocopium_ASV_Jan17.RData")

#Export host ASV matrix combined with some metadata (remove columns that aren't needed manually)
host_seq <- cbind(sample_data(psHost) %>% data.frame(), ps_otu2)
write.csv(host_seq, file = "./Dada2/PoritesOnly/PoritesHostASVmatrix_Jan17.csv")


#Now extract the dominant Porites ASVs from each sample
# Read in the matrix of values
matrix_wmeta <- read_csv("./dada2/PoritesOnly/PoritesHostASVmatrix_Jan17_2.csv")
matrix_wmeta2 <- matrix_wmeta %>% filter(matrix_wmeta$TotalSeq > 0)
matrix <- matrix_wmeta2[,c(9:215)]
matrix_wmeta %>% filter(matrix_wmeta$TotalSeq > 0)

# Find the column name of the cell with the highest value for each row
highest_col <- apply(matrix, 1, function(x) names(which.max(x)))

# Find the column name of the second highest value for each row
# Only include this value if it makes up more than 40% of the total values
second_highest_col <- apply(matrix, 1, function(x) {
  sorted_x <- sort(x, decreasing = TRUE)
  second_highest_value <- sorted_x[2]
  second_highest_value_pct <- second_highest_value / sum(x)
  if (second_highest_value_pct > 0.4) {
    names(which(x == second_highest_value))
  } else {
    ""
  }
})

# Add the new columns to the matrix
matrix$highest_col <- highest_col
matrix$second_highest_col <- second_highest_col

#
new <- data.frame(dominant1 = matrix$highest_col, dominant2 = matrix$second_highest_col, matrix_wmeta2[,c(1:8)])
write.csv(new, "./dada2/PoritesOnly/PoritesDominantASVs.csv")
