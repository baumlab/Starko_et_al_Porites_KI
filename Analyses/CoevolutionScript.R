#This script is for testing for a correlation between host phylogeny (based on RADSeq data) and coral symbiont (based on individual DIVs, ITS2 profiles and unifrac distance based on symportal outputs)

#Load library

#Read in tree file
PorTree <- read.tree("./Data/Porites_RAD_likelihoodTree.tree")
#Read in reference data
Coevol_meta <- read_csv("./Data/Porites_coevolution_meta.csv")
#Check whether tip labels and datasheet match up
PorTree$tip.label %in% Coevol_meta$filename
#Remove samples not included in tree (e.g., duplicates or failed samples)
Coevol_meta <- Coevol_meta %>% filter(Coevol_meta$filename %in% PorTree$tip.label)
dim(Coevol_meta)

Porites_DIV #The phyloseq object containing DIVs for each sample
C15_tree  <- read.tree("./Data/Trees/C15_DIVTree.newick")
#Drop non-C15 reads and re-normalize
dev <- otu_table(phyloseq(otu_table(Porites_DIV),sample_data(Porites_DIV), C15_tree)) %>% rowSums() %>% as_tibble() %>% as.data.frame()
mat <- otu_table(phyloseq(otu_table(Porites_DIV),sample_data(Porites_DIV), C15_tree)) %>% as.data.frame()
mat2 <- mapply("/",mat, dev) %>% as.data.frame()
rownames(mat2) <- rownames(sample_data(Porites_DIV))
mat2<- drop_na(mat2) 

#Create phyloseq of just C15 sequences
C15_ps<- phyloseq(otu_table(mat2, taxa_are_rows = FALSE), sample_data(Porites_DIV), C15_tree)

Coevol_meta$Sample_symbio2[which(!Coevol_meta$Sample_symbio2 %in% sample_data(C15_ps)$sample_name)]

C15_coev <- C15_ps %>% subset_samples(sample_name %in% Coevol_meta$Sample_symbio2)
otu <- otu_table(C15_coev) %>% as.data.frame()
C15_coev2 <- phyloseq(otu_table(otu[,which(colSums(otu)>0)], taxa_are_rows = FALSE), sample_data(C15_coev), C15_tree)


##Next calculate Moran's I
PorTree$tip.label
dist.mat <- cophenetic(PorTree) 

Coevol_meta
PorTree$tip.label <- Coevol_meta$Sample_symbio2[match(PorTree$tip.label, Coevol_meta$filename)] ##HERE
dist.mat <- 1/cophenetic(PorTree) 
diag(dist.mat) <- 0
sample_data(C15_coev2)$coral_tag
C15matrix <- otu_table(C15_coev2) %>% data.frame()
rownames(C15matrix) <- sample_data(C15_coev2)$sample_name
C15matrix <- C15matrix[order(match(rownames(C15matrix), PorTree$tip.label)), ]

criter <- data.frame()
for (col in colnames(C15matrix)) {
  x <- sum(C15matrix[,col] > 0)
  criter[col, 1] <- print(x)
}

C15matrix2 <- C15matrix[,criter > 6]
C15matrix3 <- C15matrix2 / rowSums(C15matrix2)

C15matrix_filtered <- C15matrix[,criter > 0] #For subsequent analysis

#Calculate for all DIVs
Mtest <-data.frame()
for (col in colnames(C15matrix)) {
  x <- C15matrix[,col]
  I <- Moran.I(x, dist.mat)
  Mtest[col,1]<-print(I$p.value)
}

#Calculate for only DIVs that were present in > 10% of samples
Mtest2 <-data.frame()
for (col in colnames(C15matrix3)) {
  x <- C15matrix3[,col]
  I <- Moran.I(x, dist.mat)
  Mtest2[col,1]<-print(I$p.value)
}

Mtest2$V2 <- p.adjust(Mtest2$V1, method = "BH")
Mtest2$V2 %>% hist(breaks = c(seq(0,1, 0.1)))
Mtest2


C15matrix_pres <- ifelse(C15matrix3 > 0, 1, 0)
Mtest_pres <-data.frame()
for (col in colnames(C15matrix_pres)) {
  x <- C15matrix3[,col]
  I <- Moran.I(x, dist.mat)
  Mtest_pres[col,1]<-print(I$p.value)
}
Mtest_pres$V2 <- Mtest_pres$V1 %>% p.adjust(method = "BH")

plot(PorTree, cex = 0.5)

ggplot(C15matrix3, )

###
Por_Profiles_filtered <- Porites_Profiles %>% subset_samples(sample_data(Porites_Profiles)$sample_name %in% PorTree$tip.label)
Por_otu <- otu_table(Por_Profiles_filtered)
rownames(Por_otu) <- sample_data(Por_Profiles_filtered)$sample_name
Por_otu_filtered <- data.frame(otu_table(Por_Profiles_filtered))[c(13,17,18,21:24,26:29,31, 33:35, 38,39,45,46,48:57,60,62,63,66:68,70:72,75, 76, 79,80, 85:89)]
rownames(Por_otu_filtered) <- sample_data(Por_Profiles_filtered)$sample_name
Profile_matrix1 <- Por_otu_filtered #for analyses further down

Por_otu_filtered <- Por_otu_filtered[,colSums(Por_otu_filtered)>5]
Por_otu_filtered2 <- Por_otu_filtered[colnames(dist.mat),]

Ptest <-data.frame()
for (col in colnames(Por_otu_filtered2)) {
  x <- Por_otu_filtered2[,col]
  I <- Moran.I(x, dist.mat)
  Ptest[col,1]<-print(I$p.value)
}
Ptest$V2 <- Ptest$V1 %>% p.adjust(method = "BH")


Por_upgma <- phangorn::upgma(PorTree %>% cophenetic())
Por_upgma$tip.label




C15matrix4 <- C15matrix3
C15matrix4$sample <- rownames(C15matrix4)
C15matrix4$order <- match(C15matrix4 %>% rownames(),tip_order )



C15.long <- gather(C15matrix4,"DIV", "Value", -sample, -order)
ggplot(C15.long,aes(y = order, x = DIV, fill =Value))+
  geom_tile()+scale_fill_gradient(high = "black", low = "white")

Por_upgma <- phangorn::upgma(PorTree %>% cophenetic())
tip_order <- jntools::get_tips_in_ape_plot_order(Por_upgma)

#heatmap(C15matrix3 %>% as.matrix(), scale = "column", margins = c(10,10),Rowv = NA, Colv = NA)



###

#phylosig(PorTree, C15matrix3[,30], method="K", test=TRUE)
Ktest <-data.frame()
for (col in colnames(Por_otu_filtered2)) {
  x <-Por_otu_filtered2[,col]
  names(x) <- rownames(Por_otu_filtered2)
  K <- phylosig(PorTree, Por_otu_filtered2[,col], method="K", test=TRUE, nsim=999)
  Ktest[col,1]<-print(K$P)
}

Ktest$V2 <- Ktest$V1 %>% p.adjust(method = "BH") 
Ktest$V2 %>% hist(seq(0,1,0.05))

Ktest



Ktest <-data.frame()
for (col in colnames(C15matrix3)) {
  x <- C15matrix3[,col]
  names(x) <- rownames(C15matrix3)
  K <- phylosig(PorTree, C15matrix3[,col], method="K", test=TRUE, nsim=999)
  Ktest[col,1]<-print(K$P)
}

Ktest$V2 <- Ktest$V1 %>% p.adjust(method = "BH") 
Ktest$V2 %>% hist(seq(0,1,0.05))


#Test for relationship between phylogenetic distance and unifrac distance of symbiont community. Note, not significant due to similarity of PKir-1 and PKir-2 symbionts.
tree_matrix <- cophenetic(PorTree)


Clad_unifrac <- read_csv("./Data/Clad_unifracdist.csv")
unifrac_dist <- Clad_unifrac[2:2125] %>% as.matrix
rownames(unifrac_dist) <- colnames(unifrac_dist)


unifrac_dist2 <- unifrac_dist[rownames(unifrac_dist) %in% PorTree$tip.label,]
unifrac_dist3 <- unifrac_dist2[,colnames(unifrac_dist2) %in% PorTree$tip.label]
unifrac_dist3 %>% dim()
unifrac_dist4 <- unifrac_dist3[PorTree$tip.label, PorTree$tip.label]
unifrac_dist4 %>% nj() 

library(ade4)
mantel.rtest(as.dist(log(unifrac_dist4)), as.dist(tree_matrix), nrepet = 9999)
ecodist::MRM(as.dist(log(unifrac_dist4)) ~ as.dist(tree_matrix))

Clad_tree <- unifrac_dist4 %>% nj()
plot(Clad_tree)

#Procrustes analysis using paco
library(paco)

Profile_matrix1 #Host-symbiont association matrix (C15 profiles)
Profile_matrix2 <- Profile_matrix1
Profile_matrix2[Profile_matrix2>0] <- 1 
Profile_matrix2

C15matrix_filtered  #Host-symbiont association matrix (C15 DIVs)

Clad_dist_all <- read_csv("./Data/Symbionts/Cladocopium_Profiles_Unifrac_dist2.csv") %>% as.data.frame()
rownames(Clad_dist_all) <- Clad_dist_all$...1
Clad_dist_all <- Clad_dist_all[,c(-1,-2)]

Clad_dist_all <- Clad_dist_all[,colnames(Clad_dist_all) %in% colnames(Profile_matrix1)]
Clad_dist_all <- Clad_dist_all[rownames(Clad_dist_all) %in% colnames(Profile_matrix1),]
Clad_dist_all <- Clad_dist_all[,!colnames(Clad_dist_all) %in% c("C1", "C42a.C1.C42.2.C1az.C1b")]
Clad_dist_all <- Clad_dist_all[!rownames(Clad_dist_all) %in% c("C1", "C42a.C1.C42.2.C1az.C1b"),]
Profile_matrix2 <- Profile_matrix2[!colnames(Profile_matrix2) %in% c("C1", "C42a.C1.C42.2.C1az.C1b")]

Clad_dist_all<- Clad_dist_all %>% as.dist() %>% nj() %>% cophenetic()

PorTree_dist <- cophenetic(PorTree)

paco_data <- prepare_paco_data(H = PorTree_dist, P = Clad_dist_all, Profile_matrix2)
paco_data <- add_pcoord(paco_data, correction = 'cailliez')
z <- PACo(paco_data, nperm = 1000, seed = 99, symmetric = FALSE)
z1 <- PACo(paco_data, nperm = 1000, seed = 99, symmetric = FALSE, shuffled = TRUE) 
z2 <- PACo(paco_data, nperm = 1000, seed = 99, symmetric = FALSE, method = "swap") 
z_final <- PACo(paco_data, nperm = 10000, seed = 99, symmetric = FALSE, shuffled = TRUE)
z_links <- paco_links(z)
residuals_paco(z_links)


hist(z_final$shuffled, xlim=c(2, 3.2), las = 1, col = "lightblue", xlab = "Sum of squared residuals", main = '')
abline(v = z_final$gof$ss, col = "blue", lwd = 3)
