
#This chunk loads in the entire Baum Lab SymPortal dataset, filters these data to remove data that are not from Porites colonies and data where sequencing had little success (< 500 Symbiodineaceae sequences for a given sample)
#Note that while data outputs from this chunk will be shared publically, non-Porites data from the Baum Lab dataset are not available at this time. This means that this chunk of code will not work for users without access to the whole Baum lab dataset

symp_matrix <- read_excel("../Data/Symbionts/BaumLabSymportalOutput_Jan2022.xlsx", sheet = "profile_sequence_matrix")

symp_meta <- read_excel("../Data/Symbionts/BaumLabSymportalOutput_Jan2022.xlsx", sheet = "Metadata", col_types = "text")

#Combine dataframes
sym_combined<-full_join(symp_matrix, symp_meta, by="sample_name")

#Remove samples with less than 500 sequences
sym_combined <- sym_combined[sym_combined$Total_sequences>500,]

#Remove non-Porites data
sym_combined <- filter(sym_combined, host_genus == "Porites")

#Separate community matrix and metadata again
mat <- sym_combined[,2:277]
dev <- sym_combined[,278]

#Convert matrix to proportion
sym_prop <- mapply("/",as.data.frame(mat), dev)
sym_prop2 <- as.data.frame(cbind(sym_combined[,"sample_name"], sym_prop))

#Drop columns with colsums of 0
sym_prop2 <- sym_prop2[,2:277]
sym_prop3<-sym_prop2[,which(colSums(sym_prop2) > 0 )]

#Create metadata dataframe
meta <- as_tibble(cbind(sym_combined[,"sample_name"], sym_combined[,278:297]))


#Create phyloseq object that contains both proportion matrix (Profiles) and associated metadata for Porites

Porites_Profiles <- phyloseq(otu_table(sym_prop3, taxa_are_rows = FALSE), sample_data(meta))

save(Porites_Profiles, file = "../Data/Symbionts/Porites_Profiles.RData")

####NOW DO THE SAME FOR THE DIV MATRIX

symp_DIVmatrix <- read_excel("../Data/Symbionts/BaumLabSymportalOutput_Jan2022.xlsx", sheet = "post_MED_sequence_matrix")


#Combine dataframes
sym_combined<-full_join(symp_DIVmatrix, symp_meta, by="sample_name")

#Remove samples with less than 500 sequences
sym_combined <- sym_combined[sym_combined$Total_sequences>500,]

#Remove non-Porites data
sym_combined <- filter(sym_combined, host_genus == "Porites")

#Separate community matrix and metadata again
mat <- sym_combined[,2:2392]
dev <- sym_combined[,2393]

#Convert matrix to proportion
sym_DIVprop <- mapply("/",as.data.frame(mat), dev)
sym_DIVprop2 <- as.data.frame(cbind(sym_combined[,"sample_name"], sym_DIVprop))

#Drop columns with colsums of 0
sym_DIVprop2 <- sym_DIVprop2[,2:2392]
sym_DIVprop3<-sym_DIVprop2[,which(colSums(sym_DIVprop2) > 0 )]

#Create metadata dataframe
meta <- as_tibble(cbind(sym_combined[,"sample_name"], sym_combined[,2393:2412]))


#Create phyloseq object that contains both proportion matrix (Profiles) and associated metadata for Porites

Porites_DIV <- phyloseq(otu_table(sym_DIVprop3, taxa_are_rows = FALSE), sample_data(meta))

#Save data file for use with subsequent analyses

#save(Porites_DIV, file = "../Data/Symbionts/Porites_DIV.RData")
