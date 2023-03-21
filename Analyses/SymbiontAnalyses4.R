#####################
###This code includes exploratory and final analyses for Cladocopium symbionts through time 
######################

#Create phyloseq objects that include only samples before and samples after the heatwave
#Samples from before heatwave
sym_before <- subset_samples(Porites_Profiles, sample_name %in% AllLins$UniqueID_symbio_before2)

#Samples from after heatwave
sym_after <- subset_samples(Porites_Profiles, sample_name %in% AllLins$UniqueID_symbio_after2)

#Extract dataframes from phyloseq object
#Before heatwave
sym_before_df<- as.data.frame(otu_table(sym_before))
sym_before_df$timeframe <- "Before"
sym_before_meta <- as.data.frame(sample_data(sym_before)) %>% as_tibble()
sym_before_df$coral_tag <- sym_before_meta$coral_tag 

#After heatwave
sym_after_df <- as.data.frame(otu_table(sym_after))
sym_after_df$timeframe <- "After"
sym_after_meta <- as.data.frame(sample_data(sym_after)) %>% as_tibble()
sym_after_df$coral_tag <- sym_after_meta$coral_tag 

#Combined
sym_combined <- rbind(sym_before_df, sym_after_df)
lineage.df <- AllLins[,c("Lineage", "coral_tag")]
sym_compare2 <- left_join(sym_combined, lineage.df, by = "coral_tag")
sym_compare2 %>% dim()

#Prep data for barplot by timeframe and lineage
sym_average <- sym_compare2 %>% 
  group_by(timeframe, Lineage) %>% 
  summarise_all(funs(mean)) 

#Convert data to long-form
sym_long<- sym_average %>% filter(Lineage!="NA") %>% gather(key="Profile", value="Abundance", 3:115)

sym_long$timeframe <- factor(sym_long$timeframe, levels = c("Before", "After"), ordered = TRUE)

#Import colour scheme
colour_scheme <- read_csv("./Analyses/colour_scheme.csv")

#Plot symbionts by lineage
ggplot(sym_long, aes(fill=Profile, y=Abundance, x=Lineage))+ 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~timeframe) +
  theme(legend.position = "none") +
  scale_fill_manual(values=colour_scheme$hex)

#Warning refers to profiles not present in any samples here

ggplot(filter(sym_long, Lineage!="Unknown"), aes(fill=Profile, y=Abundance, x=Lineage))+ 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~timeframe) +
  theme_cowplot(4) +
  scale_fill_manual(values=colour_scheme$hex2)
  #theme(legend.position = "none")


####################################
#This is for plotting symbionts through time for all colonies

#Extract dataframes from phyloseq object
sym_df<- as.data.frame(otu_table(Porites_Profiles))
sym_meta <- as.data.frame(sample_data(Porites_Profiles)) %>% as_tibble()
sym_meta$expedition <- gsub("2015a - Post", "2015a", sym_meta$expedition)
sym_meta$expedition <- gsub("2015a - Pre", "2015a", sym_meta$expedition)

sym_df$coral_tag <- sym_meta$coral_tag 
sym_df$expedition <- sym_meta$expedition
sym_df$sample_name <- sym_meta$sample_name


lineage.df <- AllLins[,c("Lineage", "coral_tag")]
sym_compare2 <- left_join(sym_df, lineage.df, by = "coral_tag")
sym_compare2 %>% dim()

sym_average <- sym_compare2 %>% 
  group_by(expedition, Lineage) %>% 
  summarise_all(funs(mean)) 

sym_long<- sym_average %>% filter(Lineage!="NA") %>% gather(key="Profile", value="Abundance", 3:115)

sym_long$expedition <- factor(sym_long$expedition, levels = c("2014", "2015a", "2015b", "2015c", "2016a", "2016b", "2017", "2018", "2019"), ordered = TRUE)

#Import colour scheme
colour_scheme <- read_csv("./Analyses/colour_scheme.csv")

#Plot symbionts by lineage
ggplot(sym_long, aes(fill=Profile, y=Abundance, x=Lineage))+ 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~expedition) +
  theme(legend.position = "none") +
  scale_fill_manual(values=colour_scheme$hex)

ggplot(sym_long, aes(fill=Profile, y=Abundance, x=expedition))+ 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~Lineage) +
  theme(legend.position = "none") +
  scale_fill_manual(values=colour_scheme$hex)+
  ylab("Relative Abundance")+
  xlab("Expedition")
#Warning describes profiles in the phyloseq object that are not present in these samples (only present in samples from other sites, not included in this paper)


###ORDINATIONS####
#This script deals with data from the DIV matrix. 

#Create phyloseq objects that include only samples before and samples after the heatwave (note that the Dominant_DIV column must be merged from AllLins dataframe)
DIV_before <- subset_samples(Porites_DIV, sample_name %in% AllLins$UniqueID_symbio_before2)

new.col <- data.frame(Dominant_DIV = AllLins$Dominant_DIV_before, sample_name = AllLins$UniqueID_symbio_before2, AllLins$Lineage)
DIV_before_meta <- sample_data(DIV_before) %>% as.data.frame() %>% as_tibble()
DIV_before_meta2 <- left_join(DIV_before_meta, new.col) %>% as.data.frame()
rownames(DIV_before_meta2) <- DIV_before_meta2$sample_name
otu_before <- otu_table(DIV_before) %>%  as.data.frame() %>% as_tibble() %>% as.data.frame()
rownames(otu_before) <- sample_data(DIV_before)$sample_name
DIV_before <- phyloseq(otu_table(otu_before, taxa_are_rows = FALSE), sample_data(DIV_before_meta2))
sample_data(DIV_before)$coral_tag <- gsub("441_1304", "441", sample_data(DIV_before)$coral_tag)



DIV_after <- subset_samples(Porites_DIV, sample_name %in% AllLins$UniqueID_symbio_after2)

new.col <- data.frame(Dominant_DIV = AllLins$Dominant_DIV_after, sample_name = AllLins$UniqueID_symbio_after2, AllLins$Lineage)
DIV_after_meta <- sample_data(DIV_after) %>% as.data.frame() %>% as_tibble()
DIV_after_meta2 <- left_join(DIV_after_meta, new.col) %>% as.data.frame()
rownames(DIV_after_meta2) <- DIV_after_meta2$sample_name
otu_after <- otu_table(DIV_after) %>% as.data.frame() %>% as_tibble() %>% as.data.frame()
rownames(otu_after) <- sample_data(DIV_after)$sample_name
DIV_after <- phyloseq(otu_table(otu_after, taxa_are_rows = FALSE), sample_data(DIV_after_meta2))
sample_data(DIV_after)$coral_tag <- gsub("441_1304", "441b_1304", sample_data(DIV_after)$coral_tag)


#Merge these phyloseq objects
DIV_comb<- merge_phyloseq(DIV_before , DIV_after)

#Remove columns that are all zero reads
DIV_comb_otu <- otu_table(DIV_comb)
DIV_sample <- sample_data(DIV_comb)
DIV_comb_otu2 <- DIV_comb_otu[,which(colSums(DIV_comb_otu) > 0 )]
DIV_comb2<- phyloseq(DIV_sample, DIV_comb_otu2)

######Exploratory ordinations
##All samples
ord <- ordinate(DIV_comb2, method = "PCoA", distance = "bray")
plot_ordination(DIV_comb2, ord, color = "expedition")+
  scale_colour_manual(values = c("Blue", "Blue", "Blue", "Blue", "Blue", "Red", "Red", "Red", "Red"))+
  ggtitle("All lineages PCoA")

plot_ordination(DIV_comb2, ord, color = "Dominant_DIV")

#Call lineages
#Lineage 3 only
DIV_lineage3 <- DIV_comb2 %>% subset_samples(coral_tag %in% Lineage3$coral_tag)
#Lineage 1 only
DIV_lineage1 <- DIV_comb2 %>% subset_samples(coral_tag %in% Lineage1$coral_tag)
#Lineage 2 only
DIV_lineage2 <- DIV_comb2 %>% subset_samples(coral_tag %in% Lineage2$coral_tag)
#Combine
DIV_combined <- merge_phyloseq(DIV_lineage1, DIV_lineage2, DIV_lineage3)


###Unifrac ordinations
library(ape)

#Import Cladocopium ITS2 tree from file
Clad_tree <- read.tree("./Data/Trees/Cladocopium_DIVTree.newick")

#Drop non-Cladocopium reads and re-normalize
dev <- otu_table(phyloseq(DIV_comb_otu2, DIV_sample, Clad_tree)) %>% rowSums() %>% as_tibble() %>% as.data.frame()
mat <- otu_table(phyloseq(DIV_comb_otu2, DIV_sample, Clad_tree)) %>% as.data.frame()
mat2 <- mapply("/",mat, dev)
rownames(mat2) <- sample_data(DIV_sample) %>% rownames()

#Create phyloseq of just Cladocopium sequences
Clad_ps<- phyloseq(otu_table(mat2, taxa_are_rows = FALSE), DIV_sample, Clad_tree)

ord <- ordinate(Clad_ps, method = "PCoA", distance = "unifrac")
plot_ordination(Clad_ps, ord, color = "expedition")+
  scale_colour_manual(values = c("Blue", "Blue", "Blue","Blue", "Blue", "Red", "Red", "Red", "Red"))+
  ggtitle("All lineages PCoA - Unifrac (Cladocopium only)")

plot_ordination(Clad_ps, ord, color = "Dominant_DIV")+
  ggtitle("All lineages PCoA - Unifrac (Cladocopium only)")

plot_ordination(Clad_ps, ord, color = "Dominant_DIV")+
  ggtitle("All lineages PCoA - Unifrac (Cladocopium only)")

plot_ordination(Clad_ps, ord, color = "AllLins.Lineage")+
  ggtitle("All lineages PCoA - Unifrac (Cladocopium only)")
# scale_colour_manual(c("Blue", "Red", "Beige", "Grey"))


#Call lineages
#Lineage 3 only
Clad_lineage3 <- Clad_ps %>% subset_samples(coral_tag %in% Lineage3$coral_tag)
#Lineage 1 only
Clad_lineage1 <- Clad_ps %>% subset_samples(coral_tag %in% Lineage1$coral_tag)
#Lineage 2 only
Clad_lineage2 <- Clad_ps %>% subset_samples(coral_tag %in% Lineage2$coral_tag)


###C15 only
C15_tree  <- read.tree("./Data/Trees/C15_DIVTree.newick")
#Drop non-C15 reads and re-normalize
dev <- otu_table(phyloseq(DIV_comb_otu2, DIV_sample, C15_tree)) %>% rowSums() %>% as_tibble() %>% as.data.frame()
mat <- otu_table(phyloseq(DIV_comb_otu2, DIV_sample, C15_tree)) %>% as.data.frame()
mat2 <- mapply("/",mat, dev) %>% as.data.frame()
rownames(mat2) <- rownames(sample_data(DIV_sample))
mat2<- drop_na(mat2) 

#Create phyloseq of just C15 sequences
C15_ps<- phyloseq(otu_table(mat2, taxa_are_rows = FALSE), DIV_sample, C15_tree)
set.seed(9)
ord <- ordinate(C15_ps, method = "PCoA", distance = "unifrac")

plot_ordination(C15_ps, ord, color = "AllLins.Lineage", shape = "AllLins.Lineage")+
  ggtitle("All lineages PCoA - Unifrac (C15 only)")+
  scale_colour_manual(values =c("#DFBE99","#729EA1","#DB5375", "grey", "grey"))+
  theme_cowplot(10)+
  scale_shape_manual(values = c(15,16,17,3))+
  stat_ellipse(level = 0.95)

#####################
##Figure 2 panel#####
#####################

#Extract each year
C15_2014 <- subset_samples(C15_ps, expedition=="2014")
C15_2015b <- subset_samples(C15_ps, expedition=="2015b")
C15_2015c <- subset_samples(C15_ps, expedition=="2015c")
C15_before <- merge_phyloseq(C15_2014,C15_2015b,C15_2015c)

#C15_2016a <- subset_samples(C15_ps, expedition=="2016a")
C15_2016b <- subset_samples(C15_ps, expedition=="2016b")
C15_2017 <- subset_samples(C15_ps, expedition=="2017")
C15_2018 <- subset_samples(C15_ps, expedition=="2018")
C15_2019 <- subset_samples(C15_ps, expedition=="2019")

#Create phyloseq of only samples from after the heatwave
C15_after <- merge_phyloseq(C15_2016b, C15_2017, C15_2018, C15_2019)

##BEFORE HEATWAVE
ord1plot <- plot_ordination(C15_before, ord, color = "AllLins.Lineage", shape = "AllLins.Lineage")+
  ggtitle("Before (DIVs)")+
  scale_colour_manual(values =c("#DFBE99","#729EA1","#DB5375", "grey", "grey"))+
  theme_cowplot(10)+
  scale_shape_manual(values = c(15,16,17,3))+
  stat_ellipse(level = 0.95)+
  ylim(-0.5,0.8)+
  xlim(-0.7,0.7)

##AFTER HEATWAVE
ord2plot <- plot_ordination(C15_after, ord, color = "AllLins.Lineage", shape = "AllLins.Lineage")+
  ggtitle("After (DIVs)")+
  scale_colour_manual(values =c("#DFBE99","#729EA1","#DB5375", "grey", "grey"))+
  theme_cowplot(10)+
  scale_shape_manual(values = c(15,16,17,3))+
  stat_ellipse()+
  ylim(-0.5,0.8)+
  xlim(-0.7,0.7)

#Plots
ord1plot
ord2plot

#Note that NA and Unknown ellipses must be removed in post-processing of these figures; error refers to convergence failures for these ellipses

#Remove colonies of unknown lineage
C15_before2 <- subset_samples(C15_before, AllLins.Lineage!="Unknown")

#Remove colonies of unknown lineage
C15_after2 <- subset_samples(C15_after, AllLins.Lineage!="Unknown")

#Fit before ordination
ord_before <- ordinate(C15_before2, method = "PCoA", distance = "unifrac")

#Fit after ordination
ord_after <- ordinate(C15_after2, method = "PCoA", distance = "unifrac")

#Plot ordinations
plot_ordination(C15_before2, ord, color = "AllLins.Lineage", shape = "AllLins.Lineage")+
  ggtitle("All lineages PCoA - Unifrac (C15 only)")+
  scale_colour_manual(values =c("#DFBE99","#729EA1","#DB5375", "grey", "grey"))+
  theme_cowplot(10)+
  scale_shape_manual(values = c(15,16,17,3))+
  stat_ellipse(level = 0.99)

plot_ordination(C15_after2, ord, color = "AllLins.Lineage", shape = "AllLins.Lineage")+
  ggtitle("All lineages PCoA - Unifrac (C15 only)")+
  scale_colour_manual(values =c("#DFBE99","#729EA1","#DB5375", "grey", "grey"))+
  theme_cowplot(10)+
  scale_shape_manual(values = c(15,16,17,3))+
  stat_ellipse()


ord2<- plot_ordination(C15_after2, ord_after, color = "AllLins.Lineage", shape = "AllLins.Lineage")+
  ggtitle("After")+
  scale_colour_manual(values =c("#DFBE99","#729EA1","#DB5375", "grey", "grey"))+
  theme_cowplot(10)+
  scale_shape_manual(values = c(15,16,17,3))+
  stat_ellipse()

ord1 <- plot_ordination(C15_before2, ord_before, color = "AllLins.Lineage", shape = "AllLins.Lineage")+
  ggtitle("Before")+
  scale_colour_manual(values =c("#DFBE99","#729EA1","#DB5375", "grey", "grey"))+
  theme_cowplot(10)+
  scale_shape_manual(values = c(15,16,17,3))+
  stat_ellipse(level = 0.99)

#PERMANOVAs
library(vegan)
sample_data(C15_before2)$time <- "Before"
sample_data(C15_after2)$time <- "After"
metadata <- as(sample_data(C15_before2), "data.frame")

#PERMANOVA of ITS2 profiles versus lineage (before heatwave)
adonis2(phyloseq::distance(C15_before2, method="bray") ~ AllLins.Lineage,
       data = metadata)

metadata <- as(sample_data(C15_after2), "data.frame")

#PERMANOVA of ITS2 profiles versus lineage (after heatwave)
adonis2(phyloseq::distance(C15_after2, method="bray") ~ AllLins.Lineage,
        data = metadata)

#Make combined phyloseq object (before and after heatwave)
C15_2 <- merge_phyloseq(C15_before2, C15_after2)


##Significant interaction (with PERMANOVA) between timeframe (before vs after) and lineage
metadata <- as(sample_data(C15_2), "data.frame")
adonis2(phyloseq::distance(C15_2, method="bray") ~ AllLins.Lineage*time,
        data = metadata)

##Call lineages
#Lineage 3 only
C15_lineage3 <- C15_ps %>% subset_samples(coral_tag %in% Lineage3$coral_tag)
#Lineage 1 only
C15_lineage1 <- C15_ps %>% subset_samples(coral_tag %in% Lineage1$coral_tag)
#Lineage 2 only
C15_lineage2 <- C15_ps %>% subset_samples(coral_tag %in% Lineage2$coral_tag)


#######Ordination plot with dada2 ASV matrix
load(file="./Dada2/PoritesOnly/Cladocopium_ASV_Jan17.RData")
Clad_ASV
sample_data(Clad_ASV)
Clad_filtered <- subset_samples(Clad_ASV, !(sample_sums(Clad_ASV) < 500))
sample_data(Clad_filtered)$coral_tag <- gsub("314c_1088", "314c_1088_1263", sample_data(Clad_filtered)$coral_tag)
sample_data(Clad_filtered)$coral_tag <- gsub("1016a_1346", "1016_1346", sample_data(Clad_filtered)$coral_tag)
#change 441_1304 from 2014, 2015 to 441 and 2016 or later to 441b_1304 (these colonies were mistakenly lumped together in the spreadsheet but represent different corals from different lineages)
sample_data(Clad_filtered)$coral_tag[225] <- "441b_1304"
sample_data(Clad_filtered)$coral_tag[246] <- "441b_1304"
sample_data(Clad_filtered)$coral_tag[280] <- "441"
sample_data(Clad_filtered)$coral_tag[416] <- "441"
sample_data(Clad_filtered)$coral_tag[470] <- "441"
sample_data(Clad_filtered)$coral_tag[504] <- "441b_1304"
sample_data(Clad_filtered)$coral_tag[591] <- "441b_1304"



#Lineage 3 only
Clad_lineage3 <- Clad_filtered %>% subset_samples(coral_tag %in% Lineage3$coral_tag)
ord <- ordinate(Clad_lineage3, method = "PCoA")
plot_ordination(Clad_lineage3, ord, color = "expedition")+
  scale_colour_manual(values = c("Blue", "Blue", "Blue", "Blue", "Blue","Blue", "Yellow", "Red", "Red", "Red", "Red"))+
  ggtitle("Lineage 3 PCoA")

Clad_lineage3_meta <- sample_data(Clad_lineage3) %>% as.data.frame()
Clad_lineage3_meta$Lineage = "Pkir-3"
Clad_lineage3 <- phyloseq(sample_data(Clad_lineage3_meta), otu_table(Clad_lineage3), tax_table(Clad_lineage3))

#Lineage 1 only
Clad_lineage1 <- Clad_filtered %>% subset_samples(coral_tag %in% Lineage1$coral_tag)
ord <- ordinate(Clad_lineage1, method = "PCoA", distance = "bray")
plot_ordination(Clad_lineage1, ord, color = "expedition")+
  scale_colour_manual(values = c("Blue", "Blue", "Blue", "Blue", "Blue","Blue", "Yellow", "Red", "Red", "Red", "Red"))+
  ggtitle("Lineage 1 PCoA - Unifrac (Cladocopium only)")

Clad_lineage1_meta <- sample_data(Clad_lineage1) %>% as.data.frame()
Clad_lineage1_meta$Lineage = "Pkir-1"
Clad_lineage1 <- phyloseq(sample_data(Clad_lineage1_meta), otu_table(Clad_lineage1), tax_table(Clad_lineage1))


#Lineage 2 only
Clad_lineage2 <- Clad_filtered %>% subset_samples(coral_tag %in% Lineage2$coral_tag)
ord <- ordinate(Clad_lineage2, method = "PCoA", distance = "bray")
plot_ordination(Clad_lineage2, ord, color = "expedition")
  #scale_colour_manual(values = c("Blue", "Blue", "Blue", "Red", "Red", "Red", "Red"))+
  #ggtitle("Lineage 1 PCoA - Unifrac (Cladocopium only)")
  
  
Clad_lineage2_meta <- sample_data(Clad_lineage2) %>% as.data.frame()
Clad_lineage2_meta$Lineage = "Pkir-2"
Clad_lineage2 <- phyloseq(sample_data(Clad_lineage2_meta), otu_table(Clad_lineage2), tax_table(Clad_lineage2))


#Unassigned only
Clad_lineageUnk <- Clad_filtered %>% subset_samples(coral_tag %in% LineageUnk$coral_tag)
ord <- ordinate(Clad_lineageUnk, method = "PCoA", distance = "bray")
plot_ordination(Clad_lineageUnk, ord, color = "expedition")
#scale_colour_manual(values = c("Blue", "Blue", "Blue", "Red", "Red", "Red", "Red"))+
#ggtitle("Lineage 1 PCoA - Unifrac (Cladocopium only)")


Clad_lineageUnk_meta <- sample_data(Clad_lineageUnk) %>% as.data.frame()
Clad_lineageUnk_meta$Lineage = "Unassigned"
Clad_lineageUnk <- phyloseq(sample_data(Clad_lineageUnk_meta), otu_table(Clad_lineageUnk), tax_table(Clad_lineageUnk))


Clad_AllLins <- merge_phyloseq(Clad_lineage1, Clad_lineage2, Clad_lineage3, Clad_lineageUnk)
ord <- ordinate(Clad_AllLins, method = "PCoA", distance = "bray")
plot_ordination(Clad_AllLins, ord, color = "Lineage")+
scale_colour_manual(values = c("Blue", "Grey", "Red", "black"))+
  ggtitle("Cladocopium ASV - Bray-curtis")+
  stat_ellipse(level = 0.95)

sample_data(Clad_AllLins)$Unique_ID <- gsub(" - Pre", "-pre", sample_data(Clad_AllLins)$Unique_ID)
sample_data(Clad_AllLins)$Unique_ID <- gsub(" - Post", "-post", sample_data(Clad_AllLins)$Unique_ID)

Clad_AllLin_before <- Clad_AllLins %>% subset_samples(Unique_ID %in% AllLins$UniqueID_symbio_before2) %>% subset_samples(Unique_ID %in% sample_data(DIV_before)$sample_name)
Clad_AllLin_after <- Clad_AllLins %>% subset_samples(Unique_ID %in% AllLins$UniqueID_symbio_after2) %>% subset_samples(Unique_ID %in% sample_data(DIV_after)$sample_name)
Clad_AllLins2 <- merge_phyloseq(Clad_AllLin_before,Clad_AllLin_after)
ord <- ordinate(Clad_AllLins2, method = "PCoA", distance = "bray")

#ord <- ordinate(Clad_AllLin_before, method = "PCoA", distance = "bray")
ord3<- plot_ordination(Clad_AllLin_before, ord, color = "Lineage", shape = "Lineage")+
  ggtitle("Before (ASVs)")+
  stat_ellipse(level = 0.95)+
  scale_colour_manual(values =c("#DFBE99","#729EA1","#DB5375", "grey", "grey"))+
  theme_cowplot(10)+
  scale_shape_manual(values = c(15,16,17, 3))+
  ylim(-0.7,0.8)+
  xlim(-0.9,0.7)

#ord <- ordinate(Clad_AllLin_after, method = "PCoA", distance = "bray")
ord4<- plot_ordination(Clad_AllLin_after, ord, color = "Lineage", shape = "Lineage")+
  ggtitle("After (ASVs)")+
  scale_colour_manual(values =c("#DFBE99","#729EA1","#DB5375", "grey", "grey"))+
  theme_cowplot(10)+
  scale_shape_manual(values = c(15,16,17,3))+
  stat_ellipse(level = 0.95)+
  ylim(-0.7,0.8)+
  xlim(-0.9,0.7)

#Compare samples to ensure that all three datasets include the same samples
#Note that the ASV dataset had more samples pass filtered, these samples were removed to ensure comparability across panels
sym_before #Profiles
sym_before_check<- sample_data(sym_before)$sample_name #161 samples

sym_after #Profiles
sym_after_check<- sample_data(sym_after)$sample_name #162 samples

DIV_before #DIVs
DIV_before_check <- sample_data(DIV_before)$sample_name #161 samples

DIV_after #DIVs
DIV_after_check <- sample_data(DIV_after)$sample_name #162 samples

Clad_AllLin_before #ASVs
Clad_before_check <- sample_data(Clad_AllLin_before)$Unique_ID #161 samples

Clad_AllLin_after #ASVs
Clad_after_check <- sample_data(Clad_AllLin_after)$Unique_ID #162 samples

plot_grid(ord1plot,ord2plot,ord3,ord4)
#Convergence fail refers to the NA and unassigned ellipse calculations which are removed manually to produce the final figure


###Site by lineage by before vs after
sample_data(sym_before)$timeframe <- "Before"
sample_data(sym_after)$timeframe <- "After"
sym_com <- merge_phyloseq(sym_before, sym_after)

#Prep data for barplots
sym_df<- as.data.frame(otu_table(sym_com))
sym_meta <- as.data.frame(sample_data(sym_com)) %>% as_tibble()
sym_meta$expedition <- gsub("2015a - Post", "2015a", sym_meta$expedition)
sym_meta$expedition <- gsub("2015a - Pre", "2015a", sym_meta$expedition)

sym_df$coral_tag <- sym_meta$coral_tag 
sym_df$expedition <- sym_meta$expedition
sym_df$timeframe <- sym_meta$timeframe %>% factor(ordered = TRUE, levels = c("Before", "After"))
sym_df$site <- sym_meta$site_name

#Note that coral tag 1198 is accidentally listed as site 24 which does not exist in the context of this study, change to M3 which is 34 (the correct site for that coral)

pubnames = data.frame(site = c("14","15", "19", "20", "25", "27", "3", "30", "32", "33", "34", "35",
                              "38", "5", "8", "9", "40", "37", "24"),
                      site2 = c("M4", "VL1", "VL2", "VL5", "M5", "VH1", "L1", "VH3", "VH2", "M12", "M3", "M2",
                              "L2", "VL3", "M1", "L4", "H1", "VL4", "M3"))

sym_df <- left_join(sym_df, pubnames)


lineage.df <- AllLins[,c("Lineage", "coral_tag")]
sym_df$coral_tag <- gsub("1016a_1346", "1016_1346", sym_df$coral_tag) #"a" is meaningless
sym_df[155,13] <- as.numeric("1") #Have to change if initial dataframe changes, standardizing this sample
sym_compare2 <- left_join(sym_df, lineage.df, by = "coral_tag")
sym_compare2 %>% dim()



sym_average <- sym_compare2 %>% 
  group_by(timeframe, Lineage, site2) %>% 
  summarise_all(funs(mean)) 

sym_long<- sym_average %>% filter(Lineage!="NA") %>% gather(key="Profile", value="Abundance", 3:115, -site2)

sym_long$expedition <- factor(sym_long$expedition, levels = c("2014", "2015a", "2015b", "2015c", "2016a", "2016b", "2017", "2018", "2019"), ordered = TRUE)
sym_long$site2 <- factor(sym_long$site2, ordered = TRUE, levels = c("VL1", "VL2", "VL3", "VL4", "VL5",
                                                                    "L1", "L2", "L4", "M1", "M2", "M3",
                                                                    "M4", "M5", "M12", "H1", "VH1", "VH2", "VH3"))

#Import colour scheme
colour_scheme <- read_csv("./Analyses/colour_scheme.csv")

#Plot symbionts by lineage
sym_long_1 <- sym_long %>% filter(Lineage == "1")
sym_long_2 <- sym_long %>% filter(Lineage == "2")
sym_long_3 <- sym_long %>% filter(Lineage == "3")

ggplot(sym_long_1, aes(fill=Profile, y=Abundance, x=timeframe))+ 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~expedition) +
  theme(legend.position = "none") +
  scale_fill_manual(values=colour_scheme$hex)+
  facet_wrap(~site2)+
  ylim(c(0,1.2))


ggplot(sym_long_2, aes(fill=Profile, y=Abundance, x=timeframe))+ 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~expedition) +
  theme(legend.position = "none") +
  scale_fill_manual(values=colour_scheme$hex)+
  facet_wrap(~site2)+
  ylim(c(0,1.2))

ggplot(sym_long_3, aes(fill=Profile, y=Abundance, x=timeframe))+ 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~expedition) +
  theme(legend.position = "none") +
  scale_fill_manual(values=colour_scheme$hex)+
  facet_wrap(~site2)+
  ylim(c(0,1.2))



####Barplot of DIVs through time
DIV_switch <- Porites_DIV %>% subset_samples(coral_tag %in% c("585_1262", "1_1425", "230", "627", "1086", "831")) %>% otu_table() %>% data.frame()
DIV_switch <- DIV_switch[,colSums(DIV_switch) > 0]
DIV_switch$coral_tag <- sample_data(Porites_DIV %>% subset_samples(coral_tag %in% c("585_1262", "1_1425", "230", "627", "1086", "831")))$coral_tag
DIV_switch$expedition <- sample_data(Porites_DIV %>% subset_samples(coral_tag %in% c("585_1262", "1_1425", "230", "627", "1086", "831")))$expedition
DIV_switch$sample_name <- sample_data(Porites_DIV %>% subset_samples(coral_tag %in% c("585_1262", "1_1425", "230", "627", "1086", "831")))$sample_name
DIV_switch <- filter(DIV_switch, !sample_name %in% c("148_2017", "177_2015c", "82_2015a-post", "318_2017", "5_2017", "374_2015c"))
rownames(DIV_switch) <- DIV_switch$sample_name
DIV_switch$expedition <- gsub("2015a - Pre", "2015a", DIV_switch$expedition)
DIV_switch$expedition <- gsub("2015a - Post", "2015a", DIV_switch$expedition)


sym_long<- DIV_switch %>% gather(key="Profile", value="Abundance", 1:145, -coral_tag, -expedition, -sample_name)
sym_long$Profile <- factor(sym_long$Profile, ordered = TRUE)
#, levels = unique(colSums(DIV_switch[,1:119]))

#sym_long$expedition <- factor(sym_long$expedition, levels = c("2014", "2015a", "2015b", "2015c", "2016a", "2016b", "2017", "2018", "2019"), ordered = TRUE)

color = read_csv("./Data/Coral_tracking/ColorScheme_DIVs2.csv")

ggplot(sym_long, aes(fill=Profile, y=Abundance, x=expedition))+ 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~expedition) +
  facet_wrap(~coral_tag)+
  theme(legend.position = "none")+
  scale_fill_manual(values = color$color)+
  theme_classic()


#Determine sample sizes 
sample_data(Porites_Profiles)$coral_tag %>% unique()




