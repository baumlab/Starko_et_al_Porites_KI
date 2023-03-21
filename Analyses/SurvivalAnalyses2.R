#This script loads in the Porites SymPortal data, combines with coral tag-level metadata and produces some key but basic exploratory plots and statistics focused on survival and analysis of raw ITS2 profiles

#Load in phyloseq objects
load("./Data/Symbionts/Porites_Profiles2.RData")
Porites_Profiles <- Porites_Profiles2
load("./Data/Symbionts/Porites_DIV.RData")

#Import coral tracking data
track <- read_excel("./Data/Coral_tracking/Workingsheet_all_fixed_version3_7.xlsx")

#Extract data from each lineage
Lineage1 <- filter(track, Clade_dom == "Red")
Lineage2 <- filter(track, Clade_dom == "Blue")
Lineage3 <- filter(track, Clade_dom == "Pale")
LineageUnk <- filter(track, Clade_dom == "Unknown")

#Create column outlining the different lineages
Lineage1$Lineage <- "1"
Lineage2$Lineage <- "2"
Lineage3$Lineage <- "3"
LineageUnk$Lineage <- "Unknown"

#Recombine dataset with new column indicating lineage
AllLins<- rbind(Lineage1, Lineage2, Lineage3, LineageUnk)
AllLins2<- rbind(Lineage1, Lineage2, Lineage3)


#####The following code are for testing for predictors of survival (binary; alive or dead)

#Filter dataset to include only colonies of known mortality status
Survival <- AllLins2 %>% filter(mortality!="NA")
Survival2 <- AllLins %>% filter(mortality!="NA")

#Create contingency table of mortality and lineage
mort_table <- table(Survival$Lineage, Survival$mortality)

#Do the same but only for colonies that were sampled with RADSeq and survival is known
Survival_RAD <- Survival %>% filter(Method == "RADSeq")
mort_table_RAD <- table(Survival_RAD$Lineage, Survival_RAD$mortality)

#Simple visualizations of mortality
balloonplot(t(mort_table), main ="Mortality", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE)
mosaicplot(mort_table, las = 2, main = "Mortality", col = c("white", "black"))

#Chi Squared test of association between mortality and lineage
chisq <- chisq.test(mort_table)
print(chisq)

#Chi Squared test of association between mortality and lineage but for colonies sequenced with RAD only
chisq_RAD <- chisq.test(mort_table_RAD) #RADSeq assigned samples only
print(chisq_RAD)

##Test and plot survival vs disturbance and lineage
#All samples
Survival$mortality_num <- Survival$mortality %>% as.numeric()
Survival$Disturbance_sqrt <- Survival$Disturbance_sqrt %>% as.numeric()

#Plot of survival versus disturbance
ggplot(Survival, aes(x = Disturbance_sqrt, y = mortality_num, color = Lineage, group = Lineage))+
  geom_jitter(cex = 3, width=3, height = 0.06)+
  theme_classic()+
  stat_smooth(formula = 'y~x', method = "bayesglm", method.args = list(family = binomial), se = TRUE, level = 0.95, aes(fill=Lineage))+
  scale_color_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  scale_fill_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  ylab("Proportion dead")+
  xlab("Human disturbance")+
  ggtitle("All samples")+
  ylim(c(-0.106,1.06))


#Binomial GLM testing for effect of disturbance and lineage on survival
z <- glm(mortality_num ~ Disturbance_sqrt*Lineage, data = Survival, family = "binomial")
summary(z)
Anova(z)

#####RAD samples only

#Plot of survival versus disturbance
ggplot(Survival_RAD, aes(x = Disturbance_sqrt, y = as.numeric(mortality), color = Lineage, group = Lineage))+
  geom_jitter(cex = 3, width=3, height = 0.06)+
  theme_classic()+
  stat_smooth(formula = 'y~x', method = "bayesglm", method.args = list(family = binomial), se = TRUE, level = 0.95, aes(fill=Lineage))+
  scale_color_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  scale_fill_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  ylab("Proportion dead")+
  xlab("Human disturbance")+
  ggtitle("RAD samples")+
  ylim(c(-0.106,1.06))

#Plot of survival versus disturbance
surv_rad <- ggplot(Survival_RAD, aes(x = Disturbance_sqrt, y = as.numeric(mortality), color = Lineage, group = Lineage))+
  geom_jitter(cex = 2, width=3, height = 0.06)+
  theme_classic()+
  stat_smooth(formula = 'y~x', method = "bayesglm", method.args = list(family = binomial), se = TRUE, level = 0.95, aes(fill=Lineage))+
  scale_color_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  scale_fill_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  ylab("Proportion dead")+
  xlab("Human disturbance")+
  ggtitle("a")+
  ylim(c(-0.106,1.06))

#Binomial GLM testing for effect of disturbance and lineage on survival
z <- glm(as.numeric(mortality) ~Disturbance_sqrt*Lineage, data = Survival_RAD, family = "binomial")
summary(z)
Anova(z)


########The following code is specifically for mortality as a continuous variable (percent tissue mortality)

#Import data
Partial <- read_excel("./Data/Coral_tracking/PartialMortality_Dec2022_2.xlsx")

#Remove colonies of unknown lineage (i.e., unassigned)
Partial <- Partial %>% filter(Lineage != "Unknown")

#Quasibinomial GLM testing for effect of disturbance and lineage on tissue mortality
z <- glm(PropDead ~ Disturbance_sqrt*Lineage, data = Partial, family = "quasibinomial")
summary(z)
Anova(z)

#Same but only colonies sequenced with RAD
Partial_RAD <- Partial %>% filter(Method=="RAD")
z <- glm(PropDead ~ Disturbance_sqrt*Lineage, data = Partial_RAD, family = "quasibinomial")
summary(z)
Anova(z)

#Plot tissue mortality versus lineage and disturbance
ggplot(Partial, aes(x = Disturbance_sqrt, y = PropDead, color =Lineage, group = Lineage))+
  geom_point()+
  theme_classic()+
  stat_smooth(formula = 'y~x', method = "bayesglm", method.args = list(family = quasibinomial), se = TRUE, level = 0.95, aes(fill=Lineage))+
  scale_color_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  scale_fill_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  ylab("Proportion dead")+
  xlab("Human disturbance")+
  ggtitle("All samples")+
  ylim(c(-0.104,1.04))

#Same but RAD samples only
ggplot(Partial_RAD, aes(x = Disturbance_sqrt, y = PropDead, color =Lineage, group = Lineage))+
  geom_point()+
  theme_classic()+
  stat_smooth(formula = 'y~x', method = "bayesglm", method.args = list(family = quasibinomial), se = TRUE, level = 0.95, aes(fill=Lineage))+
  scale_color_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  scale_fill_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  ylab("Proportion dead")+
  xlab("Human disturbance")+
  ggtitle("All samples")+
  ylim(c(-0.104,1.04))

part_RAD<- ggplot(Partial_RAD, aes(x = Disturbance_sqrt, y = PropDead, color =Lineage, group = Lineage))+
  geom_jitter(width=3, height = 0.03, cex = 2)+
  theme_classic()+
  stat_smooth(formula = 'y~x', method = "bayesglm", method.args = list(family = quasibinomial), se = TRUE, level = 0.95, aes(fill=Lineage))+
  scale_color_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  scale_fill_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  ylab("Proportion dead")+
  xlab("Human disturbance")+
  ggtitle("b")+
  ylim(c(-0.106,1.06))

part <- ggplot(Partial, aes(x = Disturbance_sqrt, y = PropDead, color =Lineage, group = Lineage))+
  geom_jitter(width=3, height = 0.03)+
  theme_classic()+
  stat_smooth(formula = 'y~x', method = "bayesglm", method.args = list(family = quasibinomial), se = TRUE, level = 0.95, aes(fill=Lineage))+
  scale_color_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  scale_fill_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  ylab("Proportion tissue death")+
  xlab("Human disturbance")+
  ggtitle("b")+
  ylim(c(-0.104,1.04))
  xlim(c(0,100))

surv <- ggplot(Survival, aes(x = Disturbance_sqrt, y = mortality_num, color = Lineage, group = Lineage))+
  geom_jitter(width=3, height = 0.03)+
  theme_classic()+
  stat_smooth(formula = 'y~x', method = "bayesglm", method.args = list(family = binomial), se = TRUE, level = 0.95, aes(fill=Lineage))+
  scale_color_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  scale_fill_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  ylab("Proportion total mortality")+
  xlab("Human disturbance")+
  ggtitle("a")+
  ylim(c(-0.104,1.04))
  xlim(c(0,110))

#Plot combined plot of binary and continuous survival metrics
plot_grid(surv, part)

######Contingency analysis from before heatwave (testing for lineage versus region on KI)
Before_col <- read_excel("./Data/Coral_tracking/Workingsheet_all_fixed_version3_7.xlsx", sheet = "Contingency_before")
Before_col <- filter(Before_col, Clade_dom!="Unknown")
Contingency_Table <- table( Before_col$Clade_dom, Before_col$region)

bf = contingencyTableBF(Contingency_Table, sampleType = "indepMulti", fixedMargin = "cols")
bf

fisher.test(Contingency_Table, simulate.p.value = TRUE)

#Now pairwise comparison of each combination
#Remember to correct p-values for multiple comparisons
Contingency_Table[,c(1,2)] %>% fisher.test(simulate.p.value = TRUE) # p = 0.26388
Contingency_Table[,c(1,3)] %>% fisher.test(simulate.p.value = TRUE) # p = 0.00899** (Bay of Wrecks vs North Shore)
Contingency_Table[,c(1,4)] %>% fisher.test(simulate.p.value = TRUE) # p = 0.00300** (Bay of Wrecks vs Vaskess/South lagoon)
Contingency_Table[,c(2,3)] %>% fisher.test(simulate.p.value = TRUE) # p = 0.64140
Contingency_Table[,c(2,4)] %>% fisher.test(simulate.p.value = TRUE) # p = 0.08394
Contingency_Table[,c(3,4)] %>% fisher.test(simulate.p.value = TRUE) # p > 0.99

#Multinomial regression for effect of disturbance
#No significant effect
Before_col$Disturbance_sqrt
Before_col$Clade_dom
model <- multinom(Clade_dom~Disturbance_sqrt, data = Before_col)
summary(model)
Anova(model)

##Test for impact of colony size (measured as "area")
Partial$full <- ifelse(Partial$PropDead==1, 1, 0)
z <- glm(PropDead ~ Size_before*Lineage, data = Partial, family = "quasibinomial")
Anova(z)

#Survival vs. size
z <- glm(full ~ Size_before*Lineage, data = Partial, family = "binomial")
Anova(z, type = "II")

#Size versus lineage
z <- lm(Size_before ~ Lineage, data = Partial)
Anova(z)

#Percentage tissue mortality versus colony size
g1 <- ggplot(Partial, aes(y = PropDead, x = Size_before, fill = Lineage, col = Lineage, group = Lineage))+
  geom_point(cex =3 )+
  stat_smooth(formula = 'y~x', method = "bayesglm", method.args = list(family = quasibinomial), se = FALSE)+
  scale_color_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  scale_fill_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))
  
#Survival (binary) versus colony size
g2 <- ggplot(Partial, aes(y = full, x = Size_before, fill = Lineage, col = Lineage, group = Lineage))+
  geom_point(cex =3 )+
  stat_smooth(formula = 'y~x', method = "bayesglm", method.args = list(family = binomial), se = FALSE)+
  scale_color_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))+
  scale_fill_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))

#Size versus lineage
g3 <- ggplot(Partial, aes(y = Size_before, x = Lineage, fill = Lineage))+
  geom_boxplot()+
  scale_fill_manual(values = c(visTree::makeTransparent("#DFBE99", 0.8), visTree::makeTransparent("#729EA1", 0.8), visTree::makeTransparent("#DB5375", 0.8)))

#Plot 3 panel figure of size results
plot_grid(g3, g1, g2, ncol = 1)

#Two panel figure of total and tissue mortality for all colonies of known survival status
plot_grid(surv, part)

#Two panel figure of total and tissue mortality for RAD-sampled colonies of known survival status
plot_grid(surv_rad, part_RAD)


