###This code is for comparing categories of bleaching between lineages and across time. Analysis is done using two different category schemes (see methods)

####McClanahan categories

#Import data
bleach_2016a <- read_excel("./Data/Coral_tracking/Mortality_and_bleaching.xlsx", sheet = "2016a")

#Remove colonies of unknown lineage
bleach_2016a <- filter(bleach_2016a, Lineage!="NA") %>% filter(Lineage!="Unknown")

#Contingency table of lineage and bleaching categories
sum2016_1<- table(bleach_2016a$Lineage, bleach_2016a$Bleaching_McClanahan5) %>% prop_row() %>% as.data.frame()

#Make lineage a factor (ordered just for the plot)
sum2016_1$Var1 <- sum2016_1$Var1 %>% factor(ordered = TRUE, levels =  c("Red","Blue",  "Pale"))

#Make "dead" colonies 6 (instead of 7)
sum2016_1$Var2 <- gsub(pattern = 7, replacement = 6, x = sum2016_1$Var2)

#Make bleaching an ordered factor
sum2016_1$Var2 <- sum2016_1$Var2 %>% factor(ordered = TRUE, levels = c(1:6))

#Plot bleaching versus lineage
plot16_1<-ggplot(sum2016_1, aes(x = Var1, y = Freq, fill = Var2, color = Var2))+
  geom_col()+
  scale_fill_manual(values = c("White", "#f0f0f0", "lightgrey", "darkgrey", "#636363", "black"))+theme_cowplot()+
  scale_color_manual(values = c("black","black", "black", "black", "black", "black" ))

#Import 2015c data
bleach_2015c <- read_excel("./Data/Coral_tracking/Mortality_and_bleaching.xlsx", sheet = "2015c")

#Remove colonies that are not assigned to lineage
bleach_2015c <- bleach_2015c %>% filter(Lineage!="NA") %>% filter(Lineage!="Unknown")

#Contingency table of lineage and bleaching categories
sum2015_1<- table(bleach_2015c$Lineage, bleach_2015c$Bleaching_McClanahan5) %>% prop_row() %>% as.data.frame()

#Make variables ordered factors
sum2015_1$Var1 <- sum2015_1$Var1 %>% factor(ordered = TRUE, levels = c("Red","Blue",  "Pale"))
sum2015_1$Var2 <- sum2015_1$Var2 %>% factor(ordered = TRUE, levels = c(1:7))

#Plot bleaching versus lineage in 2015c
plot15_1<- ggplot(sum2015_1, aes(x = Var1, y = Freq, fill = Var2, color = Var2))+
  geom_col()+
  scale_fill_manual(values = c("White", "#f0f0f0", "lightgrey", "darkgrey", "#636363", "#5A5A5A", "black"))+theme_cowplot()+
  scale_color_manual(values = c("black","black", "black", "black", "black", "black", "black" ))

#Remove colonies of unknown lineage
bleach_2016a_2 <- bleach_2016a %>% filter(Lineage!="NA") %>% filter(Lineage!="Unknown")

#Ordinal regression of bleaching versus lineage in 2016a
vglm(factor(Bleaching_McClanahan5, ordered = TRUE) ~ Lineage, data = bleach_2016a_2, family = propodds) %>% anova()

#Ordinal regression of bleaching versus lineage in 2015c
vglm(factor(bleach_2015c$Bleaching_McClanahan5, ordered = TRUE) ~ Lineage, data = bleach_2015c, family = propodds) %>% anova()

#Plot 2015c and 2016a bleaching data
plot_grid(plot15_1, plot16_1)


#####Baum categories
#Import data
bleach_2016a <- read_excel("./Data/Coral_tracking/Mortality_and_bleaching.xlsx", sheet = "2016a")

#Remove colonies of unknown lineage
bleach_2016a <- filter(bleach_2016a, Lineage!="NA") %>% filter(Lineage!="Unknown")

#Contingency table of lineage versus bleaching status
sum2016_1<- table(bleach_2016a$Lineage, bleach_2016a$Bleaching_KT_adj) %>% prop_row() %>% as.data.frame()

#Convert variables to ordered factors
sum2016_1$Var1 <- sum2016_1$Var1 %>% factor(ordered = TRUE, levels =  c("Red","Blue",  "Pale"))
sum2016_1$Var2 <- gsub(4,3, sum2016_1$Var2)
sum2016_1$Var2 <- sum2016_1$Var2 %>% factor(ordered = TRUE, levels = c(0:4))

#Plot 2016 data
plot16<-ggplot(sum2016_1, aes(x = Var1, y = Freq, fill = Var2, color = Var2))+
  geom_col()+
  scale_fill_manual(values = c("White", "#f0f0f0", "lightgrey", "black"))+theme_cowplot()+
  scale_color_manual(values = c("black","black", "black", "black" ))

#Import 2015c data
bleach_2015c <- read_excel("./Data/Coral_tracking/Mortality_and_bleaching.xlsx", sheet = "2015c")
#Remove colonies of unknown lineage
bleach_2015c <- bleach_2015c %>% filter(Lineage!="NA") %>% filter(Lineage!="Unknown")

1#Contingency table of bleaching and lineage in 2015c
sum2015_1<- table(bleach_2015c$Lineage, bleach_2015c$Bleaching_KT_adj) %>% prop_row() %>% as.data.frame()

#Make variables ordered factors
sum2015_1$Var1 <- sum2015_1$Var1 %>% factor(ordered = TRUE, levels = c("Red","Blue",  "Pale"))
sum2015_1$Var2 <- gsub(4,3, sum2015_1$Var2)
sum2015_1$Var2 <- sum2015_1$Var2 %>% factor(ordered = TRUE, levels = c(0:4))

#Plot 2015c data
plot15<- ggplot(sum2015_1, aes(x = Var1, y = Freq, fill = Var2, color = Var2))+
  geom_col()+
  scale_fill_manual(values = c("White", "#f0f0f0", "lightgrey", "black"))+theme_cowplot()+
  scale_color_manual(values = c("black","black", "black", "black"))

#Remove unknown lineages from 2016a data
bleach_2016a_2 <- bleach_2016a %>% filter(Lineage!="NA") %>% filter(Lineage!="Unknown")

#Ordinal regression comparing bleaching versus lineage in 2016a
vglm(factor(Bleaching_KT_adj, ordered = TRUE) ~ Lineage, data = bleach_2016a_2, family = propodds) %>% anova()

#Ordinal regression comparing bleaching versus lineage in 2015c
vglm(factor(bleach_2015c$Bleaching_KT_adj, ordered = TRUE) ~ Lineage, data = bleach_2015c, family = propodds) %>% anova() 

#Plot grid
plot_grid(plot15, plot16)

#Plot grid
plot_grid(plot15_1, plot16_1,plot15, plot16)

