#Bleaching comparison

#Import data
bleach_2016a <- read_excel("./Data/Coral_tracking/Mortality_and_bleaching.xlsx", sheet = "2016a")

bleach_2016a <- filter(bleach_2016a, Lineage!="NA") %>% filter(Lineage!="Unknown")

sum2016_1<- table(bleach_2016a$Lineage, bleach_2016a$bleach2) %>% prop_row() %>% as.data.frame()

sum2016_1$Var1 <- sum2016_1$Var1 %>% factor(ordered = TRUE, levels =  c("Red","Blue",  "Pale"))

plot16<-ggplot(sum2016_1, aes(x = Var1, y = Freq, fill = Var2, color = Var2))+
  geom_col()+
  scale_fill_manual(values = c("White", "lightgrey", "darkgrey", "black"))+theme_cowplot()+
  scale_color_manual(values = c("black", "black","black","black"))

bleach_2015c <- read_excel("./Data/Coral_tracking/Mortality_and_bleaching.xlsx", sheet = "2015c")

bleach_2015c <- bleach_2015c %>% filter(Lineage!="NA") %>% filter(Lineage!="Unknown")

sum2015_1<- table(bleach_2015c$Lineage, bleach_2015c$bleaching2) %>% prop_row() %>% as.data.frame()

sum2015_1$Var1 <- sum2015_1$Var1 %>% factor(ordered = TRUE, levels = c("Red","Blue",  "Pale"))

plot15<- ggplot(sum2015_1, aes(x = Var1, y = Freq, fill = Var2, color = Var2))+
  geom_col()+
  scale_fill_manual(values = c("White", "lightgrey", "darkgrey", "black"))+theme_cowplot()+
  scale_color_manual(values = c("black", "black","black","black"))


bleach_2016a_2 <- bleach_2016a %>% filter(Lineage!="NA") %>% filter(Lineage!="Unknown")

vglm(factor(bleach2, ordered = TRUE) ~ Lineage, data = bleach_2016a_2, family = propodds) %>% anova()

vglm(factor(bleaching2, ordered = TRUE) ~ Lineage, data = bleach_2015c, family = propodds) %>% anova() 

plot_grid(plot15, plot16)

##McClanahan categories

#Import data
bleach_2016a <- read_excel("./Data/Coral_tracking/Mortality_and_bleaching.xlsx", sheet = "2016a")

bleach_2016a <- filter(bleach_2016a, Lineage!="NA") %>% filter(Lineage!="Unknown")

sum2016_1<- table(bleach_2016a$Lineage, bleach_2016a$Bleaching_McClanahan5) %>% prop_row() %>% as.data.frame()

sum2016_1$Var1 <- sum2016_1$Var1 %>% factor(ordered = TRUE, levels =  c("Red","Blue",  "Pale"))
sum2016_1$Var2 <- gsub(pattern = 7, replacement = 6, x = sum2016_1$Var2)
sum2016_1$Var2 <- sum2016_1$Var2 %>% factor(ordered = TRUE, levels = c(1:6))

plot16_1<-ggplot(sum2016_1, aes(x = Var1, y = Freq, fill = Var2, color = Var2))+
  geom_col()+
  scale_fill_manual(values = c("White", "#f0f0f0", "lightgrey", "darkgrey", "#636363", "black"))+theme_cowplot()+
  scale_color_manual(values = c("black","black", "black", "black", "black", "black" ))

  
bleach_2015c <- read_excel("./Data/Coral_tracking/Mortality_and_bleaching.xlsx", sheet = "2015c")

bleach_2015c <- bleach_2015c %>% filter(Lineage!="NA") %>% filter(Lineage!="Unknown")

sum2015_1<- table(bleach_2015c$Lineage, bleach_2015c$Bleaching_McClanahan5) %>% prop_row() %>% as.data.frame()

sum2015_1$Var1 <- sum2015_1$Var1 %>% factor(ordered = TRUE, levels = c("Red","Blue",  "Pale"))
sum2015_1$Var2 <- sum2015_1$Var2 %>% factor(ordered = TRUE, levels = c(1:7))


plot15_1<- ggplot(sum2015_1, aes(x = Var1, y = Freq, fill = Var2, color = Var2))+
  geom_col()+
  scale_fill_manual(values = c("White", "#f0f0f0", "lightgrey", "darkgrey", "#636363", "#5A5A5A", "black"))+theme_cowplot()+
  scale_color_manual(values = c("black","black", "black", "black", "black", "black", "black" ))


bleach_2016a_2 <- bleach_2016a %>% filter(Lineage!="NA") %>% filter(Lineage!="Unknown")

vglm(factor(Bleaching_McClanahan5, ordered = TRUE) ~ Lineage, data = bleach_2016a_2, family = propodds) %>% anova()

vglm(factor(bleach_2015c$Bleaching_McClanahan5, ordered = TRUE) ~ Lineage, data = bleach_2015c, family = propodds) %>% anova()

plot_grid(plot15_1, plot16_1)


##kt categories

#Import data
bleach_2016a <- read_excel("./Data/Coral_tracking/Mortality_and_bleaching.xlsx", sheet = "2016a")

bleach_2016a <- filter(bleach_2016a, Lineage!="NA") %>% filter(Lineage!="Unknown")

sum2016_1<- table(bleach_2016a$Lineage, bleach_2016a$Bleaching_KTadj) %>% prop_row() %>% as.data.frame()

sum2016_1$Var1 <- sum2016_1$Var1 %>% factor(ordered = TRUE, levels =  c("Red","Blue",  "Pale"))
sum2016_1$Var2 <- gsub(4,3, sum2016_1$Var2)
sum2016_1$Var2 <- sum2016_1$Var2 %>% factor(ordered = TRUE, levels = c(0:4))

plot16<-ggplot(sum2016_1, aes(x = Var1, y = Freq, fill = Var2, color = Var2))+
  geom_col()+
  scale_fill_manual(values = c("White", "#f0f0f0", "lightgrey", "black"))+theme_cowplot()+
  scale_color_manual(values = c("black","black", "black", "black" ))


bleach_2015c <- read_excel("./Data/Coral_tracking/Mortality_and_bleaching.xlsx", sheet = "2015c")

bleach_2015c <- bleach_2015c %>% filter(Lineage!="NA") %>% filter(Lineage!="Unknown")

sum2015_1<- table(bleach_2015c$Lineage, bleach_2015c$Bleaching_KT_adj) %>% prop_row() %>% as.data.frame()

sum2015_1$Var1 <- sum2015_1$Var1 %>% factor(ordered = TRUE, levels = c("Red","Blue",  "Pale"))
sum2015_1$Var2 <- gsub(4,3, sum2015_1$Var2)
sum2015_1$Var2 <- sum2015_1$Var2 %>% factor(ordered = TRUE, levels = c(0:4))


plot15<- ggplot(sum2015_1, aes(x = Var1, y = Freq, fill = Var2, color = Var2))+
  geom_col()+
  scale_fill_manual(values = c("White", "#f0f0f0", "lightgrey", "black"))+theme_cowplot()+
  scale_color_manual(values = c("black","black", "black", "black"))


bleach_2016a_2 <- bleach_2016a %>% filter(Lineage!="NA") %>% filter(Lineage!="Unknown")

vglm(factor(Bleaching_KTadj, ordered = TRUE) ~ Lineage, data = bleach_2016a_2, family = propodds) %>% anova()

vglm(factor(bleach_2015c$Bleaching_KT_adj, ordered = TRUE) ~ Lineage, data = bleach_2015c, family = propodds) %>% anova() 

plot_grid(plot15, plot16)

plot_grid(plot15_1, plot16_1,plot15, plot16)
