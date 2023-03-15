#Script for creating heatmap of tagged corals through time
#This script makes a heatmap showing the corals of known mortality status through time

hm <- read_excel("./Data/Coral_tracking/Heatmap_Jan2023_3.xlsx", sheet = "heatmap_final3")
hm.d <- gather(hm, key = "Expedition", value = "Symbiont", 7:15, factor_key = TRUE)

hm.d$Symbiont <- hm.d$Symbiont %>% as.numeric() %>% as.factor()
hex <- hm$Color[1:35]

ggplot(hm.d, aes(Expedition, rev(order), fill= Symbiont)) +
  geom_tile() +
  scale_fill_manual(values = hex)+
  labs(x = "Time point", y = "Colony ID", fill = "Rank") +scale_y_continuous(position = "right")

Lin1<- filter(hm.d,Clade_dom=="Red")
Lin1$Symbiont <- Lin1$Symbiont %>% as.character() %>% as.factor()
hex1 <- c(hm$Color[1:35][hm$ID[1:35] %in% Lin1$Symbiont])

g1 <- ggplot(filter(hm.d,Clade_dom=="Red") , aes(Expedition, rev(order), fill= Symbiont)) +
  geom_tile() +
  scale_fill_manual(values = hex1)+
  labs(x = "Time point", y = "Colony ID", fill = "Rank") +scale_y_continuous(position = "right")+
  ggtitle("Lineage 1")+
  theme_cowplot()

Lin2<- filter(hm.d,Clade_dom=="Blue")
Lin2$Symbiont <- Lin2$Symbiont %>% as.character() %>% as.factor()
hex2 <- c(hm$Color[1:35][hm$ID[1:35] %in% Lin2$Symbiont])

g2 <- ggplot(filter(hm.d,Clade_dom=="Blue") , aes(Expedition, rev(order), fill= Symbiont)) +
  geom_tile() +
  scale_fill_manual(values = hex2)+
  labs(x = "Time point", y = "Colony ID", fill = "Rank") +scale_y_continuous(position = "right")+
  ggtitle("Lineage 2")+
  theme_cowplot()

Lin3<- filter(hm.d,Clade_dom=="Pale")
Lin3$Symbiont <- Lin3$Symbiont %>% as.character() %>% as.factor()
hex3 <- c(hm$Color[1:35][hm$ID[1:35] %in% Lin3$Symbiont])


g3 <- ggplot(filter(hm.d,Clade_dom=="Pale") , aes(Expedition, rev(order), fill= Symbiont)) +
  geom_tile() +
  scale_fill_manual(values = hex3)+
  labs(x = "Time point", y = "Colony ID", fill = "Rank") +scale_y_continuous(position = "right")+
  ggtitle("Lineage 3")+
  theme_cowplot()

g1
g2
g3

library(cowplot)
plot_grid(g1, g2, g3)


###Now ALL colonies through time (including samples with as few as 200 sequence reads)
hm <- read_excel("./Data/Coral_tracking/Heatmap_allsamples_Jan25_4.xlsx", sheet = "heatmap")
hm.d <- gather(hm, key = "Expedition", value = "Symbiont", 9:17, factor_key = TRUE)

hm.d$Symbiont <- hm.d$Symbiont %>% as.numeric() %>% as.factor()
hex <- hm$Color[1:59]

ggplot(hm.d, aes(Expedition, rev(order), fill= Symbiont)) +
  geom_tile() +
  scale_fill_manual(values = hex)+
  labs(x = "Time point", y = "Colony ID", fill = "Rank") +scale_y_continuous(position = "right")

Lin1<- filter(hm.d,Clade_dom=="Red")
Lin1$Symbiont <- Lin1$Symbiont %>% as.character() %>% as.factor()
hex1 <- c(hm$Color[1:59][hm$ID[1:59] %in% Lin1$Symbiont])

g1 <- ggplot(filter(hm.d,Clade_dom=="Red") , aes(Expedition, rev(order), fill= Symbiont)) +
  geom_tile() +
  scale_fill_manual(values = hex1)+
  labs(x = "Time point", y = "Colony ID", fill = "Rank") +scale_y_continuous(position = "right")+
  ggtitle("Lineage 1")+
  theme_cowplot()

Lin2<- filter(hm.d,Clade_dom=="Blue")
Lin2$Symbiont <- Lin2$Symbiont %>% as.character() %>% as.factor()
hex2 <- c(hm$Color[1:59][hm$ID[1:59] %in% Lin2$Symbiont])

g2 <- ggplot(filter(hm.d,Clade_dom=="Blue") , aes(Expedition, rev(order), fill= Symbiont)) +
  geom_tile() +
  scale_fill_manual(values = hex2)+
  labs(x = "Time point", y = "Colony ID", fill = "Rank") +scale_y_continuous(position = "right")+
  ggtitle("Lineage 2")+
  theme_cowplot()

Lin3<- filter(hm.d,Clade_dom=="Pale")
Lin3$Symbiont <- Lin3$Symbiont %>% as.character() %>% as.factor()
hex3 <- c(hm$Color[1:59][hm$ID[1:59] %in% Lin3$Symbiont])

g3 <- ggplot(filter(hm.d,Clade_dom=="Pale") , aes(Expedition, rev(order), fill= Symbiont)) +
  geom_tile() +
  scale_fill_manual(values = hex3)+
  labs(x = "Time point", y = "Colony ID", fill = "Rank") +scale_y_continuous(position = "right")+
  ggtitle("Lineage 3")+
  theme_cowplot()


Lin4<- filter(hm.d,Clade_dom=="Unknown")
Lin4$Symbiont <- Lin4$Symbiont %>% as.character() %>% as.factor()
hex3 <- c(hm$Color[1:59][hm$ID[1:59] %in% Lin4$Symbiont])

g4 <- ggplot(filter(hm.d,Clade_dom=="Unknown") , aes(Expedition, rev(order), fill= Symbiont)) +
  geom_tile() +
  scale_fill_manual(values = hex3)+
  labs(x = "Time point", y = "Colony ID", fill = "Rank") +scale_y_continuous(position = "right")+
  ggtitle("Unassigned")+
  theme_cowplot()

g1
g2
g3
g4

plot_grid(g1, g2, g3, g4)


