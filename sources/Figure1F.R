# scripts for making Figure 2 Panel F
# the source data table generated from the direct disorder results using script 08_meanDisorder.R
# SourceDataFig1E.csv can be retrieved from Source_Data_Fig1-4.xlsx
droplet_probabilities=read.csv("/sources/SourceDataFig1F.csv")

droplet_probabilities=droplet_probabilities[-(which(is.na(droplet_probabilities$mean))),]

# number of samples
sum(droplet_probabilities$region=="Nterminal-Extension")
sum(droplet_probabilities$region=="Cterminal-Extension")
sum(droplet_probabilities$region=="Internal")

# N-terminal extension n=809
# C-terminal extension n=335
# Conserved internal region n=809

library(ggplot2)
library(ggsignif)

droplet_probabilities$region=factor(droplet_probabilities$region, levels=c("Nterminal-Extension", "Internal","Cterminal-Extension"))

ggplot(droplet_probabilities, aes(x=region,y=as.numeric(mean), fill=region))+
  geom_violin(width=1.2) +
  geom_boxplot(width=0.05, color="gray20", alpha=0.3) +
  scale_fill_manual(values=c("#A7B7D1FF","#D5D2D1FF","#F8B58BFF"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none", axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x="", y="Mean Disorder")+
  geom_signif(data=droplet_probabilities,comparisons=list(c("Nterminal-Extension","Cterminal-Extension")),map_signif_level = TRUE,annotations="*",
              y_position = 1.0, tip_length = c(0.01, 0.01))+
  geom_signif(data=droplet_probabilities,comparisons=list(c("Nterminal-Extension","Internal")),map_signif_level = TRUE,annotations="*",
              y_position = 0.95, tip_length = c(0.01, 0.01))+
  geom_signif(data=droplet_probabilities,comparisons=list(c("Cterminal-Extension","Internal")),map_signif_level = TRUE,annotations="*",
              y_position = 0.95, tip_length = c(0.01, 0.01))+
  scale_x_discrete(labels=c("N-terminal Extension", "Internal Region", "C-terminal Extension"))+
  ylim(0,1)


# statistical test between regions
t.test(as.numeric(droplet_probabilities[droplet_probabilities$region=="Nterminal-Extension",2]),as.numeric(droplet_probabilities[droplet_probabilities$region=="Cterminal-Extension",2]))
t.test(as.numeric(droplet_probabilities[droplet_probabilities$region=="Nterminal-Extension",2]),as.numeric(droplet_probabilities[droplet_probabilities$region=="Internal",2]))
t.test(as.numeric(droplet_probabilities[droplet_probabilities$region=="Cterminal-Extension",2]),as.numeric(droplet_probabilities[droplet_probabilities$region=="Internal",2]))

