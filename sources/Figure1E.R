# scripts for making Figure 2 Panel E
# the source data table generated from the direct disorder results using script 08_meanDisorder.R
# SourceDataFig1E.csv can be retrieved from Source_Data_Fig1-4.xlsx
disorder_probabilities=read.csv("/sources/SourceDataFig1E.csv")

# number of samples
sum(disorder_probabilities$region=="Nterminal-Extension")
sum(disorder_probabilities$region=="Cterminal-Extension")
sum(disorder_probabilities$region=="Internal")

# N-terminal extension n=807
# C-terminal extension n=333
# Conserved internal region n=807

library(ggplot2)
library(ggsignif)

disorder_probabilities$region=factor(disorder_probabilities$region, levels=c("Nterminal-Extension", "Internal","Cterminal-Extension"))

ggplot(disorder_probabilities, aes(x=region,y=as.numeric(mean), fill=region))+
  geom_violin(width=1.2) +
  geom_boxplot(width=0.05, color="gray20", alpha=0.3) +
  scale_fill_manual(values=c("#A7B7D1FF","#D5D2D1FF","#F8B58BFF"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none", axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x="", y="Mean Disorder")+
  geom_signif(data=disorder_probabilities,comparisons=list(c("Nterminal-Extension","Cterminal-Extension")),map_signif_level = TRUE,annotations="*",
              y_position = 0.90, tip_length = c(0.01, 0.01))+
  geom_signif(data=disorder_probabilities,comparisons=list(c("Nterminal-Extension","Internal")),map_signif_level = TRUE,annotations="*",
              y_position = 0.85, tip_length = c(0.01, 0.01))+
  geom_signif(data=disorder_probabilities,comparisons=list(c("Cterminal-Extension","Internal")),map_signif_level = TRUE,annotations="*",
              y_position = 0.85, tip_length = c(0.01, 0.01))+
  scale_x_discrete(labels=c("N-terminal Extension", "Internal Region", "C-terminal Extension"))+
  ylim(0,1)


# statistical test between regions
t.test(as.numeric(disorder_probabilities[disorder_probabilities$region=="Nterminal-Extension",2]),as.numeric(disorder_probabilities[disorder_probabilities$region=="Cterminal-Extension",2]))
t.test(as.numeric(disorder_probabilities[disorder_probabilities$region=="Nterminal-Extension",2]),as.numeric(disorder_probabilities[disorder_probabilities$region=="Internal",2]))
t.test(as.numeric(disorder_probabilities[disorder_probabilities$region=="Cterminal-Extension",2]),as.numeric(disorder_probabilities[disorder_probabilities$region=="Internal",2]))

