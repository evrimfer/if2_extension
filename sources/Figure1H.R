# scripts for making Figure 2 Panel H
# the source data table generated from the direct disorder results using script 07_preferences.R
# SourceDataFig1E.csv can be retrieved from Source_Data_Fig1-4.xlsx

oxygen_preferences=read.csv("/sources/SourceDataFig1H.csv")

# number of samples
sum(oxygen_preferences$tolerance=="aerobic")
sum(oxygen_preferences$tolerance=="anaerobic")
sum(oxygen_preferences$tolerance=="facultative")

# aerobic n=432
# anaerobic n=236
# facultative n=127


library(ggplot2)

oxygen_preferences$tolerance=factor(oxygen_preferences$tolerance, levels=c("aerobic","facultative","anaerobic"))


ggplot(oxygen_preferences, aes(x=tolerance, y=as.numeric(L_IF2N), fill=source))+
  geom_boxplot(width=0.5, alpha=0.6) +
  geom_jitter(width=0, alpha=0.4)+
  scale_fill_manual(values=c("#fdce9d"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",
        axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x="", y="IF2 N-terminal Length")

ggplot(oxygen_preferences, aes(x=tolerance, y=as.numeric(L_IF2C), fill=source))+
  geom_boxplot(width=0.5, alpha=0.6) +
  geom_jitter(width=0, alpha=0.4)+
  scale_fill_manual(values=c("#fdce9d"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",
        axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x="", y="IF2 C-terminal Length")

# statistical test between preferences
#NTD
t.test(as.numeric(oxygen_preferences[oxygen_preferences$tolerance=="aerobic",3]),as.numeric(oxygen_preferences[oxygen_preferences$tolerance=="facultative",3]))
t.test(as.numeric(oxygen_preferences[oxygen_preferences$tolerance=="aerobic",3]),as.numeric(oxygen_preferences[oxygen_preferences$tolerance=="anaerobic",3]))
t.test(as.numeric(oxygen_preferences[oxygen_preferences$tolerance=="facultative",3]),as.numeric(oxygen_preferences[oxygen_preferences$tolerance=="anaerobic",3]))
#CTD
t.test(as.numeric(oxygen_preferences[oxygen_preferences$tolerance=="aerobic",4]),as.numeric(oxygen_preferences[oxygen_preferences$tolerance=="facultative",4]))
t.test(as.numeric(oxygen_preferences[oxygen_preferences$tolerance=="aerobic",4]),as.numeric(oxygen_preferences[oxygen_preferences$tolerance=="anaerobic",4]))
t.test(as.numeric(oxygen_preferences[oxygen_preferences$tolerance=="facultative",4]),as.numeric(oxygen_preferences[oxygen_preferences$tolerance=="anaerobic",4]))


