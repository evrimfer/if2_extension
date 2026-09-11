# scripts for making Figure 2 Panel G
# the source data table generated from the direct disorder results using script 07_preferences.R
# SourceDataFig1E.csv can be retrieved from Source_Data_Fig1-4.xlsx
temperature_preferences=read.csv("/sources/SourceDataFig1G.csv")

# number of samples
sum(temperature_preferences$temperatureType=="psychrophile")
sum(temperature_preferences$temperatureType=="mesophile")
sum(temperature_preferences$temperatureType=="thermophile")

# psychrophile n=4
# mesophile n=707
# thermophile n=84


library(ggplot2)

temperature_preferences$temperatureType=factor(temperature_preferences$temperatureType, levels=c("psychrophile","mesophile","thermophile"))


ggplot(temperature_preferences, aes(x=temperatureType, y=as.numeric(L_IF2N), fill=source))+
  geom_boxplot(width=0.5, alpha=0.6) +
  geom_jitter(width=0, alpha=0.4)+
  scale_fill_manual(values=c("#fdce9d"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",
        axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x="", y="IF2 N-terminal Length")

ggplot(temperature_preferences, aes(x=temperatureType, y=as.numeric(L_IF2C), fill=source))+
  geom_boxplot(width=0.5, alpha=0.6) +
  geom_jitter(width=0, alpha=0.4)+
  scale_fill_manual(values=c("#fdce9d"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",
        axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x="", y="IF2 C-terminal Length")

# statistical test between preferences
# NTD
t.test(as.numeric(temperature_preferences[temperature_preferences$temperatureType=="mesophile",3]),as.numeric(temperature_preferences[temperature_preferences$temperatureType=="psychrophile",3]))
t.test(as.numeric(temperature_preferences[temperature_preferences$temperatureType=="mesophile",3]),as.numeric(temperature_preferences[temperature_preferences$temperatureType=="thermophile",3]))
t.test(as.numeric(temperature_preferences[temperature_preferences$temperatureType=="psychrophile",3]),as.numeric(temperature_preferences[temperature_preferences$temperatureType=="thermophile",3]))
# CTD
t.test(as.numeric(temperature_preferences[temperature_preferences$temperatureType=="mesophile",4]),as.numeric(temperature_preferences[temperature_preferences$temperatureType=="psychrophile",4]))
t.test(as.numeric(temperature_preferences[temperature_preferences$temperatureType=="mesophile",4]),as.numeric(temperature_preferences[temperature_preferences$temperatureType=="thermophile",4]))
t.test(as.numeric(temperature_preferences[temperature_preferences$temperatureType=="psychrophile",4]),as.numeric(temperature_preferences[temperature_preferences$temperatureType=="thermophile",4]))

