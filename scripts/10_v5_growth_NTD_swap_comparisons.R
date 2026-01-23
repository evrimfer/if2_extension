# growthcurver_analysis.R
# Analyze growth curve data using R package Growthcurver
# See https://cran.r-project.org/web/packages/growthcurver/vignettes/Growthcurver-vignette.html for more information on Growthcurver package

### CHANGE THESE VARIABLES AS NEEDED ###
data <- as.data.frame(read.csv("~/data/growth_curves/NTD_Swap_Growth_curve_summary.csv", header=T, sep=","))

# Define sample groups (e.g., strains or conditions)
sample_group_names <- c("WT","dN","dG","N-Gs","NG-Gs")

# Define columns in CSV growth data file that belong to each group
# Numbering starts AFTER the "time" column
# E.g., column 1 is first sample column with measurement data
# List order must match sample group names above
sample_group_cols <- list(1:8,9:16,17:24,25:32,33:40)


### SETUP ###

#Load libraries
library(growthcurver)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(pals)
library(tidyr)
library(ggpubr)
library(pdftools)
library(grid)
library(tools)



### CALCULATE SAMPLE GROUP STATISTICS ###

# Add group labels to data summary
data_summary=data
data_summary$group <- NA

for (i in 1:(length(sample_group_names))) {
  data_summary$group[unlist(sample_group_cols[i])] <- sample_group_names[i]
}  

# Calculate summary statistics by sample group
data_summary_stats <- data_summary %>%
  group_by(group) %>%
  summarise(
    mean_t_gen=mean(t_gen),
    sd_t_gen=sd(t_gen),
    mean_t_mid=mean(t_mid),
    sd_t_mid=sd(t_mid),
    mean_r=mean(r),
    sd_r=sd(r)
  )
data_summary_stats <- data_summary_stats[match(sample_group_names, data_summary_stats$group),] # Change group ordering


### PLOT DATA ###

# Make barplots of growth parameters
data_summary_stats$group <- factor(data_summary_stats$group,  # Change group ordering
                                   levels = sample_group_names)

# Doubling time plot (t_gen)
t_gen_barplot <- ggplot(data=data_summary_stats, aes(x=group, level=sample_group_names, y=mean_t_gen, fill=group, col=group)) +
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean_t_gen-sd_t_gen, ymax=mean_t_gen+sd_t_gen), width=.05) +
  labs(x="Group", y="Doubling time (h)") +
  scale_fill_manual(name = "Group", values = c("red","black","gray30","gray60","gray80")) +
  
  scale_color_manual(name="Group",values = c("black","black","black","black","black"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

t_gen_barplot
# Midpoint time plot (t_mid)
t_mid_barplot <- ggplot(data=data_summary_stats, aes(x=group, y=mean_t_mid, fill=group, col=group)) +
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean_t_mid-sd_t_mid, ymax=mean_t_mid+sd_t_mid), width=.05) +
  labs(x="Group", y="Midpoint time (h)") +
  scale_fill_manual(name = "Group", values = c("red","black","gray30","gray60","gray80")) +
  
  scale_color_manual(name="Group",values = c("black","black","black","black","black"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
t_mid_barplot

# Rate plot (r)
rate_barplot <- ggplot(data=data_summary_stats, aes(x=group, y=mean_r, fill=group, col=group)) +
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean_r-sd_r, ymax=mean_r+sd_r), width=.05) +
  labs(x="Group", y="Growth rate") +
  scale_fill_manual(name = "Group", values = c("red","black","gray30","gray60","gray80")) +
  
  scale_color_manual(name="Group",values = c("black","black","black","black","black"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
rate_barplot


####################################
# Relative Folds
### PREP DATA ###

#Load libraries
library(growthcurver)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(multcomp)


# Import data with comma-separated vales; replace with appropriate path and file name
data <- as.data.frame(read.csv("~/data/growth_curves/NTD_Swap_Growth_curve_summary.csv", header=T, sep=","))

#####  Data Frame #####
data$group=factor(data$group, levels=c("WT","dN","dG","N-Gs","NG-Gs"))


type=rep(c("WT","dN","dG","N-Gs","NG-Gs"), each=8)

data$note=type
data$note=factor(data$note, levels=c("WT","dN","dG","N-Gs","NG-Gs"))

data$Fold.doubling=as.numeric(data$Fold.doubling)
data$Fold.growth=as.numeric(data$Fold.growth)

# for barplot - create new data frame with mean and sd values - doubling time
df_mean_dt=matrix(ncol=3)
df_mean_dt=rbind(df_mean_dt, c("WT", mean(data$Fold.doubling[1:8]), sd(data$Fold.doubling[1:8])))
df_mean_dt=rbind(df_mean_dt, c("dN", mean(data$Fold.doubling[9:16]), sd(data$Fold.doubling[9:16])))
df_mean_dt=rbind(df_mean_dt, c("dG", mean(data$Fold.doubling[17:24]), sd(data$Fold.doubling[17:24])))
df_mean_dt=rbind(df_mean_dt, c("N-Gs", mean(data$Fold.doubling[25:32]), sd(data$Fold.doubling[25:32])))
df_mean_dt=rbind(df_mean_dt, c("NG-Gs", mean(data$Fold.doubling[33:40]), sd(data$Fold.doubling[33:40])))

df_mean_dt=as.data.frame(df_mean_dt[-1,])
colnames(df_mean_dt)=c("sample","mean_dt","sd_dt")
df_mean_dt$sample=factor(df_mean_dt$sample, levels=c("WT","dN","dG","N-Gs","NG-Gs"))


# for barplot - create new data frame with mean and sd values - doubling time
df_mean_rate=matrix(ncol=3)
df_mean_rate=rbind(df_mean_rate, c("WT", mean(data$Fold.growth[1:8]), sd(data$Fold.growth[1:8])))
df_mean_rate=rbind(df_mean_rate, c("dN", mean(data$Fold.growth[9:16]), sd(data$Fold.growth[9:16])))
df_mean_rate=rbind(df_mean_rate, c("dG", mean(data$Fold.growth[17:24]), sd(data$Fold.growth[17:24])))
df_mean_rate=rbind(df_mean_rate, c("N-Gs", mean(data$Fold.growth[25:32]), sd(data$Fold.growth[25:32])))
df_mean_rate=rbind(df_mean_rate, c("NG-Gs", mean(data$Fold.growth[33:40]), sd(data$Fold.growth[33:40])))

df_mean_rate=as.data.frame(df_mean_rate[-1,])
colnames(df_mean_rate)=c("sample","mean_r","sd_r")
df_mean_rate$sample=factor(df_mean_rate$sample, levels=c("WT","dN","dG","N-Gs","NG-Gs"))


#### PLOTTING #####
# plot for doubling time
dodge=position_dodge(width=0.7)
ggplot(df_mean_dt, aes(y=as.numeric(mean_dt), x=sample, shape=sample, fill=sample, color=sample)) + 
  geom_pointrange(aes(ymin=as.numeric(mean_dt)-as.numeric(sd_dt), ymax=as.numeric(mean_dt)+as.numeric(sd_dt)), size=1.5)+
  scale_shape_manual(values=c(21, 21, 21, 21, 21))+
  scale_fill_manual(values=c("red","black","gray30","gray60","gray80"))+
  scale_color_manual(values=c("black","black","black","black","black"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none", axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x=" ", y="Doubling Time (h) Relative to WT IF2")+
  scale_y_continuous(limits = c(0.5, 1.5))

# plot for growth rate
dodge=position_dodge(width=0.7)
ggplot(df_mean_rate, aes(y=as.numeric(mean_r), x=sample, shape=sample, fill=sample, color=sample)) + 
  geom_pointrange(aes(ymin=as.numeric(mean_r)-as.numeric(sd_r), ymax=as.numeric(mean_r)+as.numeric(sd_r)), size=1.5)+
  scale_shape_manual(values=c(21, 21, 21, 21, 21))+
  scale_fill_manual(values=c("red","black","gray30","gray60","gray80"))+
  scale_color_manual(values=c("black","black","black","black","black"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none", axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x=" ", y="Growth Rate Relative to WT IF2")+
  scale_y_continuous(limits = c(0.5, 1.5))


wilcox.test(data$Fold.doubling[1:8],data$Fold.doubling[9:16])
wilcox.test(data$Fold.doubling[1:8],data$Fold.doubling[17:24])
wilcox.test(data$Fold.doubling[1:8],data$Fold.doubling[25:32])
wilcox.test(data$Fold.doubling[1:8],data$Fold.doubling[33:40])

wilcox.test(data$Fold.growth[1:8],data$Fold.growth[9:16])
wilcox.test(data$Fold.growth[1:8],data$Fold.growth[17:24])
wilcox.test(data$Fold.growth[1:8],data$Fold.growth[25:32])
wilcox.test(data$Fold.growth[1:8],data$Fold.growth[33:40])





