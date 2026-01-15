# growthcurver_analysis.R
# Analyze growth curve data using R package Growthcurver
# See https://cran.r-project.org/web/packages/growthcurver/vignettes/Growthcurver-vignette.html for more information on Growthcurver package

### CHANGE THESE VARIABLES AS NEEDED ###

# CSV growth data file path
data_fname <- "~/data/growth_curves/dN-IF2-Growth-Curves-Temperatures.csv"

# Specify output location
out_fpath <- "~/data/growth_curves/"

# Replace with "hours, seconds", or "minutes" (check spelling). 
# Will automatically convert to hours
time_units <- "minutes"

# Specify cutoff time point *in hours* to trim data.
# Cutoff is inclusive (e.g., "24" will include first 24 hours)
# Value of "0" means no trimming
# Value that is outside of data time range will result in no trimming
trim_at_time <- 0

# Define sample groups (e.g., strains or conditions)
sample_group_names <- c("WT-25","dN-25","dG1-25",
                        "WT-30","dN-30","dG1-30",
                        "WT-37","dN-37","dG1-37",
                        "WT-42","dN-42","dG1-42")

# Define columns in CSV growth data file that belong to each group
# Numbering starts AFTER the "time" column
# E.g., column 1 is first sample column with measurement data
# List order must match sample group names above
sample_group_cols <- list(1:5,6:10,11:15,
                          16:20,21:25,26:30,
                          32:35,36:40,41:45,
                          46:50,51:55,56:60)

# NOTE: Will automatically perform background correction if a "blank" column is provided.
# Values in blank column are subtracted from all other columns per row (bg_correct = "blank")
# If "blank" column is not provided, will correct background by subtracting all values in a column
# by the minimum value in the same column (bg_correct = "min").
# This method is suitable if the background is NOT expected to change during the experiment.

#############


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


### IMPORT GROWTH DATA ###

# Import data
data <- read.csv(file = data_fname, header = TRUE, sep = ",")


# Define output file basename
out_base <- paste(out_fpath, file_path_sans_ext(basename(data_fname)), sep="")


### CLEAN UP DATA ###

# Convert time units if provided in minutes or seconds
if (time_units=="minutes") {
  print("Converting time units from minutes to hours")
  data$time <- data$time / 60 # Convert minutes to hours
} else if (time_units=="seconds") {
  print("Converting time units from seconds to hours")
  data$time <- data$time / 60 / 60 # Convert seconds to hours
} else if (time_units=="hours") {
  print("Time units are in hours")
} else {
  print("Cannot recognize specified time units. Will assume in hours")
}

# Clean up na values
data[is.na(data)] <- 0 

# Trim data (optional)

if (trim_at_time != 0 & trim_at_time <= (tail(data$time, n=1)-1)) {
  data <- data %>% filter(row_number() < which(time==(trim_at_time+1)))
  print(paste("Data trimmed to", trim_at_time, "hours"))
} else {
  print("Warning: Data not trimmed because trim_at_time value is either 0 or outside data time range")
}

### SUMMARIZE DATA ###

# Summarize data and export fitted growth curves

if ("blank" %in% colnames(data)) {
  data_summary_prelim <- SummarizeGrowthByPlate(
    plate=data,
    bg_correct = "blank",
    plot_fit = TRUE,
    plot_file = paste(out_base, "_fitted_curves.pdf", sep=""))
  data_summary <- head(data_summary_prelim, -1) # Remove residual blank row from data summary
} else {
  data_summary <- SummarizeGrowthByPlate(
    plate=data,
    bg_correct = "min",
    plot_fit = TRUE,
    plot_file = paste(out_base, "_fitted_curves.pdf", sep=""))
} 

# IMPORTANT to inspect fitted growth curve PDF for samples with poor fit
# Poor fit might result from poorly resolved lag and stationary phase
# May be necessary to trim data to improve fit

### CALCULATE SAMPLE GROUP STATISTICS ###

# Add group labels to data summary
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
  scale_fill_manual(name = "Group", values = c("red","black","gray",
                                               "red","black","gray",
                                               "red","black","gray",
                                               "red","black","gray")) +
  
  scale_color_manual(name="Group",values = c("black","black","black",
                                             "black","black","black",
                                             "black","black","black",
                                             "black","black","black"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

t_gen_barplot
# Midpoint time plot (t_mid)
t_mid_barplot <- ggplot(data=data_summary_stats, aes(x=group, y=mean_t_mid, fill=group, col=group)) +
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean_t_mid-sd_t_mid, ymax=mean_t_mid+sd_t_mid), width=.05) +
  labs(x="Group", y="Midpoint time (h)") +
  scale_fill_manual(name = "Group", values = c("red","black","gray",
                                               "red","black","gray",
                                               "red","black","gray",
                                               "red","black","gray")) +
  
  scale_color_manual(name="Group",values = c("black","black","black",
                                             "black","black","black",
                                             "black","black","black",
                                             "black","black","black"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
t_mid_barplot

# Rate plot (r)
rate_barplot <- ggplot(data=data_summary_stats, aes(x=group, y=mean_r, fill=group, col=group)) +
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymin=mean_r-sd_r, ymax=mean_r+sd_r), width=.05) +
  labs(x="Group", y="Growth rate") +
  scale_fill_manual(name = "Group", values = c("red","black","gray",
                                               "red","black","gray",
                                               "red","black","gray",
                                               "red","black","gray")) +
  
  scale_color_manual(name="Group",values = c("black","black","black",
                                             "black","black","black",
                                             "black","black","black",
                                             "black","black","black"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylim(0,1.7)
rate_barplot

# Make growth curve scatter plot

# If "blank" provided, subtract blank value from all measurements
if ("blank" %in% colnames(data)) {
  data_blanked <- data
  data_blanked[2:ncol(data_blanked)] <- data_blanked[2:ncol(data_blanked)]-data_blanked$blank
  data_blanked <- data_blanked[,!names(data_blanked) %in% c("blank")]
} else {
  data_blanked <- data # do nothing
}

# Reformat data for scatter plot
## so I think averaging the 6 replicates of each biological rep to make
## total of 3 data points (1 per bio rep) per time point per set is good.
data_for_scatterplt <- data_blanked %>% gather(sample, OD, -time)

data_for_scatterplt$group <- NA

for (i in 1:(length(sample_group_names))) {
  rows = c()
  for (j in sample_group_cols[[i]]) {
    new_rows = (j*nrow(data)-(nrow(data)-1)):(j*nrow(data))
    rows = append(rows, new_rows)
  }
  data_for_scatterplt$group[rows] <- sample_group_names[i]
} 

# Make scatter plot
y_ticks <- -10:10
y_ticks_log <- sapply(y_ticks, function(x) exp(x))

growth_curve_plot <- ggplot(data_for_scatterplt, aes(x = time, y = OD)) +
  geom_point(alpha = 0.70, aes(color = group), size = 1.2) +
  scale_x_continuous(limits = c(0,15))+
  scale_y_continuous(trans = "log", breaks = y_ticks_log, labels = y_ticks) +
  scale_color_manual(name = "Group", values = c("red","black","gray",
                                                "red","black","gray",
                                                "red","black","gray",
                                                "red","black","gray"), breaks = sample_group_names) +
  labs(x="Time (h)", y="ln OD") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none", 
        axis.text.y=element_text(size=18, face="bold", color="black"), axis.text.x = element_text(size=18, face="bold", color="black")) +
  theme_classic()


# Add smoothed line per sample group
for (i in 1:(length(sample_group_names))) {
  growth_curve_plot <- growth_curve_plot + geom_smooth(
    data = data_for_scatterplt[data_for_scatterplt$group==(sample_group_names[i]),],
    se = TRUE, size=2,color = (c("red","black","gray",
                                  "red","black","gray",
                                  "red","black","gray",
                                  "red","black","gray"))[i])+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none", 
          axis.text.x=element_text(size=18, face="bold", color="black"), axis.text.y = element_text(size=18, face="bold", color="black")) +
    theme_classic()
}

growth_curve_plot

### SAVE REPORT ###

# Open PDF device
pdf(file=(paste(out_fpath, "temp_plots.pdf", sep="")), height=8.5, width=11)

# Save summary table
total_rows_per_page = 38 
start_row = 1 

if(total_rows_per_page > nrow(data_summary)){
  end_row = nrow(data_summary)
}else {
  end_row = total_rows_per_page 
}    

for(i in 1:ceiling(nrow(data_summary)/total_rows_per_page)){
  
  grid.newpage()   
  
  grid.table(data_summary[start_row:end_row, ],
             theme=ttheme_minimal(
               core=list(fg_params=list(cex=0.5)),
               colhead=list(fg_params=list(cex=0.5)),
               rowhead=list(fg_params=list(cex=0.5)))
  )             
  
  start_row = end_row + 1
  
  if((total_rows_per_page + end_row) < nrow(data_summary)){
    
    end_row = total_rows_per_page + end_row
    
  }else {
    
    end_row = nrow(data_summary)
  }    
}

# Save summary stats table
grid.newpage()
grid.table(data_summary_stats, theme=ttheme_minimal(
  core=list(fg_params=list(cex=0.5)),
  colhead=list(fg_params=list(cex=0.5)),
  rowhead=list(fg_params=list(cex=0.5))
))

# Save plots
ggarrange(growth_curve_plot)
ggarrange(t_gen_barplot, t_mid_barplot, labels = c("B", "C"), ncol = 2)

dev.off()

# Combine with fitted curve PDF into single report
pdf_combine(c(paste(out_base, "_fitted_curves.pdf", sep=""), paste(out_fpath, "temp_plots.pdf", sep="")),
            output=(paste(out_base, "_report.pdf", sep="")))
file.remove(c(paste(out_base, "_fitted_curves.pdf", sep=""), paste(out_fpath, "temp_plots.pdf", sep="")))

# Write summary and summary stats to csv files
write.csv(data_summary, file = (paste(out_base, "_summary.csv", sep="")), row.names=FALSE)
write.csv(data_summary_stats, file = (paste(out_base, "_summary_stats.csv", sep="")), row.names = FALSE, quote = FALSE)

print(paste("Saved report to ", out_base, "_report.pdf", sep=""))
print(paste("Saved data summary files to ", out_base, "_summary.csv and ",  out_base, "_summary_stats.csv", sep=""))

print("Finished with Growthcurver analysis")


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
data <- as.data.frame(read.csv("~/data/growth_curves/dN-IF2-Growth-Curves-Temperatures_summary.csv", header=T, sep=","))

#####  Data Frame #####
data$group=factor(data$group, levels=c("WT-25","dN-25","dG1-25",
                                       "WT-30","dN-30","dG1-30",
                                       "WT-37","dN-37","dG1-37",
                                       "WT-42","dN-42","dG1-42"))


type=rep(c("WT-25","dN-25","dG1-25",
           "WT-30","dN-30","dG1-30",
           "WT-37","dN-37","dG1-37",
           "WT-42","dN-42","dG1-42"), each=5)

data$note=type
data$note=factor(data$note, levels=c("WT-25","dN-25","dG1-25",
                                     "WT-30","dN-30","dG1-30",
                                     "WT-37","dN-37","dG1-37",
                                     "WT-42","dN-42","dG1-42"))

data$Fold.doubling=as.numeric(data$Fold.doubling)
data$Fold.growth=as.numeric(data$Fold.growth)

# for barplot - create new data frame with mean and sd values - doubling time
df_mean_dt=matrix(ncol=3)
df_mean_dt=rbind(df_mean_dt, c("WT-25", mean(data$Fold.doubling[1:5]), sd(data$Fold.doubling[1:5])))
df_mean_dt=rbind(df_mean_dt, c("dN-25", mean(data$Fold.doubling[6:10]), sd(data$Fold.doubling[6:10])))
df_mean_dt=rbind(df_mean_dt, c("dG1-25", mean(data$Fold.doubling[11:15]), sd(data$Fold.doubling[11:15])))
df_mean_dt=rbind(df_mean_dt, c("WT-30", mean(data$Fold.doubling[16:20]), sd(data$Fold.doubling[16:20])))
df_mean_dt=rbind(df_mean_dt, c("dN-30", mean(data$Fold.doubling[21:25]), sd(data$Fold.doubling[21:25])))
df_mean_dt=rbind(df_mean_dt, c("dG1-30", mean(data$Fold.doubling[26:30]), sd(data$Fold.doubling[26:30])))
df_mean_dt=rbind(df_mean_dt, c("WT-37", mean(data$Fold.doubling[32:35]), sd(data$Fold.doubling[32:35])))
df_mean_dt=rbind(df_mean_dt, c("dN-37", mean(data$Fold.doubling[36:40]), sd(data$Fold.doubling[36:40])))
df_mean_dt=rbind(df_mean_dt, c("dG1-37", mean(data$Fold.doubling[41:45]), sd(data$Fold.doubling[41:45])))
df_mean_dt=rbind(df_mean_dt, c("WT-42", mean(data$Fold.doubling[46:50]), sd(data$Fold.doubling[46:50])))
df_mean_dt=rbind(df_mean_dt, c("dN-42", mean(data$Fold.doubling[51:55]), sd(data$Fold.doubling[51:55])))
df_mean_dt=rbind(df_mean_dt, c("dG1-42", mean(data$Fold.doubling[56:60]), sd(data$Fold.doubling[56:60])))

df_mean_dt=as.data.frame(df_mean_dt[-1,])
colnames(df_mean_dt)=c("sample","mean_dt","sd_dt")
df_mean_dt$sample=factor(df_mean_dt$sample, levels=c("WT-25","dN-25","dG1-25",
                                                     "WT-30","dN-30","dG1-30",
                                                     "WT-37","dN-37","dG1-37",
                                                     "WT-42","dN-42","dG1-42"))


# for barplot - create new data frame with mean and sd values - doubling time
df_mean_rate=matrix(ncol=3)
df_mean_rate=rbind(df_mean_rate, c("WT-25", mean(data$Fold.doubling[1:5]), sd(data$Fold.doubling[1:5])))
df_mean_rate=rbind(df_mean_rate, c("dN-25", mean(data$Fold.doubling[6:10]), sd(data$Fold.doubling[6:10])))
df_mean_rate=rbind(df_mean_rate, c("dG1-25", mean(data$Fold.doubling[11:15]), sd(data$Fold.doubling[11:15])))
df_mean_rate=rbind(df_mean_rate, c("WT-30", mean(data$Fold.doubling[16:20]), sd(data$Fold.doubling[16:20])))
df_mean_rate=rbind(df_mean_rate, c("dN-30", mean(data$Fold.doubling[21:25]), sd(data$Fold.doubling[21:25])))
df_mean_rate=rbind(df_mean_rate, c("dG1-30", mean(data$Fold.doubling[26:30]), sd(data$Fold.doubling[26:30])))
df_mean_rate=rbind(df_mean_rate, c("WT-37", mean(data$Fold.doubling[32:35]), sd(data$Fold.doubling[32:35])))
df_mean_rate=rbind(df_mean_rate, c("dN-37", mean(data$Fold.doubling[36:40]), sd(data$Fold.doubling[36:40])))
df_mean_rate=rbind(df_mean_rate, c("dG1-37", mean(data$Fold.doubling[41:45]), sd(data$Fold.doubling[41:45])))
df_mean_rate=rbind(df_mean_rate, c("WT-42", mean(data$Fold.doubling[46:50]), sd(data$Fold.doubling[46:50])))
df_mean_rate=rbind(df_mean_rate, c("dN-42", mean(data$Fold.doubling[51:55]), sd(data$Fold.doubling[51:55])))
df_mean_rate=rbind(df_mean_rate, c("dG1-42", mean(data$Fold.doubling[56:60]), sd(data$Fold.doubling[56:60])))

df_mean_rate=as.data.frame(df_mean_rate[-1,])
colnames(df_mean_rate)=c("sample","mean_r","sd_r")
df_mean_rate$sample=factor(df_mean_rate$sample, levels=c("WT-25","dN-25","dG1-25",
                                                         "WT-30","dN-30","dG1-30",
                                                         "WT-37","dN-37","dG1-37",
                                                         "WT-42","dN-42","dG1-42"))


#### PLOTTING #####
# plot for doubling time
dodge=position_dodge(width=0.7)
ggplot(df_mean_dt, aes(y=as.numeric(mean_dt), x=sample, shape=sample, fill=sample, color=sample)) + 
  geom_pointrange(aes(ymin=as.numeric(mean_dt)-as.numeric(sd_dt), ymax=as.numeric(mean_dt)+as.numeric(sd_dt)), size=1.5)+
  scale_shape_manual(values=c(21, 21, 21,21, 21, 21,21, 21, 21,21, 21, 21))+
  scale_fill_manual(values=c("red","black","gray","red","black","gray","red","black","gray","red","black","gray"))+
  scale_color_manual(values=c("black","black","black","black","black","black","black","black","black","black","black","black"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none", axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x=" ", y="Doubling Time (h) Relative to WT IF2")+
  scale_y_continuous(limits = c(0.5, 2.2))

# plot for growth rate
dodge=position_dodge(width=0.7)
ggplot(df_mean_rate, aes(y=as.numeric(mean_r), x=sample, shape=sample, fill=sample, color=sample)) + 
  geom_pointrange(aes(ymin=as.numeric(mean_r)-as.numeric(sd_r), ymax=as.numeric(mean_r)+as.numeric(sd_r)), size=1.5)+
  scale_shape_manual(values=c(21, 21, 21,21, 21, 21,21, 21, 21,21, 21, 21))+
  scale_fill_manual(values=c("red","black","gray","red","black","gray","red","black","gray","red","black","gray"))+
  scale_color_manual(values=c("black","black","black","black","black","black","black","black","black","black","black","black"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none", axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x=" ", y="Growth Rate Relative to WT IF2")+
  scale_y_continuous(limits = c(0.5, 2.2))



wilcox.test(data$Fold.doubling[1:5],data$Fold.doubling[6:10])
wilcox.test(data$Fold.doubling[1:5],data$Fold.doubling[11:15])

wilcox.test(data$Fold.doubling[16:20],data$Fold.doubling[21:25])
wilcox.test(data$Fold.doubling[16:20],data$Fold.doubling[26:30])

wilcox.test(data$Fold.doubling[32:35],data$Fold.doubling[36:40])
wilcox.test(data$Fold.doubling[32:35],data$Fold.doubling[41:45])

wilcox.test(data$Fold.doubling[46:50],data$Fold.doubling[51:55])
wilcox.test(data$Fold.doubling[46:50],data$Fold.doubling[56:60])


wilcox.test(data$Fold.growth[1:5],data$Fold.growth[6:10])
wilcox.test(data$Fold.growth[1:5],data$Fold.growth[11:15])

wilcox.test(data$Fold.growth[16:20],data$Fold.growth[21:25])
wilcox.test(data$Fold.growth[16:20],data$Fold.growth[26:30])

wilcox.test(data$Fold.growth[32:35],data$Fold.growth[36:40])
wilcox.test(data$Fold.growth[32:35],data$Fold.growth[41:45])

wilcox.test(data$Fold.growth[46:50],data$Fold.growth[51:55])
wilcox.test(data$Fold.growth[46:50],data$Fold.growth[56:60])




