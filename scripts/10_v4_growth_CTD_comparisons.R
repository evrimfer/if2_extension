# growthcurver_analysis.R
# Analyze growth curve data using R package Growthcurver
# See https://cran.r-project.org/web/packages/growthcurver/vignettes/Growthcurver-vignette.html for more information on Growthcurver package


### CHANGE THESE VARIABLES AS NEEDED ###

# CSV growth data file path, change based on temperature or pH, or anaerobic
data_fname <- "~/data/growth_curves/CTD_IF2_Growth_curves_30_summary.csv"

# Specify output location
out_fpath <- "~/data/growth_curves/"

# Replace with "hours, seconds", or "minutes" (check spelling). 
# Will automatically convert to hours
time_units <- "minutes"

# Specify cutoff time point *in hours* to trim data.
# Cutoff is inclusive (e.g., "24" will include first 24 hours)
# Value of "0" means no trimming
# Value that is outside of data time range will result in no trimming
trim_at_time <-0

# Define sample groups (e.g., strains or conditions)
sample_group_names <- c("IF2-WT","IF2-DSM5456","IF2-DSM15908","IF2-H-salinum","IF2-G-aggregans","IF2-T-modestius","IF2-P-aerophilum","IF2-S-osmophilus")

# Define columns in CSV growth data file that belong to each group
# Numbering starts AFTER the "time" column
# E.g., column 1 is first sample column with measurement data
# List order must match sample group names above
#sample_group_cols <- list(1:5,6:10,11:15,16:20,21:25,26:30,31:35,36:40)

# FOR ANAEROBIC DATA USE THIS INSTEAD
sample_group_cols <- list(1:3,4:6,7:9,10:12,13:15,16:18,19:21,22:24)

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
    sd_r=sd(r),
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
  scale_fill_manual(name = "Group", values = c("red","black","gray15",
                                               "gray35","gray55","gray75",
                                               "gray90","white")) +
  
  scale_color_manual(name="Group",values = c("black","black","black",
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
  scale_fill_manual(name = "Group", values = c("red","black","gray15",
                                               "gray35","gray55","gray75",
                                               "gray90","white")) +
  
  scale_color_manual(name="Group",values = c("black","black","black",
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
  scale_fill_manual(name = "Group", values = c("red","black","gray15",
                                               "gray35","gray55","gray75",
                                               "gray90","white")) +
  
  scale_color_manual(name="Group",values = c("black","black","black",
                                             "black","black","black",
                                             "black","black","black"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
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
data_for_scatterplt$group=factor(data_for_scatterplt$group, levels=c("IF2-WT","IF2-DSM5456","IF2-DSM15908","IF2-H-salinum","IF2-G-aggregans","IF2-T-modestius","IF2-P-aerophilum","IF2-S-osmophilus"))

growth_curve_plot <- ggplot(data_for_scatterplt, aes(x = time, y = OD)) +
  geom_point(alpha = 1, aes(color=group, fill = group, shape=group), size = 1) +
  scale_y_continuous(trans = "log", breaks = y_ticks_log, labels = y_ticks) +
  scale_color_manual(name = "Group", values = c("red","black","gray10",
                                                 "gray20","gray35","gray45",
                                                 "gray55","gray65"), breaks = sample_group_names) +
  scale_fill_manual(name = "Group", values = c("red","black","gray10",
                                                "gray20","gray35","gray45",
                                                "gray55","gray65"), breaks = sample_group_names) +
  scale_shape_manual(values=c(21, 21, 21, 22, 22, 24, 24, 25))+
  
  labs(x="Time (h)", y="ln OD") +
  theme_classic()

# Add smoothed line per sample group
for (i in 1:(length(sample_group_names))) {
  growth_curve_plot <- growth_curve_plot + geom_smooth(
    data = data_for_scatterplt[data_for_scatterplt$group==(sample_group_names[i]),],
    se = TRUE, size=1,color = (c("red","black","gray10",
                                 "gray20","gray35","gray45",
                                 "gray55","gray65"))[i])+
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
ggarrange(growth_curve_plot, 
          ggarrange(t_gen_barplot, t_mid_barplot, labels = c("B", "C"), ncol = 2),
          labels = "A",
          nrow = 2
)

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

### PREP DATA ###

#Load libraries
library(growthcurver)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(multcomp)


# Import data with comma-separated vales; replace with appropriate path and file name
data <- as.data.frame(read.csv("~/data/growth_curves/CTD_IF2_Growth_curves_37_anaerobic_summary.csv", header=T, sep=","))

#####  Data Frame #####
data$group=factor(data$group, levels=c("IF2-WT","IF2-DSM5456","IF2-DSM15908","IF2-H-salinum","IF2-G-aggregans","IF2-T-modestius","IF2-P-aerophilum","IF2-S-osmophilus"))

# for anaerobic analysis use this line 
# type=rep(c("IF2-WT","IF2-DSM5456","IF2-DSM15908","IF2-H-salinum","IF2-G-aggregans","IF2-T-modestius","IF2-P-aerophilum","IF2-S-osmophilus"), each=3)
type=rep(c("IF2-WT","IF2-DSM5456","IF2-DSM15908","IF2-H-salinum","IF2-G-aggregans","IF2-T-modestius","IF2-P-aerophilum","IF2-S-osmophilus"), each=5)

data$note=type
data$note=factor(data$note, levels=c("IF2-WT","IF2-DSM5456","IF2-DSM15908","IF2-H-salinum","IF2-G-aggregans","IF2-T-modestius","IF2-P-aerophilum","IF2-S-osmophilus"))


# for barplot - create new data frame with mean and sd values - doubling time - FOR ANAEROBIC
df_mean=matrix(ncol=3)
df_mean=rbind(df_mean, c("IF2-WT", mean(data$Fold.doubling[1:3]), sd(data$Fold.doubling[1:3])))
df_mean=rbind(df_mean, c("IF2-DSM5456", mean(data$Fold.doubling[4:6]), sd(data$Fold.doubling[4:6])))
df_mean=rbind(df_mean, c("IF2-DSM15908", mean(data$Fold.doubling[7:9]), sd(data$Fold.doubling[7:9])))
df_mean=rbind(df_mean, c("IF2-H-salinum", mean(data$Fold.doubling[10:12]), sd(data$Fold.doubling[10:12])))
df_mean=rbind(df_mean, c("IF2-G-aggregans", mean(data$Fold.doubling[13:15]), sd(data$Fold.doubling[13:15])))
df_mean=rbind(df_mean, c("IF2-T-modestius", mean(data$Fold.doubling[16:18]), sd(data$Fold.doubling[16:18])))
df_mean=rbind(df_mean, c("IF2-P-aerophilum", mean(data$Fold.doubling[19:21]), sd(data$Fold.doubling[19:21])))
df_mean=rbind(df_mean, c("IF2-S-osmophilus", mean(data$Fold.doubling[22:24]), sd(data$Fold.doubling[22:24])))

# for barplot - create new data frame with mean and sd values - doubling time
df_mean=matrix(ncol=3)
df_mean=rbind(df_mean, c("IF2-WT", mean(data$Fold.doubling[1:5]), sd(data$Fold.doubling[1:5])))
df_mean=rbind(df_mean, c("IF2-DSM5456", mean(data$Fold.doubling[6:10]), sd(data$Fold.doubling[6:10])))
df_mean=rbind(df_mean, c("IF2-DSM15908", mean(data$Fold.doubling[11:15]), sd(data$Fold.doubling[11:15])))
df_mean=rbind(df_mean, c("IF2-H-salinum", mean(data$Fold.doubling[16:20]), sd(data$Fold.doubling[16:20])))
df_mean=rbind(df_mean, c("IF2-G-aggregans", mean(data$Fold.doubling[21:25]), sd(data$Fold.doubling[21:25])))
df_mean=rbind(df_mean, c("IF2-T-modestius", mean(data$Fold.doubling[26:30]), sd(data$Fold.doubling[26:30])))
df_mean=rbind(df_mean, c("IF2-P-aerophilum", mean(data$Fold.doubling[31:35]), sd(data$Fold.doubling[31:35])))
df_mean=rbind(df_mean, c("IF2-S-osmophilus", mean(data$Fold.doubling[36:40]), sd(data$Fold.doubling[36:40])))

df_mean=as.data.frame(df_mean[-1,])
colnames(df_mean)=c("sample","mean_dt","sd_dt")
df_mean$sample=factor(df_mean$sample, levels=c("IF2-WT","IF2-DSM5456","IF2-DSM15908","IF2-H-salinum","IF2-G-aggregans","IF2-T-modestius","IF2-P-aerophilum","IF2-S-osmophilus"))


#### PLOTTING #####
# barplot
dodge=position_dodge(width=0.7)
ggplot(df_mean, aes(y=as.numeric(mean_dt), x=sample, shape=sample, fill=sample, color=sample)) + 
  geom_pointrange(aes(ymin=as.numeric(mean_dt)-as.numeric(sd_dt), ymax=as.numeric(mean_dt)+as.numeric(sd_dt)), size=1.5)+
  scale_shape_manual(values=c(21, 21, 21, 22, 22, 24, 24, 25))+
  scale_fill_manual(values=c("red","black","white","black","white","black","white","black"))+
  scale_color_manual(values=c("red","black","black","black","black","black","black","black"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none", axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x=" ", y="Doubling Time (h) Relative to WT IF2")+
  scale_y_continuous(breaks = c(0.4,0.6,0.8,1.0,1.2,1.4,1.6), limits = c(0.4, 1.6))

## use this for anaerobic test
# t.test(data$Fold.doubling[1:3],data$Fold.doubling[4:6])
# t.test(data$Fold.doubling[1:3],data$Fold.doubling[7:9])
# t.test(data$Fold.doubling[1:3],data$Fold.doubling[10:12])
# t.test(data$Fold.doubling[1:3],data$Fold.doubling[13:15])
# t.test(data$Fold.doubling[1:3],data$Fold.doubling[16:18])
# t.test(data$Fold.doubling[1:3],data$Fold.doubling[19:21])
# t.test(data$Fold.doubling[1:3],data$Fold.doubling[22:24])

## use this for anaerobic test
# wilcox.test(data$Fold.doubling[1:3],data$Fold.doubling[4:6])
# wilcox.test(data$Fold.doubling[1:3],data$Fold.doubling[7:9])
# wilcox.test(data$Fold.doubling[1:3],data$Fold.doubling[10:12])
# wilcox.test(data$Fold.doubling[1:3],data$Fold.doubling[13:15])
# wilcox.test(data$Fold.doubling[1:3],data$Fold.doubling[16:18])
# wilcox.test(data$Fold.doubling[1:3],data$Fold.doubling[19:21])
# wilcox.test(data$Fold.doubling[1:3],data$Fold.doubling[22:24])

wilcox.test(data$Fold.doubling[1:5],data$Fold.doubling[6:10])
wilcox.test(data$Fold.doubling[1:5],data$Fold.doubling[11:15])
wilcox.test(data$Fold.doubling[1:5],data$Fold.doubling[16:20])
wilcox.test(data$Fold.doubling[1:5],data$Fold.doubling[21:25])
wilcox.test(data$Fold.doubling[1:5],data$Fold.doubling[26:30])
wilcox.test(data$Fold.doubling[1:5],data$Fold.doubling[31:35])
wilcox.test(data$Fold.doubling[1:5],data$Fold.doubling[36:40])
