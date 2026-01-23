# 12_invitro.R
# Script for PURExpress luminescence analysis for measuring N-terminal truncation effect on translation

library(ggplot2)

# Import data with comma-separated vales; replace with appropriate path and file name
data <- read.table(file = "~/data/in_vitro/IF2-WT-deltaN-PURE-Express-LUMINESCENCE-cleaned-averages-for-wt-dN-average-rep-fold.csv", header = TRUE, sep = ",", row.names = 1)

# if there is any, Cut out all non-deltaN-IF2 
# data_delta_IF2 <- data[-c(1:3),]
# data_delta_NEB <- data_delta_IF2[-c(4:4),]

# transpose data 
data_transpose <- t(data)
data_dataframe <- as.data.frame(data_transpose)

values <- c(data_dataframe$`noluc-noif2`, data_dataframe$`yesluc-noif2`,data_dataframe$`yesluc-wtif2`,data_dataframe$`yesluc-dN`)

groups <- rep(colnames(data_dataframe), each=3)

data_for_plot <- as.data.frame(cbind(values, groups))

# box plot
data_for_plot$groups=factor(data_for_plot$groups, levels=c("noluc-noif2", "yesluc-noif2","yesluc-wtif2","yesluc-dN"))

ggplot(data_for_plot, aes(x = groups, y = as.numeric(values), fill=groups))+ 
  geom_boxplot()+
  geom_point()+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = "" , y = 'Luminescence (AU)', size=16)+
  ggtitle("PURExpress∆IF2")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(angle=0))

Neg.vs.WT <- t.test(data_dataframe$`yesluc-wtif2`, data_dataframe$`yesluc-dN`)
minus.Luc.vs.plus.Luc <- t.test(data_dataframe$`∆IF2 -Luc`, data_dataframe$`∆IF2 +Luc`)
minus.Luc.vs.Neg <- t.test(data_dataframe$`∆IF2 -Luc`, data_dataframe$`∆IF2 -Luc +NEB IF2`)
plus.Luc.vs.Neg <- t.test(data_dataframe$`∆IF2 -Luc +NEB IF2`, data_dataframe$`∆IF2 +Luc`)

