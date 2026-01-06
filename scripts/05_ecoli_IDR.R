# 05_ecoliIDR.R
# plot the IDR and pLDDT values for E. coli IF2

# Load libraries
library(ggplot2)
library(tidyverse)

# read the predicted data
idr_predict=read.csv("~/intrinsic_disorder/Ecoli_IF2_metapredict.csv", header=T)

# modify the data frame for plotting
idr_df=as.data.frame(cbind(rep(idr_predict$Index,2),c(rep("Disorder",890), rep("pLDDT", 890)),c(idr_predict$Disorder,(idr_predict$ppLDDT)/100)))
colnames(idr_df)=c("Index","Prediction","Values")

# plot
ggplot(idr_df, aes(x=as.numeric(Index), y=as.numeric(Values), col=Prediction))+
  scale_color_manual(values=c("red","blue"))+
  geom_line()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(title=" ", x="Residue Position", y="Disorder (red) / pLDDT (blue)")+
  geom_vline(xintercept = c(157,294,387,561,658,793), linetype="dotted", 
                color = "blue", size=1.5)

# to visualize pLDDT on the structure 
a=as.matrix(paste("\t",paste("/A:",idr_predict$Index, sep=""), "\t",idr_predict$ppLDDT, sep=""))
write.table(a[,1],"plddt.defattr", col.names=F, row.names=F, quote=F)
