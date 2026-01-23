# 08_meanDisorder.R
# script for IDR

# Load libraries
library(ggtree)
library(ggplot2)
library(ape)
library(seqinr)
library(see)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(ggsignif)

# Load the sequence lengths for N-terminal and C-terminal extensions
nterminal=read.csv("~/data/lengths/IF2_NTD_lengths.tsv", sep="\t", header=T)
cterminal=read.csv("~/data/lengths/IF2_CTD_lengths.tsv", sep="\t", header=T)

# Load the data frame coming from Metapredict v3
idr_predictions=read.csv("~/data/intrinsic_disorder/IF2_disorder_predictions.csv", header=F, sep=",")

# find mean IDRs in each region of each sequence
df_mean_idr=matrix(ncol=3)
colnames(df_mean_idr)=c("sequence_name","mean","region")

for(sequence_name in unique(idr_predictions$V1)){
  probabilities=idr_predictions[idr_predictions$V1==sequence_name,3:ncol(idr_predictions)]
  
  # remove NAs if there is any
  if(sum(is.na(probabilities))>0){
    probabilities=probabilities[-(which(is.na(probabilities)))]
    
    sequence_length=length(probabilities)
    
    nterminal_length=nterminal[nterminal$Name==sequence_name,2]
    if(length(nterminal_length)>0){
      meanIDR_N_ext=mean(as.numeric(probabilities[1:nterminal_length]))
    }else{
      nterminal_length=0
      meanIDR_N_ext=NA
    }
    
    cterminal_length=cterminal[cterminal$Name==sequence_name,2]
    if(length(cterminal_length)>0){
      meanIDR_C_ext=mean(as.numeric(probabilities[((sequence_length-cterminal_length)+1):sequence_length]))
    }else{
      cterminal_length=0
      meanIDR_C_ext=NA
    }
    
    meanIDR_internal=mean(as.numeric(probabilities[(nterminal_length+1):(sequence_length-cterminal_length)]))
    
    df_mean_idr=rbind(df_mean_idr, c(sequence_name, meanIDR_N_ext, "Nterminal-Extension"), c(sequence_name, meanIDR_C_ext, "Cterminal-Extension"), c(sequence_name, meanIDR_internal, "Internal"))
  } 
  
  else {
    sequence_length=length(probabilities)
    nterminal_length=nterminal[nterminal$Name==sequence_name,2]
    if(length(nterminal_length)>0){
      meanIDR_N_ext=mean(as.numeric(probabilities[1:nterminal_length]))
    }else{
      nterminal_length=0
      meanIDR_N_ext=NA
    }
    
    cterminal_length=cterminal[cterminal$Name==sequence_name,2]
    if(length(cterminal_length)>0){
      meanIDR_C_ext=mean(as.numeric(probabilities[((sequence_length-cterminal_length)+1):sequence_length]))
    }else{
      cterminal_length=0
      meanIDR_C_ext=NA
    }
    
    meanIDR_internal=mean(as.numeric(probabilities[(nterminal_length+1):(sequence_length-cterminal_length)]))
    
    df_mean_idr=rbind(df_mean_idr, c(sequence_name, meanIDR_N_ext, "Nterminal-Extension"), c(sequence_name, meanIDR_C_ext, "Cterminal-Extension"), c(sequence_name, meanIDR_internal, "Internal"))
  }
}

df_mean_idr=as.data.frame(df_mean_idr[-1,])


# if a sequence do not have N-terminal and C-terminal extension remove that sequence from the data
# if a sequence has one of the extensions, remove the absent one from the data
to_remove=c()
for( i in seq(1,nrow(df_mean_idr),3)){
  nstatus=df_mean_idr[i,2]
  cstatus=df_mean_idr[(i+1),2]
  istatus=df_mean_idr[(i+2),2]
  if ((is.na(nstatus)) & (is.na(cstatus))){to_remove=c(to_remove, i, i+1, i+2)}
  else if ((!is.na(nstatus)) & (is.na(cstatus))){to_remove=c(to_remove, i+1)}
  else if ((is.na(nstatus)) & (!is.na(cstatus))){to_remove=c(to_remove, i)}
}

df_mean_idr_no_na=df_mean_idr[(-(to_remove)),]


#plot
#png("/Volumes/bkacar/kacarlab/data/Translation/05_project_extensions/figures/idr_phase_separation/2025-02-05-IF2-meanDisorder-violinboxplot.png", height=1900, width=2000, res=300)
ggplot(df_mean_idr_no_na, aes(x=region,y=as.numeric(mean), fill=region))+
  geom_violin(width=1.2) +
  geom_boxplot(width=0.05, color="gray20", alpha=0.3) +
  scale_fill_manual(values=c("#A7B7D1FF","#D5D2D1FF","#F8B58BFF"))+
  #scale_fill_manual(values=c("#D4C2AD","#669F85","#D7A184"))+
  #scale_fill_manual(values=c("#F7D6D2","#D7D7D9","#4F6F8C"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none", axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x="", y="Mean Disorder")+
  geom_signif(data=df_mean_idr_no_na,comparisons=list(c("Nterminal-Extension","Cterminal-Extension")),map_signif_level = TRUE,annotations="*",
              y_position = 1.0, tip_length = c(0.01, 0.01))+
  geom_signif(data=df_mean_idr_no_na,comparisons=list(c("Nterminal-Extension","Internal")),map_signif_level = TRUE,annotations="*",
              y_position = 0.96, tip_length = c(0.01, 0.01))+
  geom_signif(data=df_mean_idr_no_na,comparisons=list(c("Cterminal-Extension","Internal")),map_signif_level = TRUE,annotations="*",
              y_position = 0.95, tip_length = c(0.01, 0.01))+
  scale_x_discrete(labels=c("N-terminal Extension", "Internal Region", "C-terminal Extension"))+
  ylim(0,1)
#dev.off()

# statistical test between regions
t.test(as.numeric(df_mean_idr_no_na[df_mean_idr_no_na$region=="Nterminal-Extension",2]),as.numeric(df_mean_idr_no_na[df_mean_idr_no_na$region=="Cterminal-Extension",2]))
t.test(as.numeric(df_mean_idr_no_na[df_mean_idr_no_na$region=="Nterminal-Extension",2]),as.numeric(df_mean_idr_no_na[df_mean_idr_no_na$region=="Internal",2]))
t.test(as.numeric(df_mean_idr_no_na[df_mean_idr_no_na$region=="Cterminal-Extension",2]),as.numeric(df_mean_idr_no_na[df_mean_idr_no_na$region=="Internal",2]))

