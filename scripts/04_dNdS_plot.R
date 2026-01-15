# 04_dNdSplot.R
# this script is used to visualize dN/dS results

# load libraries
library(ggplot2)
library(reshape2)
library(seqinr)

# load the data - change the dNdS result file name as necessary
rst=as.matrix(readLines("~/selection/bacterial/bacteria_if2_acidobacteria_m1.rst"))
rst=rst[12:(nrow(rst)-4),1]

# convert the data into matrix 
rates=matrix(ncol=5)
colnames(rates)=c("position","residue","pp1","pp2","meanw")
for(i in 1:length(rst)){
  if(i<10){
    position=strsplit(rst[i]," ")[[1]][4]
    residue=strsplit(rst[i]," ")[[1]][5]
    pp1=strsplit(rst[i]," ")[[1]][8]
    pp2=strsplit(rst[i]," ")[[1]][9]
    meanw=strsplit(rst[i]," ")[[1]][13]
    rates=rbind(rates, c(position, residue, pp1, pp2, meanw))
  }else if(i>=10 & i<100){
    position=strsplit(rst[i]," ")[[1]][3]
    residue=strsplit(rst[i]," ")[[1]][4]
    pp1=strsplit(rst[i]," ")[[1]][7]
    pp2=strsplit(rst[i]," ")[[1]][8]
    meanw=strsplit(rst[i]," ")[[1]][12]
    rates=rbind(rates, c(position, residue, pp1, pp2, meanw))
  }else if(i>=100 & i<1000){
    position=strsplit(rst[i]," ")[[1]][2]
    residue=strsplit(rst[i]," ")[[1]][3]
    pp1=strsplit(rst[i]," ")[[1]][6]
    pp2=strsplit(rst[i]," ")[[1]][7]
    meanw=strsplit(rst[i]," ")[[1]][11]
    rates=rbind(rates, c(position, residue, pp1, pp2, meanw))
  }else if(i>=1000){
    position=strsplit(rst[i]," ")[[1]][1]
    residue=strsplit(rst[i]," ")[[1]][2]
    pp1=strsplit(rst[i]," ")[[1]][5]
    pp2=strsplit(rst[i]," ")[[1]][6]
    meanw=strsplit(rst[i]," ")[[1]][10]
    rates=rbind(rates, c(position, residue, pp1, pp2, meanw))
  }
}
rates=as.data.frame(rates[-1,])

# remove the gap regions
rates=rates[-which(rates$residue=="-"),]

# load the protein alignment extract for the E. coli and selection reference sequence - change the alignment file name as necessary
subset_aln=read.fasta("~/selection/dn-ds-representatives-and-ecoli.fasta",seqtype="AA")

# go over the alignment, when you see a residue check rate
# white for gaps, gray under negative selection, red for extension under neutral evolution
# black for conserved regions (does not matter if under selection or neutral)

# change the representative species name in the "subset_aln" each time based on the dNds-phyla result you loaded

rate_color=c()
count=1
for(i in 1:1349){
  if(subset_aln$GCA_026000775.1_Archaea_Thermoproteota_Vulcanisaeta_souniana_JCM_11219[i]=="-"){rate_color=c(rate_color,"white")}
  else if(subset_aln$GCA_026000775.1_Archaea_Thermoproteota_Vulcanisaeta_souniana_JCM_11219[i]!="-" & (count <= 9 | count >= 569)){
    if (rates$meanw[count]<1){
      rate_color=c(rate_color,"gray")
      count=count+1
    }else{
      rate_color=c(rate_color,"red")
      count=count+1
      }
  }
  else if(subset_aln$GCA_026000775.1_Archaea_Thermoproteota_Vulcanisaeta_souniana_JCM_11219[i]!="-" & (count > 9 | count < 569)){
    if (rates$meanw[count]<1){
      rate_color=c(rate_color,"black")
      count=count+1
    }else{
      rate_color=c(rate_color,"black")
      count=count+1
    }
  }
}

# Plot
df <- data.frame(x = seq_along(rate_color), y = 1, col = rate_color)

ggplot(df, aes(x, y, fill = col)) +
  geom_raster() +
  scale_fill_identity() +
  theme_void() +
  theme(aspect.ratio = 0.02)

# 
# arc_thermo_cterm=rates[569:598,]
# arc_halocterm=rates[562:596,]
# arc_methano_cterm=rates[557:596,]
# euk_cterm=rates[1017:1045,]
# arc=rbind(arc_thermo_cterm,arc_methano_cterm, arc_halocterm,euk_cterm)
# 
# mean(as.numeric(arc$meanw))
# sd(as.numeric(arc$meanw))
# 
# 
# bac_acido_nterm=rates[1:354,]
# bac_actino_nterm=rates[1:480,]
# bac_aqui_nterm=rates[1:244,]
# bac_bacteroi_nterm=rates[1:652,]
# bac_plancto_nterm=rates[1:262,]
# bac_pseudo_nterm=rates[1:461,]
# bac_bacillo_nterm=rates[1:421,]
# bac_cyano_nterm=rates[1:540,]
# 
# bac=rbind(bac_acido_nterm,bac_actino_nterm,bac_aqui_nterm,bac_bacteroi_nterm,bac_plancto_nterm,bac_pseudo_nterm,bac_bacillo_nterm,bac_cyano_nterm)
# mean(as.numeric(bac$meanw))
# sd(as.numeric(bac$meanw))
# hist(as.numeric(bac$meanw))
