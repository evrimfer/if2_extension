# 06_idr_Composition.R
# order promoting and disorder promoting amino acid composition in IF2 extensions

# Load libraries

library(ggplot2)
library(seqinr)

# Read alignment file
# for IF2 the boundaries of extensions are N:1-3531, internal:3532-4763, C:4764-4815 on the alignment
alignment=read.fasta("~/phase_separation/IF2-GCF-aligned.fasta", seqtype="AA", set.attributes = F)

df=matrix(ncol=3)
colnames(df)=c("residue","region","mean_freq")

aminoacid_MeanFreq=function(alignment, aminoacid, df){
  all_freqN=c()
  all_freqI=c()
  all_freqC=c()
  for(i in 1:length(names(alignment))){
    nterminal_extension=alignment[[i]][1:3531]
    internal=alignment[[i]][3532:4763]
    cterminal_extension=alignment[[i]][4764:4815]
    
    nterminal_extension_nogaps=nterminal_extension[-(which(nterminal_extension=="-"))]
    internal_nogaps=internal[-(which(internal=="-"))]
    cterminal_extension_nogaps=cterminal_extension[-(which(cterminal_extension=="-"))]
    
    freqN=((sum(nterminal_extension_nogaps==aminoacid))/(length(nterminal_extension_nogaps)))*100
    freqI=((sum(internal_nogaps==aminoacid))/(length(internal_nogaps)))*100
    freqC=((sum(cterminal_extension_nogaps==aminoacid))/(length(cterminal_extension_nogaps)))*100
    
    all_freqN=c(all_freqN, freqN)
    all_freqI=c(all_freqI, freqI)
    all_freqC=c(all_freqC, freqC)
  }
  
  residue=aminoacid
  region="NTD"
  mean_freq=mean(all_freqN, na.rm=T)
  df=rbind(df, c(residue,region, mean_freq))
  
  residue=aminoacid
  region="Internal"
  mean_freq=mean(all_freqI, na.rm=T)
  df=rbind(df, c(residue,region, mean_freq))
  
  residue=aminoacid
  region="CTD"
  mean_freq=mean(all_freqC, na.rm=T)
  df=rbind(df, c(residue,region, mean_freq))
  
  return(df)
}

df=aminoacid_MeanFreq(alignment, "C", df)
df=aminoacid_MeanFreq(alignment, "W", df)
df=aminoacid_MeanFreq(alignment, "Y", df)
df=aminoacid_MeanFreq(alignment, "F", df)
df=aminoacid_MeanFreq(alignment, "I", df)
df=aminoacid_MeanFreq(alignment, "L", df)
df=aminoacid_MeanFreq(alignment, "H", df)
df=aminoacid_MeanFreq(alignment, "V", df)
df=aminoacid_MeanFreq(alignment, "N", df)
df=aminoacid_MeanFreq(alignment, "M", df)
df=aminoacid_MeanFreq(alignment, "T", df)
df=aminoacid_MeanFreq(alignment, "D", df)
df=aminoacid_MeanFreq(alignment, "K", df)
df=aminoacid_MeanFreq(alignment, "E", df)
df=aminoacid_MeanFreq(alignment, "Q", df)
df=aminoacid_MeanFreq(alignment, "S", df)
df=aminoacid_MeanFreq(alignment, "P", df)
df=aminoacid_MeanFreq(alignment, "R", df)
df=aminoacid_MeanFreq(alignment, "G", df)
df=aminoacid_MeanFreq(alignment, "A", df)

df=as.data.frame(df[-1,])

df$residue=factor(df$residue, levels=c("C","W","Y","F","I","L","H","V","N","M","K","E","Q","S","P","R","G","A","T","D"))
df$region=factor(df$region,levels=c("NTD","Internal","CTD"))

df[df$residue=="C",]

NTD=c()
for(i in unique(df$residue)){
  NTD=c(NTD, (as.numeric(df[df$residue==i,][1,3])-as.numeric(df[df$residue==i,][2,3]))/as.numeric(df[df$residue==i,][2,3]))
}

CTD=c()
for(i in unique(df$residue)){
  CTD=c(CTD,(as.numeric(df[df$residue==i,][3,3])-as.numeric(df[df$residue==i,][2,3]))/as.numeric(df[df$residue==i,][2,3]))
}


NTD_df=as.data.frame(cbind(NTD,levels(df$residue)))
colnames(NTD_df)=c("Fold","Residue")
NTD_df$Residue=factor(NTD_df$Residue, levels=c("C","W","Y","F","I","L","H","V","N","M","K","E","Q","S","P","R","G","A","T","D"))
CTD_df=as.data.frame(cbind(CTD,levels(df$residue)))
colnames(CTD_df)=c("Fold","Residue")
CTD_df$Residue=factor(CTD_df$Residue, levels=c("C","W","Y","F","I","L","H","V","N","M","K","E","Q","S","P","R","G","A","T","D"))

ggplot(NTD_df, aes(x=Residue, y=as.numeric(Fold))) + 
  geom_bar(position="dodge", stat="identity",  fill="#acb8d0")+
  geom_hline(aes(yintercept = 0))+
  theme_bw()+
  ylim(-1,8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none", axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x=" ", y="Fold Mean Frequency")

ggplot(CTD_df, aes(x=Residue, y=as.numeric(Fold))) + 
  geom_bar(position="dodge", stat="identity", fill="#efb990")+
  geom_hline(aes(yintercept = 0))+
  theme_bw()+
  ylim(-1,8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none", axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x=" ", y="Fold Mean Frequency")
  