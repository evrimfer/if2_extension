# 09_meanDP.R
# script for Phase Separation

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

# Step 1: Read phylogeny tree
tree=read.tree("~/data/species_tree/Species_tree_1003_universal_rproteins_aligned.fasta.treefile.rooted") 

# Step 2: Read alignment file
alignment=read.fasta("~/data/phase_separation/IF2-GCF-aligned.fasta", seqtype="AA",  set.attributes = F)

# Step 3: Load the data frame coming from the python script
df=read.csv("~/data/phase_separation/IF2_combined_probabilities_with_gaps.tsv", header=F, sep="\t")
colnames(df)=c("sequence","position","residue","probability","sbind")

# for IF1 the boundaries of extensions are N:1-125, internal:126-201, C:202-295 on the alignment
# for IF2 the boundaries of extensions are N:1-3531, internal:3532-4763, C:4764-4815 on the alignment
# for EFTU the boundaries of extensions are N:1-162, internal:163-692, C:693-711 on the alignment
# for EFG the boundaries of extensions are N:1-444, internal:445-1820, C:1821-1842 on the alignment
df_mean_region=matrix(ncol=3)
colnames(df_mean_region)=c("sequence_name","mean","region")
for(sequence_name in unique(df$sequence)){
  subdata=df[df$sequence==sequence_name,]
  nterminal_extension=subdata[1:3531,]
  internal=subdata[3532:4763,]
  cterminal_extension=subdata[4764:4815,]
  
  nterminal_extension_nogaps=nterminal_extension[-(which(nterminal_extension$residue=="-")),]
  internal_nogaps=internal[-(which(internal$residue=="-")),]
  cterminal_extension_nogaps=cterminal_extension[-(which(cterminal_extension$residue=="-")),]
  
  meanDP_N_ext=mean(as.numeric(nterminal_extension_nogaps[,4]))
  meanDP_internal=mean(as.numeric(internal_nogaps[,4]))
  meanDP_C_ext=mean(as.numeric(cterminal_extension_nogaps[,4]))
  
  df_mean_region=rbind(df_mean_region,c(sequence_name, meanDP_N_ext, "Nterminal-Extension"),
                       c(sequence_name, meanDP_internal, "Internal"),
                       c(sequence_name, meanDP_C_ext, "Cterminal-Extension"))
}

df_mean_region=as.data.frame(df_mean_region[-1,])
df_mean_region$region=factor(df_mean_region$region, levels=c("Nterminal-Extension", "Internal", "Cterminal-Extension"))

# Plot
#png("/Volumes/bkacar/kacarlab/data/Translation/05_project_extensions/figures/idr_phase_separation/2025-02-05-IF1-meanDP-violinboxplot.png", height=1900, width=2000, res=300)
ggplot(df_mean_region, aes(x=region,y=as.numeric(mean), fill=region))+
  geom_violin(width=1.2) +
  geom_boxplot(width=0.05, color="gray20", alpha=0.3) +
  scale_fill_manual(values=c("#A7B7D1FF","#D5D2D1FF","#F8B58BFF"))+
  #scale_fill_manual(values=c("#D4C2AD","#669F85","#D7A184"))+
  #scale_fill_manual(values=c("#F7D6D2","#D7D7D9","#4F6F8C"))+
  scale_x_discrete(labels=c("N-terminal Extension", "Internal Region", "C-terminal Extension"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none", axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x="", y="MeanDP")+
  ylim(0,1)
#dev.off()

# statistical test between regions
t.test(as.numeric(df_mean_region[df_mean_region$region=="Nterminal-Extension",2]),as.numeric(df_mean_region[df_mean_region$region=="Cterminal-Extension",2]))
t.test(as.numeric(df_mean_region[df_mean_region$region=="Nterminal-Extension",2]),as.numeric(df_mean_region[df_mean_region$region=="Internal",2]))
t.test(as.numeric(df_mean_region[df_mean_region$region=="Cterminal-Extension",2]),as.numeric(df_mean_region[df_mean_region$region=="Internal",2]))

