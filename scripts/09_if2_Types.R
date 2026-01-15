# 09_if2Types.R
# map different types of IF2 proteins on the tree

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
library(treeio)
library(tidytree)
library(TDbook)
library(ggnewscale)
library(ggtreeExtra)
library(ggstar)
library(cowplot)
library(ggplotify)
library(FSA)

# Step 1: Read phylogeny tree
load("~/rdata/domain-genome-gs-gc-ogt-ox-tree.Rdata")

groups_if2=read.csv("~/structures/if2_clustering.tsv", header=T, sep="\t")
colnames(groups_if2)[2]="Groups_IF2"

# change the order of phyla groups
dat$Phyla=factor(dat$Phyla, levels=c("Acidobacteriota","Actinomycetota","Aquificota","Bacillota","Bacteroidota",
                                     "Cyanobacteriota","Chloroflexota","Deinococcota","Desulfobacterota","Hydrogenedentota", 
                                     "Myxococcota", "Planctomycetota","Pseudomonadota", "Spirochaetota","Verrucomicrobiota",
                                     "Candidatus archaeon","Halobacteriota","Methanobacteriota", "Thermoplasmatota", "Thermoproteota",
                                     "Amoebozoa","Fungi","Metazoa","Sar","Viridiplantae"))

# assign colors for phyla and domain groups
dat2=unique(dat)
dat2$Phyla=factor(dat2$Phyla, levels=c("Actinomycetota","Cyanobacteriota","Chloroflexota","Bacillota","Deinococcota",
                                       "Bacteroidota","Planctomycetota","Verrucomicrobiota","Aquificota","Pseudomonadota",
                                       "Myxococcota","Desulfobacterota","Acidobacteriota","Hydrogenedentota","Spirochaetota",
                                       "Halobacteriota","Thermoplasmatota","Thermoproteota","Amoebozoa","Fungi",
                                       "Metazoa","Sar","Viridiplantae","Methanobacteriota","Candidatus archaeon"))

# plot tree and IF2 types
p <- ggtree(tree, layout="circular", size=0.4)+  xlim(-1,NA)+ ylim(-1000, NA)+
  geom_treescale(fontsize=1, linesize=0.1, offset=1, x=0, y=45)

p1=p %<+% dat+
  geom_highlight(data=dat2, 
                 aes(node=node, fill=Phyla),
                 alpha=0.5,
                 align="right",
                 extend=0.8,
                 show.legend=FALSE) +
  scale_fill_manual(values=c("#f5f5f5","#d1d3d4","#f5f5f5","#d1d3d4","#f5f5f5","#d1d3d4",
                             "#f5f5f5","#d1d3d4","#f5f5f5","#d1d3d4","#f5f5f5","#d1d3d4",
                             "#f5f5f5","#d1d3d4","#f5f5f5","#d1d3d4","#f5f5f5","#d1d3d4",
                             "#f5f5f5","#f5f5f5","#f5f5f5","#f5f5f5","#f5f5f5","#d1d3d4","#f5f5f5"))+
  new_scale_fill() +
  geom_fruit(data=groups_if2, geom=geom_tile, mapping=aes(x=Types, y=Name, fill=Groups_IF2), offset=0.1,  pwidth=0.5, size=1)+
  scale_fill_manual(values=c("#729CB2","#B2E0E1","#A89ECA","#573B8B","#E39E85","#D8AD4E","#293477"))

rotate_tree(p1,180)   


# plot tree and properties
p <- ggtree(tree, size=0.4)+
  geom_treescale(fontsize=1, linesize=0.1, offset=1, x=0, y=45)

p1=p %<+% dat+
  geom_highlight(data=dat2, 
                 aes(node=node, fill=Phyla),
                 alpha=0.5,
                 align="right",
                 extend=0.5,
                 show.legend=FALSE) +
  scale_fill_manual(values=c("#f5f5f5","#d1d3d4","#f5f5f5","#d1d3d4","#f5f5f5","#d1d3d4",
                             "#f5f5f5","#d1d3d4","#f5f5f5","#d1d3d4","#f5f5f5","#d1d3d4",
                             "#f5f5f5","#d1d3d4","#f5f5f5","#d1d3d4","#f5f5f5","#d1d3d4",
                             "#f5f5f5","#f5f5f5","#f5f5f5","#f5f5f5","#f5f5f5","#d1d3d4","#f5f5f5"))+
  # if2 whole sequence
  new_scale_fill() +
  geom_fruit(data=finalData,geom=geom_tile,mapping=aes(x=IF2, y=genome_name, fill=as.numeric(L_IF2)), offset=0.05, pwidth=0.2, size=0.1)+
  scale_fill_gradient(low="#E9EDC9",high="#6B705C", na.value="white")+
  
  # if2 nterminal extension
  new_scale_fill() +
  geom_fruit(data=finalData,geom=geom_tile,mapping=aes(x=IF2_N, y=genome_name, fill=as.numeric(L_IF2N)), offset=0.025, pwidth=0.2, size=0.1)+
  scale_fill_gradient(low="#E9EDC9",high="#6B705C", na.value="white")+
  
  # if2 cterminal extension
  new_scale_fill() +
  geom_fruit(data=finalData,geom=geom_tile,mapping=aes(x=IF2_C, y=genome_name, fill=as.numeric(L_IF2C)), offset=0.025, pwidth=0.2, size=0.1)+
  scale_fill_gradient(low="#E9EDC9",high="#6B705C", na.value="white")+
  
  new_scale_fill() +
  geom_fruit(data=finalData,geom=geom_tile,mapping=aes(x=Temperature, y=genome_name, fill=as.numeric(predicted_ogt)), offset=0.06,  pwidth=0.3, size=0.1)+
  #scale_fill_gradient(low="#D98594",high="#863339", na.value="white")+
  scale_fill_gradientn(colors=c("#FFFDD0","#F4A464FF","#65000B"), na.value="white",breaks=c(15,58.5,100),labels=c(15,60,100), limits=c(15,100))+

  # oxygen tolerance
  new_scale_fill() +
  geom_fruit(data=finalData, geom=geom_tile,mapping=aes(x=Ox_Tolerance, y=genome_name, fill=factor(tolerance)), offset=0.05,  pwidth=0.3, size=0.2)+
  scale_fill_manual(values=c("aerobic"="#C0C0C0","anaerobic"="#3C4C81", "facultative"="#848482"))+
  theme(legend.position="none")


p1
 
  


