# 03_tree_length_properties.R
# plot tree with sequence length and predicted features

# load libraries
library(ggtree)
library(treeio)
library(tidytree)
library(ggplot2)
library(TDbook)
library(ggnewscale)
library(dplyr)
library(ggtreeExtra)
library(ggstar)
library(cowplot)
library(ggplotify)
library(FSA)

### Plot
# can start from here 

load("~/rdata/domain-genome-gs-gc-ogt-ox-tree.Rdata")

# change the order of phyla groups
dat$Phyla=factor(dat$Phyla, levels=c("Acidobacteriota","Actinomycetota","Aquificota","Bacillota","Bacteroidota",
                                     "Cyanobacteriota","Chloroflexota","Deinococcota","Desulfobacterota","Hydrogenedentota", 
                                     "Myxococcota", "Planctomycetota","Pseudomonadota", "Spirochaetota","Verrucomicrobiota",
                                     "Candidatus archaeon","Halobacteriota","Methanobacteriota", "Thermoplasmatota", "Thermoproteota",
                                     "Amoebozoa","Fungi","Metazoa","Sar","Viridiplantae"))

# assign colors for phyla groups

phyla_colors=c("#F05B43FF","#E78429FF","#F78462FF","#F9D14AFF","#F7DEA3FF",
               "#C7C45EFF","#B1A866FF","#FEAC81FF","#FFC6C4FF","#F4A3A8FF",
               "#E38191FF","#AD466CFF","#CC607DFF","#8B3058FF","#672044FF",
               "#B4D9CCFF","#89C0B6FF","#63A6A0FF","#448C8AFF","#287274FF",
               "#404058FF","#484868FF","#5870A0FF","#88B0E0FF","#3868A0FF")

# plot tree and properties
p <- ggtree(tree, layout="rectangular", size=0.3) 
p %<+% dat+geom_tippoint(aes(color=Phyla))+
  scale_color_manual(values=phyla_colors)+
  
  # if1 whole sequence
  new_scale_fill() +
  geom_fruit(data=finalData, geom=geom_tile, mapping=aes(x=IF1, y=genome_name, fill=as.numeric(L_IF1)), offset=0.05,  pwidth=0.2, size=0.1)+
  scale_fill_gradient(low="#DEC5DAFF",high="#574571FF", na.value="white")+
  
  # if1 nterminal extension
  new_scale_fill() +
  geom_fruit(data=finalData,geom=geom_tile,mapping=aes(x=IF1_N, y=genome_name, fill=as.numeric(L_IF1N)), offset=0.025,  pwidth=0.2, size=0.1)+
  scale_fill_gradient(low="#DEC5DAFF",high="#574571FF", na.value="white")+
  
  # if1 cterminal extension
  new_scale_fill() +
  geom_fruit(data=finalData,geom=geom_tile,mapping=aes(x=IF1_C, y=genome_name, fill=as.numeric(L_IF1C)), offset=0.025,  pwidth=0.2, size=0.1)+
  scale_fill_gradient(low="#DEC5DAFF",high="#574571FF", na.value="white")+
  
  # if2 whole sequence
  new_scale_fill() +
  geom_fruit(data=finalData,geom=geom_tile,mapping=aes(x=IF2, y=genome_name, fill=as.numeric(L_IF2)), offset=0.05, pwidth=0.2, size=0.1)+
  scale_fill_gradient(low="#85C4C9FF",high="#3B738FFF", na.value="white")+
  
  # if2 nterminal extension
  new_scale_fill() +
  geom_fruit(data=finalData,geom=geom_tile,mapping=aes(x=IF2_N, y=genome_name, fill=as.numeric(L_IF2N)), offset=0.025, pwidth=0.2, size=0.1)+
  scale_fill_gradient(low="#85C4C9FF",high="#3B738FFF", na.value="white")+
  
  # if2 cterminal extension
  new_scale_fill() +
  geom_fruit(data=finalData,geom=geom_tile,mapping=aes(x=IF2_C, y=genome_name, fill=as.numeric(L_IF2C)), offset=0.025, pwidth=0.2, size=0.1)+
  scale_fill_gradient(low="#85C4C9FF",high="#3B738FFF", na.value="white")+
  
  # genome size
  new_scale_fill() +
  geom_fruit(data=finalData,geom=geom_tile,mapping=aes(x=Size, y=genome_name, fill=log(as.numeric(genome_size))), offset=0.05,  pwidth=0.2, size=0.1)+
  scale_fill_gradient(low="#E0CEBE", high="#BF8A7C", na.value="white")+
  
  # genomic GC content
  new_scale_fill() +
  geom_fruit(data=finalData,geom=geom_tile,mapping=aes(x=`GC%`, y=genome_name, fill=as.numeric(gc_percentage)), offset=0.025,  pwidth=0.2, size=0.1)+
  scale_fill_gradient(low="#D0E8F8FF",high="#6090C8FF", na.value="white")+ 
  
  # predicted temperature
  new_scale_fill() +
  geom_fruit(data=finalData,geom=geom_tile,mapping=aes(x=Temperature, y=genome_name, fill=as.numeric(predicted_ogt)), offset=0.035,  pwidth=0.2, size=0.1)+
  scale_fill_gradient(low="#F2C88FFF",high="#B9563FFF", na.value="white")+ 
  
  # oxygen tolerance
  new_scale_fill() +
  geom_fruit(data=finalData, geom=geom_tile,mapping=aes(x=Ox_Tolerance, y=genome_name, fill=factor(tolerance)), offset=0.025,  pwidth=0.2, size=0.2)+
  scale_fill_manual(values=c("aerobic"="#C0C0C0","anaerobic"="#4A5D74", "facultative"="#949FA6"))+
  
  theme(legend.position="none")
  

# load the data having domain, genome, genome size, gc-content, predicted OGT, oxygen enzyme number, oxygen tolerance and tree data
load("~/rdata/domain-genome-gs-gc-ogt-ox-tree.Rdata")

# read lengths of whole sequence length, N-terminal extension and C-terminal extension of IF2
IF2_WH=as.data.frame(read.csv("~/lengths/IF2_properties.tsv", sep="\t", header=T))
IF2_NTD=as.data.frame(read.csv("~/lengths/IF2_Nterm_properties.tsv", sep="\t", header=T))
IF2_CTD=as.data.frame(read.csv("~/lengths/IF2_Cterm_properties.tsv", sep="\t", header=T))

# function to add the sequence length to the properties data
add_sequence_lengths=function(subData, factor_lengths, protein_name){
  to_add=c()
  
  for(genome in subData[,2]){
    genome_id=paste(strsplit(genome,"_")[[1]][1],strsplit(genome,"_")[[1]][2], sep="_")
    
    if(length(grep(genome_id, factor_lengths[,1]))>0){to_add=c(to_add, factor_lengths[grep(genome_id, factor_lengths[,1]),2])}
    else if(length(grep(genome_id, factor_lengths[,1]))==0){to_add=c(to_add, 0)}
  }
  
  subData2=cbind(subData, to_add, rep(protein_name, nrow(subData)))
  return(subData2)
}

# run the function for IF2 (whole, N-term extensions, C-term extensions)
all_data4=add_sequence_lengths(all_data3, IF2_WH$Sequence_Length, "IF2")
all_data5=add_sequence_lengths(all_data4, IF2_NTD$Sequence_Length, "IF2_N")
all_data6=add_sequence_lengths(all_data5, IF2_CTD$Sequence_Length, "IF2_C")

# rename column names
finalData=as.data.frame(all_data6)
colnames(finalData)=c(colnames(subData), "L_IF2","IF2","L_IF2N","IF2_N","L_IF2C","IF2_C")

# add new columns for property names
finalData=cbind(finalData, rep("genome_size", nrow(finalData)), rep("gc_percentage", nrow(finalData)),rep("Temperature", nrow(finalData)), rep("oxygen_tolerance", nrow(finalData)))
colnames(finalData)[14:17]=c("Size", "GC%","Temperature","Ox_Tolerance")



                 



