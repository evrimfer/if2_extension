# 14_nusA-if2.R
# Script for NusA and IF2 extension correlations

load("/Users/evrimfer/Documents/kacar-lab/04_project_extensions/data/Rdata/2025-02-04-domain-genome-gs-gc-ogt-ox-tree.Rdata")

nusa_blast=read.csv("/Users/evrimfer/Documents/kacar-lab/04_project_extensions/results/BLAST/nusA/nusa.blastp-results.csv", header=T, sep="\t")
nusa_archaea_blast=read.csv("/Users/evrimfer/Documents/kacar-lab/04_project_extensions/results/BLAST/nusA/nusA_archaea.blastp-results.csv", header=T, sep="\t")
finalData2=cbind(finalData, rep(NA, nrow(finalData)))

for (i in 1:nrow(finalData2)){
  genome_id=paste(strsplit(finalData2[i,2], "_")[[1]][1],strsplit(finalData2[i,2], "_")[[1]][2], sep="_")
  if(length(grep(genome_id, nusa_blast[,3]))>0){
    print(grep(genome_id, nusa_blast[,3]))
    finalData2[i,24]=paste0(nusa_blast[grep(genome_id, nusa_blast[,3]),3],collapse=",") 
  }else{finalData2[i,24]=NA}
}

for (i in 1:nrow(finalData2)){
  genome_id=paste(strsplit(finalData2[i,2], "_")[[1]][1],strsplit(finalData2[i,2], "_")[[1]][2], sep="_")
  if(is.na(finalData2[i,24])){
    if(length(grep(genome_id, nusa_archaea_blast[,3]))>0){
      print(grep(genome_id, nusa_archaea_blast[,3]))
      finalData2[i,24]=paste0(nusa_archaea_blast[grep(genome_id, nusa_archaea_blast[,3]),3],collapse=",") 
    }
  }
}

finalData3=cbind(finalData2, rep(NA, nrow(finalData2)))

for (i in 1:nrow(finalData3)){
  genome_id=paste(strsplit(finalData3[i,2], "_")[[1]][1],strsplit(finalData3[i,2], "_")[[1]][2], sep="_")
  genome_id=gsub("GCA","GCF", genome_id)
  if(length(grep(genome_id, nusa_blast[,3]))>0){
    print(grep(genome_id, nusa_blast[,3]))
    finalData3[i,25]=paste0(nusa_blast[grep(genome_id, nusa_blast[,3]),4],collapse=",") 
  }else{finalData3[i,25]=NA}
}

for (i in 1:nrow(finalData3)){
  genome_id=paste(strsplit(finalData3[i,2], "_")[[1]][1],strsplit(finalData3[i,2], "_")[[1]][2], sep="_")
  if(is.na(finalData3[i,25])){
    if(length(grep(genome_id, nusa_blast[,3]))>0){
      print(grep(genome_id, nusa_blast[,3]))
      finalData3[i,25]=paste0(nusa_blast[grep(genome_id, nusa_blast[,3]),4],collapse=",") 
    }
  }
}

colnames(finalData3)[24]="nusA"

finalData4=finalData3[(-(which(is.na(finalData3$nusA)))),]

nusA_list=c()
for(i in 1:nrow(finalData4)){
  if(length(grep(",",finalData4$nusA[i]))>0){
    first=strsplit(finalData4$nusA[i],",")[[1]][1]
    second=strsplit(first,"\\|")[[1]][2]
    nusA_list=c(nusA_list, second)
  }
  else{
    second=strsplit(finalData4$nusA[i],"\\|")[[1]][2]
    nusA_list=c(nusA_list, second)
  }
}
nusA_list=as.matrix(nusA_list)


install.packages("rentrez")
library(rentrez)

fasta_seqs <- entrez_fetch(
  db = "protein",
  id = nusA_list[,1],
  rettype = "fasta",
  retmode = "text"
)

nusa_ctd=read.csv("/Users/evrimfer/Documents/kacar-lab/publications/2026-Fer-Extensions/data/nusa/nusa-ctd.csv", header=T)

for(i in 1:nrow(nusa_ctd)){
  if(length(grep(nusa_ctd[i,1],finalData4$nusA))>0){
    finalData4[grep(nusa_ctd[i,1],finalData4$nusA),25]=nusa_ctd[i,2]
  }
}

finalData4B=finalData4[finalData4$domain=="Bacteria",]
cor.test(as.numeric(finalData4B[,25]), as.numeric(finalData4B$L_IF2N))

plot(as.numeric(finalData4B[,25]), as.numeric(finalData4B$L_IF2N))
