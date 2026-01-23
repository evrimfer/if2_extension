# 01_gtdbFiltering.R
# this script is used for curating a dataset of genomes from GTDB metadata 

# Read the metadata
gtdbBacteria=read.csv("~/genomes/gtdb/bac120_metadata.tsv", sep="\t", header=T)
gtdbArchaea=read.csv("~/genomes/gtdb/ar53_metadata.tsv", sep="\t", header=T)
gtdbAll=rbind(gtdbArchaea,gtdbBacteria)

# filter by completeness and contamination
gtdb_filtered=gtdbAll[((gtdbAll$checkm_completeness>95)&(gtdbAll$checkm_contamination<5)&(gtdbAll$checkm2_completeness>95)&(gtdbAll$checkm2_contamination<5)),]
# keep complete and scaffold genomes only
gtdb_filtered=gtdb_filtered[(gtdb_filtered$ncbi_assembly_level=="Complete Genome")|(gtdb_filtered$ncbi_assembly_level=="Scaffold"),]

# keep useful columns
colnames(gtdb_filtered)
gtdb_filtered=gtdb_filtered[,c(18,1,58,20,82,83,66,77,49,17,15,16,67,6,7,3,4,8,10,11,14)]

# separate taxonomy 
gtdb_superkingdom=rep(NA, nrow(gtdb_filtered))
gtdb_phylum=rep(NA, nrow(gtdb_filtered))
gtdb_class=rep(NA, nrow(gtdb_filtered))
gtdb_order=rep(NA, nrow(gtdb_filtered))
gtdb_family=rep(NA, nrow(gtdb_filtered))
gtdb_genus=rep(NA, nrow(gtdb_filtered))
gtdb_species=rep(NA, nrow(gtdb_filtered))

gtdb_filtered=cbind(gtdb_filtered[,1:4],gtdb_superkingdom,gtdb_phylum,gtdb_class,gtdb_order,gtdb_family,gtdb_genus,gtdb_species, gtdb_filtered[,5:21])

for(i in 1:nrow(gtdb_filtered)){
  taxonomy=gtdb_filtered$gtdb_taxonomy[i]
  if(length(grep("d__",taxonomy))>0){
    gtdb_filtered$gtdb_superkingdom[i]=strsplit(strsplit(taxonomy, "d__")[[1]][2],"\\;")[[1]][1]
  }
  if(length(grep("p__",taxonomy))>0){
    gtdb_filtered$gtdb_phylum[i]=strsplit(strsplit(taxonomy, "p__")[[1]][2],"\\;")[[1]][1]
  }
  if(length(grep("c__",taxonomy))>0){
    gtdb_filtered$gtdb_class[i]=strsplit(strsplit(taxonomy, "c__")[[1]][2],"\\;")[[1]][1]
  }
  if(length(grep("o__",taxonomy))>0){
    gtdb_filtered$gtdb_order[i]=strsplit(strsplit(taxonomy, "o__")[[1]][2],"\\;")[[1]][1]
  }
  if(length(grep("f__",taxonomy))>0){
    gtdb_filtered$gtdb_family[i]=strsplit(strsplit(taxonomy, "f__")[[1]][2],"\\;")[[1]][1]
  }
  if(length(grep("g__",taxonomy))>0){
    gtdb_filtered$gtdb_genus[i]=strsplit(strsplit(taxonomy, "g__")[[1]][2],"\\;")[[1]][1]
  }
  if(length(grep("s__",taxonomy))>0){
    gtdb_filtered$gtdb_species[i]=strsplit(taxonomy, "s__")[[1]][2]
  }
}

# remove duplicated genome accessions
gtdb_filtered=gtdb_filtered[-(which(duplicated(gtdb_filtered$gtdb_genome_representative))),]
for (i in 1:nrow(gtdb_filtered)){
  gtdb_filtered$gtdb_genome_representative[i]=paste(strsplit(gtdb_filtered$gtdb_genome_representative[i],"_")[[1]][2],strsplit(gtdb_filtered$gtdb_genome_representative[i],"_")[[1]][3],sep="_")
}

accessions=as.matrix(gtdb_filtered[,1])
GCA=c()
GCF=c()
for(i in accessions){
  if(length(grep("GCA_",i))>0){
    GCA=c(GCA,i)
  }else{GCF=c(GCF,i)}
}

GCA=as.matrix(GCA)
GCF=as.matrix(GCF)

# write the final data
write.table(gtdb_filtered,"~/genomes/gtdb/GTDB_filtered_metadata.csv", col.names=T, row.names=F, quote=F, sep=",")
write.table(gtdb_filtered[,1],"~/genomes/gtdb/GTDB_filtered_accessions.txt", col.names=F, row.names=F, quote=F)


sum(gtdb_filtered$gtdb_superkingdom=="Archaea")
sum(gtdb_filtered$gtdb_superkingdom=="Bacteria")

