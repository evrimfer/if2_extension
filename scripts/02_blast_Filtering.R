# 02_blastFiltering.R
# this script is used for filtering the BLAST results

# read the blast results
proteins=c("IF2","L1","L2","L3","L4","L5","L6","L10","L11","L13","L14","L15","L16","L18","L22","L23","L24","L29","L30",
           "S2","S4","S5","S7","S8","S9","S10","S11","S12","S13","S14","S15","S17","S19")

for (protein in proteins) {
  df=read.csv(paste("~/blastp/",protein,".blastp-results.csv", sep=""), header=T, sep="\t")
  assign(protein, df)
}

# separate genome accession and protein accession
separateColumns=function(rp){
  for(i in 1:nrow(rp)){
    genome=strsplit(rp[i,3],"\\|")[[1]][1]
    protein=strsplit(rp[i,3],"\\|")[[1]][2]
    rp[i,1]=genome
    rp[i,3]=protein
    
    description=strsplit(rp[i,5],"\\[")[[1]][1]
    organism=strsplit(rp[i,5],"\\[")[[1]][2]
    rp[i,5]=description
    rp[i,6]=organism
  }
  return(rp)
}

for (protein in proteins) {
  df <- get(protein)
  df <- separateColumns(df)
  assign(protein, df)
}

# only keep target protein names
unique(IF2$description)
IF2=IF2[c(grep("IF-2",IF2$description),grep("IF-5B",IF2$description),grep("actor 2",IF2$description)),]
IF2=IF2[-(c(grep("partial",IF2$description),grep("domain-containing",IF2$description),grep("probable",IF2$description),grep("putative",IF2$description),grep("unclassified",IF2$description),grep("UNVERIFIED",IF2$description))),]
IF2=IF2[which(IF2$length>=500),]

unique(L1$description)
L1=L1[c(grep("L1",L1$description),grep("plA",L1$description),grep("pl1",L1$description)),]
L1=L1[-(c(grep("partial",L1$description),grep("domain-containing",L1$description),grep("probable",L1$description),grep("putative",L1$description),grep("unclassified",L1$description),grep("UNVERIFIED",L1$description))),]

unique(L2$description)
L2=L2[c(grep("L2",L2$description),grep("plB",L2$description)),]
L2=L2[-(c(grep("partial",L2$description),grep("UNVERIFIED",L2$description))),]
          
unique(L3$description)
L3=L3[c(grep("L3",L3$description),grep("plC",L3$description)),]
L3=L3[-(c(grep("partial",L3$description),grep("UNVERIFIED",L3$description),grep("hypothetical",L3$description))),]

unique(L4$description)
L4=L4[c(grep("L4",L4$description),grep("plD",L4$description)),]
L4=L4[-(c(grep("partial",L4$description),grep("UNVERIFIED",L4$description),grep("regulates expression",L4$description))),]

unique(L5$description)
L5=L5[c(grep("L5",L5$description),grep("plE",L5$description)),]
L5=L5[-(c(grep("partial",L5$description),grep("UNVERIFIED",L5$description),grep("C-terminal",L5$description))),]

unique(L6$description)
L6=L6[c(grep("L6",L6$description),grep("plF",L6$description)),]
L6=L6[-(c(grep("partial",L6$description),grep("UNVERIFIED",L6$description))),]

unique(L10$description)
L10=L10[c(grep("L10",L10$description),grep("plJ",L10$description)),]
L10=L10[-(c(grep("partial",L10$description),grep("UNVERIFIED",L10$description),grep("Fused",L10$description))),]

unique(L11$description)
L11=L11[c(grep("L11",L11$description),grep("plK",L11$description)),]
L11=L11[-(c(grep("partial",L11$description),grep("UNVERIFIED",L11$description))),]

unique(L13$description)
L13=L13[c(grep("L13",L13$description),grep("plM",L13$description)),]
L13=L13[-(c(grep("partial",L13$description),grep("UNVERIFIED",L13$description),grep("modular",L13$description),grep("fragment",L13$description))),]

unique(L14$description)
L14=L14[c(grep("L14",L14$description),grep("plN",L14$description)),]
L14=L14[-(c(grep("partial",L14$description),grep("UNVERIFIED",L14$description))),]

unique(L15$description)
L15=L15[c(grep("L15",L15$description),grep("plO",L15$description)),]
L15=L15[-(c(grep("partial",L15$description),grep("UNVERIFIED",L15$description))),]

unique(L16$description)
L16=L16[c(grep("L16",L16$description),grep("plP",L16$description)),]
L16=L16[-(c(grep("partial",L16$description),grep("UNVERIFIED",L16$description))),]

unique(L18$description)
L18=L18[c(grep("L18",L18$description),grep("plR",L18$description)),]
L18=L18[-(c(grep("partial",L18$description),grep("UNVERIFIED",L18$description))),]

unique(L22$description)
L22=L22[c(grep("L22",L22$description),grep("plV",L22$description)),]
L22=L22[-(c(grep("partial",L22$description),grep("UNVERIFIED",L22$description))),]

unique(L23$description)
L23=L23[c(grep("L23",L23$description),grep("plW",L23$description)),]
L23=L23[-(c(grep("partial",L23$description),grep("UNVERIFIED",L23$description))),]

unique(L24$description)
L24=L24[c(grep("L24",L24$description),grep("plX",L24$description)),]
L24=L24[-(c(grep("partial",L24$description),grep("UNVERIFIED",L24$description))),]

unique(L29$description)
L29=L29[c(grep("L29",L29$description),grep("pmC",L29$description)),]
L29=L29[-(c(grep("partial",L29$description),grep("UNVERIFIED",L29$description))),]

unique(L30$description)
L30=L30[c(grep("L30",L30$description),grep("pmD",L30$description)),]
L30=L30[-(c(grep("partial",L30$description),grep("UNVERIFIED",L30$description))),]

unique(S2$description)
S2=S2[c(grep("S2",S2$description),grep("psB",S2$description)),]
S2=S2[-(c(grep("partial",S2$description),grep("UNVERIFIED",S2$description),grep("BS2a",S2$description))),]

unique(S4$description)
S4=S4[c(grep("S4",S4$description),grep("psD",S4$description)),]
S4=S4[-(c(grep("partial",S4$description),grep("UNVERIFIED",S4$description),grep("domain-containing",S4$description),grep("zinc-independent",S4$description))),]

unique(S5$description)
S5=S5[c(grep("S5",S5$description),grep("psE",S5$description)),]
S5=S5[-(c(grep("partial",S5$description),grep("UNVERIFIED",S5$description),grep("domain-containing",S5$description),grep("putative",S5$description))),]

unique(S7$description)
S7=S7[c(grep("S7",S7$description),grep("psG",S7$description)),]
S7=S7[-(c(grep("partial",S7$description),grep("UNVERIFIED",S7$description),grep("domain-containing",S7$description),grep("putative",S7$description))),]

unique(S8$description)
S8=S8[c(grep("S8",S8$description),grep("psH",S8$description)),]
S8=S8[-(c(grep("partial",S8$description),grep("UNVERIFIED",S8$description),grep("regulator",S8$description),grep("putative",S8$description))),]

unique(S9$description)
S9=S9[c(grep("S9",S9$description),grep("psI",S9$description)),]
S9=S9[-(c(grep("partial",S9$description),grep("UNVERIFIED",S9$description),grep("like",S9$description),grep("probable",S9$description))),]

unique(S10$description)
S10=S10[c(grep("S10",S10$description),grep("psJ",S10$description)),]
S10=S10[-(c(grep("partial",S10$description),grep("UNVERIFIED",S10$description),grep("transcription",S10$description),grep("fragment",S10$description))),]

unique(S11$description)
S11=S11[c(grep("S11",S11$description),grep("psK",S11$description)),]
S11=S11[-(c(grep("partial",S11$description),grep("UNVERIFIED",S11$description),grep("putative",S11$description),grep("HmaS11",S11$description))),]

unique(S12$description)
S12=S12[c(grep("S12",S12$description),grep("psL",S12$description)),]
S12=S12[-(c(grep("partial",S12$description),grep("UNVERIFIED",S12$description))),]

unique(S13$description)
S13=S13[c(grep("S13",S13$description),grep("psM",S13$description)),]
S13=S13[-(c(grep("partial",S13$description),grep("UNVERIFIED",S13$description))),]

unique(S14$description)
S14=S14[c(grep("S14",S14$description),grep("psN",S14$description)),]
S14=S14[-(c(grep("partial",S14$description),grep("UNVERIFIED",S14$description),grep("lternat",S14$description),grep("putative",S14$description),grep("Zinc",S14$description))),]

unique(S15$description)
S15=S15[c(grep("S15",S15$description),grep("psO",S15$description)),]
S15=S15[-(c(grep("partial",S15$description),grep("UNVERIFIED",S15$description),grep("putative",S15$description))),]

unique(S17$description)
S17=S17[c(grep("S17",S17$description),grep("psQ",S17$description)),]
S17=S17[-(c(grep("partial",S17$description),grep("UNVERIFIED",S17$description),grep("putative",S17$description))),]

unique(S19$description)
S19=S19[c(grep("S19",S19$description),grep("psS",S19$description)),]
S19=S19[-(c(grep("partial",S19$description),grep("UNVERIFIED",S19$description),grep("putative",S19$description))),]

# remove duplicates based on protein accession
for (protein in proteins) {
  df <- get(protein)
  df <- df[-(which(duplicated(df$accession))), ]
  assign(protein, df)
}

# remove duplicates based on genome accession
for (protein in proteins) {
  df <- get(protein)
  df <- df[-(which(duplicated(df$query_accession))), ]
  assign(protein, df)
}

# find the taxonomy of remained genomes from GTDB
for (protein in proteins) {
  df <- get(protein)
  #df_superkingdom=rep(NA,nrow(df))
  #df_phylum=rep(NA,nrow(df))
  #df_class=rep(NA,nrow(df))
  #df_order=rep(NA,nrow(df))
  #df_family=rep(NA,nrow(df))
  #df_genus=rep(NA,nrow(df))
  #df_species=rep(NA, nrow(df))
  #df=cbind(df, df_superkingdom, df_phylum, df_class, df_order, df_family, df_genus, df_species)
  for(i in 1:nrow(df)){
    print(i)
    n=grep(df$query_accession[i], gtdb_filtered$gtdb_genome_representative)
    if(length(n)>0){
      df$df_superkingdom[i]=gtdb_filtered$gtdb_superkingdom[n]
      df$df_phylum[i]=gtdb_filtered$gtdb_phylum[n]
      df$df_class[i]=gtdb_filtered$gtdb_class[n]
      df$df_order[i]=gtdb_filtered$gtdb_order[n]
      df$df_family[i]=gtdb_filtered$gtdb_family[n]
      df$df_genus[i]=gtdb_filtered$gtdb_genus[n]
      df$df_species[i]=gtdb_filtered$gtdb_species[n]
    }else{
      df$df_superkingdom[i]=NA
      df$df_phylum[i]=NA
      df$df_class[i]=NA
      df$df_order[i]=NA
      df$df_family[i]=NA
      df$df_genus[i]=NA
      df$df_species[i]=NA
    }
  }
  assign(protein, df)
}

which(is.na(S19$df_superkingdom))
S19=S19[-(which(is.na(S19$df_superkingdom))),]

sum(S19$df_superkingdom=="Bacteria")
unique(S19[S19$df_superkingdom=="Bacteria",]$df_phylum)
length(unique(S19[S19$df_superkingdom=="Bacteria",]$df_phylum))

sum(S19$df_superkingdom=="Archaea")
unique(S19[S19$df_superkingdom=="Archaea",]$df_phylum)
length(unique(S19[S19$df_superkingdom=="Archaea",]$df_phylum))

# take intersects based on remained genomes
a=intersect(IF2$query_accession, L2$query_accession)
a=intersect(a,S9 $query_accession) # repeat this for each

common=c(L1,L2,L5,L11,L14,S9,S10,S11,S19,IF2)

keep_row=c()
for (i in a){
  keep_row=c(keep_row, grep(i,IF2$query_accession))
}

IF2_intersect=IF2[keep_row,]
IF2_intersect=IF2_intersect[-(grep("sp0", IF2_intersect$df_species)),]

IF2_intersect_archaea=IF2_intersect[IF2_intersect$df_superkingdom=="Archaea",]
IF2_intersect_bacteria=IF2_intersect[IF2_intersect$df_superkingdom=="Bacteria",]
IF2_intersect_bacteria=IF2_intersect_bacteria[-(grep("sp0", IF2_intersect_bacteria$df_species)),]
IF2_intersect_bacteria=IF2_intersect_bacteria[-(grep("sp9", IF2_intersect_bacteria$df_species)),]
 
unique(IF2_intersect_bacteria$df_phylum)
unique(IF2_intersect_archaea$df_phylum)
  
# remove strain names from genus and species
for(i in 1:nrow(IF2_intersect_bacteria)){
  IF2_intersect_bacteria$df_species[i]=gsub("_AA","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AB","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AC","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AD","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AE","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AF","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AG","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AH","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AI","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AJ","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AK","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AL","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AM","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AN","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AO","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AP","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AQ","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AR","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AS","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AT","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AU","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AV","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AW","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AY","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AZ","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_AX","",IF2_intersect_bacteria$df_species[i])
  
  IF2_intersect_bacteria$df_species[i]=gsub("_BA","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BB","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BC","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BD","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BE","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BF","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BG","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BH","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BI","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BJ","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BK","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BL","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BM","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BN","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BO","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BP","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BQ","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BR","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BS","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BT","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BU","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BV","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BW","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BY","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BZ","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_BX","",IF2_intersect_bacteria$df_species[i])
  
  IF2_intersect_bacteria$df_species[i]=gsub("_A","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_B","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_C","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_D","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_E","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_F","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_G","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_H","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_I","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_J","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_K","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_L","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_M","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_N","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_O","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_P","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_Q","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_R","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_S","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_T","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_U","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_V","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_W","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_Y","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_Z","",IF2_intersect_bacteria$df_species[i])
  IF2_intersect_bacteria$df_species[i]=gsub("_X","",IF2_intersect_bacteria$df_species[i])
}

for(i in 1:nrow(IF2_intersect_bacteria)){
  IF2_intersect_bacteria$df_genus[i]=gsub("_AA","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AB","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AC","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AD","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AE","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AF","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AG","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AH","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AI","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AJ","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AK","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AL","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AM","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AN","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AO","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AP","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AQ","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AR","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AS","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AT","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AU","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AV","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AW","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AY","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AZ","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_AX","",IF2_intersect_bacteria$df_genus[i])
  
  IF2_intersect_bacteria$df_genus[i]=gsub("_BA","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BB","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BC","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BD","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BE","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BF","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BG","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BH","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BI","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BJ","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BK","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BL","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BM","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BN","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BO","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BP","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BQ","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BR","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BS","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BT","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BU","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BV","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BW","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BY","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BZ","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_BX","",IF2_intersect_bacteria$df_genus[i])
  
  IF2_intersect_bacteria$df_genus[i]=gsub("_A","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_B","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_C","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_D","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_E","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_F","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_G","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_H","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_I","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_J","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_K","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_L","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_M","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_N","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_O","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_P","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_Q","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_R","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_S","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_T","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_U","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_V","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_W","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_Y","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_Z","",IF2_intersect_bacteria$df_genus[i])
  IF2_intersect_bacteria$df_genus[i]=gsub("_X","",IF2_intersect_bacteria$df_genus[i])
}

# remove species duplicates
sum(duplicated(IF2_intersect_bacteria$df_species))
IF2_intersect_bacteria=IF2_intersect_bacteria[-which(duplicated(IF2_intersect_bacteria$df_species)),]

# keep one species per genus
sum(duplicated(IF2_intersect_bacteria$df_genus))
IF2_intersect_bacteria=IF2_intersect_bacteria[-which(duplicated(IF2_intersect_bacteria$df_genus)),]

# re-combine filtered archaea and bacteria IF2 again
IF2_filtered_archaea_bacteria=rbind(IF2_intersect_archaea,IF2_intersect_bacteria)
OGT=read.csv("~/genomes/gtdb/gtdb_ssu_with_OGT.csv", header=T)

ogt_predicted=c()
for(i in IF2_filtered_archaea_bacteria$query_accession){
  if(length(grep(i,OGT$genome))>0){
    ogt_predicted=c(ogt_predicted,i)
  }
}

ogt_keep=c()
for(i in ogt_predicted){
  if(length(grep(i, IF2_filtered_archaea_bacteria$query_accession))>0){
    ogt_keep=c(ogt_keep, grep(i, IF2_filtered_archaea_bacteria$query_accession))
  }
}

IF2_ogt=IF2_filtered_archaea_bacteria[ogt_keep,]
ogt=c()
for(i in IF2_ogt$query_accession){
  ogt=c(ogt,OGT[grep(i, OGT$genome),5])
}
IF2_ogt=cbind(IF2_ogt, ogt)

psychrophile=IF2_ogt[IF2_ogt$ogt<20,]
thermophile=IF2_ogt[IF2_ogt$ogt>45,]
mesophile=IF2_ogt[((IF2_ogt$ogt>20)&(IF2_ogt$ogt<45)),]

plot(IF2_ogt$length, IF2_ogt$ogt)
plot(psychrophile$length, psychrophile$ogt)
plot(mesophile$length, mesophile$ogt)
plot(thermophile$length, thermophile$ogt)

cor.test(IF2_ogt$length, IF2_ogt$ogt)
summary(lm(IF2_ogt$ogt~IF2_ogt$length))
summary(lm(IF2_ogt$length~IF2_ogt$ogt))

summary(lm(psychrophile$length~psychrophile$ogt))
summary(lm(mesophile$length~mesophile$ogt))
summary(lm(thermophile$length~thermophile$ogt))
