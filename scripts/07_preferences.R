# 07_preferences.R
# this script compares the distribution of extension lengths with respect to organisms' environmental preferences
library(ggplot2)

load("~/data/rdata/domain-genome-gs-gc-ogt-ox-tree.Rdata")
doubleCheck=finalData[,c(1,2,3,5,7,16,17,18,19)]

# temperature
bd_psychrophile=read.csv("~/data/traits/BacDive_psychrophile.csv")
bd_mesophile=read.csv("~/data/traits/BacDive_mesophile.csv")
bd_thermophile=read.csv("~/data/traits/BacDive_thermophile.csv")

# oxygen tolerance
bd_aerobe=read.csv("~/data/traits/BacDive_aerobe.csv")
bd_aerotolerant=read.csv("~/data/traits/BacDive_aerotolerant.csv")
bd_anaerobe=read.csv("~/data/traits/BacDive_anaerobe.csv")
bd_facultative_aerobe=read.csv("~/data/traits/BacDive_facultative_aerobe.csv")
bd_facultative_anaerobe=read.csv("~/data/traits/BacDive_facultative_anaerobe.csv")
bd_microaerophile=read.csv("~/data/traits/BacDive_microaerophile.csv")
bd_microaerotolerant=read.csv("~/data/traits/BacDive_microaerotolerant.csv")
bd_obligate_aerobe=read.csv("~/data/traits/BacDive_obligate_aerobe.csv")
bd_obligate_anaerobe=read.csv("~/data/traits/BacDive_obligate_anaerobe.csv")

# habitat
bd_habitat=read.csv("~/traits/BacDive_habitat.csv")

temperatureType=c()
for(i in 1:nrow(doubleCheck)){
  id=strsplit(doubleCheck$genome_name[i],"\\.")[[1]][1]
  temp_type=c()
  if (id %in% bd_psychrophile$Genome.seq..accession.number){temp_type=c(temp_type,"psychrophile")}
  if (id %in% bd_mesophile$Genome.seq..accession.number){temp_type=c(temp_type,"mesophile")}
  if (id %in% bd_thermophile$Genome.seq..accession.number){temp_type=c(temp_type,"thermophile")}
  else if((!id %in% bd_psychrophile$Genome.seq..accession.number)&(!id %in% bd_mesophile$Genome.seq..accession.number)&(!id %in% bd_thermophile$Genome.seq..accession.number)) {temp_type=c(temp_type,"NA")}
  temperatureType=c(temperatureType, paste(temp_type,collapse="_"))
}

doubleCheck=cbind(doubleCheck,temperatureType)


# plot known temperatures and temperature preference from BacDive
temperatureType=as.matrix(temperatureType)
# for having two categories
doubleCheck[c(72,705,769,823),10]="thermophile"
doubleCheck[doubleCheck$temperatureType=="psychrophile_mesophile_thermophile",10]="mesophile"

doubleCheck[c(727,741,750,767,785,796,806,810,812,814,819,831,834),10]="mesophile"
doubleCheck[doubleCheck$temperatureType=="mesophile_thermophile",10]="thermophile"

doubleCheck[c(569,606),10]="psychrophile"
doubleCheck[doubleCheck$temperatureType=="psychrophile_mesophile",10]="mesophile"

doubleCheck_noNA=doubleCheck[-(which(doubleCheck$temperatureType=="NA")),]

doubleCheck_noNA$temperatureType=factor(doubleCheck_noNA$temperatureType,levels=c("psychrophile","mesophile","thermophile"))

doubleCheck_noNA=doubleCheck_noNA[-(which((doubleCheck_noNA$L_IF2N==0)&(doubleCheck_noNA$L_IF2C==0))),]


ggplot(doubleCheck_noNA, aes(x=temperatureType,y=as.numeric(L_IF2N), fill=temperatureType))+
  geom_boxplot(width=0.5, alpha=0.6) +
  scale_fill_manual(values=c("#6387B5","#4DAA8C","#D63A31"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none", axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x="", y="IF2 N-terminal Length")

# plot predicted temperatures and temperature preference
predicted_doubleCheck=finalData[,c(1,2,3,5,7,16,17,18,19)]
temperatureTypePredicted=rep(NA,nrow(predicted_doubleCheck))

for(i in 1:nrow(predicted_doubleCheck)){
  if(!is.na(predicted_doubleCheck$predicted_ogt[i])){
    if(predicted_doubleCheck$predicted_ogt[i]<=20){
      temperatureTypePredicted[i]="psychrophile"
    }else if(predicted_doubleCheck$predicted_ogt[i]>45){
      temperatureTypePredicted[i]="thermophile" 
    }else if((predicted_doubleCheck$predicted_ogt[i]>20)&(predicted_doubleCheck$predicted_ogt[i]<=45)){
      temperatureTypePredicted[i]="mesophile" 
    }
  }else{temperatureTypePredicted[i]="NA"}
}

predicted_doubleCheck=cbind(predicted_doubleCheck,temperatureTypePredicted)
predicted_doubleCheck_noNA=predicted_doubleCheck[-(which(predicted_doubleCheck$temperatureTypePredicted=="NA")),]
predicted_doubleCheck_noNA$temperatureTypePredicted=factor(predicted_doubleCheck_noNA$temperatureTypePredicted,levels=c("psychrophile","mesophile","thermophile"))
predicted_doubleCheck_noNA=predicted_doubleCheck_noNA[-(which((predicted_doubleCheck_noNA$L_IF2N==0)&(predicted_doubleCheck_noNA$L_IF2C==0))),]

ggplot(predicted_doubleCheck_noNA, aes(x=temperatureTypePredicted,y=as.numeric(L_IF2N), fill=temperatureTypePredicted))+
  geom_boxplot(width=0.5, alpha=0.6) +
  scale_fill_manual(values=c("#6387B5","#4DAA8C","#D63A31"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none", axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x="", y="IF2 N-terminal Length")

# inferred and known temperatures together
predicted_doubleCheck_noNA=cbind(predicted_doubleCheck_noNA, rep("inferred", nrow(predicted_doubleCheck_noNA)))
doubleCheck_noNA=cbind(doubleCheck_noNA, rep("BacDive", nrow(doubleCheck_noNA)))
colnames(predicted_doubleCheck_noNA)[10]="temperatureType"
colnames(predicted_doubleCheck_noNA)[11]="source"
colnames(doubleCheck_noNA)[11]="source"

source_combined=rbind(predicted_doubleCheck_noNA,doubleCheck_noNA)

ggplot(source_combined, aes(x=temperatureType, y=as.numeric(L_IF2N), fill=source))+
  geom_boxplot(width=0.5, alpha=0.6) +
  scale_fill_manual(values=c("#A87BB6","#E09C87"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",
        axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x="", y="IF2 N-terminal Length")


# oxygen
oxygen_tolerance=c()
for(i in 1:nrow(doubleCheck)){
  id=strsplit(doubleCheck$genome_name[i],"\\.")[[1]][1]
  oxy_type=c()
  if (id %in% bd_aerobe$Genome.seq..accession.number){oxy_type=c(oxy_type,"aerobe")}
  if (id %in% bd_aerotolerant$Genome.seq..accession.number){oxy_type=c(oxy_type,"aerotolerant")}
  if (id %in% bd_anaerobe$Genome.seq..accession.number){oxy_type=c(oxy_type,"anaerobe")}
  if (id %in% bd_facultative_aerobe$Genome.seq..accession.number){oxy_type=c(oxy_type,"facultativeAerobe")}
  if (id %in% bd_facultative_anaerobe$Genome.seq..accession.number){oxy_type=c(oxy_type,"facultativeAnaerobe")}
  if (id %in% bd_microaerophile$Genome.seq..accession.number){oxy_type=c(oxy_type,"microaerophile")}
  if (id %in% bd_microaerotolerant$Genome.seq..accession.number){oxy_type=c(oxy_type,"microaerotolerant")}
  if (id %in% bd_obligate_aerobe$Genome.seq..accession.number){oxy_type=c(oxy_type,"obligateAerobe")}
  if (id %in% bd_obligate_anaerobe$Genome.seq..accession.number){oxy_type=c(oxy_type,"obligateAnaerobe")}
  
  else if((!id %in% bd_aerobe$Genome.seq..accession.number)&
          (!id %in% bd_aerotolerant$Genome.seq..accession.number)&
          (!id %in% bd_anaerobe$Genome.seq..accession.number)&
          (!id %in% bd_facultative_aerobe$Genome.seq..accession.number)&
          (!id %in% bd_facultative_anaerobe$Genome.seq..accession.number)&
          (!id %in% bd_microaerophile$Genome.seq..accession.number)&
          (!id %in% bd_microaerotolerant$Genome.seq..accession.number)&
          (!id %in% bd_obligate_aerobe$Genome.seq..accession.number)&
          (!id %in% bd_obligate_anaerobe$Genome.seq..accession.number)) {oxy_type=c(oxy_type,"NA")}
  oxygen_tolerance=c(oxygen_tolerance, paste(oxy_type,collapse="_"))
}

doubleCheck=cbind(doubleCheck,oxygen_tolerance)
doubleCheck_noNA=doubleCheck[-(which(doubleCheck$oxygen_tolerance=="NA")),]

doubleCheck_noNA$oxygen_tolerance[doubleCheck_noNA$oxygen_tolerance=="aerobe_obligateAerobe"]="aerobe"
doubleCheck_noNA$oxygen_tolerance[doubleCheck_noNA$oxygen_tolerance=="obligateAerobe"]="aerobe"
doubleCheck_noNA$oxygen_tolerance[doubleCheck_noNA$oxygen_tolerance=="obligateAnaerobe"]="anaerobe"
doubleCheck_noNA$oxygen_tolerance[doubleCheck_noNA$oxygen_tolerance=="anaerobe_obligateAnaerobe"]="anaerobe"

doubleCheck_noNA$oxygen_tolerance[doubleCheck_noNA$oxygen_tolerance=="facultativeAerobe"]="facultative"
doubleCheck_noNA$oxygen_tolerance[doubleCheck_noNA$oxygen_tolerance=="aerobe_anaerobe"]="facultative"
doubleCheck_noNA$oxygen_tolerance[doubleCheck_noNA$oxygen_tolerance=="aerobe_facultativeAnaerobe"]="facultative"
doubleCheck_noNA$oxygen_tolerance[doubleCheck_noNA$oxygen_tolerance=="facultativeAerobe_microaerophile"]="facultative"
doubleCheck_noNA$oxygen_tolerance[doubleCheck_noNA$oxygen_tolerance=="facultativeAnaerobe"]="facultative"
doubleCheck_noNA$oxygen_tolerance[doubleCheck_noNA$oxygen_tolerance=="facultativeAnaerobe_obligateAerobe"]="facultative"
doubleCheck_noNA$oxygen_tolerance[doubleCheck_noNA$oxygen_tolerance=="microaerophile"]="facultative"
doubleCheck_noNA$oxygen_tolerance[doubleCheck_noNA$oxygen_tolerance=="microaerophile_obligateAerobe"]="facultative"

doubleCheck_noNA=doubleCheck_noNA[-(which((doubleCheck_noNA$L_IF2N==0)&(doubleCheck_noNA$L_IF2C==0))),]

doubleCheck_noNA_predicted=doubleCheck
doubleCheck_noNA_predicted$oxygen_tolerance=doubleCheck_noNA_predicted$tolerance
doubleCheck_noNA_predicted=doubleCheck_noNA_predicted[-(which((doubleCheck_noNA_predicted$L_IF2N==0)&(doubleCheck_noNA_predicted$L_IF2C==0))),]
doubleCheck_noNA_predicted=cbind(doubleCheck_noNA_predicted, rep("predicted",nrow(doubleCheck_noNA_predicted)))
doubleCheck_noNA=cbind(doubleCheck_noNA, rep("BacDive",nrow(doubleCheck_noNA)))

colnames(doubleCheck_noNA)[11]="source"
colnames(doubleCheck_noNA_predicted)[11]="source"

source_combined=rbind(doubleCheck_noNA, doubleCheck_noNA_predicted)
source_combined$oxygen_tolerance[source_combined$oxygen_tolerance=="anaerobic"]="anaerobe"
source_combined$oxygen_tolerance[source_combined$oxygen_tolerance=="aerobic"]="aerobe"

source_combined$oxygen_tolerance=factor(source_combined$oxygen_tolerance, levels=c("aerobe","facultative","anaerobe"))

ggplot(source_combined, aes(x=oxygen_tolerance, y=as.numeric(L_IF2N), fill=source))+
  geom_boxplot(width=0.5, alpha=0.6) +
  scale_fill_manual(values=c("#A87BB6","#E09C87"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",
        axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x="", y="IF2 N-terminal Length")


ggplot(source_combined, aes(x=oxygen_tolerance, y=as.numeric(L_IF2C), fill=source))+
  geom_boxplot(width=0.5, alpha=0.6) +
  scale_fill_manual(values=c("#A87BB6","#E09C87"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",
        axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"),
        plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
  labs(title=" ", x="", y="IF2 C-terminal Length")


source_combined=source_combined[-(which((source_combined$L_IF2N==0)&(source_combined$L_IF2C==0))),]

#statistical tests
bd_aerobe_if2=doubleCheck_noNA[doubleCheck_noNA$oxygen_tolerance=="aerobe",]
bd_anaerobe_if2=doubleCheck_noNA[doubleCheck_noNA$oxygen_tolerance=="anaerobe",]
bd_facultative_if2=doubleCheck_noNA[doubleCheck_noNA$oxygen_tolerance=="facultative",]

pd_aerobe_if2=doubleCheck_noNA_predicted[doubleCheck_noNA_predicted$oxygen_tolerance=="aerobic",]
pd_anaerobe_if2=doubleCheck_noNA_predicted[doubleCheck_noNA_predicted$oxygen_tolerance=="anaerobic",]
pd_facultative_if2=doubleCheck_noNA_predicted[doubleCheck_noNA_predicted$oxygen_tolerance=="facultative",]

t.test(bd_aerobe_if2$L_IF2C,pd_aerobe_if2$L_IF2C, p.adjust="bonferroni")
t.test(bd_anaerobe_if2$L_IF2C,pd_anaerobe_if2$L_IF2C, p.adjust="bonferroni")
t.test(bd_facultative_if2$L_IF2C,pd_facultative_if2$L_IF2C, p.adjust="bonferroni")

wilcox.test(bd_aerobe_if2$L_IF2C,pd_aerobe_if2$L_IF2C, p.adjust="bonferroni")
wilcox.test(bd_anaerobe_if2$L_IF2C,pd_anaerobe_if2$L_IF2C, p.adjust="bonferroni")
wilcox.test(bd_facultative_if2$L_IF2C,pd_facultative_if2$L_IF2C, p.adjust="bonferroni")

t.test(pd_aerobe_if2$L_IF2C,pd_anaerobe_if2$L_IF2C, p.adjust="bonferroni")
t.test(pd_aerobe_if2$L_IF2C,pd_facultative_if2$L_IF2C, p.adjust="bonferroni")
t.test(pd_anaerobe_if2$L_IF2C,pd_facultative_if2$L_IF2C, p.adjust="bonferroni")

wilcox.test(pd_aerobe_if2$L_IF2C,pd_anaerobe_if2$L_IF2C, p.adjust="bonferroni")
wilcox.test(pd_aerobe_if2$L_IF2C,pd_facultative_if2$L_IF2C, p.adjust="bonferroni")
wilcox.test(pd_anaerobe_if2$L_IF2C,pd_facultative_if2$L_IF2C, p.adjust="bonferroni")
