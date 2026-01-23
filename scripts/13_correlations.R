# 13_correlations.R
# script for correlation plots

#load libraries
library(phylolm)
library(ggplot2)
library(cowplot)

# load the latest data with properties
load("~/data/rdata/domain-genome-gs-gc-ogt-ox-tree.Rdata")

# change the IF2 extension between N and C for following lines when necessary

##### TEMPERATURE
correlation3=as.data.frame(cbind(finalData$genome_name, finalData$predicted_ogt, finalData$L_IF2N))
colnames(correlation3)=c("genome_name", "predicted_ogt","L_IF2N")
correlation3$predicted_ogt=as.numeric(correlation3$predicted_ogt)
correlation3$L_IF2N=as.numeric(correlation3$L_IF2N)
rownames(correlation3)=correlation3[,1]
correlation3=correlation3[,-1]

Domain=c()
for(i in 1:nrow(correlation3)){Domain=c(Domain,strsplit(rownames(correlation3)[i],"_")[[1]][3])}
correlation3=as.data.frame(cbind(correlation3,Domain))

fit0_temp=lm(correlation3$predicted_ogt~correlation3$L_IF2N)
coeff0_temp<-coefficients(fit0_temp)           
intercept0_temp<-coeff0_temp[1] 
slope0_temp<- coeff0_temp[2]

fit1_temp= phylolm(predicted_ogt~L_IF2N,data=correlation3,phy=tree,model="lambda",boot=100)
coeff1_temp<-coefficients(fit1_temp)           
intercept1_temp<-coeff1_temp[1] 
slope1_temp <- coeff1_temp[2]

cor.test(correlation3$predicted_ogt,correlation3$L_IF2N)

##### OXYGEN
correlation4=as.data.frame(cbind(finalData$genome_name, finalData$oxygen_enzyme, finalData$L_IF2N))
colnames(correlation4)=c("genome_name","oxygen_enzyme","L_IF2N")
correlation4$oxygen_enzyme=as.numeric(correlation4$oxygen_enzyme)
correlation4$L_IF2N=as.numeric(correlation4$L_IF2N)
rownames(correlation4)=correlation4[,1]
correlation4=correlation4[,-1]

Domain=c()
for(i in 1:nrow(correlation4)){Domain=c(Domain,strsplit(rownames(correlation4)[i],"_")[[1]][3])}
correlation4=as.data.frame(cbind(correlation4,Domain))

fit0_oxy=lm(correlation4$oxygen_enzyme~correlation4$L_IF2N)
coeff0_oxy<-coefficients(fit0_oxy)           
intercept0_oxy<-coeff0_oxy[1] 
slope0_oxy<- coeff0_oxy[2]

fit1_oxy=phylolm(oxygen_enzyme~L_IF2N,data=correlation4,phy=tree,model="lambda",boot=100)
coeff1_oxy<-coefficients(fit1_oxy)           
intercept1_oxy<-coeff1_oxy[1] 
slope1_oxy <- coeff1_oxy[2]

cor.test(correlation4$oxygen_enzyme,correlation4$L_IF2)

#plot
colors=c("#749E89FF","#C385A2FF","#7D87B2FF") 


p3=ggplot(correlation3,aes(x=L_IF2N, y=predicted_ogt))+geom_point(color="black",shape=21, size=4, aes(fill=factor(Domain)))+scale_fill_manual(values=colors)+
  geom_abline(intercept = intercept0_temp, slope = slope0_temp, color="gray20",linetype="solid", size=1)+
  geom_abline(intercept = intercept1_temp, slope = slope1_temp, color="gray20",linetype="dashed", size=1)+
  theme_bw(base_size=12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none", axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"))+
  xlab("IF2 N-Terminal Extension Length") +
  ylab("Optimal Growth Temperature")
p3

p4=ggplot(correlation4,aes(x=L_IF2N, y=oxygen_enzyme))+geom_point(color="black",shape=21, size=4, aes(fill=factor(Domain)))+scale_fill_manual(values=colors)+
  geom_abline(intercept = intercept0_oxy, slope = slope0_oxy, color="gray20",linetype="solid", size=1)+
  geom_abline(intercept = intercept1_oxy, slope = slope1_oxy, color="gray20",linetype="dashed", size=1)+
  theme_bw(base_size=12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none", axis.text=element_text(size=12, face="bold", color="black"), axis.title = element_text(size=14, face="bold", color="black"))+
  xlab("IF2 N-Terminal Extension Length") +
  ylab("# of O2 Utilizing Enzymes")
p4

# for within domains
correlation3A=correlation3[correlation3$Domain=="Archaea",]
correlation3B=correlation3[correlation3$Domain=="Bacteria",]


correlation4A=correlation4[correlation4$Domain=="Archaea",]
correlation4B=correlation4[correlation4$Domain=="Bacteria",]
correlation4E=correlation4[correlation4$Domain=="Eukaryote",]

# Pearson correlation tests
# for temperature (eukaryote not included)
cor.test(correlation3$predicted_ogt,correlation3$L_IF2N)
cor.test(correlation3A$predicted_ogt,correlation3A$L_IF2N)
cor.test(correlation3B$predicted_ogt,correlation3B$L_IF2N)

# for oxygen
cor.test(correlation4$oxygen_enzyme,correlation4$L_IF2N)
cor.test(correlation4A$oxygen_enzyme,correlation4A$L_IF2N)
cor.test(correlation4B$oxygen_enzyme,correlation4B$L_IF2N)
cor.test(correlation4E$oxygen_enzyme,correlation4E$L_IF2N)

