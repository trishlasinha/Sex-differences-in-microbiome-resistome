############################ Science Male Female ###############################

rm(list = ls())

rd<-read.delim("RD_Age_Gender_1135", row.names = 1, header = T, sep = "\t")


# Taxa file 
taxa<-read.delim("Taxa_filterd_present_more_than_5%", row.names = 1, header = T, sep = "\t")

# Pathways 

path<-read.delim("Pathways_filterd_present_more_than_5%", row.names = 1, header = T, sep = "\t")

# Phenotypes 
phenos <- read.delim("all_phenotpes_imputed_TS", header=T, row.names=1)


# Loading all the packages 
library(vegan)
library(foreach)
library(doMC)
library(reshape2)
registerDoMC(2)


# Taxa 
#meta.data=taxa


# Phenotypes 
#phenotypes=phenos

# Matching rownames
#meta.data<-meta.data[rownames(phenotypes),]

# Making a distance matrix 
#dist.matrix <- vegdist(meta.data,method = "bray")



#adon<-foreach(i=1:ncol(phenotypes),.combine=rbind)%do%{
  
  #ad1<-adonis(dist.matrix ~ phenotypes[,i],permutations=1000,parallel=2)
  #ad1$aov.tab[1,]
#}

#adon
#results<-adon
#rownames(adon) = colnames(phenotypes)
#adon$qFDR=NULL
#adon2<-adon[-c(1,123,136,146),]
#row.names(adon)


adon<-read.delim("Adonis_1135_204_factors_Science.txt", row.names = 1, header = T, sep = "\t")
adon[,7] = p.adjust(adon[,6],method = "BH")
colnames(adon)[7] = "qFDR"


onlysigpheno<-adon[adon$qFDR<0.05,]

phenossig<-intersect(row.names(onlysigpheno), names(phenos))
phenos3<-phenos[,phenossig]
phenos2<-phenos3


spearman <- function(x,y) {
  matchID <- intersect(rownames(x), rownames(y))
  x1 <- x[matchID,]
  y1 <- y[matchID,]
  result_cor <- matrix(nrow=ncol(x1), ncol=ncol(y1))
  rownames(result_cor) <- colnames(x1)
  colnames(result_cor) <- colnames(y1)
  result_pvalue <- matrix(nrow=ncol(x1), ncol=ncol(y1))
  for (i in 1:ncol(y1)) {
    for (j in 1:ncol(x1)) {
      cor1<-cor.test(x1[,j], y1[,i], method = "spearman")
      result_cor[j,i]= cor1$estimate 
      result_pvalue[j,i]= cor1$p.value
    }}
  result = list()
  result$p.val= result_pvalue
  result$cor= result_cor
  return(result)
}


result <- function(x,y){
  correlation <- spearman(x,y)
  a<- melt(correlation$cor) 
  a<- cbind(a, melt(correlation$p.val)[,"value"])
  result= a[order(a[,4]),]
  colnames(result)=c("factor1", "factor2", "CorCoefficient","pvalue")
  return(result)
}

P<-result(phenos2,phenos2)
warnings()
P$FDR<-p.adjust(P$pvalue, method = "BH")
sig<-P[P$FDR<0.05,]
sig1<-sig[complete.cases(sig),]

trial<-sig1[sig1$CorCoefficient>0.8,]

### remove all the factors with a strong correlation with each other (>0.8)

phenos2$how_often_tea=NULL
phenos2$protein.total_log=NULL
phenos2$protein.plant_log=NULL
phenos2$how_often_alcohol=NULL
phenos2$BlCells_Granulo_log=NULL
phenos2$how_often_pasta=NULL
phenos2$how_often_rice=NULL
phenos2$fat.total_log=NULL
phenos2$how_often_muesli=NULL
phenos2$how_often_soda=NULL

write.table(phenos2, "90_factors_influencing_Bray_Curtis_Distance_1135_after_removing_10_highly_correlated_0.8 ", sep="\t", quote = F, row.names = T)
phenos<-read.delim("90_factors_influencing_Bray_Curtis_Distance_1135_after_removing_10_highly_correlated_0.8 ", row.names = 1, header = T, sep = "\t")
phenos$RD<-rd$RD
phenos<- phenos[,c(1,91,2:90)]
colnames(phenos) [1]<- "Sex"
phenos$Sex<-rd$antrop_gender.F1M2

########## Maaslin for species ##################

#Creating a my input file 
my_input<-merge(phenos, taxa, by="row.names")
colnames(my_input)[1] <- "SID"


# Important to do this as a tsv file
write.table(my_input, "Phenotypes_for Maaslin.tsv", sep = "\t", quote = F, row.names = F)
names(my_input) # to get the name of the first species 

# Open text editor and make output file "input.read.config"
library(Maaslin)
#Univariate
Maaslin("Phenotypes_for Maaslin.tsv", "Males_versus_female_Force_Read_Depth_Age", strInputConfig = "INPUT_input.read.config", strForcedPredictors = c("RD","antrop_age"),fAllvAll = T, dSignificanceLevel = 0.05)
#Multivariate
Maaslin("Phenotypes_for Maaslin.tsv", "Males_versus_female_Force_Read_Depth_Age", strInputConfig = "INPUT_input.read.config", strForcedPredictors = c("RD","antrop_age"),fAllvAll = T, dSignificanceLevel = 0.05)


########## Maaslin for pathways ##################

# Making the pathway into percentages

Path <- read.delim("Pathways_filterd_present_more_than_5%",sep="\t", header=T, row.names = 1)

# Creating a my input file 
my_input<-merge(phenos, path, by="row.names")
colnames(my_input)[1] <- "SID"
names(my_input)

# Important to do this as a tsv file
write.table(my_input, "Pathway_input.tsv", sep = "\t", quote = F, row.names = F)

# Open text editor and make output file "input.read.config
#Univariate
Maaslin("Pathway_input.tsv", "Pathway_Males_versus_females_force_RD_Age", strInputConfig = "INPUT_input.read.config", strForcedPredictors = c("RD", "antrop_age"), fAllvAll = T,dSignificanceLevel = 0.05)
#Multivariate
Maaslin("Pathway_input.tsv", "Pathway_Males_versus_females_force_RD_Age", strInputConfig = "INPUT_input.read.config", strForcedPredictors = c("RD", "antrop_age"), fAllvAll = T,dSignificanceLevel = 0.05)

#### Extracting only the females 

females<-all[all$antrop_gender.F1M2=="Female",]
fh<-read.table("Female_health_Nov_2017.txt", sep = "\t", header = T, row.names = 1)

femaleall<-merge(all,fh,by="row.names")
femaleall$antrop_gender.F1M2=NULL
row.names(femaleall)<-femaleall$Row.names
femaleall$Row.names=NULL

#Creating a my input file 
my_input_female<-merge(femaleall, Path, by="row.names")
colnames(my_input_female)[1] <- "SID"

# Important to do this as a tsv file
write.table(my_input_female, "Pathway_female.tsv", sep = "\t", quote = F, row.names = F)
names(my_input_female) # to get the name of the first species 

# Open text editor and make output file "input.read.config"
library(Maaslin)
Maaslin("Pathway_female.tsv", "Pathway_Female_factors_force_RD_AGe", strInputConfig = "Pathway_female_input.read.config", strForcedPredictors = c("RD","antrop_age"),fAllvAll = T,dSignificanceLevel = 0.05)

summary(all)


ames(my_input)[1:95]

#### Extracting only the females 

females<-all[all$antrop_gender.F1M2=="Female",]
fh<-read.table("Female_health_Nov_2017.txt", sep = "\t", header = T, row.names = 1)

femaleall<-merge(all,fh,by="row.names")
femaleall$antrop_gender.F1M2=NULL
row.names(femaleall)<-femaleall$Row.names
femaleall$Row.names=NULL

#Creating a my input file 
my_input_female<-merge(femaleall, taxa, by="row.names")
colnames(my_input_female)[1] <- "SID"


# Important to do this as a tsv file
write.table(my_input_female, "my_input_female.tsv", sep = "\t", quote = F, row.names = F)
names(my_input_female) # to get the name of the first species 

# Open text editor and make output file "input.read.config"
library(Maaslin)
Maaslin("my_input_female.tsv", "Female_factors_force_RD_AGe", strInputConfig = "FEMALE_input.read.config", strForcedPredictors = c("RD","antrop_age"),fAllvAll = T,dSignificanceLevel = 0.05)

Maaslin("my_input_female.tsv", "Female_factors_exploratative", strInputConfig = "FEMALE_input.read.config",dSignificanceLevel = 0.25)

