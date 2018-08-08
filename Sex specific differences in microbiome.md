Sex specific differences in microbiome
-------------
Author of code: Trishla Sinha

Year: 2018

Aim: 
1) Using package vegan to detect differences in overall gut microbiome between males and females and factors influencing overall gut microbiome in men and women
2) Using MaAsLin to detect differences in microbial species and pathways in gut microbiome between males and females 


Vegan to detect sex differences
-------------

#Importing Taxa file

```
taxa<-read.delim("Taxa_filterd_present_more_than_5%", row.names = 1, header = T, sep = "\t")
```
#Importing Pathways file
```
path<-read.delim("Pathways_filterd_present_more_than_5%", row.names = 1, header = T, sep = "\t")
```
#Importing Phenotypes file
```
phenos <- read.delim("all_phenotpes_imputed_TS", header=T, row.names=1)
```

#Loading all the required packages
```
library(vegan)
library(foreach)
library(doMC)
library(reshape2)
registerDoMC(2)
```

#Naming files
```
# Taxa 
meta.data=taxa

# Phenotypes 
#phenotypes=phenos

# Matching rownames
#meta.data<-meta.data[rownames(phenotypes),]

# Making a distance matrix 
#dist.matrix <- vegdist(meta.data,method = "bray")
```
**1 a) Beta diversity: Adonis test**
```
adon<-foreach(i=1:ncol(phenotypes),.combine=rbind)%do%{
  
  ad1<-adonis(dist.matrix ~ phenotypes[,i],permutations=1000,parallel=2)
  ad1$aov.tab[1,]
}

adon
results<-adon
rownames(adon) = colnames(phenotypes)
adon[,7] = p.adjust(adon[,6],method = "BH")
colnames(adon)[7] = "qFDR"
onlysigpheno<-adon[adon$qFDR<0.05,]
phenossig<-intersect(row.names(onlysigpheno), names(phenos))
phenos3<-phenos[,phenossig]
phenos2<-phenos3
```
#Spearman's correlation between phenotypes
```
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
```
#Removing one pair of the factors displaying strong correlation with each other (>0.8)(factor removed chosen at random)
```
phenos2$protein.plant_log=NULL
phenos2$protein.total_log=NULL
phenos2$protein.animal_log=NULL
phenos2$BlCells_Granulo_log=NULL
phenos2$how_often_pasta=NULL
phenos2$how_often_rice=NULL
phenos2$carbohydrates.total_log=NULL
phenos2$fat.total_log=NULL
phenos2$how_often_muesli=NULL
phenos2$how_often_alcohol=NULL
phenos2$how_often_soda=NULL

```
**1 b) Alpha diversity**
```
alpha<-diversity(taxa, index = "shannon")
library(stats)
phenos$Shannon<-alpha
names(phenos2)
IQR(phenos$Shannon)
quantile(phenos$Shannon) # remove outliers greater then Q3+(IQRx3) and smaller than (Q1-(IQRx3)) # In my case only removed those with a value less than (Q1-(IQRx3)
phenos2 <- phenos[ which(phenos$Shannon > 1.51), ]
WRSLLD <- wilcox.test(phenos2$Shannon ~ phenos2$antrop_gender.F1M2, conf.int=TRUE)
WRSLLD
```

**2. a)  Maaslin for species**

```
library(Maaslin)
#Creating a my input file 
my_input<-merge(phenos, taxa, by="row.names")
colnames(my_input)[1] <- "SID"
```

#Important to do this as a tsv file
```
write.table(my_input, "Phenotypes_for Maaslin.tsv", sep = "\t", quote = F, row.names = F)
names(my_input) # to get the name of the first species 

```

#Open text editor and make output file "input.read.config"

MaAsLin requires the tsv/csv file, the name of the output file in which MaAsLin puts all the results, and the R-script which says which columns/rows he should analyze (i.e. file.read.config). Furthermore, I forced a zero inflated model (fZeroInlfated = T) and force correction of phenotypes.

#Univariate Maaslin analysis
```
Maaslin("Phenotypes_for Maaslin.tsv", "Males_versus_female_Force_Read_Depth_Age", strInputConfig = "INPUT_input.read.config", strForcedPredictors = c("RD","antrop_age"),fAllvAll = T, dSignificanceLevel = 0.05)
```
#Multivariate Maaslin analysis
```
Maaslin("Phenotypes_for Maaslin.tsv", "Males_versus_female_Force_Read_Depth_Age", strInputConfig = "INPUT_input.read.config", strForcedPredictors = c("RD","antrop_age"),fAllvAll = T, dSignificanceLevel = 0.05)
```

**2. b)Maaslin for pathways**

```
Path <- read.delim("Pathways_filterd_present_more_than_5%",sep="\t", header=T, row.names = 1)
# Creating a my input file 
my_input<-merge(phenos, path, by="row.names")
colnames(my_input)[1] <- "SID"
names(my_input)

# Important to do this as a tsv file
write.table(my_input, "Pathway_input.tsv", sep = "\t", quote = F, row.names = F)
```
#Open text editor and make output file "input.read.config

#Univariate Maaslin analysis
```
Maaslin("Pathway_input.tsv", "Pathway_Males_versus_females_force_RD_Age", strInputConfig = "INPUT_input.read.config", strForcedPredictors = c("RD", "antrop_age"), fAllvAll = T,dSignificanceLevel = 0.05)
```
#Multivariate Maaslin analysis
```
Maaslin("Pathway_input.tsv", "Pathway_Males_versus_females_force_RD_Age", strInputConfig = "INPUT_input.read.config", strForcedPredictors = c("RD", "antrop_age"), fAllvAll = T,dSignificanceLevel = 0.05)
```



