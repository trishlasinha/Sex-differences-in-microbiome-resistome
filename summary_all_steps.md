Sex specific differences in microbiome
-------------
Author of code: Trishla Sinha

Year: 2018

Steps:
1) Microbiome characterization using MetaPhlAn and HUMAnN2

2) Antibiotic characterization using ShortBRED

3) Richness and compositional calculations

4) Phenotypes correlations

5) Univariate analyses linear models (MaAsLin)

6) Multivariate analyses linea models (MaAsLin)

7) Age effect in males and females

8) Antibiotic prevalences

9) Antibiotic correlations to relative abundances


1.Microbiome characterization using MetaPhlAn and HUMAnN2
------------------------------------------------------------

Bash script for SLURM environment  

```{shell}

#!/bin/bash
#SBATCH --job-name=$SAMPLE_ID_metagenomes
#SBATCH --error=$SAMPLE_ID.err
#SBATCH --output=$SAMPLE_ID.out
#SBATCH --mem=70gb
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task=6

#BAM to FASTQ
java -jar ${EBROOTPICARD}/picard.jar SamToFastq I=$SAMPLE_ID.bam F=./$SAMPLE_ID/filtering_data/$SAMPLE_ID.fastq1 F2=./$SAMPLE_ID/filtering_data/$SAMPLE_ID.fastq2

#QC and trimming
kneaddata --input ./$SAMPLE_ID/filtering_data/$SAMPLE_ID_1.fq -t 6 -p 7 --input ./$SAMPLE_ID/filtering_data/$SAMPLE_ID_2.fq -db /groups/umcg-gastrocol/tmp03/metagenomic_tools/kneaddata-0.5.4/Homo_sapiens_Bowtie2_v0.1/ --output ./$SAMPLE_ID/filtering_data/ --log ./$SAMPLE_ID/clean_reads/$SAMPLE_ID.log

#Merge files (prepare input)
cat ./$SAMPLE_ID/filtering_data/$SAMPLE_ID_1_kneaddata_paired_1.fastq > ./$SAMPLE_ID/filtering_data/$SAMPLE_ID_kneaddata_merged.fastq
cat ./$SAMPLE_ID/filtering_data/$SAMPLE_ID_1_kneaddata_paired_2.fastq >> ./$SAMPLE_ID/filtering_data/$SAMPLE_ID_kneaddata_merged.fastq
mv ./$SAMPLE_ID/filtering_data/*kneaddata_paired_1.fastq ./$SAMPLE_ID/clean_reads/
mv ./$SAMPLE_ID/filtering_data/*kneaddata_paired_2.fastq ./$SAMPLE_ID/clean_reads/
mv ./$SAMPLE_ID/filtering_data/*kneaddata_merged.fastq ./$SAMPLE_ID/clean_reads/
rm -r ./$SAMPLE_ID/filtering_data/

#MetaPhlAn2 2.7.2 (taxonomy characterization)

metaphlan2.py ./$SAMPLE_ID/clean_reads/$SAMPLE_ID_kneaddata_merged.fastq  --input_type multifastq --mpa_pkl /groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_2/db_v20/mpa_v20_m200.pkl --nproc 6 -o ./$SAMPLE_ID/metaphlan/$SAMPLE_ID_metaphlan.txt --tmp_dir ./$SAMPLE_ID/clean_reads/

#HUMAnN2 v.0.10.0 (pathways characterization, we used Uniref90 & Chochophlan db)

humann2 --input ./$SAMPLE_ID/clean_reads/$SAMPLE_ID_kneaddata_merged.fastq --output ./$SAMPLE_ID/humann2/ --taxonomic-profile ./$SAMPLE_ID/metaphlan/$SAMPLE_ID_metaphlan.txt --threads 6 --o-log ./$SAMPLE_ID/clean_reads/$SAMPLE_ID.full.humann2.log --remove-temp-output
```

2.Antibiotic characterization using ShortBRED (0.9.5)
---------------------------------------------------

Bash script for SLURM environment

```{shell}

#!/bin/bash
#SBATCH --job-name=$SAMPLE_ID_metagenomes
#SBATCH --error=$SAMPLE_ID.err
#SBATCH --output=$SAMPLE_ID.out
#SBATCH --mem=25gb
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task=6

module load Biopython/1.65-foss-2015b-Python-2.7.11
python ./ShortBred/shortbred_quantify.py --markers ./ShortBred/ShortBRED_CARD_2017_markers.faa --wgs "$sample_id".tar.gz --results "$sample_id"_out --tmp "$sample_id"_tmp --usearch ./usearch6.1.544 --threads 6
```

3.Richness and compositional calculations
---------------------------------------------

 Calculations using R

```{R}
# Loading packages

library(vegan)
library(foreach)
library(doMC)
library(reshape2)

registerDoMC(2)

#Importing Taxa file
taxa<-read.delim("Taxa_filterd_present_more_than_5%", row.names = 1, header = T, sep = "\t")

#Importing Pathways file

path<-read.delim("Pathways_filterd_present_more_than_5%", row.names = 1, header = T, sep = "\t")

#Importing Phenotypes file

phenos <- read.delim("all_phenotpes_imputed_TS", header=T, row.names=1)

#Changing names to variables

# Taxa
meta.data=taxa

# Phenotypes

phenotypes=phenos

# Matching rownames

meta.data<-meta.data[rownames(phenotypes),]

# Calculating distance matrix using Bray-Curtis distances

dist.matrix <- vegdist(meta.data,method = "bray")

```
**3a. Variance explain Bray-Curtis vs phenotypes**

```{R}
dist.species.bray = vegdist(species, method = "bray")
dist.species.jaccard = vegdist(species, method = "jaccard")

dist.pathways.bray = vegdist(pathways, method = "bray")
dist.pathways.jaccard = vegdist(pathways, method = "jaccard")

#adonis

adonis.species.bray = adonis(dist.species.bray ~ RD + Age + Sex, data = phenotypes,permutations = 1000)
adonis.species.jaccard = adonis(dist.species.jaccard ~ RD + Age + Sex, data = phenotypes,permutations = 1000)

adonis.pathways.bray = adonis(dist.pathways.bray ~ RD + Age + Sex, data = phenotypes,permutations = 1000)
adonis.pathways.jaccard = adonis(dist.pathways.jaccard ~ RD + Age + Sex, data = phenotypes,permutations = 1000)

#betadisp
betadisp.species.bray = betadisper(dist.species.bray, phenotypes$Sex, type = "centroid")
betadisp.species.jaccard = betadisper(dist.species.jaccard, phenotypes$Sex, type = "centroid")
betadisp.pathways.bray = betadisper(dist.pathways.bray, phenotypes$Sex, type = "centroid")
betadisp.pathways.jaccard = betadisper(dist.pathways.jaccard, phenotypes$Sex, type = "centroid")

#betadisp.test
permutest(betadisp.species.bray,permutations = 1000)
permutest(betadisp.species.jaccard,permutations = 1000)
permutest(betadisp.pathways.bray,permutations = 1000)
permutest(betadisp.pathways.jaccard,permutations = 1000)
```

**3b. Shannon Index**

```{R}
alpha<-diversity(taxa, index = "shannon")
library(stats)
phenos$Shannon<-alpha
names(phenos2)
IQR(phenos$Shannon)
quantile(phenos$Shannon)                                  # remove outliers greater then Q3+(IQRx3) and smaller than (Q1-(IQRx3)) # In my case only removed those with a value less than (Q1-(IQRx3)
phenos2 <- phenos[ which(phenos$Shannon > 1.51), ]
WRSLLD <- wilcox.test(phenos2$Shannon ~ phenos2$antrop_gender.F1M2, conf.int=TRUE)
WRSLLD

```

4.Phenotypes correlations
---------------------------------------------

```{R}
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
 **4a. Removing one pair of the factors displaying strong correlation with each other (>0.8)(factor removed chosen at random)**

```{R}
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
 **4b. Finally deriving the file that shows the influence of all the significant factors to be used for correction in Maaslin on the Bray-Curtis distance**
```{R}

write.table(phenos2, "Phenotypes_84_factors_removing_11_highly_correlated_0.8_including_sex ", sep="\t", quote = F, row.names = T)
adosissup1<-intersect(row.names(onlysigpheno), names(phenos2))
adosissup2<-onlysigpheno[adosissup1,]
write.table(adosissup2, "84_factors_influencing_Bray_Curtis_Distance_1135_after_removing_11_highly_correlated_0.8_including_sex ", sep="\t", quote = F, row.names = T)

```

5.Univariate analyses linear models (MaAsLin)
-----------------------------------------------
**5a. MaAsLin for species**

```{R}
library(Maaslin)
#Creating a my input file
my_input<-merge(phenos, taxa, by="row.names")
colnames(my_input)[1] <- "SID"
```

Open text editor and make output file "input.read.config"

MaAsLin requires the tsv/csv file, the name of the output file in which MaAsLin puts all the results, and the R-script which says which columns/rows he should analyze (i.e. file.read.config).

```{R}
Maaslin("Phenotypes_for Maaslin.tsv", "Univariable_Males_versus_female", strInputConfig = "INPUT_input.read.config")
```

**5b. MaAsLin for pathways**

```{R}
Path <- read.delim("Pathways_filterd_present_more_than_5%",sep="\t", header=T, row.names = 1)
# Creating a my input file
my_input<-merge(phenos, path, by="row.names")
colnames(my_input)[1] <- "SID"
names(my_input)
# Important to do this as a tsv file
write.table(my_input, "Pathway_input.tsv", sep = "\t", quote = F, row.names = F)
```

Open text editor and make output file "input.read.config

```{R}
Maaslin("Pathway_input.tsv", "Pathway_Males_versus_females_force_RD_Age", strInputConfig = "INPUT_input.read.config")
```


6.Multivariate analyses linear models (MaAsLin)
-----------------------------------------------
**6a.MaAsLin for species**


```{R}
Maaslin("Phenotypes_for Maaslin.tsv", "Multivariable_Males_versus_female_", strInputConfig = "INPUT_input.read.config",,fAllvAll = F)
```

**6b. MaAsLin for pathways**

```{R}
Maaslin("Pathway_input.tsv", "Pathway_Males_versus_females_force_RD_Age", strInputConfig = "INPUT_input.read.config")
```

7.Age effect in males and females
---------------------------------


```{R}
mydata<-read.table("Phenotypes_for_Maaslin_comparing_male_versus_female_08_08_2018", sep = "\t", header = T, row.names = 1)
femalespecific<-read.table("Female_factors_09_08_2018.txt", sep = "\t", header = T, row.names = 1)

### Randomly choose 450 females and males from the sample
females<-mydata[mydata$Sex=="Female",]
males<-mydata[mydata$Sex=="Male",]
female450<-females[sample(nrow(females), 450), ]
male450<-males[sample(nrow(males), 450), ]
histogram(female450$antrop_age)
histogram(male450$antrop_age)
summary(female450$antrop_age)
summary(male450$antrop_age)

# Females 450 samples species
my_input<-merge(female450,taxa,by="row.names")
my_input$Sex=NULL
write.table(my_input, "Female_450_input.tsv", sep = "\t", quote = F, row.names = F)

# Males 450 samples species
my_input<-merge(male450,taxa,by="row.names")
my_input$Sex=NULL
write.table(my_input, "Male_450_input.tsv", sep = "\t", quote = F, row.names = F)
names(my_input)

# Effect of age in Males
Maaslin("Male_450_input.tsv", "Effect_Age_in_Males", strInputConfig = "INPUT_input.read.config")
# Effect of age in Females
Maaslin("Female_450_input.tsv", "Effect_Age_in_Females", strInputConfig = "INPUT_input.read.config")

#### Maaslin Pathways
# Females 450 samples pathways
my_input<-merge(female450,path,by="row.names")
my_input$Sex=NULL
write.table(my_input, "Female_450_input_pathways.tsv", sep = "\t", quote = F, row.names = F)

# Males 450 samples pathways
my_input<-merge(male450,path,by="row.names")
my_input$Sex=NULL
write.table(my_input, "Male_450_input_pathways.tsv", sep = "\t", quote = F, row.names = F)
names(my_input)

# Effect of age in Males
Maaslin("Male_450_input_pathways.tsv", "Effect_Age_in_Males_PATHWAYS", strInputConfig = "INPUT_input.read.config")
# Effect of age in Females
Maaslin("Female_450_input_pathways.tsv", "EffectAge_in_Females_PATHWAYS", strInputConfig = "INPUT_input.read.config")

# Post and pre menopausal age groups in men and women
prepost<-merge(femalespecific,females,by="row.names")
row.names(prepost)<-prepost$Row.names
prepost$Row.names=NULL
prepost$agegroup <- as.numeric(cut(prepost$antrop_age, breaks=c(17,35,58,65,Inf),labels=FALSE))
prepost$agegroup[prepost$agegroup == 1] <- "prememopausal"
prepost$agegroup[prepost$agegroup == 2] <- NA
prepost$agegroup[prepost$agegroup == 3] <- "postmenopausal"
prepost$agegroup[prepost$agegroup == 4] <- NA
prepost[is.na(prepost)] <- 0
femalesmen <- prepost[ which(prepost$Oral_contraceptive==0 & prepost$Hormonal_contraception_last_month == 0 & prepost$agegroup !=0), ]
femalesmen$agegroup<-as.factor(femalesmen$agegroup)
summary(femalesmen$agegroup)
males<-mydata[mydata$Sex=="Male",]
maleselder<-males
maleselder$agegroup <- as.numeric(cut(maleselder$antrop_age, breaks=c(17,30,55,65,Inf),labels=FALSE))
maleselder$agegroup[maleselder$agegroup == 1] <- "18-30"
maleselder$agegroup[maleselder$agegroup == 2] <- NA
maleselder$agegroup[maleselder$agegroup == 3] <- "55-65"
maleselder$agegroup[maleselder$agegroup == 4] <- NA
maleselder[is.na(maleselder)] <- 0
maleselder <- maleselder[ which(maleselder$agegroup!=0), ]
maleselder$agegroup<-as.factor(maleselder$agegroup)
summary(maleselder$agegroup)

# Females prepost menopausal
my_input<-merge(femalesmen,taxa,by="row.names")
my_input$Sex=NULL
write.table(my_input, "Females_pre_post_OR_removed_input.tsv", sep = "\t", quote = F, row.names = F)
names(my_input)

# Males prepost/above 50
my_input<-merge(maleselder,taxa,by="row.names")
my_input$Sex=NULL
write.table(my_input, "Male_age_groups_input.tsv", sep = "\t", quote = F, row.names = F)
names(my_input)

# Effect of pre-post menopausal in females (OC removed)
Maaslin("Females_pre_post_OR_removed_input.tsv", "Female_pre_post_OR_removed_Species", strInputConfig = "INPUT_input.read.config")
# Effect of age in males
Maaslin("Male_age_groups_input.tsv", "Male_age_groups_Species", strInputConfig = "INPUT_input.read.config")


```


8.Antibiotic prevalences
---------------------------------

**8a. Male versus females resistome**

Note: custom made logistic regression function is provided as a separate script

```{R}

library (psych)
library(plyr)

phenotypes <- read.table("./phenotypes_resistome.txt", sep="\t", header=T, row.names = 1)
aro_index <- read.csv("./annotation_all_markers.txt", sep='\t', row.names = 1, header = T)

#Remove antibiotic users

phenos=phenotypes[phenotypes$Antibiotics_merged>0,]

#Subset, Sex, Age, Seq.depth phenotypes

phenotypes=phenos[,c(2,1,3)]
pheno_male=phenos[phenos$Sex=="Male",]
pheno_female=phenos[phenos$Sex=="Female",]
pheno_female$Sex=NULL
pheno_male$Sex=NULL
pheno_female$Antibiotics_merged=NULL
pheno_male$Antibiotics_merged=NULL
pheno_male$Menstruations_present=NULL
pheno_female2=pheno_female[,c(3,1,2,4)]
pheno_male2=pheno_male[,c(3,1,2)]

#Analysis at AR gene level

markers_RAW <- read.csv("../Female_health/ARG_all_recode_LLD.tsv", sep='\t', header=T)

#Consolidate based on ARO number (CARD id / genes)

GF_RAW=ddply(markers_RAW,"ARO.Accession",numcolwise(sum)) #Collapse to 877 gene families

row.names(GF_RAW) <- GF_RAW$ARO.Accession

GF_RAW$ARO.Accession <- NULL

#Select samples
GF_samples=GF_RAW[colnames(GF_RAW)%in%rownames(phenotypes)]

GF_LLD <- GF_samples[rowSums(GF_samples>0) > round(ncol (GF_samples) * 0.05), ] #152 Genes

LLD_GF <- as.data.frame(t(GF_LLD))

logistic_regression(phenotypes, LLD_GF, 1)
```

```{R}

#Repeat for different AR gene-families and class levels
#Consolidate based on Gene families
RAW_annotation=merge(aro_index,markers_RAW,by = "ARO.Accession")
GF2_RAW=ddply(RAW_annotation,"AMR.Gene.Family",numcolwise(sum))
#217 Gene families present
GF2_RAW=GF2_RAW[-218,]
rownames(GF2_RAW)=GF2_RAW$AMR.Gene.Family
GF2_RAW$AMR.Gene.Family=NULL
GF2_RAW=GF2_RAW[colnames(GF2_RAW)%in%rownames(phenotypes)]
GF2_RAW2<- GF2_RAW[rowSums(GF2_RAW>0) > round(ncol (GF2_RAW) * 0.05), ] #69 Gene families
GF2_RAW2=as.data.frame(t(GF2_RAW2))
logistic_regression(phenotypes,GF2_RAW2, 1)

## Analysis AB class
#RAW_annotation=merge(classification,markers_RAW,by = "ARO.Accession")
GF2_RAW=ddply(RAW_annotation,"Drug.Class",numcolwise(sum))
#120 AB classes
GF2_RAW=GF2_RAW[-121,]
rownames(GF2_RAW)=GF2_RAW$Drug.Class
GF2_RAW$Drug.Class=NULL
GF2_RAW=GF2_RAW[colnames(GF2_RAW)%in%rownames(phenotypes)]
GF2_RAW2<- GF2_RAW[rowSums(GF2_RAW>0) > round(ncol (GF2_RAW) * 0.05), ] #53 AB classes
GF2_RAW2=as.data.frame(t(GF2_RAW2))
logistic_regression(phenotypes,GF2_RAW2, 1)
```

**8b. IBS: male resistome case controls**

```{R}

# Male gene IBS vs AB
GF_RAW=ddply(markers_RAW,"ARO.Accession",numcolwise(sum))
row.names(GF_RAW) <- GF_RAW$ARO.Accession
GF_RAW$ARO.Accession <- NULL
GF_samples=GF_RAW[colnames(GF_RAW)%in%rownames(pheno_male2)]
GF_LLD <- GF_samples[rowSums(GF_samples>0) > round(ncol (GF_samples) * 0.05), ]
LLD_GF <- as.data.frame(t(GF_LLD))
logistic_regression(pheno_male2, LLD_GF, 1)

# Male gene-family IBS vs AB

GF2_RAW=ddply(RAW_annotation,"AMR.Gene.Family",numcolwise(sum))
GF2_RAW=GF2_RAW[-218,]
rownames(GF2_RAW)=GF2_RAW$AMR.Gene.Family
GF2_RAW$AMR.Gene.Family=NULL
GF2_RAW=GF2_RAW[colnames(GF2_RAW)%in%rownames(pheno_male2)]
GF2_RAW2<- GF2_RAW[rowSums(GF2_RAW>0) > round(ncol (GF2_RAW) * 0.05), ]
GF2_RAW2=as.data.frame(t(GF2_RAW2))
logistic_regression(pheno_male2,GF2_RAW2, 1)

# Male class IBS vs AB
GF2_RAW=ddply(RAW_annotation,"Drug.Class",numcolwise(sum))
GF2_RAW=GF2_RAW[-121,]
rownames(GF2_RAW)=GF2_RAW$Drug.Class
GF2_RAW$Drug.Class=NULL
GF2_RAW=GF2_RAW[colnames(GF2_RAW)%in%rownames(pheno_male2)]
GF2_RAW2<- GF2_RAW[rowSums(GF2_RAW>0) > round(ncol (GF2_RAW) * 0.05), ]
GF2_RAW2=as.data.frame(t(GF2_RAW2))
logistic_regression(pheno_male2,GF2_RAW2, 1)
```

**8c. IBS: female resistome case controls**

```{R}

#Genes
GF_RAW=ddply(markers_RAW,"ARO.Accession",numcolwise(sum))
row.names(GF_RAW) <- GF_RAW$ARO.Accession
GF_RAW$ARO.Accession <- NULL
GF_samples=GF_RAW[colnames(GF_RAW)%in%rownames(pheno_female3)]
GF_LLD <- GF_samples[rowSums(GF_samples>0) > round(ncol (GF_samples) * 0.05), ]
LLD_GF <- as.data.frame(t(GF_LLD))
logistic_regression(pheno_female3, LLD_GF, 1)

#Gene families

GF2_RAW=ddply(RAW_annotation,"AMR.Gene.Family",numcolwise(sum))
GF2_RAW=GF2_RAW[-218,]
rownames(GF2_RAW)=GF2_RAW$AMR.Gene.Family
GF2_RAW$AMR.Gene.Family=NULL
GF2_RAW=GF2_RAW[colnames(GF2_RAW)%in%rownames(pheno_female3)]
GF2_RAW2<- GF2_RAW[rowSums(GF2_RAW>0) > round(ncol (GF2_RAW) * 0.05), ]
GF2_RAW2=as.data.frame(t(GF2_RAW2))
logistic_regression(pheno_female3,GF2_RAW2, 1)

# AB class

GF2_RAW=ddply(RAW_annotation,"Drug.Class",numcolwise(sum))
GF2_RAW=GF2_RAW[-121,]
rownames(GF2_RAW)=GF2_RAW$Drug.Class
GF2_RAW$Drug.Class=NULL
GF2_RAW=GF2_RAW[colnames(GF2_RAW)%in%rownames(pheno_female3)]
GF2_RAW2<- GF2_RAW[rowSums(GF2_RAW>0) > round(ncol (GF2_RAW) * 0.05), ]
GF2_RAW2=as.data.frame(t(GF2_RAW2))
logistic_regression(pheno_female3,GF2_RAW2, 1)
```

9.Antibiotic correlations to relative abundances
---------------------------------

**9a. Genes vs Species**

```{R}
genes_LLD_rab <- ddply(markers_LLD,"ARO.Accession",numcolwise(mean))
row.names(genes_LLD_rab) <- genes_LLD_rab$ARO.Accession
genes_LLD_rab$ARO.Accession <- NULL
genes_LLD_rab <- genes_LLD_rab[rowSums(genes_LLD_rab) > 0, ]
LLD_genes_rab <- as.data.frame(t(genes_LLD_rab))
LLD_genes_rab_phenos <- merge(phenotypes_ab_excl, LLD_genes_rab, by="row.names")
row.names(LLD_genes_rab_phenos) <- LLD_genes_rab_phenos$Row.names
LLD_genes_rab_phenos$Row.names <- NULL
LLD_genes_rab <- LLD_genes_rab_phenos[,-c(1:6)]
genes_LLD_rab <- as.data.frame(t(LLD_genes_rab))
#choose only those genes that passed 5% filtering
list_LLD_genes <- genes_LLD[,c(1,2)]
genes_LLD_rab <- merge(list_LLD_genes, genes_LLD_rab, by="row.names")
row.names(genes_LLD_rab) <- genes_LLD_rab$Row.names
genes_LLD_rab$Row.names <- NULL
genes_LLD_rab <- genes_LLD_rab[,-c(1:2)]
names(genes_LLD_rab)[names(genes_LLD_rab)=="LLDeep_0001.y"] <- "LLDeep_0001"
names(genes_LLD_rab)[names(genes_LLD_rab)=="LLDeep_0002.y"] <- "LLDeep_0002"
LLD_genes_rab <- as.data.frame(t(genes_LLD_rab))
rnp <- rcorr(as.matrix(taxa), as.matrix(LLD_genes_rab), type="spearman")
temp.matrix <- rnp$r
correlated_species <- temp.matrix[223:363,1:222]
correlated_species <- as.data.frame(correlated_species)
temp.matrix <- rnp$P
correlated_species_pval <- temp.matrix[223:363,1:222]
correlated_species_pval_bonf <- as.data.frame(matrix(p.adjust(correlated_species_pval, method='bonferroni'), ncol = 222))
row.names(correlated_species_pval_bonf) <- row.names(correlated_species_pval)
colnames(correlated_species_pval_bonf) <- colnames(correlated_species_pval)
colnames(correlated_species) = sub(".*s__","",colnames(correlated_species))
colnames(correlated_species_pval_bonf) = sub(".*s__","",colnames(correlated_species_pval_bonf))
genes_species_corr <- merge(list_genes_significant, correlated_species, by="row.names")
row.names(genes_species_corr) <- genes_species_corr$Name
genes_species_corr$Name <- NULL
genes_species_corr <- genes_species_corr[,-c(1:4)]
species_genes_corr <- as.data.frame(t(genes_species_corr))
genes_species_pval <- merge(list_genes_significant, correlated_species_pval_bonf, by="row.names")
row.names(genes_species_pval) <- genes_species_pval$Name
genes_species_pval$Name <- NULL
genes_species_pval <- genes_species_pval[,-c(1:4)]
species_genes_pval <- as.data.frame(t(genes_species_pval))
correlation_species_genes <- merge(species_genes_corr, species_genes_pval, by="row.names")
row.names(correlation_species_genes) <- correlation_species_genes$Row.names
correlation_species_genes$Row.names <- NULL
colnames(correlation_species_genes) = sub(".x","",colnames(correlation_species_genes))
colnames(correlation_species_genes) = sub("y","p-val",colnames(correlation_species_genes))
correlation_species_genes <- correlation_species_genes[ , order(names(correlation_species_genes))]
```
**9b. Gene families vs Species**

```{R}
gene_to_families <- genes_LLD_rab
gene_to_families[,1123] <- row.names(gene_to_families)
gene_to_families <- gene_to_families[,c(1123, 1:1122)]
names(gene_to_families)[names(gene_to_families)=="V1123"] <- "ARO.Accession"
gene_to_families <- merge(classification, gene_to_families, by="ARO.Accession")
gene_families_LLD_rab <- ddply(gene_to_families,"AMR.Gene.Family",numcolwise(mean))
row.names(gene_families_LLD_rab) <- gene_families_LLD_rab$AMR.Gene.Family
gene_families_LLD_rab$AMR.Gene.Family <- NULL
LLD_gene_families_rab <- as.data.frame(t(gene_families_LLD_rab))
rnp <- rcorr(as.matrix(taxa), as.matrix(LLD_gene_families_rab), type="spearman")
temp.matrix <- rnp$r
correlated_species_gf <- temp.matrix[223:284,1:222]
correlated_species_gf <- as.data.frame(correlated_species_gf)
temp.matrix <- rnp$P
correlated_species_gf_pval <- temp.matrix[223:284,1:222]
correlated_species_gf_pval_bonf <- as.data.frame(matrix(p.adjust(correlated_species_gf_pval, method='bonferroni'), ncol = 222))
row.names(correlated_species_gf_pval_bonf) <- row.names(correlated_species_gf_pval)
colnames(correlated_species_gf_pval_bonf) <- colnames(correlated_species_gf_pval)
colnames(correlated_species_gf) = sub(".*s__","",colnames(correlated_species_gf))
colnames(correlated_species_gf_pval_bonf) = sub(".*s__","",colnames(correlated_species_gf_pval_bonf))
gf_species_corr <- merge(list_gf_significant, correlated_species_gf, by="row.names")
row.names(gf_species_corr) <- gf_species_corr$Row.names
gf_species_corr$Row.names <- NULL
gf_species_corr <- gf_species_corr[,-c(1:2)]
species_gf_corr <- as.data.frame(t(gf_species_corr))
gf_species_pval <- merge(list_gf_significant, correlated_species_gf_pval_bonf, by="row.names")
row.names(gf_species_pval) <- gf_species_pval$Row.names
gf_species_pval$Row.names <- NULL
gf_species_pval <- gf_species_pval[,-c(1:2)]
species_gf_pval <- as.data.frame(t(gf_species_pval))
correlation_species_gf <- merge(species_gf_corr, species_gf_pval, by="row.names")
row.names(correlation_species_gf) <- correlation_species_gf$Row.names
correlation_species_gf$Row.names <- NULL
colnames(correlation_species_gf) = sub(".x","",colnames(correlation_species_gf))
colnames(correlation_species_gf) = sub("y","p-val",colnames(correlation_species_gf))
correlation_species_gf <- correlation_species_gf[ , order(names(correlation_species_gf))]
```

**9c. Genes vs Genus**

```{R}
##genes & genus
rnp <- rcorr(as.matrix(taxa_genus), as.matrix(LLD_genes_rab), type="spearman")
temp.matrix <- rnp$r
correlated_genus <- temp.matrix[89:229,1:88]
correlated_genus <- as.data.frame(correlated_genus)
temp.matrix <- rnp$P
correlated_genus_pval <- temp.matrix[89:229,1:88]
correlated_genus_pval_bonf <- as.data.frame(matrix(p.adjust(correlated_genus_pval, method='bonferroni'), ncol = 88))
row.names(correlated_genus_pval_bonf) <- row.names(correlated_genus_pval)
colnames(correlated_genus_pval_bonf) <- colnames(correlated_genus_pval)
colnames(correlated_genus) = sub(".*g__","",colnames(correlated_genus))
colnames(correlated_genus_pval_bonf) = sub(".*g__","",colnames(correlated_genus_pval_bonf))
genes_genus_corr <- merge(list_genes_significant, correlated_genus, by="row.names")
row.names(genes_genus_corr) <- genes_genus_corr$Name
genes_genus_corr$Name <- NULL
genes_genus_corr <- genes_genus_corr[,-c(1:4)]
genus_genes_corr <- as.data.frame(t(genes_genus_corr))
genes_genus_pval <- merge(list_genes_significant, correlated_genus_pval_bonf, by="row.names")
row.names(genes_genus_pval) <- genes_genus_pval$Name
genes_genus_pval$Name <- NULL
genes_genus_pval <- genes_genus_pval[,-c(1:4)]
genus_genes_pval <- as.data.frame(t(genes_genus_pval))
correlation_genus_genes <- merge(genus_genes_corr, genus_genes_pval, by="row.names")
row.names(correlation_genus_genes) <- correlation_genus_genes$Row.names
correlation_genus_genes$Row.names <- NULL
colnames(correlation_genus_genes) = sub(".x","",colnames(correlation_genus_genes))
colnames(correlation_genus_genes) = sub("y","p-val",colnames(correlation_genus_genes))
correlation_genus_genes <- correlation_genus_genes[ , order(names(correlation_genus_genes))]
```

**9d. Gene families vs Genus**

```{R}
##gene families & genus
rnp <- rcorr(as.matrix(taxa_genus), as.matrix(LLD_gene_families_rab), type="spearman")
temp.matrix <- rnp$r
correlated_genus_gf <- temp.matrix[89:150,1:88]
correlated_genus_gf <- as.data.frame(correlated_genus_gf)
temp.matrix <- rnp$P
correlated_genus_pval_gf <- temp.matrix[89:150,1:88]
correlated_genus_pval_gf_bonf <- as.data.frame(matrix(p.adjust(correlated_genus_pval_gf, method='bonferroni'), ncol = 88))
row.names(correlated_genus_pval_gf_bonf) <- row.names(correlated_genus_pval_gf)
colnames(correlated_genus_pval_gf_bonf) <- colnames(correlated_genus_pval_gf)
colnames(correlated_genus_gf) = sub(".*g__","",colnames(correlated_genus_gf))
colnames(correlated_genus_pval_gf_bonf) = sub(".*g__","",colnames(correlated_genus_pval_gf_bonf))
gf_genus_corr <- merge(list_gf_significant, correlated_genus_gf, by="row.names")
row.names(gf_genus_corr) <- gf_genus_corr$Row.names
gf_genus_corr$Row.names <- NULL
gf_genus_corr <- gf_genus_corr[,-c(1:2)]
genus_gf_corr <- as.data.frame(t(gf_genus_corr))
gf_genus_pval <- merge(list_gf_significant, correlated_genus_pval_gf_bonf, by="row.names")
row.names(gf_genus_pval) <- gf_genus_pval$Row.names
gf_genus_pval$Row.names <- NULL
gf_genus_pval <- gf_genus_pval[,-c(1:2)]
genus_gf_pval <- as.data.frame(t(gf_genus_pval))
correlation_genus_gf <- merge(genus_gf_corr, genus_gf_pval, by="row.names")
row.names(correlation_genus_gf) <- correlation_genus_gf$Row.names
correlation_genus_gf$Row.names <- NULL
colnames(correlation_genus_gf) = sub(".x","",colnames(correlation_genus_gf))
colnames(correlation_genus_gf) = sub("y","p-val",colnames(correlation_genus_gf))
correlation_genus_gf <- correlation_genus_gf[ , order(names(correlation_genus_gf))]
```
