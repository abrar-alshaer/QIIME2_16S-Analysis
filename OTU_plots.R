setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/16S Microbiome Analysis HF+EPA Study/Week 13-14 10.23.18 Analysis/")
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(gplots)
library("dplyr")
library("ggpubr")
library(car)

qiime <- read.csv("n3_qiime_phylum_Log10_norm_metaMerged.csv", header = TRUE, row.names = 1) #load in normalized OTU file
qiime <- qiime[c(1:5,7:18,20:28),] #removing mouse #832 & #823 because they were singly housed
qiime <- qiime[order(qiime$Diet),] #sort by diet group

###Making OTU Boxplots

#Actinobacteria

#checking for assumption of normality for each diet within a microbe
shapiro.test(qiime[1:10,2]) #1:10 are the control diets, column 2 is actinobacteria
hist(qiime[1:10,2])
qqnorm(qiime[1:10,2])
qqline(qiime[1:10,2])
#Checking for the assumption of equal variance with normal data
bartlett.test(qiime[1:10,2], qiime$Diet.group)
#Checking for the assumption of equal variance with non-normal data
leveneTest(qiime[1:10,2] ~ qiime$Diet.group, data = qiime)
#if the data has unequal variance but is normal, you can use a Welch T-test in place of an ANOVA

shapiro.test(qiime[11:17,2]) #11:17 are the HF diets, column 2 is actinobacteria
hist(qiime[11:17,2])
qqnorm(qiime[11:17,2])
qqline(qiime[11:17,2])
bartlett.test(qiime[11:17,2], qiime$Diet.group)
leveneTest(qiime[11:17,2] ~ qiime$Diet.group, data = qiime)

shapiro.test(qiime[18:26,2]) #18:26 are the HF+EPA diets, column 2 is actinobacteria
hist(qiime[18:26,2])
qqnorm(qiime[18:26,2])
qqline(qiime[18:26,2])
bartlett.test(qiime[18:26,2], qiime$Diet.group)
leveneTest(qiime[18:26,2] ~ qiime$Diet.group, data = qiime)

compare_means(Actinobacteria ~ Diet.group,  data = qiime, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("Control", "High Fat"), c("Control", "High Fat + EPA"))
p1 <- ggboxplot(qiime, x = "Diet.group", y = "Actinobacteria",
          title = "Actinobacteria OTU Plot", xlab = "Diet", ylab = "Log10 Normalized Absolute Abundance",
          color = "Diet.group", add = c("point"), ylim = c(min(qiime$Actinobacteria), max(qiime$Actinobacteria)+1)) + 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Bacteroidetes

#checking for assumption of normality
shapiro.test(qiime[1:10,3])
hist(qiime[1:10,3])
qqnorm(qiime[1:10,3])
qqline(qiime[1:10,3])
#Checking for the assumption of equal variance with normal data
bartlett.test(qiime[1:10,3], qiime$Diet.group)
#Checking for the assumption of equal variance with non-normal data
leveneTest(qiime[1:10,3] ~ qiime$Diet.group, data = qiime)
#if the data has unequal variance but is normal, you can use a Welch T-test in place of an ANOVA

shapiro.test(qiime[11:17,3]) #11:17 are the HF diets
hist(qiime[11:17,3])
qqnorm(qiime[11:17,3])
qqline(qiime[11:17,3])
bartlett.test(qiime[11:17,3], qiime$Diet.group)
leveneTest(qiime[11:17,3] ~ qiime$Diet.group, data = qiime)

shapiro.test(qiime[18:26,3]) #18:26 are the HF+EPA diets
hist(qiime[18:26,3])
qqnorm(qiime[18:26,3])
qqline(qiime[18:26,3])
bartlett.test(qiime[18:26,3], qiime$Diet.group)
leveneTest(qiime[18:26,3] ~ qiime$Diet.group, data = qiime)

compare_means(Bacteroidetes ~ Diet.group,  data = qiime, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("Control", "High Fat"), c("Control", "High Fat + EPA"))
p2 <- ggboxplot(qiime, x = "Diet.group", y = "Bacteroidetes",
                title = "Bacteroidetes OTU Plot", xlab = "Diet", ylab = "Log10 Normalized Absolute Abundance",
                color = "Diet.group", add = c("point"), ylim = c(min(qiime$Bacteroidetes), max(qiime$Bacteroidetes)+0.4)) + 
  theme(legend.position='none') #palette = "jco"
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Firmicutes

#checking for assumption of normality
shapiro.test(qiime[1:10,7])
hist(qiime[1:10,7])
qqnorm(qiime[1:10,7])
qqline(qiime[1:10,7])
#Checking for the assumption of equal variance with normal data
bartlett.test(qiime[1:10,7], qiime$Diet.group)
#Checking for the assumption of equal variance with non-normal data
leveneTest(qiime[1:10,7] ~ qiime$Diet.group, data = qiime)
#if the data has unequal variance but is normal, you can use a Welch T-test in place of an ANOVA

shapiro.test(qiime[11:17,7]) #11:17 are the HF diets
hist(qiime[11:17,7])
qqnorm(qiime[11:17,7])
qqline(qiime[11:17,7])
bartlett.test(qiime[11:17,7], qiime$Diet.group)
leveneTest(qiime[11:17,7] ~ qiime$Diet.group, data = qiime)

shapiro.test(qiime[18:26,7]) #18:26 are the HF+EPA diets
hist(qiime[18:26,7])
qqnorm(qiime[18:26,7])
qqline(qiime[18:26,7])
bartlett.test(qiime[18:26,7], qiime$Diet.group)
leveneTest(qiime[18:26,7] ~ qiime$Diet.group, data = qiime)

compare_means(Firmicutes ~ Diet.group,  data = qiime, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
res.aov <- aov(Firmicutes ~ Diet.group,  data = qiime)
TukeyHSD(res.aov)
#my_comparisons <- list(c("Control", "High Fat"), c("Control", "High Fat + EPA"))
p3 <- ggboxplot(qiime, x = "Diet.group", y = "Firmicutes",
                title = "Firmicutes OTU Plot", xlab = "Diet", ylab = "Log10 Normalized Absolute Abundance",
                color = "Diet.group", add = c("point"), ylim = c(min(qiime$Firmicutes), max(qiime$Firmicutes)+0.4)) + 
  theme(legend.position='none') #palette = "jco"
                #label = "Sample_Diet", repel = TRUE
#  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Proteobacteria

#checking for assumption of normality
shapiro.test(qiime[1:10,9])
hist(qiime[1:10,9])
qqnorm(qiime[1:10,9])
qqline(qiime[1:10,9])
#Checking for the assumption of equal variance with normal data
bartlett.test(qiime[1:10,9], qiime$Diet.group)
#Checking for the assumption of equal variance with non-normal data
leveneTest(qiime[1:10,9] ~ qiime$Diet.group, data = qiime)
#if the data has unequal variance but is normal, you can use a Welch T-test in place of an ANOVA

shapiro.test(qiime[11:17,9]) #11:17 are the HF diets
hist(qiime[11:17,9])
qqnorm(qiime[11:17,9])
qqline(qiime[11:17,9])
bartlett.test(qiime[11:17,9], qiime$Diet.group)
leveneTest(qiime[11:17,9] ~ qiime$Diet.group, data = qiime)

shapiro.test(qiime[18:26,9]) #18:26 are the HF+EPA diets
hist(qiime[18:26,9])
qqnorm(qiime[18:26,9])
qqline(qiime[18:26,9])
bartlett.test(qiime[18:26,9], qiime$Diet.group)
leveneTest(qiime[18:26,9] ~ qiime$Diet.group, data = qiime)

compare_means(Proteobacteria ~ Diet.group,  data = qiime, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
res.aov <- aov(Proteobacteria ~ Diet.group,  data = qiime)
TukeyHSD(res.aov) #more conservative than T-test

my_comparisons <- list(c("Control", "High Fat"), c("Control", "High Fat + EPA"))
p4 <- ggboxplot(qiime, x = "Diet.group", y = "Proteobacteria",
                title = "Proteobacteria OTU Plot", xlab = "Diet", ylab = "Log10 Normalized Absolute Abundance",
                color = "Diet.group", add = c("point"), ylim = c(min(qiime$Proteobacteria), max(qiime$Proteobacteria)+0.4)) + 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Tenericutes

#checking for assumption of normality
shapiro.test(qiime[1:10,10])
hist(qiime[1:10,10])
qqnorm(qiime[1:10,10])
qqline(qiime[1:10,10])
#Checking for the assumption of equal variance with normal data
bartlett.test(qiime[1:10,10], qiime$Diet.group)
#Checking for the assumption of equal variance with non-normal data
leveneTest(qiime[1:10,10] ~ qiime$Diet.group, data = qiime)
#if the data has unequal variance but is normal, you can use a Welch T-test in place of an ANOVA

shapiro.test(qiime[11:17,10]) #11:17 are the HF diets
hist(qiime[11:17,10])
qqnorm(qiime[11:17,10])
qqline(qiime[11:17,10])
bartlett.test(qiime[11:17,10], qiime$Diet.group)
leveneTest(qiime[11:17,10] ~ qiime$Diet.group, data = qiime)

shapiro.test(qiime[18:26,10]) #18:26 are the HF+EPA diets
hist(qiime[18:26,10])
qqnorm(qiime[18:26,10])
qqline(qiime[18:26,10])
bartlett.test(qiime[18:26,10], qiime$Diet.group)
leveneTest(qiime[18:26,10] ~ qiime$Diet.group, data = qiime)

compare_means(Tenericutes ~ Diet.group,  data = qiime, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
res.aov <- aov(Tenericutes ~ Diet.group,  data = qiime)
TukeyHSD(res.aov) #more conservative than T-test

#my_comparisons <- list(c("Control", "High Fat"), c("Control", "High Fat + EPA"))
p5 <- ggboxplot(qiime, x = "Diet.group", y = "Tenericutes",
                title = "Tenericutes OTU Plot", xlab = "Diet", ylab = "Log10 Normalized Absolute Abundance",
                color = "Diet.group", add = c("point"), ylim = c(min(qiime$Tenericutes), max(qiime$Tenericutes)+0.4))+ 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Verrucomicrobia

#checking for assumption of normality
shapiro.test(qiime[1:10,11])
hist(qiime[1:10,11])
qqnorm(qiime[1:10,11])
qqline(qiime[1:10,11])
#Checking for the assumption of equal variance with normal data
bartlett.test(qiime[1:10,11], qiime$Diet.group)
#Checking for the assumption of equal variance with non-normal data
leveneTest(qiime[1:10,11] ~ qiime$Diet.group, data = qiime)
#if the data has unequal variance but is normal, you can use a Welch T-test in place of an ANOVA

shapiro.test(qiime[11:17,11]) #11:17 are the HF diets, column 11 is Verrucomicrobia
hist(qiime[11:17,11])
qqnorm(qiime[11:17,11])
qqline(qiime[11:17,11])
bartlett.test(qiime[11:17,11], qiime$Diet.group)
leveneTest(qiime[11:17,11] ~ qiime$Diet.group, data = qiime)

shapiro.test(qiime[18:26,11]) #18:26 are the HF+EPA diets, column 11 is Verrucomicrobia
hist(qiime[18:26,11])
qqnorm(qiime[18:26,11])
qqline(qiime[18:26,11])
bartlett.test(qiime[18:26,11], qiime$Diet.group)
leveneTest(qiime[18:26,11] ~ qiime$Diet.group, data = qiime)

compare_means(Verrucomicrobia ~ Diet.group,  data = qiime, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
res.aov <- aov(Verrucomicrobia ~ Diet.group,  data = qiime)
TukeyHSD(res.aov) #more conservative than T-test

my_comparisons <- list(c("Control", "High Fat"), c("Control", "High Fat + EPA"))
p6 <- ggboxplot(qiime, x = "Diet.group", y = "Verrucomicrobia",
                title = "Verrucomicrobia OTU Plot", xlab = "Diet", ylab = "Log10 Normalized Absolute Abundance",
                color = "Diet.group", add = c("point"), ylim = c(min(qiime$Verrucomicrobia), max(qiime$Verrucomicrobia)+0.4))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

ggarrange(p1,p2,p3,p4,p5,p6, ncol = 2, nrow = 3, legend = "none")
