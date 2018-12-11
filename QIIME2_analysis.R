rm(list=ls())
setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/16S Microbiome Analysis HF+EPA Study/DEMUX_2018/DEMUX_2018/ALL_RUNS/Week13_14/Paired End Analysis/")
library(phyloseq) 
library(ggplot2)
library(dplyr)
library(car)
library(ggpubr)

meta <- read.csv("metadata_W13-14_10.23.18.csv", header = TRUE, row.names = 9) #load in metadata
meta <- meta["Diet"]

###Making PHYLUM Boxplots#######

#qiime <- read.csv("table.from_biom_phylum_filter.csv", header = FALSE, row.names = 1)
#Normalized_DESeq2_counts_biomTableComboData
qiime <- read.csv("Normalized_DESeq2_SequenceCountsComboData.csv", header = FALSE, row.names = 1)
qiime <- qiime[c(1:7),]
qiime <- data.frame(t(qiime)) #transpose the dataframe
qiime$V1 <- substring(qiime$V1, 1) #remove the first character of the first column strings
#assign row names
qiime2 <- qiime[,-1]
qiime2 <- mutate_all(qiime2, function(x) as.numeric(as.character(x)))
rownames(qiime2) <- qiime[,1]
#merge dataframe by sample IDs to append diet information 
qiime2 <- merge(qiime2, meta, by="row.names", all=TRUE)
#sort by diet group
qiime2 <- qiime2[order(qiime2$Diet),] 
row.names(qiime2) <- qiime2[,1]
qiime2 <- qiime2[,-1] #remove the first column
qiime2_log <- log10(qiime2[,1:6]) #only Log10 first 6 columns (don't include the Diet column)
qiime2_log[qiime2_log < 0] <- 0 #replacing all negative infinaties with 0 (-Inf happens when you Log(0))
qiime2_log$Diet <- qiime2$Diet #add back the diet column

#Actinobacteria

#checking for assumption of normality for each diet within a microbe
shapiro.test(qiime2_log[1:11,1]) #1:11 are the control diets, column 1 is actinobacteria
hist(qiime2_log[1:11,1])
qqnorm(qiime2_log[1:11,1])
qqline(qiime2_log[1:11,1])

shapiro.test(qiime2_log[12:19,1]) #12:19 are the HF diets, column 1 is actinobacteria
hist(qiime2_log[12:19,1])
qqnorm(qiime2_log[12:19,1])
qqline(qiime2_log[12:19,1])

shapiro.test(qiime2_log[20:28,1]) #20:28 are the HF+EPA diets, column 1 is actinobacteria
hist(qiime2_log[20:28,1])
qqnorm(qiime2_log[20:28,1])
qqline(qiime2_log[20:28,1])

compare_means(Actinobacteria ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF"), c("C", "HF_EPA"))
p1 <- ggboxplot(qiime2_log, x = "Diet", y = "Actinobacteria",
                title = "Actinobacteria", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Actinobacteria), max(qiime2_log$Actinobacteria)+1)) + 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Bacteroidetes

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,2])
hist(qiime2_log[1:11,2])
qqnorm(qiime2_log[1:11,2])
qqline(qiime2_log[1:11,2])

shapiro.test(qiime2_log[12:19,2]) 
hist(qiime2_log[12:19,2])
qqnorm(qiime2_log[12:19,2])
qqline(qiime2_log[12:19,2])

shapiro.test(qiime2_log[20:28,2]) 
hist(qiime2_log[20:28,2])
qqnorm(qiime2_log[20:28,2])
qqline(qiime2_log[20:28,2])

compare_means(Bacteroidetes ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF"))
p2 <- ggboxplot(qiime2_log, x = "Diet", y = "Bacteroidetes",
                title = "Bacteroidetes", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Bacteroidetes), max(qiime2_log$Bacteroidetes)+0.4)) + 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Firmicutes

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,3])
hist(qiime2_log[1:11,3])
qqnorm(qiime2_log[1:11,3])
qqline(qiime2_log[1:11,3])

shapiro.test(qiime2_log[12:19,3]) 
hist(qiime2_log[12:19,3])
qqnorm(qiime2_log[12:19,3])
qqline(qiime2_log[12:19,3])

shapiro.test(qiime2_log[20:28,3]) 
hist(qiime2_log[20:28,3])
qqnorm(qiime2_log[20:28,3])
qqline(qiime2_log[20:28,3])

compare_means(Firmicutes ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
p3 <- ggboxplot(qiime2_log, x = "Diet", y = "Firmicutes",
                title = "Firmicutes", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Firmicutes), max(qiime2_log$Firmicutes)+0.4)) + 
  theme(legend.position='none') #palette = "jco"
#label = "Sample_Diet", repel = TRUE

#Proteobacteria

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,4])
hist(qiime2_log[1:11,4])
qqnorm(qiime2_log[1:11,4])
qqline(qiime2_log[1:11,4])

shapiro.test(qiime2_log[12:19,4]) 
hist(qiime2_log[12:19,4])
qqnorm(qiime2_log[12:19,4])
qqline(qiime2_log[12:19,4])

shapiro.test(qiime2_log[20:28,4]) 
hist(qiime2_log[20:28,4])
qqnorm(qiime2_log[20:28,4])
qqline(qiime2_log[20:28,4])

compare_means(Proteobacteria ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("HF", "HF_EPA"), c("C", "HF_EPA"))
#I am excluding the 0 value (row number 5)
p4 <- ggboxplot(qiime2_log[c(1:4,6:28),], x = "Diet", y = "Proteobacteria",
                title = "Proteobacteria", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log[c(1:4,6:28),4]), max(qiime2_log$Proteobacteria)+0.5)) + 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Tenericutes

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,5])
hist(qiime2_log[1:11,5])
qqnorm(qiime2_log[1:11,5])
qqline(qiime2_log[1:11,5])

shapiro.test(qiime2_log[12:19,5]) 
hist(qiime2_log[12:19,5])
qqnorm(qiime2_log[12:19,5])
qqline(qiime2_log[12:19,5])

shapiro.test(qiime2_log[20:28,5]) 
hist(qiime2_log[20:28,5])
qqnorm(qiime2_log[20:28,5])
qqline(qiime2_log[20:28,5])

compare_means(Tenericutes ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF_EPA"))
p5 <- ggboxplot(qiime2_log, x = "Diet", y = "Tenericutes",
                title = "Tenericutes", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Tenericutes), max(qiime2_log$Tenericutes)+0.4))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Verrucomicrobia

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,6])
hist(qiime2_log[1:11,6])
qqnorm(qiime2_log[1:11,6])
qqline(qiime2_log[1:11,6])

shapiro.test(qiime2_log[12:19,6]) 
hist(qiime2_log[12:19,6])
qqnorm(qiime2_log[12:19,6])
qqline(qiime2_log[12:19,6])

shapiro.test(qiime2_log[20:28,6]) 
hist(qiime2_log[20:28,6])
qqnorm(qiime2_log[20:28,6])
qqline(qiime2_log[20:28,6])

compare_means(Verrucomicrobia ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF_EPA"))
p6 <- ggboxplot(qiime2_log, x = "Diet", y = "Verrucomicrobia",
                title = "Verrucomicrobia", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Verrucomicrobia), max(qiime2_log$Verrucomicrobia)+0.4))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

ggarrange(p1,p2,p3,p4,p5,p6, ncol = 2, nrow = 3, legend = "none")

###Making CLASS Boxplots####

qiime <- read.csv("Normalized_DESeq2_SequenceCountsComboData.csv", header = FALSE, row.names = 1)
qiime <- qiime[c(1,8:13),]
qiime <- data.frame(t(qiime)) #transpose the dataframe
qiime$V1 <- substring(qiime$V1, 1) #remove the first character of the first column strings
#assign row names
qiime2 <- qiime[,-1]
qiime2 <- mutate_all(qiime2, function(x) as.numeric(as.character(x)))
rownames(qiime2) <- qiime[,1]
#merge dataframe by sample IDs to append diet information 
qiime2 <- merge(qiime2, meta, by="row.names", all=TRUE)
#sort by diet group
qiime2 <- qiime2[order(qiime2$Diet),] 
row.names(qiime2) <- qiime2[,1]
qiime2 <- qiime2[,-1]
qiime2_log <- log10(qiime2[,1:6]) #include just numbers, not diet group
qiime2_log[qiime2_log < 0] <- 0
qiime2_log$Diet <- qiime2$Diet

#Bacteroidia

#checking for assumption of normality for each diet within a microbe
shapiro.test(qiime2_log[1:11,1]) #1:11 are the control diets, column 1 is Bacteroidia
hist(qiime2_log[1:11,1])
qqnorm(qiime2_log[1:11,1])
qqline(qiime2_log[1:11,1])

shapiro.test(qiime2_log[12:19,1]) #12:19 are the HF diets, column 1 is Bacteroidia
hist(qiime2_log[12:19,1])
qqnorm(qiime2_log[12:19,1])
qqline(qiime2_log[12:19,1])

shapiro.test(qiime2_log[20:28,1]) #20:28 are the HF+EPA diets, column 1 is Bacteroidia
hist(qiime2_log[20:28,1])
qqnorm(qiime2_log[20:28,1])
qqline(qiime2_log[20:28,1])

compare_means(Bacteroidia ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF"), c("C", "HF+EPA"))
p1 <- ggboxplot(qiime2_log, x = "Diet", y = "Bacteroidia",
                title = "Bacteroidia", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Bacteroidia), max(qiime2_log$Bacteroidia)+1)) + 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Bacilli

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,2])
hist(qiime2_log[1:11,2])
qqnorm(qiime2_log[1:11,2])
qqline(qiime2_log[1:11,2])

shapiro.test(qiime2_log[12:19,2]) 
hist(qiime2_log[12:19,2])
qqnorm(qiime2_log[12:19,2])
qqline(qiime2_log[12:19,2])

shapiro.test(qiime2_log[20:28,2]) 
hist(qiime2_log[20:28,2])
qqnorm(qiime2_log[20:28,2])
qqline(qiime2_log[20:28,2])

#qiime2_log$Coriobacteriia[qiime2_log$Coriobacteriia == 0] <- NA
#qiime2_log_new <- na.omit(qiime2_log)

compare_means(Bacilli ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("HF", "HF+EPA"))
p2 <- ggboxplot(qiime2_log, x = "Diet", y = "Bacilli",
                title = "Bacilli", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Bacilli), max(qiime2_log$Bacilli)+0.4)) + 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Clostridia

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,3])
hist(qiime2_log[1:11,3])
qqnorm(qiime2_log[1:11,3])
qqline(qiime2_log[1:11,3])

shapiro.test(qiime2_log[12:19,3]) 
hist(qiime2_log[12:19,3])
qqnorm(qiime2_log[12:19,3])
qqline(qiime2_log[12:19,3])

shapiro.test(qiime2_log[20:28,3]) 
hist(qiime2_log[20:28,3])
qqnorm(qiime2_log[20:28,3])
qqline(qiime2_log[20:28,3])

compare_means(Clostridia ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
p3 <- ggboxplot(qiime2_log, x = "Diet", y = "Clostridia",
                title = "Clostridia", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Clostridia), max(qiime2_log$Clostridia)+0.4)) + 
  theme(legend.position='none') #palette = "jco"

#Gammaproteobacteria

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,4])
hist(qiime2_log[1:11,4])
qqnorm(qiime2_log[1:11,4])
qqline(qiime2_log[1:11,4])

shapiro.test(qiime2_log[12:19,4]) 
hist(qiime2_log[12:19,4])
qqnorm(qiime2_log[12:19,4])
qqline(qiime2_log[12:19,4])

shapiro.test(qiime2_log[20:28,4]) 
hist(qiime2_log[20:28,4])
qqnorm(qiime2_log[20:28,4])
qqline(qiime2_log[20:28,4])

#The 5th control mouse is weird - he is void of any proteobacteria
compare_means(Gammaproteobacteria ~ Diet,  data = qiime2_log[c(1:4,6:28),], method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("HF", "HF_EPA"), c("C", "HF_EPA"))
p4 <- ggboxplot(qiime2_log[c(1:4,6:28),], x = "Diet", y = "Gammaproteobacteria",
                title = "Gammaproteobacteria", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log[c(1:4,6:28),]$Gammaproteobacteria), max(qiime2_log[c(1:4,6:28),]$Gammaproteobacteria)+0.4)) + 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Mollicutes

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,5])
hist(qiime2_log[1:11,5])
qqnorm(qiime2_log[1:11,5])
qqline(qiime2_log[1:11,5])

shapiro.test(qiime2_log[12:19,5]) 
hist(qiime2_log[12:19,5])
qqnorm(qiime2_log[12:19,5])
qqline(qiime2_log[12:19,5])

shapiro.test(qiime2_log[20:28,5]) 
hist(qiime2_log[20:28,5])
qqnorm(qiime2_log[20:28,5])
qqline(qiime2_log[20:28,5])

compare_means(Mollicutes ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF_EPA"))
p5 <- ggboxplot(qiime2_log, x = "Diet", y = "Mollicutes",
                title = "Mollicutes", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Mollicutes), max(qiime2_log$Mollicutes)+0.2))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Verrucomicrobiae

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,6])
hist(qiime2_log[1:11,6])
qqnorm(qiime2_log[1:11,6])
qqline(qiime2_log[1:11,6])

shapiro.test(qiime2_log[12:19,6]) 
hist(qiime2_log[12:19,6])
qqnorm(qiime2_log[12:19,6])
qqline(qiime2_log[12:19,6])

shapiro.test(qiime2_log[20:28,6]) 
hist(qiime2_log[20:28,6])
qqnorm(qiime2_log[20:28,6])
qqline(qiime2_log[20:28,6])

compare_means(Verrucomicrobiae ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF_EPA"))
p6 <- ggboxplot(qiime2_log, x = "Diet", y = "Verrucomicrobiae",
                title = "Verrucomicrobiae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Verrucomicrobiae), max(qiime2_log$Verrucomicrobiae)+0.7))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

ggarrange(p1,p2,p3,p4,p5,p6, ncol = 2, nrow = 3, legend = "none")

###Making ORDER Boxplots####

qiime <- read.csv("Normalized_DESeq2_SequenceCountsComboData.csv", header = FALSE, row.names = 1)
qiime <- qiime[c(1,14:20),]
qiime <- data.frame(t(qiime)) #transpose the dataframe
qiime$V1 <- substring(qiime$V1, 1) #remove the first character of the first column strings
#assign row names
qiime2 <- qiime[,-1]
qiime2 <- mutate_all(qiime2, function(x) as.numeric(as.character(x)))
rownames(qiime2) <- qiime[,1]
#merge dataframe by sample IDs to append diet information 
qiime2 <- merge(qiime2, meta, by="row.names", all=TRUE)
#sort by diet group
qiime2 <- qiime2[order(qiime2$Diet),] 
row.names(qiime2) <- qiime2[,1]
qiime2 <- qiime2[,-1]
qiime2_log <- log10(qiime2[,1:7]) #include just numbers, not diet group
qiime2_log[qiime2_log < 0] <- 0
qiime2_log$Diet <- qiime2$Diet

#Bacteroidales

#checking for assumption of normality for each diet within a microbe
shapiro.test(qiime2_log[1:11,1]) 
hist(qiime2_log[1:11,1])
qqnorm(qiime2_log[1:11,1])
qqline(qiime2_log[1:11,1])

shapiro.test(qiime2_log[12:19,1])
hist(qiime2_log[12:19,1])
qqnorm(qiime2_log[12:19,1])
qqline(qiime2_log[12:19,1])

shapiro.test(qiime2_log[20:28,1]) 
hist(qiime2_log[20:28,1])
qqnorm(qiime2_log[20:28,1])
qqline(qiime2_log[20:28,1])

compare_means(Bacteroidales ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("HF", "HF+EPA"))
p1 <- ggboxplot(qiime2_log, x = "Diet", y = "Bacteroidales",
                title = "Bacteroidales", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Bacteroidales), max(qiime2_log$Bacteroidales)+1)) + 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Lactobacillales

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,2])
hist(qiime2_log[1:11,2])
qqnorm(qiime2_log[1:11,2])
qqline(qiime2_log[1:11,2])

shapiro.test(qiime2_log[12:19,2]) 
hist(qiime2_log[12:19,2])
qqnorm(qiime2_log[12:19,2])
qqline(qiime2_log[12:19,2])

shapiro.test(qiime2_log[20:28,2]) 
hist(qiime2_log[20:28,2])
qqnorm(qiime2_log[20:28,2])
qqline(qiime2_log[20:28,2])

compare_means(Lactobacillales ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("Control", "High Fat"), c("Control", "High Fat + EPA"))
p2 <- ggboxplot(qiime2_log, x = "Diet", y = "Lactobacillales",
                title = "Lactobacillales", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Lactobacillales), max(qiime2_log$Lactobacillales)+0.4)) + 
  theme(legend.position='none') #palette = "jco"
#stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Clostridiales

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,3])
hist(qiime2_log[1:11,3])
qqnorm(qiime2_log[1:11,3])
qqline(qiime2_log[1:11,3])

shapiro.test(qiime2_log[12:19,3]) 
hist(qiime2_log[12:19,3])
qqnorm(qiime2_log[12:19,3])
qqline(qiime2_log[12:19,3])

shapiro.test(qiime2_log[20:28,3]) 
hist(qiime2_log[20:28,3])
qqnorm(qiime2_log[20:28,3])
qqline(qiime2_log[20:28,3])

compare_means(Clostridiales ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
p3 <- ggboxplot(qiime2_log, x = "Diet", y = "Clostridiales",
                title = "Clostridiales", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Clostridiales), max(qiime2_log$Clostridiales)+0.4)) + 
  theme(legend.position='none') #palette = "jco"

#Erysipelotrichales

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,4])
hist(qiime2_log[1:11,4])
qqnorm(qiime2_log[1:11,4])
qqline(qiime2_log[1:11,4])

shapiro.test(qiime2_log[12:19,4]) 
hist(qiime2_log[12:19,4])
qqnorm(qiime2_log[12:19,4])
qqline(qiime2_log[12:19,4])

shapiro.test(qiime2_log[20:28,4]) 
hist(qiime2_log[20:28,4])
qqnorm(qiime2_log[20:28,4])
qqline(qiime2_log[20:28,4])

compare_means(Erysipelotrichales ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("C", "HF_EPA"))
p4 <- ggboxplot(qiime2_log, x = "Diet", y = "Erysipelotrichales",
                title = "Erysipelotrichales", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Erysipelotrichales), max(qiime2_log$Erysipelotrichales)+0.7)) + 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Betaproteobacteriales

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,5])
hist(qiime2_log[1:11,5])
qqnorm(qiime2_log[1:11,5])
qqline(qiime2_log[1:11,5])

shapiro.test(qiime2_log[12:19,5]) 
hist(qiime2_log[12:19,5])
qqnorm(qiime2_log[12:19,5])
qqline(qiime2_log[12:19,5])

shapiro.test(qiime2_log[20:28,5]) 
hist(qiime2_log[20:28,5])
qqnorm(qiime2_log[20:28,5])
qqline(qiime2_log[20:28,5])

compare_means(Betaproteobacteriales ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("HF", "HF_EPA"), c("C", "HF_EPA"))
p5 <- ggboxplot(qiime2_log, x = "Diet", y = "Betaproteobacteriales",
                title = "Betaproteobacteriales", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Betaproteobacteriales), max(qiime2_log$Betaproteobacteriales)+1))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Mollicutes_RF39

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,6])
hist(qiime2_log[1:11,6])
qqnorm(qiime2_log[1:11,6])
qqline(qiime2_log[1:11,6])

shapiro.test(qiime2_log[12:19,6]) 
hist(qiime2_log[12:19,6])
qqnorm(qiime2_log[12:19,6])
qqline(qiime2_log[12:19,6])

shapiro.test(qiime2_log[20:28,6]) 
hist(qiime2_log[20:28,6])
qqnorm(qiime2_log[20:28,6])
qqline(qiime2_log[20:28,6])

compare_means(Mollicutes_RF39 ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF"), c("C", "HF+EPA"))
p6 <- ggboxplot(qiime2_log, x = "Diet", y = "Mollicutes_RF39",
                title = "Mollicutes_RF39", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Mollicutes_RF39), max(qiime2_log$Mollicutes_RF39)+0.2))+ 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Verrucomicrobiales

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,7])
hist(qiime2_log[1:11,7])
qqnorm(qiime2_log[1:11,7])
qqline(qiime2_log[1:11,7])

shapiro.test(qiime2_log[12:19,7]) 
hist(qiime2_log[12:19,7])
qqnorm(qiime2_log[12:19,7])
qqline(qiime2_log[12:19,7])

shapiro.test(qiime2_log[20:28,7]) 
hist(qiime2_log[20:28,7])
qqnorm(qiime2_log[20:28,7])
qqline(qiime2_log[20:28,7])

compare_means(Verrucomicrobiales ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF_EPA"))
p7 <- ggboxplot(qiime2_log, x = "Diet", y = "Verrucomicrobiales",
                title = "Verrucomicrobiales", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Verrucomicrobiales), max(qiime2_log$Verrucomicrobiales)+0.2))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

# #Anaeroplasmatales
# 
# #checking for assumption of normality
# shapiro.test(qiime2_log[1:11,8])
# hist(qiime2_log[1:11,8])
# qqnorm(qiime2_log[1:11,8])
# qqline(qiime2_log[1:11,8])
# 
# shapiro.test(qiime2_log[12:19,8]) 
# hist(qiime2_log[12:19,8])
# qqnorm(qiime2_log[12:19,8])
# qqline(qiime2_log[12:19,8])
# 
# shapiro.test(qiime2_log[20:28,8]) 
# hist(qiime2_log[20:28,8])
# qqnorm(qiime2_log[20:28,8])
# qqline(qiime2_log[20:28,8])
# 
# compare_means(Anaeroplasmatales ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
# my_comparisons <- list(c("C", "HF"), c("C", "HF+EPA"), c("HF", "HF+EPA"))
# p8 <- ggboxplot(qiime2_log, x = "Diet", y = "Anaeroplasmatales",
#                 title = "Anaeroplasmatales", subtitle ="*42% of samples are zero-inflated", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
#                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Anaeroplasmatales), max(qiime2_log$Anaeroplasmatales)+1.2))+ 
#   theme(legend.position='none')+ #palette = "jco"
#   stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
# 
# #RF39
# 
# #checking for assumption of normality
# shapiro.test(qiime2_log[1:11,9])
# hist(qiime2_log[1:11,9])
# qqnorm(qiime2_log[1:11,9])
# qqline(qiime2_log[1:11,9])
# 
# shapiro.test(qiime2_log[12:19,9]) 
# hist(qiime2_log[12:19,9])
# qqnorm(qiime2_log[12:19,9])
# qqline(qiime2_log[12:19,9])
# 
# shapiro.test(qiime2_log[20:28,9]) 
# hist(qiime2_log[20:28,9])
# qqnorm(qiime2_log[20:28,9])
# qqline(qiime2_log[20:28,9])
# 
# compare_means(RF39 ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
# #my_comparisons <- list(c("C", "HF"), c("C", "HF+EPA"))
# p9 <- ggboxplot(qiime2_log, x = "Diet", y = "RF39",
#                 title = "RF39", subtitle = "Phylum: Tenericutes; Class: Mollicutes", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
#                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$RF39), max(qiime2_log$RF39)+0.2))+ 
#   theme(legend.position='none') #palette = "jco"
# #  stat_compare_means(comparisons = my_comparisons, method = "t.test")
# 
# #Verrucomicrobiales
# 
# #checking for assumption of normality
# shapiro.test(qiime2_log[1:11,10])
# hist(qiime2_log[1:11,10])
# qqnorm(qiime2_log[1:11,10])
# qqline(qiime2_log[1:11,10])
# 
# shapiro.test(qiime2_log[12:19,10]) 
# hist(qiime2_log[12:19,10])
# qqnorm(qiime2_log[12:19,10])
# qqline(qiime2_log[12:19,10])
# 
# shapiro.test(qiime2_log[20:28,10]) 
# hist(qiime2_log[20:28,10])
# qqnorm(qiime2_log[20:28,10])
# qqline(qiime2_log[20:28,10])
# 
# compare_means(Verrucomicrobiales ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
# my_comparisons <- list(c("C", "HF"), c("C", "HF+EPA"))
# p10 <- ggboxplot(qiime2_log, x = "Diet", y = "Verrucomicrobiales",
#                 title = "Verrucomicrobiales", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
#                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Verrucomicrobiales), max(qiime2_log$Verrucomicrobiales)+0.4))+ 
#   theme(legend.position='none')+ #palette = "jco"
#   stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Combo plots
ggarrange(p1,p2,p3,p4,p5,p6, ncol = 2, nrow = 3, legend = "none")
ggarrange(p7, ncol = 2, nrow = 3, legend = "none")

###Making FAMILY Boxplots####

qiime <- read.csv("Normalized_DESeq2_SequenceCountsComboData.csv", header = FALSE, row.names = 1)
qiime <- qiime[c(1,21:34),]
qiime <- data.frame(t(qiime)) #transpose the dataframe
qiime$V1 <- substring(qiime$V1, 1) #remove the first character of the first column strings
#assign row names
qiime2 <- qiime[,-1]
qiime2 <- mutate_all(qiime2, function(x) as.numeric(as.character(x)))
rownames(qiime2) <- qiime[,1]
#merge dataframe by sample IDs to append diet information 
qiime2 <- merge(qiime2, meta, by="row.names", all=TRUE)
#sort by diet group
qiime2 <- qiime2[order(qiime2$Diet),] 
row.names(qiime2) <- qiime2[,1]
qiime2 <- qiime2[,-1]
qiime2_log <- log10(qiime2[,1:14]) #include just numbers, not diet group
qiime2_log[qiime2_log < 0] <- 0
qiime2_log$Diet <- qiime2$Diet

#Bacteroidaceae

#checking for assumption of normality for each diet within a microbe
shapiro.test(qiime2_log[1:11,1]) 
hist(qiime2_log[1:11,1])
qqnorm(qiime2_log[1:11,1])
qqline(qiime2_log[1:11,1])

shapiro.test(qiime2_log[12:19,1])
hist(qiime2_log[12:19,1])
qqnorm(qiime2_log[12:19,1])
qqline(qiime2_log[12:19,1])

shapiro.test(qiime2_log[20:28,1]) 
hist(qiime2_log[20:28,1])
qqnorm(qiime2_log[20:28,1])
qqline(qiime2_log[20:28,1])

compare_means(Bacteroidaceae ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("C", "HF_EPA"))
p1 <- ggboxplot(qiime2_log, x = "Diet", y = "Bacteroidaceae",
                title = "Bacteroidaceae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Bacteroidaceae), max(qiime2_log$Bacteroidaceae)+0.5)) + 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Muribaculaceae

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,2])
hist(qiime2_log[1:11,2])
qqnorm(qiime2_log[1:11,2])
qqline(qiime2_log[1:11,2])

shapiro.test(qiime2_log[12:19,2]) 
hist(qiime2_log[12:19,2])
qqnorm(qiime2_log[12:19,2])
qqline(qiime2_log[12:19,2])

shapiro.test(qiime2_log[20:28,2]) 
hist(qiime2_log[20:28,2])
qqnorm(qiime2_log[20:28,2])
qqline(qiime2_log[20:28,2])

compare_means(Muribaculaceae ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF"), c("C", "HF+EPA"))
p2 <- ggboxplot(qiime2_log, x = "Diet", y = "Muribaculaceae",
                title = "Muribaculaceae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Muribaculaceae), max(qiime2_log$Muribaculaceae)+0.4)) + 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Lactobacillaceae

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,3])
hist(qiime2_log[1:11,3])
qqnorm(qiime2_log[1:11,3])
qqline(qiime2_log[1:11,3])

shapiro.test(qiime2_log[12:19,3]) 
hist(qiime2_log[12:19,3])
qqnorm(qiime2_log[12:19,3])
qqline(qiime2_log[12:19,3])

shapiro.test(qiime2_log[20:28,3]) 
hist(qiime2_log[20:28,3])
qqnorm(qiime2_log[20:28,3])
qqline(qiime2_log[20:28,3])

compare_means(Lactobacillaceae ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"))
p3 <- ggboxplot(qiime2_log, x = "Diet", y = "Lactobacillaceae",
                title = "Lactobacillaceae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Lactobacillaceae), max(qiime2_log$Lactobacillaceae)+0.4)) + 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Streptococcaceae

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,4])
hist(qiime2_log[1:11,4])
qqnorm(qiime2_log[1:11,4])
qqline(qiime2_log[1:11,4])

shapiro.test(qiime2_log[12:19,4]) 
hist(qiime2_log[12:19,4])
qqnorm(qiime2_log[12:19,4])
qqline(qiime2_log[12:19,4])

shapiro.test(qiime2_log[20:28,4]) 
hist(qiime2_log[20:28,4])
qqnorm(qiime2_log[20:28,4])
qqline(qiime2_log[20:28,4])

compare_means(Streptococcaceae ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF"))
p4 <- ggboxplot(qiime2_log, x = "Diet", y = "Streptococcaceae",
                title = "Streptococcaceae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Streptococcaceae), max(qiime2_log$Streptococcaceae)+0.7)) + 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Christensenellaceae

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,5])
hist(qiime2_log[1:11,5])
qqnorm(qiime2_log[1:11,5])
qqline(qiime2_log[1:11,5])

shapiro.test(qiime2_log[12:19,5]) 
hist(qiime2_log[12:19,5])
qqnorm(qiime2_log[12:19,5])
qqline(qiime2_log[12:19,5])

shapiro.test(qiime2_log[20:28,5]) 
hist(qiime2_log[20:28,5])
qqnorm(qiime2_log[20:28,5])
qqline(qiime2_log[20:28,5])

compare_means(Christensenellaceae ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF+EPA"))
p5 <- ggboxplot(qiime2_log, x = "Diet", y = "Christensenellaceae",
                title = "Christensenellaceae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Christensenellaceae), max(qiime2_log$Christensenellaceae)+0.2))+ 
  theme(legend.position='none') #palette = "jco"
#stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Clostridiaceae

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,6])
hist(qiime2_log[1:11,6])
qqnorm(qiime2_log[1:11,6])
qqline(qiime2_log[1:11,6])

shapiro.test(qiime2_log[12:19,6]) 
hist(qiime2_log[12:19,6])
qqnorm(qiime2_log[12:19,6])
qqline(qiime2_log[12:19,6])

shapiro.test(qiime2_log[20:28,6]) 
hist(qiime2_log[20:28,6])
qqnorm(qiime2_log[20:28,6])
qqline(qiime2_log[20:28,6])

compare_means(Clostridiaceae ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("C", "HF_EPA"))
p6 <- ggboxplot(qiime2_log, x = "Diet", y = "Clostridiaceae",
                title = "Clostridiaceae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Clostridiaceae), max(qiime2_log$Clostridiaceae)+0.8))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Clostridiales_vadinBB60

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,7])
hist(qiime2_log[1:11,7])
qqnorm(qiime2_log[1:11,7])
qqline(qiime2_log[1:11,7])

shapiro.test(qiime2_log[12:19,7]) 
hist(qiime2_log[12:19,7])
qqnorm(qiime2_log[12:19,7])
qqline(qiime2_log[12:19,7])

shapiro.test(qiime2_log[20:28,7]) 
hist(qiime2_log[20:28,7])
qqnorm(qiime2_log[20:28,7])
qqline(qiime2_log[20:28,7])

compare_means(Clostridiales_vadinBB60 ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("HF", "HF_EPA"), c("C", "HF_EPA"))
p7 <- ggboxplot(qiime2_log, x = "Diet", y = "Clostridiales_vadinBB60",
                title = "Clostridiales_vadinBB60", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Clostridiales_vadinBB60), max(qiime2_log$Clostridiales_vadinBB60)+1))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Lachnospiraceae

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,8])
hist(qiime2_log[1:11,8])
qqnorm(qiime2_log[1:11,8])
qqline(qiime2_log[1:11,8])

shapiro.test(qiime2_log[12:19,8]) 
hist(qiime2_log[12:19,8])
qqnorm(qiime2_log[12:19,8])
qqline(qiime2_log[12:19,8])

shapiro.test(qiime2_log[20:28,8]) 
hist(qiime2_log[20:28,8])
qqnorm(qiime2_log[20:28,8])
qqline(qiime2_log[20:28,8])

compare_means(Lachnospiraceae ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF"), c("C", "HF+EPA"))
p8 <- ggboxplot(qiime2_log, x = "Diet", y = "Lachnospiraceae",
                title = "Lachnospiraceae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Lachnospiraceae), max(qiime2_log$Lachnospiraceae)))+ 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Peptococcaceae

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,9])
hist(qiime2_log[1:11,9])
qqnorm(qiime2_log[1:11,9])
qqline(qiime2_log[1:11,9])

shapiro.test(qiime2_log[12:19,9]) 
hist(qiime2_log[12:19,9])
qqnorm(qiime2_log[12:19,9])
qqline(qiime2_log[12:19,9])

shapiro.test(qiime2_log[20:28,9]) 
hist(qiime2_log[20:28,9])
qqnorm(qiime2_log[20:28,9])
qqline(qiime2_log[20:28,9])

compare_means(Peptococcaceae ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF_EPA"))
p9 <- ggboxplot(qiime2_log, x = "Diet", y = "Peptococcaceae",
                title = "Peptococcaceae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Peptococcaceae), max(qiime2_log$Peptococcaceae)+0.2))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Peptostreptococcaceae

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,10])
hist(qiime2_log[1:11,10])
qqnorm(qiime2_log[1:11,10])
qqline(qiime2_log[1:11,10])

shapiro.test(qiime2_log[12:19,10]) 
hist(qiime2_log[12:19,10])
qqnorm(qiime2_log[12:19,10])
qqline(qiime2_log[12:19,10])

shapiro.test(qiime2_log[20:28,10]) 
hist(qiime2_log[20:28,10])
qqnorm(qiime2_log[20:28,10])
qqline(qiime2_log[20:28,10])

compare_means(Peptostreptococcaceae ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF_EPA"))
p10 <- ggboxplot(qiime2_log, x = "Diet", y = "Peptostreptococcaceae",
                 title = "Peptostreptococcaceae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Peptostreptococcaceae), max(qiime2_log$Peptostreptococcaceae)+0.2))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Ruminococcaceae

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,11])
hist(qiime2_log[1:11,11])
qqnorm(qiime2_log[1:11,11])
qqline(qiime2_log[1:11,11])

shapiro.test(qiime2_log[12:19,11]) 
hist(qiime2_log[12:19,11])
qqnorm(qiime2_log[12:19,11])
qqline(qiime2_log[12:19,11])

shapiro.test(qiime2_log[20:28,11]) 
hist(qiime2_log[20:28,11])
qqnorm(qiime2_log[20:28,11])
qqline(qiime2_log[20:28,11])

compare_means(Ruminococcaceae ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"))
p11 <- ggboxplot(qiime2_log, x = "Diet", y = "Ruminococcaceae",
                 title = "Ruminococcaceae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Ruminococcaceae), max(qiime2_log$Ruminococcaceae)+0.4))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Erysipelotrichaceae

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,12])
hist(qiime2_log[1:11,12])
qqnorm(qiime2_log[1:11,12])
qqline(qiime2_log[1:11,12])

shapiro.test(qiime2_log[12:19,12]) 
hist(qiime2_log[12:19,12])
qqnorm(qiime2_log[12:19,12])
qqline(qiime2_log[12:19,12])

shapiro.test(qiime2_log[20:28,12]) 
hist(qiime2_log[20:28,12])
qqnorm(qiime2_log[20:28,12])
qqline(qiime2_log[20:28,12])

compare_means(Erysipelotrichaceae ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("C", "HF_EPA"))
p12 <- ggboxplot(qiime2_log, x = "Diet", y = "Erysipelotrichaceae",
                 title = "Erysipelotrichaceae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Erysipelotrichaceae), max(qiime2_log$Erysipelotrichaceae)+0.8))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Burkholderiaceae

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,13])
hist(qiime2_log[1:11,13])
qqnorm(qiime2_log[1:11,13])
qqline(qiime2_log[1:11,13])

shapiro.test(qiime2_log[12:19,13]) 
hist(qiime2_log[12:19,13])
qqnorm(qiime2_log[12:19,13])
qqline(qiime2_log[12:19,13])

shapiro.test(qiime2_log[20:28,13]) 
hist(qiime2_log[20:28,13])
qqnorm(qiime2_log[20:28,13])
qqline(qiime2_log[20:28,13])

compare_means(Burkholderiaceae ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("HF","HF_EPA"), c("C", "HF_EPA"))
p13 <- ggboxplot(qiime2_log, x = "Diet", y = "Burkholderiaceae",
                 title = "Burkholderiaceae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Burkholderiaceae), max(qiime2_log$Burkholderiaceae)+1))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Akkermansiaceae

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,14])
hist(qiime2_log[1:11,14])
qqnorm(qiime2_log[1:11,14])
qqline(qiime2_log[1:11,14])

shapiro.test(qiime2_log[12:19,14]) 
hist(qiime2_log[12:19,14])
qqnorm(qiime2_log[12:19,14])
qqline(qiime2_log[12:19,14])

shapiro.test(qiime2_log[20:28,14]) 
hist(qiime2_log[20:28,14])
qqnorm(qiime2_log[20:28,14])
qqline(qiime2_log[20:28,14])

compare_means(Akkermansiaceae ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF_EPA"))
p14 <- ggboxplot(qiime2_log, x = "Diet", y = "Akkermansiaceae",
                 title = "Akkermansiaceae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Akkermansiaceae), max(qiime2_log$Akkermansiaceae)+0.2))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

# #Erysipelotrichaceae
# 
# #checking for assumption of normality
# shapiro.test(qiime2_log[1:11,15])
# hist(qiime2_log[1:11,15])
# qqnorm(qiime2_log[1:11,15])
# qqline(qiime2_log[1:11,15])
# 
# shapiro.test(qiime2_log[12:19,15]) 
# hist(qiime2_log[12:19,15])
# qqnorm(qiime2_log[12:19,15])
# qqline(qiime2_log[12:19,15])
# 
# shapiro.test(qiime2_log[20:28,15]) 
# hist(qiime2_log[20:28,15])
# qqnorm(qiime2_log[20:28,15])
# qqline(qiime2_log[20:28,15])
# 
# compare_means(Erysipelotrichaceae ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
# my_comparisons <- list(c("C", "HF"), c("C", "HF+EPA"))
# p15 <- ggboxplot(qiime2_log, x = "Diet", y = "Erysipelotrichaceae",
#                  title = "Erysipelotrichaceae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
#                  color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Erysipelotrichaceae), max(qiime2_log$Erysipelotrichaceae)+0.8))+ 
#   theme(legend.position='none')+ #palette = "jco"
#   stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
# 
# #Alcaligenaceae
# 
# #checking for assumption of normality
# shapiro.test(qiime2_log[1:11,16])
# hist(qiime2_log[1:11,16])
# qqnorm(qiime2_log[1:11,16])
# qqline(qiime2_log[1:11,16])
# 
# shapiro.test(qiime2_log[12:19,16]) 
# hist(qiime2_log[12:19,16])
# qqnorm(qiime2_log[12:19,16])
# qqline(qiime2_log[12:19,16])
# 
# shapiro.test(qiime2_log[20:28,16]) 
# hist(qiime2_log[20:28,16])
# qqnorm(qiime2_log[20:28,16])
# qqline(qiime2_log[20:28,16])
# 
# compare_means(Alcaligenaceae ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
# my_comparisons <- list(c("C", "HF"), c("C", "HF+EPA"))
# p16 <- ggboxplot(qiime2_log, x = "Diet", y = "Alcaligenaceae",
#                  title = "Alcaligenaceae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
#                  color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Alcaligenaceae), max(qiime2_log$Alcaligenaceae)+0.8))+ 
#   theme(legend.position='none')+ #palette = "jco"
#   stat_compare_means(comparisons = my_comparisons, method = "t.test")
# 
# #Verrucomicrobiaceae
# 
# #checking for assumption of normality
# shapiro.test(qiime2_log[1:11,17])
# hist(qiime2_log[1:11,17])
# qqnorm(qiime2_log[1:11,17])
# qqline(qiime2_log[1:11,17])
# 
# shapiro.test(qiime2_log[12:19,17]) 
# hist(qiime2_log[12:19,17])
# qqnorm(qiime2_log[12:19,17])
# qqline(qiime2_log[12:19,17])
# 
# shapiro.test(qiime2_log[20:28,17]) 
# hist(qiime2_log[20:28,17])
# qqnorm(qiime2_log[20:28,17])
# qqline(qiime2_log[20:28,17])
# 
# compare_means(Verrucomicrobiaceae ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
# my_comparisons <- list(c("C", "HF"), c("C", "HF+EPA"))
# p17 <- ggboxplot(qiime2_log, x = "Diet", y = "Verrucomicrobiaceae",
#                  title = "Verrucomicrobiaceae", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
#                  color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Verrucomicrobiaceae), max(qiime2_log$Verrucomicrobiaceae)+0.6))+ 
#   theme(legend.position='none')+ #palette = "jco"
#   stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Combo plots
ggarrange(p1,p2,p3,p4,p5,p6, ncol = 2, nrow = 3, legend = "none")
ggarrange(p7,p8,p9,p10,p11,p12, ncol = 2, nrow = 3, legend = "none")
ggarrange(p13,p14, ncol = 2, nrow = 3, legend = "none")

###Making GENUS Boxplots####

qiime <- read.csv("Normalized_DESeq2_SequenceCountsComboData.csv", header = FALSE, row.names = 1)
qiime <- qiime[c(1,35:55),]
qiime <- data.frame(t(qiime)) #transpose the dataframe
qiime$V1 <- substring(qiime$V1, 1) #remove the first character of the first column strings
#assign row names
qiime2 <- qiime[,-1]
qiime2 <- mutate_all(qiime2, function(x) as.numeric(as.character(x)))
rownames(qiime2) <- qiime[,1]
#merge dataframe by sample IDs to append diet information 
qiime2 <- merge(qiime2, meta, by="row.names", all=TRUE)
#sort by diet group
qiime2 <- qiime2[order(qiime2$Diet),] 
row.names(qiime2) <- qiime2[,1]
qiime2 <- qiime2[,-1]
qiime2_log <- log10(qiime2[,1:21]) #include just numbers, not diet group
qiime2_log[qiime2_log < 0] <- 0
qiime2_log$Diet <- qiime2$Diet

#Bacteroides

#checking for assumption of normality for each diet within a microbe
shapiro.test(qiime2_log[1:11,1]) 
hist(qiime2_log[1:11,1])
qqnorm(qiime2_log[1:11,1])
qqline(qiime2_log[1:11,1])

shapiro.test(qiime2_log[12:19,1])
hist(qiime2_log[12:19,1])
qqnorm(qiime2_log[12:19,1])
qqline(qiime2_log[12:19,1])

shapiro.test(qiime2_log[20:28,1]) 
hist(qiime2_log[20:28,1])
qqnorm(qiime2_log[20:28,1])
qqline(qiime2_log[20:28,1])

compare_means(Bacteroides ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("C", "HF_EPA"))
p1 <- ggboxplot(qiime2_log, x = "Diet", y = "Bacteroides",
                title = "Bacteroides", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Bacteroides), max(qiime2_log$Bacteroides)+0.4)) + 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Lactobacillus

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,2])
hist(qiime2_log[1:11,2])
qqnorm(qiime2_log[1:11,2])
qqline(qiime2_log[1:11,2])

shapiro.test(qiime2_log[12:19,2]) 
hist(qiime2_log[12:19,2])
qqnorm(qiime2_log[12:19,2])
qqline(qiime2_log[12:19,2])

shapiro.test(qiime2_log[20:28,2]) 
hist(qiime2_log[20:28,2])
qqnorm(qiime2_log[20:28,2])
qqline(qiime2_log[20:28,2])

compare_means(Lactobacillus ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"))
p2 <- ggboxplot(qiime2_log, x = "Diet", y = "Lactobacillus",
                title = "Lactobacillus", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Lactobacillus), max(qiime2_log$Lactobacillus)+0.4)) + 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Lactococcus

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,3])
hist(qiime2_log[1:11,3])
qqnorm(qiime2_log[1:11,3])
qqline(qiime2_log[1:11,3])

shapiro.test(qiime2_log[12:19,3]) 
hist(qiime2_log[12:19,3])
qqnorm(qiime2_log[12:19,3])
qqline(qiime2_log[12:19,3])

shapiro.test(qiime2_log[20:28,3]) 
hist(qiime2_log[20:28,3])
qqnorm(qiime2_log[20:28,3])
qqline(qiime2_log[20:28,3])

compare_means(Lactococcus ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF"))
p3 <- ggboxplot(qiime2_log, x = "Diet", y = "Lactococcus",
                title = "Lactococcus", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Lactococcus), max(qiime2_log$Lactococcus))) + 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Eubacterium_brachy

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,4])
hist(qiime2_log[1:11,4])
qqnorm(qiime2_log[1:11,4])
qqline(qiime2_log[1:11,4])

shapiro.test(qiime2_log[12:19,4]) 
hist(qiime2_log[12:19,4])
qqnorm(qiime2_log[12:19,4])
qqline(qiime2_log[12:19,4])

shapiro.test(qiime2_log[20:28,4]) 
hist(qiime2_log[20:28,4])
qqnorm(qiime2_log[20:28,4])
qqline(qiime2_log[20:28,4])

compare_means(Eubacterium_brachy ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("C", "HF_EPA"))
p4 <- ggboxplot(qiime2_log, x = "Diet", y = "Eubacterium_brachy",
                title = "Eubacterium_brachy", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Eubacterium_brachy), max(qiime2_log$Eubacterium_brachy)+0.5)) + 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Eubacterium_nodatum

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,5])
hist(qiime2_log[1:11,5])
qqnorm(qiime2_log[1:11,5])
qqline(qiime2_log[1:11,5])

shapiro.test(qiime2_log[12:19,5]) 
hist(qiime2_log[12:19,5])
qqnorm(qiime2_log[12:19,5])
qqline(qiime2_log[12:19,5])

shapiro.test(qiime2_log[20:28,5]) 
hist(qiime2_log[20:28,5])
qqnorm(qiime2_log[20:28,5])
qqline(qiime2_log[20:28,5])

#removed one control with zero count (outlier that skewed the plot)
compare_means(Eubacterium_nodatum ~ Diet,  data = qiime2_log[c(1:2,4:28),], method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("C", "HF_EPA"))
p5 <- ggboxplot(qiime2_log[c(1:2,4:28),], x = "Diet", y = "Eubacterium_nodatum",
                title = "Eubacterium_nodatum", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log[c(1:2,4:28),]$Eubacterium_nodatum), max(qiime2_log[c(1:2,4:28),]$Eubacterium_nodatum)+0.5))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Acetatifactor

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,6])
hist(qiime2_log[1:11,6])
qqnorm(qiime2_log[1:11,6])
qqline(qiime2_log[1:11,6])

shapiro.test(qiime2_log[12:19,6]) 
hist(qiime2_log[12:19,6])
qqnorm(qiime2_log[12:19,6])
qqline(qiime2_log[12:19,6])

shapiro.test(qiime2_log[20:28,6]) 
hist(qiime2_log[20:28,6])
qqnorm(qiime2_log[20:28,6])
qqline(qiime2_log[20:28,6])

compare_means(Acetatifactor ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("HF", "HF_EPA"))
p6 <- ggboxplot(qiime2_log, x = "Diet", y = "Acetatifactor",
                title = "Acetatifactor", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Acetatifactor), max(qiime2_log$Acetatifactor)+0.2))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Blautia

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,7])
hist(qiime2_log[1:11,7])
qqnorm(qiime2_log[1:11,7])
qqline(qiime2_log[1:11,7])

shapiro.test(qiime2_log[12:19,7]) 
hist(qiime2_log[12:19,7])
qqnorm(qiime2_log[12:19,7])
qqline(qiime2_log[12:19,7])

shapiro.test(qiime2_log[20:28,7]) 
hist(qiime2_log[20:28,7])
qqnorm(qiime2_log[20:28,7])
qqline(qiime2_log[20:28,7])

compare_means(Blautia ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF"), c("HF", "HF+EPA"))
p7 <- ggboxplot(qiime2_log, x = "Diet", y = "Blautia",
                title = "Blautia", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Blautia), max(qiime2_log$Blautia)))+ 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Lachnoclostridium

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,8])
hist(qiime2_log[1:11,8])
qqnorm(qiime2_log[1:11,8])
qqline(qiime2_log[1:11,8])

shapiro.test(qiime2_log[12:19,8]) 
hist(qiime2_log[12:19,8])
qqnorm(qiime2_log[12:19,8])
qqline(qiime2_log[12:19,8])

shapiro.test(qiime2_log[20:28,8]) 
hist(qiime2_log[20:28,8])
qqnorm(qiime2_log[20:28,8])
qqline(qiime2_log[20:28,8])

compare_means(Lachnoclostridium ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("C", "HF_EPA"))
p8 <- ggboxplot(qiime2_log, x = "Diet", y = "Lachnoclostridium",
                title = "Lachnoclostridium", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Lachnoclostridium), max(qiime2_log$Lachnoclostridium)+0.5))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Lachnospiraceae_FCS020

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,9])
hist(qiime2_log[1:11,9])
qqnorm(qiime2_log[1:11,9])
qqline(qiime2_log[1:11,9])

shapiro.test(qiime2_log[12:19,9]) 
hist(qiime2_log[12:19,9])
qqnorm(qiime2_log[12:19,9])
qqline(qiime2_log[12:19,9])

shapiro.test(qiime2_log[20:28,9]) 
hist(qiime2_log[20:28,9])
qqnorm(qiime2_log[20:28,9])
qqline(qiime2_log[20:28,9])

compare_means(Lachnospiraceae_FCS020 ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("HF", "HF_EPA"))
p9 <- ggboxplot(qiime2_log, x = "Diet", y = "Lachnospiraceae_FCS020",
                title = "Lachnospiraceae_FCS020", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Lachnospiraceae_FCS020), max(qiime2_log$Lachnospiraceae_FCS020)+0.5))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Lachnospiraceae_NK4A136

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,10])
hist(qiime2_log[1:11,10])
qqnorm(qiime2_log[1:11,10])
qqline(qiime2_log[1:11,10])

shapiro.test(qiime2_log[12:19,10]) 
hist(qiime2_log[12:19,10])
qqnorm(qiime2_log[12:19,10])
qqline(qiime2_log[12:19,10])

shapiro.test(qiime2_log[20:28,10]) 
hist(qiime2_log[20:28,10])
qqnorm(qiime2_log[20:28,10])
qqline(qiime2_log[20:28,10])

compare_means(Lachnospiraceae_NK4A136 ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF"), c("C", "HF+EPA"))
p10 <- ggboxplot(qiime2_log, x = "Diet", y = "Lachnospiraceae_NK4A136",
                 title = "Lachnospiraceae_NK4A136", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Lachnospiraceae_NK4A136), max(qiime2_log$Lachnospiraceae_NK4A136)+0.2))+ 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Tyzzerella

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,11])
hist(qiime2_log[1:11,11])
qqnorm(qiime2_log[1:11,11])
qqline(qiime2_log[1:11,11])

shapiro.test(qiime2_log[12:19,11]) 
hist(qiime2_log[12:19,11])
qqnorm(qiime2_log[12:19,11])
qqline(qiime2_log[12:19,11])

shapiro.test(qiime2_log[20:28,11]) 
hist(qiime2_log[20:28,11])
qqnorm(qiime2_log[20:28,11])
qqline(qiime2_log[20:28,11])

compare_means(Tyzzerella ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF+EPA"))
p11 <- ggboxplot(qiime2_log, x = "Diet", y = "Tyzzerella",
                 title = "Tyzzerella", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Tyzzerella), max(qiime2_log$Tyzzerella)+0.2))+ 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Romboutsia

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,12])
hist(qiime2_log[1:11,12])
qqnorm(qiime2_log[1:11,12])
qqline(qiime2_log[1:11,12])

shapiro.test(qiime2_log[12:19,12]) 
hist(qiime2_log[12:19,12])
qqnorm(qiime2_log[12:19,12])
qqline(qiime2_log[12:19,12])

shapiro.test(qiime2_log[20:28,12]) 
hist(qiime2_log[20:28,12])
qqnorm(qiime2_log[20:28,12])
qqline(qiime2_log[20:28,12])

compare_means(Romboutsia ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF_EPA"))
p12 <- ggboxplot(qiime2_log, x = "Diet", y = "Romboutsia",
                 title = "Romboutsia", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Romboutsia), max(qiime2_log$Romboutsia)+0.2))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Anaerotruncus

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,13])
hist(qiime2_log[1:11,13])
qqnorm(qiime2_log[1:11,13])
qqline(qiime2_log[1:11,13])

shapiro.test(qiime2_log[12:19,13]) 
hist(qiime2_log[12:19,13])
qqnorm(qiime2_log[12:19,13])
qqline(qiime2_log[12:19,13])

shapiro.test(qiime2_log[20:28,13]) 
hist(qiime2_log[20:28,13])
qqnorm(qiime2_log[20:28,13])
qqline(qiime2_log[20:28,13])

compare_means(Anaerotruncus ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF+EPA"))
p13 <- ggboxplot(qiime2_log, x = "Diet", y = "Anaerotruncus",
                 title = "Anaerotruncus", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Anaerotruncus), max(qiime2_log$Anaerotruncus)+0.2))+ 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Intestinimonas

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,14])
hist(qiime2_log[1:11,14])
qqnorm(qiime2_log[1:11,14])
qqline(qiime2_log[1:11,14])

shapiro.test(qiime2_log[12:19,14]) 
hist(qiime2_log[12:19,14])
qqnorm(qiime2_log[12:19,14])
qqline(qiime2_log[12:19,14])

shapiro.test(qiime2_log[20:28,14]) 
hist(qiime2_log[20:28,14])
qqnorm(qiime2_log[20:28,14])
qqline(qiime2_log[20:28,14])

compare_means(Intestinimonas ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF+EPA"))
p14 <- ggboxplot(qiime2_log, x = "Diet", y = "Intestinimonas",
                 title = "Intestinimonas", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Intestinimonas), max(qiime2_log$Intestinimonas)+0.2))+ 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Oscillibacter

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,15])
hist(qiime2_log[1:11,15])
qqnorm(qiime2_log[1:11,15])
qqline(qiime2_log[1:11,15])

shapiro.test(qiime2_log[12:19,15]) 
hist(qiime2_log[12:19,15])
qqnorm(qiime2_log[12:19,15])
qqline(qiime2_log[12:19,15])

shapiro.test(qiime2_log[20:28,15]) 
hist(qiime2_log[20:28,15])
qqnorm(qiime2_log[20:28,15])
qqline(qiime2_log[20:28,15])

compare_means(Oscillibacter ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF_EPA"))
p15 <- ggboxplot(qiime2_log, x = "Diet", y = "Oscillibacter",
                 title = "Oscillibacter", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Oscillibacter), max(qiime2_log$Oscillibacter)+0.2))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Ruminiclostridium

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,16])
hist(qiime2_log[1:11,16])
qqnorm(qiime2_log[1:11,16])
qqline(qiime2_log[1:11,16])

shapiro.test(qiime2_log[12:19,16]) 
hist(qiime2_log[12:19,16])
qqnorm(qiime2_log[12:19,16])
qqline(qiime2_log[12:19,16])

shapiro.test(qiime2_log[20:28,16]) 
hist(qiime2_log[20:28,16])
qqnorm(qiime2_log[20:28,16])
qqline(qiime2_log[20:28,16])

compare_means(Ruminiclostridium ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF+EPA"))
p16 <- ggboxplot(qiime2_log, x = "Diet", y = "Ruminiclostridium",
                 title = "Ruminiclostridium", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Ruminiclostridium), max(qiime2_log$Ruminiclostridium)+0.6))+ 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Ruminiclostridium_5

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,17])
hist(qiime2_log[1:11,17])
qqnorm(qiime2_log[1:11,17])
qqline(qiime2_log[1:11,17])

shapiro.test(qiime2_log[12:19,17]) 
hist(qiime2_log[12:19,17])
qqnorm(qiime2_log[12:19,17])
qqline(qiime2_log[12:19,17])

shapiro.test(qiime2_log[20:28,17]) 
hist(qiime2_log[20:28,17])
qqnorm(qiime2_log[20:28,17])
qqline(qiime2_log[20:28,17])

compare_means(Ruminiclostridium_5 ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF+EPA"))
p17 <- ggboxplot(qiime2_log, x = "Diet", y = "Ruminiclostridium_5",
                 title = "Ruminiclostridium_5", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Ruminiclostridium_5), max(qiime2_log$Ruminiclostridium_5)+0.4))+ 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Ruminiclostridium_9

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,18])
hist(qiime2_log[1:11,18])
qqnorm(qiime2_log[1:11,18])
qqline(qiime2_log[1:11,18])

shapiro.test(qiime2_log[12:19,18]) 
hist(qiime2_log[12:19,18])
qqnorm(qiime2_log[12:19,18])
qqline(qiime2_log[12:19,18])

shapiro.test(qiime2_log[20:28,18]) 
hist(qiime2_log[20:28,18])
qqnorm(qiime2_log[20:28,18])
qqline(qiime2_log[20:28,18])

compare_means(Ruminiclostridium_9 ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
#my_comparisons <- list(c("C", "HF"), c("C", "HF+EPA"))
p18 <- ggboxplot(qiime2_log, x = "Diet", y = "Ruminiclostridium_9",
                 title = "Ruminiclostridium_9", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Ruminiclostridium_9), max(qiime2_log$Ruminiclostridium_9)+0.8))+ 
  theme(legend.position='none') #palette = "jco"
#  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Dubosiella

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,19])
hist(qiime2_log[1:11,19])
qqnorm(qiime2_log[1:11,19])
qqline(qiime2_log[1:11,19])

shapiro.test(qiime2_log[12:19,19]) 
hist(qiime2_log[12:19,19])
qqnorm(qiime2_log[12:19,19])
qqline(qiime2_log[12:19,19])

shapiro.test(qiime2_log[20:28,19]) 
hist(qiime2_log[20:28,19])
qqnorm(qiime2_log[20:28,19])
qqline(qiime2_log[20:28,19])

compare_means(Dubosiella ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("C", "HF_EPA"))
p19 <- ggboxplot(qiime2_log, x = "Diet", y = "Dubosiella",
                 title = "Dubosiella", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Dubosiella), max(qiime2_log$Dubosiella)+0.8))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Parasutterella

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,20])
hist(qiime2_log[1:11,20])
qqnorm(qiime2_log[1:11,20])
qqline(qiime2_log[1:11,20])

shapiro.test(qiime2_log[12:19,20]) 
hist(qiime2_log[12:19,20])
qqnorm(qiime2_log[12:19,20])
qqline(qiime2_log[12:19,20])

shapiro.test(qiime2_log[20:28,20]) 
hist(qiime2_log[20:28,20])
qqnorm(qiime2_log[20:28,20])
qqline(qiime2_log[20:28,20])

compare_means(Parasutterella ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("HF", "HF_EPA"), c("C", "HF_EPA"))
p20 <- ggboxplot(qiime2_log, x = "Diet", y = "Parasutterella",
                 title = "Parasutterella", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Parasutterella), max(qiime2_log$Parasutterella)+1))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Akkermansia

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,21])
hist(qiime2_log[1:11,21])
qqnorm(qiime2_log[1:11,21])
qqline(qiime2_log[1:11,21])

shapiro.test(qiime2_log[12:19,21]) 
hist(qiime2_log[12:19,21])
qqnorm(qiime2_log[12:19,21])
qqline(qiime2_log[12:19,21])

shapiro.test(qiime2_log[20:28,21]) 
hist(qiime2_log[20:28,21])
qqnorm(qiime2_log[20:28,21])
qqline(qiime2_log[20:28,21])

compare_means(Akkermansia ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF_EPA"))
p21 <- ggboxplot(qiime2_log, x = "Diet", y = "Akkermansia",
                 title = "Akkermansia", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Akkermansia), max(qiime2_log$Akkermansia)+0.2))+ 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#Combo plots
ggarrange(p1,p2,p3,p4,p5,p6, ncol = 2, nrow = 3, legend = "none")
ggarrange(p7,p8,p9,p10,p11,p12, ncol = 2, nrow = 3, legend = "none")
ggarrange(p13,p14,p15,p16,p17,p18, ncol = 2, nrow = 3, legend = "none")
ggarrange(p19,p20,p21, ncol = 2, nrow = 3, legend = "none")

###Making SPECIES Boxplots####

qiime <- read.csv("Normalized_DESeq2_SequenceCountsComboData.csv", header = FALSE, row.names = 1)
qiime <- qiime[c(1,56:57),]
qiime <- data.frame(t(qiime)) #transpose the dataframe
qiime$V1 <- substring(qiime$V1, 1) #remove the first character of the first column strings
#assign row names
qiime2 <- qiime[,-1]
qiime2 <- mutate_all(qiime2, function(x) as.numeric(as.character(x)))
rownames(qiime2) <- qiime[,1]
#merge dataframe by sample IDs to append diet information 
qiime2 <- merge(qiime2, meta, by="row.names", all=TRUE)
#sort by diet group
qiime2 <- qiime2[order(qiime2$Diet),] 
row.names(qiime2) <- qiime2[,1]
qiime2 <- qiime2[,-1]
qiime2_log <- log10(qiime2[,1:2]) #include just numbers, not diet group
qiime2_log[qiime2_log < 0] <- 0
qiime2_log$Diet <- qiime2$Diet

#Clostridium_sensu_stricto_1

#checking for assumption of normality for each diet within a microbe
shapiro.test(qiime2_log[1:11,1]) 
hist(qiime2_log[1:11,1])
qqnorm(qiime2_log[1:11,1])
qqline(qiime2_log[1:11,1])

shapiro.test(qiime2_log[12:19,1])
hist(qiime2_log[12:19,1])
qqnorm(qiime2_log[12:19,1])
qqline(qiime2_log[12:19,1])

shapiro.test(qiime2_log[20:28,1]) 
hist(qiime2_log[20:28,1])
qqnorm(qiime2_log[20:28,1])
qqline(qiime2_log[20:28,1])

compare_means(Clostridium_sensu_stricto_1 ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("C", "HF_EPA"))
p1 <- ggboxplot(qiime2_log, x = "Diet", y = "Clostridium_sensu_stricto_1",
                title = "Clostridium_sensu_stricto_1", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Clostridium_sensu_stricto_1), max(qiime2_log$Clostridium_sensu_stricto_1)+0.7)) + 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#Dubosiella_newyorkensis

#checking for assumption of normality
shapiro.test(qiime2_log[1:11,2])
hist(qiime2_log[1:11,2])
qqnorm(qiime2_log[1:11,2])
qqline(qiime2_log[1:11,2])

shapiro.test(qiime2_log[12:19,2]) 
hist(qiime2_log[12:19,2])
qqnorm(qiime2_log[12:19,2])
qqline(qiime2_log[12:19,2])

shapiro.test(qiime2_log[20:28,2]) 
hist(qiime2_log[20:28,2])
qqnorm(qiime2_log[20:28,2])
qqline(qiime2_log[20:28,2])

compare_means(Dubosiella_newyorkensis ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
my_comparisons <- list(c("C", "HF"), c("C", "HF_EPA"))
p2 <- ggboxplot(qiime2_log, x = "Diet", y = "Dubosiella_newyorkensis",
                title = "Dubosiella_newyorkensis", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
                color = "Diet", add = c("point"), ylim = c(min(qiime2_log$Dubosiella_newyorkensis), max(qiime2_log$Dubosiella_newyorkensis)+0.8)) + 
  theme(legend.position='none')+ #palette = "jco"
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

# #g_Clostridium_s_cocleatum
# 
# #checking for assumption of normality
# shapiro.test(qiime2_log[1:11,3])
# hist(qiime2_log[1:11,3])
# qqnorm(qiime2_log[1:11,3])
# qqline(qiime2_log[1:11,3])
# 
# shapiro.test(qiime2_log[12:19,3]) 
# hist(qiime2_log[12:19,3])
# qqnorm(qiime2_log[12:19,3])
# qqline(qiime2_log[12:19,3])
# 
# shapiro.test(qiime2_log[20:28,3]) 
# hist(qiime2_log[20:28,3])
# qqnorm(qiime2_log[20:28,3])
# qqline(qiime2_log[20:28,3])
# 
# compare_means(g_Clostridium_s_cocleatum ~ Diet,  data = qiime2_log, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
# my_comparisons <- list(c("C", "HF+EPA"))
# p3 <- ggboxplot(qiime2_log, x = "Diet", y = "g_Clostridium_s_cocleatum",
#                 title = "g_Clostridium_s_cocleatum", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
#                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$g_Clostridium_s_cocleatum), max(qiime2_log$g_Clostridium_s_cocleatum)+0.4)) + 
#   theme(legend.position='none')+ #palette = "jco"
#   stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
# 
# #g_Akkermansia_s_muciniphila
# 
# #checking for assumption of normality
# shapiro.test(qiime2_log[1:11,4])
# hist(qiime2_log[1:11,4])
# qqnorm(qiime2_log[1:11,4])
# qqline(qiime2_log[1:11,4])
# 
# shapiro.test(qiime2_log[12:19,4]) 
# hist(qiime2_log[12:19,4])
# qqnorm(qiime2_log[12:19,4])
# qqline(qiime2_log[12:19,4])
# 
# shapiro.test(qiime2_log[20:28,4]) 
# hist(qiime2_log[20:28,4])
# qqnorm(qiime2_log[20:28,4])
# qqline(qiime2_log[20:28,4])
# 
# compare_means(g_Akkermansia_s_muciniphila ~ Diet,  data = qiime2_log, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 
# my_comparisons <- list(c("C", "HF"), c("C", "HF+EPA"))
# p4 <- ggboxplot(qiime2_log, x = "Diet", y = "g_Akkermansia_s_muciniphila",
#                 title = "g_Akkermansia_s_muciniphila", xlab = "Diet", ylab = "Log10 Normalized Sequence Abundance",
#                 color = "Diet", add = c("point"), ylim = c(min(qiime2_log$g_Akkermansia_s_muciniphila), max(qiime2_log$g_Akkermansia_s_muciniphila)+0.6)) + 
#   theme(legend.position='none')+ #palette = "jco"
#   stat_compare_means(comparisons = my_comparisons, method = "t.test")

ggarrange(p1,p2, ncol = 2, nrow = 3, legend = "none")
