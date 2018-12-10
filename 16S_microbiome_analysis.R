setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/16S Microbiome Analysis HF+EPA Study/Week 13-14 10.23.18 Analysis/")
rm(list=ls())
library(ggplot2)
library(reshape2)
library(ggpubr)
library(vegan)
library(phyloseq)
library(cluster)
require("igraph")
library(gplots)
library("dplyr")
library("ggpubr")
library(tidyverse)
library(car)

meta <- read.csv("metadata_W13-14_10.23.18.csv", header = TRUE, row.names = 1) #load in metadata
#meta <- meta[c(1:18, 20:28),] #removing mouse #832 because he was singly housed
qiime <- read.csv("n3_qiime_phylum_Log10_norm_metaMerged.csv", header = TRUE, row.names = 1) #load in normalized OTU file
qiime <- qiime[c(1:5,7:18,20:28),] #removing mouse #832 & #823 because they were singly housed
#n3_qiime_phylum_Log10_norm_metaMerged
#n3_qiime_phylum_Log10_median_impute
qiime <- qiime[order(qiime$Diet),] #sort by diet group
#qiime <- qiime[,c(1:4,8:length(qiime))] #remove bacteria with zero inflated OTUs

#reading in relative abundances
rel <- read.csv("n3_qiime_phylum_RelativeAbundances.csv", header = TRUE) #load in normalized relative OTU abundances file

#plotting absolute abundances (IDs)
qiime_phyla_only <- qiime[,c(1:11)] #only selecting bacteria & sample IDs
qiime_phyla_only <- qiime_phyla_only[,c(1,3:length(qiime_phyla_only))]
df_long <- melt(qiime_phyla_only, id.vars = "X.SampleID", variable.name = "Phyla")
ggplot(df_long, aes(x = X.SampleID, y = value, fill = Phyla)) + geom_bar(stat = "identity")

#plotting aboslute abundances (diets)
qiime_phyla_only_2 <- qiime[,c(2:11,15)]
df_long2 <- melt(qiime_phyla_only_2, id.vars = "Diet_Num", variable.name = "Phyla")
ggplot(df_long2, aes(x = Diet_Num, y = value, fill = Phyla)) + geom_bar(stat = "identity")

#plotting absolute abundances (diets averaged)
qiime_phyla_only_3 <- qiime[,c(2:11,16)]
df_long3 <- melt(qiime_phyla_only_3, id.vars = "Diet", variable.name = "Phyla")
ggplot(df_long3, aes(x = Diet, y = value, fill = Phyla)) + geom_bar(stat = "identity")

#plotting RELATIVE abundances (diet_IDs) #THIS IS THE CORRECT WAY TO PLOT RELATIVE ABUNDANCES
qiime_phyla_only_4 <- rel[,c(2:11,15)] #2:11 are the bacteria, 15 is the column with the diets
df_long4 <- melt(qiime_phyla_only_4, id.vars = "Diet_Num", variable.name = "Phyla")
ggplot(df_long4, aes(x = Diet_Num, y = value, fill = Phyla)) + geom_bar(stat = "identity")

qiime_phyla_only_5 <- rel[,c(2:11,17)] #2:11 are the bacteria, 17 is the column with the diet groups
agg_phyla <- aggregate(qiime_phyla_only_5[,1:10], list(qiime_phyla_only_5$Diet.group), mean) #1-10 are the bacteria, we are aggregating by Diet.group
#plotting RELATIVE abundances grouping by diets
#Use this if you want to take the means of each phyla per group of inidividuals/samples
agg_phyla_1 <- melt(agg_phyla, id.vars = "Group.1", variable.name = "Phyla")
#to make a high resolution image at 600 dpi
png("RelativeAbundances-Phyla_W13-14.png", units="in", width=9, height=7, res=600)
ggplot(agg_phyla_1, aes(x = Group.1, y = value, fill = Phyla)) + geom_bar(stat = "identity")+
  ggtitle("Phylum Level Relative Abundances") +
  xlab("Diet Group") + ylab("Relative Abundance")
dev.off()

#qiime_phyla_only_k <- setDT(qiime_phyla_only_k, keep.rownames = TRUE)[] #making rows into first column as well

#plotting OTU abundances heatmap 1 (IDs)
transpose_phyla1 <- t(qiime_phyla_only)
colnames(transpose_phyla1) <- transpose_phyla1[1, ] # the first row will be the header
transpose_phyla1 = transpose_phyla1[-1, ] 
#transpose_phyla1 = transpose_phyla1[,c(1:3,6,8:22)] #remove bad controls (samples 829,830,832)
class(transpose_phyla1) <- "numeric"
#creating an OTU table object
otu <- otu_table(transpose_phyla1, taxa_are_rows = TRUE)
plot_heatmap(otu, method="PCoA", distance="bray")

#Plotting heatmap 2 (diets)
qiime_phyla_only_2 <- qiime[,c(15,2:11)]
transpose_phyla2 <- t(qiime_phyla_only_2)
colnames(transpose_phyla2) = transpose_phyla2[1, ] # the first row will be the header
transpose_phyla2 = transpose_phyla2[-1, ] 
class(transpose_phyla2) <- "numeric"
#creating an OTU table object
otu2 <- otu_table(transpose_phyla2, taxa_are_rows = TRUE)
plot_heatmap(otu2, method="PCoA", distance="bray", sample.order = c("C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "HF_1", "HF_2", "HF_3", "HF_4", "HF_5", "HF_6", "HF_7", "HF+EPA_1", "HF+EPA_2", 'HF+EPA_3', 'HF+EPA_4', 'HF+EPA_5', 'HF+EPA_6', 'HF+EPA_7', 'HF+EPA_8'))

#creating a sampleData object for the meta data
sampleData <- sample_data(meta)
#making a phyloseq object with the metadata annotations combined with the otus
phylo <- phyloseq(otu, sampleData) 

#PCoA & MDS Plots are the same!!
#PCoA Plot with Bray-Curtis distance
ord = ordinate(phylo, method="PCoA", distance="bray")
png("PCoA_Phyla_W13-14_noLabel.png", units="in", width=9, height=7, res=600)
plot_ordination(phylo, ord, color = "Diet", shape = "Diet", label = "X.SampleID") + #label = "BarcodeSequence"
  ggplot2::ggtitle("Bray-Curtis Principal Coordinates Analysis")#+ stat_ellipse(type = "t", linetype = 2, level = 0.90)
dev.off()

ord2 = ordinate(phylo, method = "NMDS", distance="bray")
plot_ordination(phylo, ord2, color = "Diet", shape = "Diet", label = "X.SampleID") + #label = "BarcodeSequence"
  ggplot2::ggtitle("NMDS Plot")
#An NMDS plot is the same as a PCoA of relative abundances

#make a histogram of the OTU absolute abundances
#only selecting columns 2:end because my first column contains sample IDs
qiime_phyla_only_nums <- qiime_phyla_only[,c(2:length(qiime_phyla_only))]
hist(as.numeric(unlist(qiime_phyla_only_nums)), main = "Histogram of all OTU abundance values \n non-imputed", xlab = "OTU Abundances")

#clustering with vegan
#only selecting columns 2:end because my first column contains sample IDs
dis <- vegdist(qiime_phyla_only_nums)
clua <- hclust(dis, "average") #average linkage clustering
plot(clua)
rect.hclust(clua, 2)
ord <- metaMDS(qiime_phyla_only) #NMDS (based on rankings of values and not actual numbers) - like a non-parametric test
plot(ord, display = "sites")
ordicluster(ord, clua, col="blue")

#getting beta diversity using vegan package
betad <- betadiver(qiime_phyla_only_nums, "z") #getting beta diversity (see vegan documentation for more info)
adonis(betad ~ Diet.group, qiime, perm=200) #multivariate ANOVA
#betadisper can be used to analyze beta diversities with respect to diet groups
mod <- with(qiime, betadisper(betad, Diet.group))
plot(mod)
# Box plot below of beta diversity measured as the average steepness (z) of the species area curve
# in the Arrhenius model S = cXz
png("Phylum_beta_diversity.png", units="in", width=9, height=7, res=600)
boxplot(mod, main="Phylum Level Beta Diversity")
dev.off()
#run statistical tests to see if there is a significant difference in beta diversity. 
anova(mod)
permutest(mod)
TukeyHSD(mod)

#K-means clustering
qiime_phyla_only_k <- qiime_phyla_only[,2:length(qiime_phyla_only)] #removing second column to include only phyla info
#qiime_phyla_only_k <- qiime_phyla_only_k[!rownames(qiime_phyla_only_k) %in% influential, ] #removing outliers
#using PCoA ordination vectors as input for kmeans
fviz_nbclust(ord$vectors, kmeans, method = "wss") + geom_vline(xintercept = 3, linetype = 2) + labs(subtitle = "Elbow Method - Kmeans Clustering_Imputed") #k = 4
fviz_nbclust(ord$vectors, kmeans, method = "silhouette") + labs(subtitle = "Silhouette Method Imputed") #k=2
fviz_nbclust(ord$vectors, kmeans, method = "gap_stat", nboot = 500) + labs(subtitle = "Gap Statistic Method - K-means Clustering, nboot = 500") #k=4 with nboot=50, k=5 with nboot=500

km <- kmeans(ord$vectors, 3, nstart = 50)
fviz_cluster(km, data = qiime_phyla_only_k, ellipse.type = "euclid", star.plot = TRUE, repel = TRUE, ggtheme = theme_minimal(), main = "K-means Clustering", subtitle = "16S Microbiome Dataset - Imputed, PCoA used as input")

#CLARA clustering
#optimal clusters = 3
fviz_nbclust(qiime_phyla_only_k, clara, method = "wss") + geom_vline(xintercept = 3, linetype = 2) #+ theme_classic()
clara.res <- clara(qiime_phyla_only_k,3,samples = 100, pamLike = TRUE)
autoplot(clara.res, frame = TRUE, label = TRUE, label.repel = TRUE)

#PAM clustering #did not work well
fviz_nbclust(qiime_phyla_only_k, pam, method = "silhouette") + labs(subtitle = "Silhouette Method - Imputed") #k=3
pam.res <- pam(qiime_phyla_only_k, 2)
fviz_cluster(pam.res, ellipse.type = "t", repel = TRUE) #ellipse.type = norm, euclide

#DBSCAN Clustering #did not work well
set.seed(123)
dbscan::kNNdistplot(qiime_phyla_only_k, k=5)
abline(h=0.47, lty=2)
db <- fpc::dbscan(qiime_phyla_only_k, eps=0.47, MinPts = 5)
fviz_cluster(db, data = qiime_phyla_only_k, stand = FALSE, ellipse = FALSE, show.clust.cent = FALSE, geom = "point", palette = "jco", repel = TRUE)

#hierarchal agglomerative clustering
res.dist <- dist(qiime_phyla_only_nums, method = "euclidean")
res.hc <- hclust(d = res.dist, method = "ward.D")
grp <- cutree(res.hc, k=3)
png("Phylum_dendogram_W13-14.png", units="in", width=9, height=7, res=600)
fviz_dend(res.hc, cex=0.5, k=3, color_labels_by_k = TRUE, rect = TRUE, main = "Phylum Level Agglomerative Hierarchical Clustering")
dev.off()
fviz_dend(res.hc, cex=0.7, k=3, type = "circular")
#assess cluster accuracy (correlated clustered groups with actual data groups)
res.coph <- cophenetic(res.hc)
cor(res.dist, res.coph) #0.62 correlation - ugly
#making cluster plot from hierarchal clustering
#fviz_cluster(list(data=qiime_phyla_only_k, cluster=grp), ellipse.type = "convex", repel = TRUE, show.clust.cent = FALSE)

#phylogenetic dendogram
fviz_dend(res.hc, k=3, type = "phylogenic", repel = TRUE)

#computing multivariate outliers with cook's distance method
qiime_phyla_only <- qiime_phyla_only[,c(2:length(qiime_phyla_only))]
model <- lm(Encoding ~ ., data=qiime_phyla_only)
cooksd <- cooks.distance(model)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 2.5*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>2.5*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels
influential <- names(cooksd)[(cooksd > 2.5*mean(cooksd, na.rm=T))]  # influential row numbers
head(qiime_phyla_only[influential, ])  # influential observations.
car::outlierTest(model)

#####PCoA Plot using relative abundances and Bray-Curtis distance#####
# qiime_phyla_only <- rel[,c(1:11)] #only selecting bacteria & samples
# transpose_phyla1 <- t(qiime_phyla_only)
# colnames(transpose_phyla1) <- transpose_phyla1[1, ] # the first row will be the header
# transpose_phyla1 = transpose_phyla1[-1, ] 
# #transpose_phyla1 = transpose_phyla1[,c(1:3,6,8:22)] #remove bad controls (samples 829,830,832)
# class(transpose_phyla1) <- "numeric"
# #creating an OTU table object
# otu <- otu_table(transpose_phyla1, taxa_are_rows = TRUE)
# sampleData <- sample_data(meta)
# phylo <- phyloseq(otu, sampleData) 
# ord = ordinate(phylo, method="PCoA", distance="bray")
# plot_ordination(phylo, ord, color = "Diet", shape = "Diet", label = "Samples") + #label = "BarcodeSequence"
#   ggplot2::ggtitle("Bray-Curtis Principal Coordinates Analysis")
######

#Plotting alpha diversity - Box plot with dot plot
p <- ggplot(qiime, aes(x=Diet.group, y=shannon_normalized_alpha, color=Diet.group)) + geom_boxplot()
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + labs(title="Phylum Level Shannon Alpha Diversity")

#stats for shannon alpha diversity
compare_means(shannon_normalized_alpha ~ Diet.group, data = qiime, method = "t.test") #method = anova

#better looking boxplot
png("Phylum_alpha_diversity.png", units="in", width=9, height=7, res=600)
ggplot(qiime, aes(x=Diet.group, y=shannon_normalized_alpha, color=Diet.group, group=factor(Diet.group))) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(height = .2, width = .2))+ labs(title="Phylum Level Shannon Alpha Diversity")
dev.off()

