rm(list=ls())
setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/16S Microbiome Analysis HF+EPA Study/DEMUX_2018/DEMUX_2018/ALL_RUNS/Week13_14/Paired End Analysis/")
library(phyloseq) 
library(ggplot2)
library(gridExtra)
library("ggpubr")
library("reshape2")
library("dplyr")
library("DESeq2")
library(plyr)

meta <- read.csv("metadata_W13-14_10.23.18.csv", header = TRUE, row.names = 9) #load in metadata
meta2 <- meta["Diet"]

#table.from_biom_COMBO_allTaxa_filter
qiime <- read.csv("sequence abundances all taxa_COMBO.csv", header = FALSE, row.names = 1)
qiime <- data.frame(t(qiime)) #transpose the dataframe
qiime$Microbe <- substring(qiime$Microbe, 1) #remove the first character of the first column strings
#assign row names
qiime2 <- qiime[,-1] #remove the sample column to apply the numeric funciton to all the coutns first
qiime2 <- mutate_all(qiime2, function(x) as.numeric(as.character(x)))
rownames(qiime2) <- qiime[,1] #add back the samples as row names
qiime2 <- merge(qiime2, meta2, by="row.names", all=TRUE)
#sort by diet group
qiime2 <- qiime2[order(qiime2$Diet),] 
row.names(qiime2) <- qiime2[,1]
qiime2 <- qiime2[,-1]

#creating OTU object for phyloseq
transpose_phyla1 <- t(qiime2[,c(1:length(qiime2)-1)])
class(transpose_phyla1) <- "numeric"
#creating an OTU table object
otu <- otu_table(transpose_phyla1, taxa_are_rows = TRUE)

#creating dummy taxonomy table
taxmat <- as.matrix(read.csv("taxa.csv", header = TRUE, row.names = 1)) #must be a matrix
TAX = tax_table(taxmat)

#creating a sampleData object for the meta data
sampleData <- sample_data(meta)
#making a phyloseq object with the metadata annotations combined with the otus
phylo <- phyloseq(otu, sampleData, TAX) 

#non-normalized plots
#normalized and non-normalized NMDS plots look the same because it's non-parametric and utilizes ranks 
ord = ordinate(phylo, method="NMDS", distance="bray") #PCoA & MDS are the same thing
plot_ordination(phylo, ord, color = "Diet", shape = "Diet") + #label = "BarcodeSequence"
  ggplot2::ggtitle("Bray-Curtis NMDS Plot \n norm") + stat_ellipse(type = "norm", size = 1, linetype = 2, level = 0.8)
  #Principal Coordinates Analysis

####Plotting before/after transformation into relative abundances####
plot_abundance = function(phylo,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(phylo, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Diet",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}

# Transform to relative abundance. Save as new object.
phylo_rel = transform_sample_counts(phylo, function(x){x / sum(x)})

#plot before and after transformation
plotBefore = plot_abundance(phylo,"")
plotAfter = plot_abundance(phylo_rel,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2,  plotBefore, plotAfter)


head(sample_data(phylo)$Diet, 10)
#remove control samples because we want to see differences between HF & HF+EPA
phylo <- subset_samples(phylo, Diet != "C") 


####PhyloSeq to DeSeq2 conversion####
diagdds = phyloseq_to_deseq2(phylo, ~ Diet)
#Wald is a method for testing significance with a negative bionial distribution
diagdds = DESeq(diagdds, test="Wald", fitType="local") #diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds)
#alpha = 0.01
#sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(res, "data.frame"), as(tax_table(phylo)[rownames(res), ], "matrix"))
head(sigtab)
write.csv(sigtab, "HF-EPA_vs_HF_DESeq2_normalization_stats.csv")

#MAKE SURE TO DO THE STEPS BELOW WITHOUT RUNNING THE SUBSET_SAMPLES COMMAND THAT EXCLUDES SAMPLES
#So that you get normalized data for all your samples
counts <- counts(diagdds, normalized=TRUE) #retrieving normalized DEXSeq counts
results_counts <- as.data.frame(counts) #converting the counts to a dataframe
write.csv(results_counts, "Normalized_DESeq2_SequenceCountsComboData.csv")
#counts_rmNA <- na.omit(results_counts) #remove genes with NA values
rld <- rlogTransformation(diagdds, fitType="local") 
plotPCA(rld, intgroup = "Diet")
#Sparsity plot
plotSparsity(diagdds, normalized = TRUE)

#all ordination methods: ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
otu_norm <- otu_table(results_counts, taxa_are_rows = TRUE)
phylo_norm <- phyloseq(otu_norm, sampleData, TAX) 
ord_norm = ordinate(phylo_norm, method="RDA", distance="bray") #PCoA & MDS are the same thing
plot_ordination(phylo_norm, ord_norm, color = "Diet", shape = "Diet") + #label = "BarcodeSequence"
  ggplot2::ggtitle("Bray-Curtis Normalized RDA Plot") 
