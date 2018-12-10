rm(list=ls())
setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/16S Microbiome Analysis HF+EPA Study/DEMUX_2018/DEMUX_2018/ALL_RUNS/Week13_14/")
library(dada2); packageVersion("dada2")
library(phyloseq) #packageVersion("phyloseq")
library(ggplot2) #packageVersion("ggplot2")

path <- "C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/16S Microbiome Analysis HF+EPA Study/DEMUX_2018/DEMUX_2018/ALL_RUNS/Week13_14/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:10])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out) 

#Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#Inspecting the returned dada-class object:
dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#construct sequencing table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
#DADA2 says that is more than 25% of sequences (100-sum above) are chimeras that means you may have to do
#primer trimming 

#track # of reads throughout pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Assigning Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa_s <- addSpecies(taxa, "silva_species_assignment_v132.fa.gz")

#Plotting output in phyloseq
meta <- read.csv("metadata_W13-14_10.23.18.csv", header = TRUE, row.names = 9) #load in metadata
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa_s))
ps

#Visualize alpha-diversity:
plot_richness(ps, x="Diet", measures=c("Shannon", "Simpson"), color="Diet")

#Ordination
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray") #or PCoA
plot_ordination(ps.prop, ord.nmds.bray, color="Diet", title="Bray-Curtis NMDS", shape = "Diet")

ord.pcoa <- ordinate(ps.prop, method="PCoA", distance="bray")
plot_ordination(ps.prop, ord.pcoa, color="Diet", title="Bray-Curtis Principal Coordinates Analysis", shape = "Diet")
#plot_ordination(ps, ord.pcoa, color="Diet", title="Bray PCoA") #yielded same result as above

#abundances
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:500] #pick out top 20 sequences
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Diet", fill="Phylum") + facet_wrap(~Diet, scales="free_x")
