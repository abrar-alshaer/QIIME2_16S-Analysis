rm(list=ls())
setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/16S Microbiome Analysis HF+EPA Study/Week 13-14 10.23.18 Analysis/")
library("nlme")
library(ggplot2)
library(ggpubr)
library(gplots)

file <- read.csv("n3_qiime_phylum_Log10_norm_metaMerged_cages.csv", header = TRUE)
file <- file[order(file$Diet),] #sort by diet group
set <- c("2_C", "5_C", "7_C", "1_HF", "4_HF", "11_HF", "3_HF+EPA", "6_HF+EPA", "8_HF+EPA", "9_HF+EPA", "Single_HF", "Single_C")

par(mfrow=c(3,2))
plot(file[,4] ~ file$Cages)
plot(file[,5] ~ file$Cages)
plot(file[,6] ~ file$Cages)
plot(file[,7] ~ file$Cages)
plot(file[,8] ~ file$Cages)
plot(file[,9] ~ file$Cages)
par(mfrow=c(3,2))
plot(file[,10] ~ file$Cages)
plot(file[,11] ~ file$Cages)
plot(file[,12] ~ file$Cages)
plot(file[,13] ~ file$Cages)

pdf("OTU_cage_plots.pdf")
for (col in colnames(file)[4:13]) {
print(ggboxplot(file, x = "Cages", y = col,
          title = paste(col, " OTU Plot by Cages"), xlab = "Cage", ylab = "Log10 Normalized Absolute Abundance",
          color = "Diet.group", add = c("point"), order = set)+ font("x.text", size = 7)+
  theme(legend.position='top')) #palette = "jco"
}
dev.off()

#running cage effects linear model
pvalue = vector()
rho = vector()
for(i in 4:13) #column indeces 
{
  print(colnames(file)[i])
  microbe <- file[,i] #inidividual microbe column
  cage <- file$Cages
  diet <- file$Diet
  myFrame <- data.frame(microbe, cage, diet)
  #check for assumption of normality for each diet group, skipping zero inflated microbes
  if(colnames(file)[i] != "Cyanobacteria" && colnames(file)[i] != "Deferribacteres" && colnames(file)[i] != "Euryarchaeota" && colnames(file)[i] != "Patescibacteria")
  {
    print(shapiro.test(myFrame[1:11,1])) #controls are 1:10
    print(shapiro.test(myFrame[12:19,1])) #HF are 12:19
    print(shapiro.test(myFrame[20:28,1])) #HF+EPA are 20:28
  }
  #Create linear model (generalized least squares with mixed effects & REML - restricted maximum likelihood function)
  #corCompSym is to calculate the correlation (Rho) between the cage and diet group under a specific microbe...
  #...and you are only testing 1 correlation (cage effect).
  M.gls <- gls(microbe ~ diet , method = "REML", correlation = corCompSymm( form = ~ 1 | cage),data=myFrame)
  print(summary(M.gls))
  rho[i] <- coef(M.gls$modelStruct[1]$corStruct, unconstrained = FALSE)[[1]]
  pvalue[i] <- anova(M.gls)$'p-value'[2]
  #we are running an ANOVA to test between the various diet groups and see if there are differences...
  #...within a diet group between different cages. -> An ANOVA takes random effects into account! 
}

adjusted_pvals <- p.adjust(pvalue[c(4:13)], method="BH")
rho <- rho[c(4:13)] #make sure 4:13 match your column numbers that you are using
