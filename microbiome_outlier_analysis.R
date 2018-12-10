setwd("C:/Users/Abrar/OneDrive - University of North Carolina at Chapel Hill/Shaikh Lab/16S Microbiome Analysis HF+EPA Study/Week 13-14 10.23.18 Analysis/")
rm(list=ls())

meta <- read.csv("metadata_W13-14_10.23.18.csv", header = TRUE, row.names = 1) #load in metadata
qiime <- read.csv("n3_qiime_phylum_Log10_norm_metaMerged.csv", header = TRUE) #load in normalized OTU file
qiime <- qiime[order(qiime$Diet),] 

groups = levels(qiime$Diet) #create a variable called groups for all levels (row groups) of the dataframe
columns = colnames(qiime)[2:11] #store columns in variable

datalist = list()
c = 1

#loop through each group (rows) for microbiome data
for(i in 1:length(groups)) {
  print(groups[i]) #print group
  sub <- subset(qiime, qiime$Diet %in% groups[i]) #subset dataframe by the biological group
  #loop through each column within each group
  for(j in 1:length(columns)){
    #print(paste("COLUMN (Time Point):",columns[j]))
    outliers <- boxplot.stats(sub[,columns[j]])$out #output of outlier values to a list, calculate outlier values based on the IQR formula (Q1-1.5*IQR, Q3+1.5*IQR)
    #if there are outliers (list isn't empty)
    if(length(outliers)>0)
      {
        datalist[[c]] <- c(groups[i],columns[j],boxplot.stats(sub[,columns[j]])$out)
        c = c+1
      }
  }
}

list_data = data.frame(datalist = unlist(datalist))
matrix_data = as.matrix(datalist)
#big_data = do.call(rbind, datalist) #combines all previous dataframes from for loop
#Check to see if the list is empty (no outliers) - if not, add column headers & send dataframe to a file
if(length(matrix_data) < 1){
  print("There are no outliers in your dataset")
} else {
  #colnames(big_data) <- c("Sample", "Microbe", "Abundance")
  write.csv(matrix_data, "microbiomeW13-14_outliers.csv")
}

#calculate outlier values based on the IQR formula MANUALLY
#quantile(sub$X1.413664)[2]-1.5*IQR(sub$X1.413664) Q1-1.5*IQR
#quantile(sub$X1.413664)[4]+1.5*IQR(sub$X1.413664) Q3+1.5*IQR