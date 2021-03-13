Differential Expression Analysis using DESeq2
================

``` r
#Create a folder for the input files and results of the analysis
#Enter path of the folder that was created by you.
#at the end of path add "/".

path<-"C:/Users/ckpav/Documents/RNA-Seq_DEA_Pathway/"


#Copy the input file in to the path
#Enter input file name of the raw counts file
Input<-"GSE150392.csv"


#Total number of samples
n=6
```

\#Loading reuired
packages

``` r
if(!require('BiocManager',quietly = T)) install.packages("BiocManager");library(BiocManager)
if(!require('DESeq2',quietly = T))
  BiocManager::install("DESeq2");library(DESeq2)
if(!require('ggplot2',quietly = T)) install.packages("ggplot2");library(ggplot2)
if(!require('dplyr',quietly = T)) install.packages("dplyr");library(dplyr)
if(!require('stringr',quietly = T)) install.packages("stringr");library(stringr)
if(!require('pheatmap',quietly = T)) install.packages("pheatmap");library(pheatmap)
if(!require('EnhancedVolcano',quietly = T))
  BiocManager::install("EnhancedVolcano");library(EnhancedVolcano)
```

\#Reading the raw counts(Counts\_matrix) and creating metadata file

``` r
# Reading the raw counts
counts_matrix <- read.table(paste0(path,Input),header=T,sep = ",")
dimnames(counts_matrix)[[1]] <- counts_matrix[,1]
counts_matrix <- counts_matrix[,-1]
dim(counts_matrix)
```

    ## [1] 36941     6

``` r
# Creating a metadata file 
meta_data<-data.frame(colnames(counts_matrix))
colnames(meta_data)[1] <- "Samples"

meta_data <- meta_data %>%
  mutate(Condition = if_else(str_detect(meta_data$Samples,"\\Cov"),"COVID","Control"))

Condition <- factor(meta_data[,2])
```

\#Creating DESeqDataSet object from counts\_matrix

``` r
#Creating a DESeq Dataset object
dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = meta_data, design = ~Condition)

#Filtering the Dataset object based on chosen Rowsum filter
dds <- dds[rowSums(counts(dds)) >= 10,]
dim(dds)
```

    ## [1] 28941     6

``` r
#Differential expression analysis
dds <- DESeq(dds)
```

# Calculating Library sizes

``` r
par(mfrow=c(1,2))

barplot(colSums(counts_matrix), main = "Before Normalization", las =2, cex.names = 0.75)

dds_norm <- counts(dds, normalized=T)
barplot(colSums(dds_norm), main = "After Normalization", las=2, cex.names=0.75)
```

![](DEA_DESeq2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# PCA plot

``` r
#Number of samples
n=6
if (n<50){
rld <- rlog(dds, blind = T)
#png(paste0(path,"PCA_DESeq2",Input,".png"), 700, 500, pointsize=20)
pca <- plotPCA(rld, intgroup="Condition")
print(pca)
#dev.off()

}else{
vsd <- varianceStabilizingTransformation(dds, blind=T)
#png(paste0(path,"PCA_DESeq2",Input,".png"), 700, 500, pointsize=20)
pca <- plotPCA(vsd, intgroup="Condition")
print(pca)
#dev.off()
}
```

![](DEA_DESeq2_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

# Boxplot for Identifying Outliers

``` r
#png(paste0(path,"Boxplot_DESeq2",Input,".png"), 700, 500, pointsize=20)
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
```

![](DEA_DESeq2_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
#dev.off()
```

# Heatmap for understanding the clustering

``` r
if (n<30){
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$Condition
colnames(sampleDistMatrix) <- rld$Condition
#png(paste0(path,"Heatmap_DESeq2",Input,".png"), 700, 500, pointsize=20)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)

}else{
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Labels
colnames(sampleDistMatrix) <- vsd$Condition
# png(paste0(path,"Heatmap_DESeq2",Input,".png"), 700, 500, pointsize=20)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
}
```

![](DEA_DESeq2_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

# Comparing between the groups

``` r
# Enter comparison
# If Control vs COVID --c("Condition","Control","COVID")
# If COVID vs Control --c("Condition","COVID","Control")
compare<-c("Condition","COVID","Control")

# Extract results from a DESeq analysis
res <- results(dds,compare)

head(res)
```

    ## log2 fold change (MLE): Condition COVID vs Control 
    ## Wald test p-value: Condition COVID vs Control 
    ## DataFrame with 6 rows and 6 columns
    ##                              baseMean log2FoldChange     lfcSE      stat
    ##                             <numeric>      <numeric> <numeric> <numeric>
    ## ENSG00000000003.15_TSPAN6   2364.2595     -0.9325238  0.283927 -3.284377
    ## ENSG00000000005.6_TNMD      1421.7657     -4.1624703  0.467955 -8.895026
    ## ENSG00000000419.12_DPM1     1675.1150     -0.4458774  0.445809 -1.000153
    ## ENSG00000000457.14_SCYL3     860.9035      1.5166261  0.387429  3.914594
    ## ENSG00000000460.17_C1orf112  242.8127     -0.0964369  0.442937 -0.217721
    ## ENSG00000000938.13_FGR         6.3513      1.6077389  1.264198  1.271746
    ##                                  pvalue        padj
    ##                               <numeric>   <numeric>
    ## ENSG00000000003.15_TSPAN6   1.02208e-03 7.05295e-03
    ## ENSG00000000005.6_TNMD      5.84057e-19 1.08100e-16
    ## ENSG00000000419.12_DPM1     3.17236e-01 5.30771e-01
    ## ENSG00000000457.14_SCYL3    9.05565e-05 9.47885e-04
    ## ENSG00000000460.17_C1orf112 8.27646e-01 9.13671e-01
    ## ENSG00000000938.13_FGR      2.03463e-01 3.93838e-01

# Visuvalization of results through Volcano plot

``` r
#png(paste0(path,"VolcanoPlot_DESeq2",Input,".png"), 700, 500, pointsize=20)
volcano_plot<-EnhancedVolcano(res,lab = NA,x = 'log2FoldChange',y = 'pvalue',pCutoff = (10e-2)/2,FCcutoff = 2.0,xlim = c(-10, 10),ylim = c(0, -log10(0.01e-12)))
print(volcano_plot)
```

![](DEA_DESeq2_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
#dev.off()
```

# Visuvalization of results through MA plot

``` r
#png(paste0(path,"MAPlot_DESeq2",Input,".png"), 700, 500, pointsize=20)
DESeq2::plotMA(res,alpha=0.01, ylim=c(-6,6),colSig = "red3")
```

![](DEA_DESeq2_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
#dev.off()
```

# Writing the Significantly differntially expressed results in to a CSV file

``` r
# Passing the results in to Dataframe
A<-data.frame(res)

# Passing the rownames to a new column
A$Ensembl_ID<-row.names(A)

# Removing the remaining charcters after "."
A$Ensembl_ID <- gsub("\\..*","", row.names(A))

# REmoving the duplicates
A<-unique(A)

# Changing the rownames to clean EnsemblIDs
row.names(A)<- A$Ensembl_ID

# Statistical filtering of the results
B<-A%>%
  filter(abs(log2FoldChange)>1)%>%
  arrange(desc(log2FoldChange))%>%
  filter(pvalue<0.05)%>%
  filter(padj<0.01)

# Creating the regulation column
B$Regulation <- ifelse(B$log2FoldChange>0, "Up", "Down")

# Passing the required columns to Dataframe
B1<-B%>%
  select(pvalue,log2FoldChange,padj,Regulation)

# writing the results to a csv
write.csv(B1,paste0(path,"DiffExp_DESeq2_",Input,"samples.csv"))
```