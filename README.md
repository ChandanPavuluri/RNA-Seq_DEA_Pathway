# RNA-Seq_DEA_Pathway
Differntial Expression Analysis of the RNA-Seq data

I have performed Differential expression anlysis on publicly available GEO dataset(GSE150392). 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150392

I have done the analysis using DESeq2 and edgeR.

# DESeq2_Pathway:
Differntial expression anlysis is ro identify the genes whose expression differs under different conditions, this cannot be done directly by using the raw counts.
The raw counts needs to processed.
1. Normalization
2. Clustering Analysis
3. Modelling
4. Statistical testing

**Normalization**: Raw data needs to be adjusted to account factors that prevent direct comparison of expression measures. Library sizes totally differ for each and every sample even if equal there can some are highly expressed and lowly expressed genes may be present. To over come these factors like Library size, RNA composition and Sequencing depth. Raw data is adjusted according to these factors.

**Clustering Analysis**: RNA-Seq data needs to be checked for data quality by assesing the variance between replicates through visualization, so that we can remove the outliers.

**Modelling**: Read counts cannot be distributed in poisson distribution as the mean and variance is not equal between samples  because the samples are biological replicates, they may share the same condition, but the RNA originate from different samples. So that they Have larger variance compared to mean. SO Negative Binomial Distribution is suitable for the counts.

**Statistical testing**: After modelling, Testing for differntial expression to determine which genes shoe evidence for differnce in expression levels between groups. So we can get the upregulated and downregulated genes.

