# Chandan Pavuluri

# RNA-Seq_DEA_Pathway
Differntial Expression Analysis of the RNA-Seq data

I have performed Differential expression anlysis on publicly available GEO dataset(GSE150392). 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150392


# Differntial Expression Analysis:
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

I have done the analysis using DESeq2 and edgeR. 

**Comapring DESeq2 and edgeR**
Both the methods are almost similar but they differ at Normalization, Dispersion estimation, Statistical testing.
Normalization:- DESeq2 uses the Median of ratios Normalization, where as in edgeR I have used the Trimmed Mean of Mvalues(TMM) normalization.

Dispersion estimation:-  DESeq2 uses maximum likelihood estimation, in edgeR I have used the Quasi likelihood estimation.

Statistical testing:- DESeq2 uses Wald test and edgeR uses the F test.

Both the methods uses the Negative Binomial distribution and FDR calculation by Benjamini Hochberg.

I tried to identify the percentage of similarilty in significantly regulated genes of both the methods.
In this case while filtering the genes only both the methods have a lot of difference in thier total count after filtering.
Total genes in dataset 36941

DESeq2 filtered gene count 28941
edger filtered gene count 19809

Genes passing the pvalue<0.05,abs(lfc)>1,FDR<0.01
DESeq2 significant count - 3575
edgeR significant count - 1302
Common on both count is 1302

All the genes from the edgeR are matched with DESeq2, Here there is a 100% match, but generally the datasets that I worked earlier when the filtered gene count is almost similar then they have 85% similarity but it is obvious that they differ with pvalues and lfc in both methods as the testing is different. It is not good idea to compare both the methods, both has their own pros and cons.


# Pathway Mapping

Need to convert the EnsemblID to entrez ID so that we can mapthese genes to pathways.

Clusterprofiler reads the Kegg annotation online. KEGG Enrichment Analysis of a gene set. Input is a vector of genes, the enrichkegg will return the enrichment KEGG categories. 
Both upregulated and downregulated can be known by doing it separately. 

Visualization of the path ways can be done by using the pathview but the logfoldchanges are required and genes are coloured based on their expression levels in the data. The output is generated as a png file.

# GO enrichment 

By using GO.db and GOstats need to give parmeters for both upraguated and downregulated entrezIDs and need to select the ontology for the GO like BP(Biological Process), CC(Cellular Components), MF(Molecular Functions). we can change it easily with using all parmeters just by using ontology(params_up) = "MF" so that params_up is changed to MF from BP.
Then these are subjected to the Hypergeometric tests. then by looking at the summary of the results we can understand the enrichment of the results.

