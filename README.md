# Signature Gene optimization
We propose a method for identifying a set of de novo representative genes, termed signature genes (SGs), which can be used to measure the relative abundance and as phylogenetic markers of each metagenomic species with high precision. An initial set of the 100 genes that correlate with the median gene abundance profile of the metagenomic species (MGS) is selected. However, even in samples with high sequencing depth and species abundances, some genes in the initial set may be undetected, leading to inconsistencies in the estimation of metagenomic species abundance. A variant of the coupon collector’s problem was utilized to evaluate the probability of identifying a certain number of genes in a sample, given their presence, and score the performance of a gene set. This allows us to reject the abundance measurements that are significantly deviating from the expected number of detected genes from the set. Within each sample the expected read counts per gene can be approximated by the discrete negative binomial (NB) distribution, as the reads are assumed to map in proportion to the gene length and show biological variability. A rank-based negative binomial model is used to assess the performance of different gene sets across a large set of samples, facilitating identification of an optimal signature gene set for the MGS


In order to perform the SG identification, a read count matrix and a record of the species with its respective genes and information of gene lengths is needed. 
The data used for the analysis provided by Borderes M, et al* can be found at https://zenodo.org/record/4306051#.Yg5xgy8w2u5

* Functions.R: Contains the functions used for SG optimization
* SG_refinement_SCG.R: Main for optimizing the SGs
* abundance_profiles.R: Takes the optimized SGs from SG_refinement_SCG.R and creates relative abundance estimates.
* mspminer_format_conversion.R: Takes the data fra Borderes M, et al* and converts it to the format used for SG_refinement_SCG.R
* tax_df.RDS: Created using the information in the supplementaty file, Table S3 provided by Borderes M, et al*

*Borderes M, et al. A comprehensive evaluation of binning methods to recover human gut microbial species from a non-redundant reference gene catalog. NAR Genom Bioinform. 2021 Mar 1;3(1):lqab009.
