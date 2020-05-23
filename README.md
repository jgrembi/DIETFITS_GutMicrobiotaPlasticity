# DIETFITS_GutMicrobiotaPlasticity
Gut microbiota analysis from the DIETFITS study; found that plasticity (variability) was important for sustained weight loss. Manuscript can be found [here](https://doi.org/10.1038/s41598-020-58000-y)

This github project contains the analysis files for the DIETFITS gut microbiota study by Grembi, Nguyen, Haggerty, Gardner, Holmes, and Parsonnet.

The files include:

- **dada2_pipeline.R** - This file contains the script to go from raw FASTQ files to an amplicon sequence variant (ASV) table (analgous to the common OTU table, but with individual, ungrouped sequences).  This file executes quality filtering, learning error rates, dereplication, sample inference, merging paired-end reads, chimera removal, and taxonomic assignment.  

- **generate_phyloseq_object.rmd** - In this notebook we process the outputs from the DADA2 amplicon read counting 
pipeline. We provide the read coverage statistics, filter samples and sequences
which have zero counts, or do not meet required criteria. 
We apply decontamination procedure to filter out sequences suspected
to be contaminants. We also estimate a phylogenetic tree for
the sequences. All data is then combined into a `phyloseq` object.

- **run_phangorm.rmd** - This file contains code to make a phylogenetic tree from multiple sequence alignment using the _phangorn_ package.

- **read_stats.rmd** - This file contains code to report the sequence read statistics.

- **alpha_div_and_pb_ratio.rmd** - In this notebook we calculate alpha diversity for samples and also do analysis on the _Prevotella_/_Bacteroides_ ratio (recently published as important in long-term weight loss by Hjorth et. al., 2017)

- **Ordination.rmd** - In this notebook, we do a Principal Coordinates Analysis of the dataset and show there is no clustering of the data by weight-loss success, by age, gender, or any of the other measured covariates.

- **plasticity.Rmd** - In this notebook, we include all of the analysis for the plasticity findings.  This includes pre-diet daily plasticity, daily plasticity while one the diet (taken after participants have been on the diet for 10 weeks), and also plasticity between baseline and 10 weeks. Also included is the analysis of dietary adherence and it's relationship with plasticity.

- **differential_abundance_limma.Rmd** - This notebook contains code for the differnetial abundance analysis for the discovery cohort.

- **differential_abundance_validation_cohort.Rmd** - This notebook contains code for the differnetial abundance analysis for the validation cohort.

- **voom_ihs.R** - This script contains modifications to the \code{voom} command from \code{limma} package by using arcsinh transformations.

[DOI](https://zenodo.org/badge/182343382.svg)
