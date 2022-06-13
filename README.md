# nanobiome
Tools to convert EPI2ME output for 16S sequencing and analyze it in R.

epi2meR is an R script that reads .csv files created with the EPI2ME 16S pipeline from ONT, and creates a clean merged table of OTU counts.
Such a table is sparse and similar to counts tables from RNAseq experiments, so it could be treated in a similar way. We use edgeR to obtain differentially represented OTUs (DROs) in two groups of samples.
