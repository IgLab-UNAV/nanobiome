# nanobiome
Tools to convert EPI2ME output for 16S sequencing and analyze it in R.

epi2meR is an R script that reads .csv files created with the EPI2ME 16S pipeline from ONT, and creates a clean merged table of OTU counts.
Such a table is sparse and similar to counts tables from RNAseq experiments, so we use edgeR to obtain differentially represented OTUs (DROs) in two groups of samples.
Also, since it has been shown that the Wilcoxon ranks test controls better for the FDR without losing power, when n >= 8, we implement this test. Users can compare the results from both strategies, check heatmaps and interactive volcanoplots created with Glimma.

The script can be run line-by-line. THE ONLY CHANGE USERS SHOULD MAKE is in line 72, where conditions (groups to be compared) should be specified.
