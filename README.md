# nanobiome
Tools to convert EPI2ME output for 16S sequencing and analyze it in R.

epi2meR is an R script that reads .csv files created with the EPI2ME 16S pipeline from ONT, and creates a clean merged table of OTU counts for quick preliminary analysis, creating an interactive MDS plot.
Assuming that the OTU counts table (very sparse) could be treated like to counts tables from RNAseq experiments (a heavily contested assumption, so beware!), the script applies edgeR to obtain differentially represented OTUs (DROs) in two groups of samples.
It has been shown that the Wilcoxon ranks test controls better for the FDR than tests based on the negative binomial (like edgeR or DSeq2) without losing power, when n >= 8, so the script also implements this test. Users can compare the results from both strategies, check heatmaps and interactive volcanoplots created with Glimma.

The script can be run line-by-line. THE ONLY CHANGE USERS SHOULD MAKE is in line 72, where conditions (groups to be compared) should be specified.

For more in-depth analyses, users are advised to download the wimp .csv files from EPI2ME and run the epi2me2r script from https://github.com/cran/epi2me2r which will create an object of the Phyloseq class that can be directly used in R using https://github.com/joey711/phyloseq
