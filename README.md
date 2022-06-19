# nanobiome
A pipeline to convert EPI2ME output from Oxford Nanopore Technologies (ONT) 16S sequencing protocols and analyze it in R.

**nanobiome** is an R script that reads .csv files created with the EPI2ME 16S pipeline from ONT.
After sequencing and uploading FASTQ files to EPI2ME for alignment, ONT users can download a 16S .csv file. This script first reads all .csv files in the folder and creates a clean merged table of OTU counts. We filter barcodes with < 500 counts and OTUs with < 2 total counts in all samples in order to decrease sparsity, but we do not rarefy the table (users can change these thresholds). This is followed by the creation of an interactive MDS plot using the Glimma library, to check for unsupervised clustering of samples and potential batch effects.

For Differential Abundance (DA) analysis, we follow the recommendations from https://doi.org/10.1038/s41467-022-28034-z and implement two different methods: limma-voom (with TMMwsp normalization) and ALDEx2 (which has shown very good performance in benchmarking studies).
Users can compare the results from both tools, and check heatmaps for Differentially Abundant OTUs (DAOs).

The script can be run line-by-line, changing the metadata with groups of samples to be compared, both for limma and for ALDEx2. Be careful to ensure that you csv files are comma-delimited (check comments in the script about this issue). 

For further analyses, like diversity, Sankey plots, etc, users can download the wimp .csv files from EPI2ME and run the epi2me2r script from https://github.com/cran/epi2me2r which will create an object of the Phyloseq class to be directly used in R with https://github.com/joey711/phyloseq
