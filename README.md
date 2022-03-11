# Discovering Dominant Tumor Immune Archetypes in A Pan-Cancer Census

This repository contains code used to create figures for this paper
All RNAseq and single cell data is available using the GEO accession GSE184398. Some data that is not available on GEO or [UCSF Data Library](https://datalibrary.ucsf.edu/node/121/) has been added to this folder `files_used_for_plots`

# Analysis

### **Figure 1**

1. **Figure 1C** `get_flow_score.py` input files are in `files_used_for_plots`
2. **Figure S1D**

		python get_score.py -h
		
		Usage: get_score.py -o <output file name> -g <gene signature file name> -f <gene expresson tsv>  [options] ALL OPTIONS ARE REQUIRED
		
		Get gene signature score
		
		Options:
		  -h, --help      show this help message and exit
		  -t <title>      title used for all output files
		  -g <gene sig>   path to a file with a list of genes (HUGO Names) in gene
		                  signature.
		  -f <gene file>  path to normalized gene expression file (TPM or logCPM) with
		                  samples in the columns and genes in the rows 
	
1.   **Figure 1D** `correlation_score_from_flow.py `  

	Input files should have this format (example file for Tcell flow vs Score correlation `files_used_for_plots/tcell_flow_and_Score.tsv`). 
		

					|flow_Tcell |Score
		IPIBLAD032.T1.rna.live	|  31.3	    |  46.98
		IPIBLAD033.T1.rna.live	|   0.66    |   2.68
	
		Input files can be made from calculating scores with get_score.py and joining the corresponding flow data

	`python correlation_score_from_flow.py ../files_used_for_plots/tcell_flow_and_Score.tsv`
1. **Figure 1E** `cross_whisker_plots_flow.py` 

	input files needed 
	
			1. feature scores tsv file output of `get_score.py` (example file in `/files_used_for_plots/tcell_score_percentile_of_percentiles.tsv`
			2. flow scores tsv file output of `get_flow_score.py` (example file in `/files_used_for_plots/Flow_score_Sept10_2020_Tcell_percentile_of_percent.tsv`

	


# Contact Info

Bushra Samad (email: Bushra[dot]Samad[at]ucsf[dot]edu)
Alexis Combes (email: Alexis[dot]Combes[at]ucsf[dot]edu)
Matthew Krummel (email: Matthew[dot]Krummel[at]ucsf[dot]edu)
