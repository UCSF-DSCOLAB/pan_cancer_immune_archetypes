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
	
3.   **Figure 1D** `correlation_score_from_flow.py `  

	Input files should have this format (example file for Tcell flow vs Score correlation `files_used_for_plots/tcell_flow_and_Score.tsv`). 
		

					|flow_Tcell |Score
		IPIBLAD032.T1.rna.live	|  31.3	    |  46.98
<<<<<<< HEAD
		IPIBLAD033.T1.rna.live	|  	0.66    |   2.68
	
		Input files can be made from calculating scores with get_score.py and joining the corresponding flow data

	`python correlation_score_from_flow.py ../files_used_for_plots/tcell_flow_and_Score.tsv`
1. **Figure 1E** `cross_whisker_plots_flow.py` 

	input files needed 
	
			1. feature scores tsv file output of `get_score.py` (example file in `/files_used_for_plots/tcell_score_percentile_of_percentiles.tsv`
			2. flow scores tsv file output of `get_flow_score.py` (example file in `/files_used_for_plots/Flow_score_Sept10_2020_Tcell_percentile_of_percent.tsv`

### **Figure 2**
1. **Figures 2A-C, 2E-F, 2J, 2L** `3_feature_UMAP_IPI.py`

	input file needed 
		
				1. feature matrix tsv file output (example file in `/files_used_for_plots/3_feature_matrix.txt for IPI. 
				2. The feature scores for TCGA using not included but can be generated using TCGA data and `get_score.py``
				
1.  **Figures 2H-I, 2K** boxplots made using `seaborn.boxplot` with data in `/files_used_for_plots/3_feature_matrix.txt`
2.  **Figure 2M** heatmap was made using `seaborn.clustermap` with medians calculated from TPMs in `/files_used_for_plots/3_feature_matrix.txt`
3. **Figure 2N** barplot and bubble plots were made using `seaborn.barplot` and `seaborn.scatter` respectively with medians calculated from TPMs in ``/files_used_for_plots/3_feature_matrix.txt` for IPI the chemokine TPMs for TCGA are not included
 
### **Figure 3**
1. **Figures 3A** `3_feature_UMAP_IPI.py`

	input file needed 
		
				1. feature matrix tsv file output (example file in `/files_used_for_plots/3_feature_matrix.txt for IPI. 
1. **Figures 3B-C** `TCGA_survival_3_features.py`

input files needed

Can be found here
 [XENA Browser](https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)		
	
				
			   1. A folder with TCGA name mappings files for each indication in the format below
				   
					  new_name	| tcga
					0	BLCA_0		| TCGA-2F-A9KO-01
					1	BLCA_1		| TCGA-2F-A9KP-01
					2	BLCA_2		| TCGA-2F-A9KQ-01
					
 			   2. A folder with TCGA clinical data files for each indication 
					


2. **Figures 3D-E** boxplots made using `seaborn.boxplot` with data in `/files_used_for_plots/3_feature_matrix.txt`

### **Figure 4**
1. **Figure 4A** `get_flow_score.py`
2. **Figure 4B** heatmap was made using `seaborn.clustermap`
3. **Figures 4C** `6_feature_UMAP_IPI.py`

	input file needed 
		
				1. feature matrix tsv file output (example file in `/files_used_for_plots/6_feature_matrix.txt for IPI.
1. **Figures 4D,4G** boxplots made using `seaborn.boxplot` with data in `/files_used_for_plots/6_feature_matrix.txt`
2. **Figures 4E** the alluvial plot was made using [RAWgraphs 2.0](https://www.rawgraphs.io/)
3. **Figure 4F** heatmap was made using `seaborn.clustermap` with medians calculated from TPMs in `/files_used_for_plots/6_feature_matrix.txt`

### **Figure 5**
1. **Figure 5A** `get_flow_score.py`
2. **Figures 5A, 5O, 5R**
`10_feature_UMAP_IPI.py`

	input file needed 
		
				1. feature matrix tsv file output (example file in `/files_used_for_plots/10_feature_matrix.txt for IPI.
				
1. **Figures 5B** the alluvial plot was made using [RAWgraphs 2.0](https://www.rawgraphs.io/)
2. **Figures 5D-N** boxplots and scatter plots were made using `seaborn.boxplot` and `seaborn.scatter` respectively.
3. **Figures 5P-Q, 5S-T** boxplots made using `seaborn.boxplot` with data in `/files_used_for_plots/10_feature_matrix.txt`

### **Figure 6**
1. **Figures 6A-D** heatmap and bubble plots were made using `seaborn.clustermap` and `seaborn.scatter` respectively with medians of chemokines calculated from TPMs in ``/files_used_for_plots/10_feature_matrix.txt`. TPMS for other genes can be 	found on using the GEO accession GSE184398

### **Figure 7**
1. **Figures 7A-B**
`10_feature_UMAP_IPI.py`

	input file needed 
		
				1. feature matrix tsv file output (example file in `/files_used_for_plots/10_feature_matrix.txt for IPI.
1. **Figures 7B-E** boxplot, heatmaps and bubble plot were made using `seaborn.boxplot` `seaborn.clustermap` and `seaborn.scatter` respectively.TPMS for genes can be 	found on using the GEO accession GSE184398
2. **Figures 7F** `cox_regression.py` 				


##Supplementary Figures

### **Figure 1S**

1. **Figure 1SA** 
	
	Kmeans clustering using `sklearn.cluster.KMeans` of and dataframe of TPMS of all compartments. These can be 	found on using the GEO accession GSE184398

1. **Figure 1SB**

		makes a volcano plot from a standard Limma DGE output
		e.g.
		           logFC    |   AveExpr	    |       t	    |  P.Value	|  adj.P.Val	|      B
		 PLK1	7.374675473	| 1.415522735	| 6.886987802	| 3.92E-06	| 0.050802876	| 3.971754014
		 ZWINT	5.772193613	| 1.340166868	| 5.872941797	| 2.48E-05	| 0.065626979	| 2.450364751
				
		
		python volcano_plots.py -h
		
		Usage: volcano_plots.py -i <input file name> [options]
		
		Make Volcano Plots
		
		Options:
		  -h, --help      show this help message and exit
		  -i <input file>   Path and name of input file
		  -o <output file>  Path and names of output files e.g DGE_treg_vs_tcell (a .png,.pdf & .svg file extension will automatically be added)
          -p <p-value>      P-value cutoff DEFAULT 0.005           
		  -t <title>      Title in quotes e.g. "Tregs vs Tcells"
		 
1. 	**Figure 1SC**	 `correlation_score_from_flow.py ` see figure 1D
2. **Figure 1SD** `get_flow_score.py` see figure 1C
3. **Figure 1SE** `cross_whisker_plots_flow.py`  see figure 1E


### **Figure 2S**	         	

1. **Figure 2SA & 2SB** `DBI_frequency_of_clusters.py`
	
	input file needed 
		
				1. feature matrix tsv file output (example file in `/files_used_for_plots/3_feature_matrix.txt for IPI.
			
1. **Figures 2SD** `3_feature_UMAP_IPI.py`

	input file needed 
		
				1. feature matrix tsv file output (example file in `/files_used_for_plots/3_feature_matrix.txt for IPI.

1. **Figures 2SE-F** boxplots and scatter plots were made using `seaborn.boxplot` and `seaborn.scatter` respectively.

### **Figure 3S**	  

1. **Figures 3SA** `3_feature_UMAP_IPI.py`
2. **Figures 3SB-D** `cox_regression.py` 	

### **Figure 4S**	 

1. **Figure 4SB**

		makes a volcano plot from a standard Limma DGE output
		e.g.
		           logFC    |   AveExpr	    |       t	    |  P.Value	|  adj.P.Val	|      B
		 PLK1	7.374675473	| 1.415522735	| 6.886987802	| 3.92E-06	| 0.050802876	| 3.971754014
		 ZWINT	5.772193613	| 1.340166868	| 5.872941797	| 2.48E-05	| 0.065626979	| 2.450364751
				
		
		python volcano_plots.py -h
		
		Usage: volcano_plots.py -i <input file name> [options]
		
		Make Volcano Plots
		
		Options:
		  -h, --help      show this help message and exit
		  -i <input file>   Path and name of input file
		  -o <output file>  Path and names of output files e.g DGE_treg_vs_tcell (a .png,.pdf & .svg file extension will automatically be added)
          -p <p-value>      P-value cutoff DEFAULT 0.005           
		  -t <title>      Title in quotes e.g. "Tregs vs Tcells"
1. **Figure 4SC, 4Sk**	 `correlation_score_from_flow.py ` see figure 1D
2. **Figures 4SD, 4SF-G, 4SL, 4SP** `6_feature_UMAP_IPI.py`

	input file needed 
		
				1. feature matrix tsv file output (example file in `/files_used_for_plots/6_feature_matrix.txt

1.  **Figure 4SE** `DBI_frequency_of_clusters.py`
	
	input file needed 
		
				1. feature matrix tsv file output (example file in `/files_used_for_plots/6_feature_matrix.txt.
1. **Figures 4SG-H, 4SI, 4SN-N** boxplots and scatter plots were made using `seaborn.boxplot` and `seaborn.scatter` respectively.

### **Figure 5S**	

1.  **Figure 5SF**	 `correlation_score_from_flow.py ` see figure 1D
2. **Figure 5SG** `get_flow_score.py` see figure 1C
3. **Figure 4SE** `DBI_frequency_of_clusters.py`
	
	input file needed 
		
				1. feature matrix tsv file output (example file in `/files_used_for_plots/10_feature_matrix.txt.

1. **Figures 5SI-J** `10_feature_UMAP_IPI.py`

	input file needed 
		
				1. feature matrix tsv file output (example file in `/files_used_for_plots/10_feature_matrix.txt


### **Figure 6S**

1. **Figures 6SA, 6SE-F** `10_feature_UMAP_IPI.py`

	input file needed 
		
				1. feature matrix tsv file output (example file in `/files_used_for_plots/10_feature_matrix.txt
	
1. **Figures 6SC-E, 6SG** boxplots made using `seaborn.boxplot` with data in `/files_used_for_plots/10_feature_matrix.txt`

### **Figure 7S**
1. **Figures 7SA-K** `cox_regression.py`

# Contact Info

Bushra Samad (email: Bushra[dot]Samad[at]ucsf[dot]edu)
Alexis Combes (email: Alexis[dot]Combes[at]ucsf[dot]edu)
Matthew Krummel (email: Matthew[dot]Krummel[at]ucsf[dot]edu)
=======
		IPIBLAD033.T1.rna.live	|   0.66    |   2.68
	
		Input files can be made from calculating scores with get_score.py and joining the corresponding flow data

	`python correlation_score_from_flow.py ../files_used_for_plots/tcell_flow_and_Score.tsv`
4. **Figure 1E** `cross_whisker_plots_flow.py` 

	input files needed 
	
			1. feature scores tsv file output of `get_score.py` (example file in `/files_used_for_plots/tcell_score_percentile_of_percentiles.tsv`
			2. flow scores tsv file output of `get_flow_score.py` (example file in `/files_used_for_plots/Flow_score_Sept10_2020_Tcell_percentile_of_percent.tsv`

	


# Contact Info

Bushra Samad (email: Bushra[dot]Samad[at]ucsf[dot]edu)
Alexis Combes (email: Alexis[dot]Combes[at]ucsf[dot]edu)
Matthew Krummel (email: Matthew[dot]Krummel[at]ucsf[dot]edu)
>>>>>>> a2fe0a0d50729d1c8ab85fa1eee7e9e44d81622b
