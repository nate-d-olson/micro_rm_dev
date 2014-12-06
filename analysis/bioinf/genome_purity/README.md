## Genomic Purity Analysis
Comparison of genomic sequence data using pathoscope2 to the NCBI nt database

## Pathoscope and dependencies setup 
* Pathoscope version 2.0 downloaded from http://sourceforge.net/projects/pathoscope/
* Pathoqc version 1.0 downloaded from http://sourceforge.net/projects/pathoscope/
	* followed README for install, had to include pp python module in pathoqc directory to run
* Commands used in pathoscope_pipeline.py based on Pathoscope tutorial
http://www.microbiomejournal.com/content/supplementary/2049-2618-2-33-s5.pdf
* Downloaded nt sequence data from ftp://oligomer.
bumc.bu.edu/data/nt_ti.fa.gz on 11/19/2014
* Using bowtie 2 version 2.2.3 for mapping, version installed as part of biolinux8

## Running pathoscope pipeline
1. create a run parameters file 
	* Parameters and values seperated by '=', with each one individual lines in a file
		* qc_loc - path to pathoqc
		* path_loc - path to pathoscope.py
		* ref - location of reference sequence file
		* index_dir - directory with indexed reference files
		* analysis_out_dir - output directory for analysis
		* datasets - comma seperated list of SRR:plat for each dataset to run
		* n - number of threads to use
		* root_dir - project root directors - full system path
2. Commands for pathoqc and the pathoscope MAP and ID modules are in the pathocope_pipeline.py script
3. can run from command line run_pathoscope_pipeline.py parameters.txt