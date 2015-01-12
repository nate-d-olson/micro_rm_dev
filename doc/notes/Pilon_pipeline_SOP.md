# Pilon pipeline SOP
Requirements mapped sequence files

## Steps
1. Starting docker
	* for mac open boot2docker (double click application), otherwise run docker from commandline in linux
2. Mapping 
	2.1. MiSeq reads
		2.1.1 generate param file for run
			* bwa_mem_pe.py has functions for running bwa for Illumina PE data, run_bwa_mem_pe.py runs datasets based on parameters provided input file
			* see bwa_mem_pipeline_params.txt for example
		2.2.2 start docker container `docker run -i -t -v ~/Desktop/micro_rm_dev:/notebooks varcall_notebook sh`
		2.2.3 inside the container run `python run_bwa_mem_pe.py params_file.txt` from the `analysis/bioinf/sequence_purity` directory
	2.2 PGM reads
		* need to workout, either generate a new pipeline or using TorrentSuite

3. Running pilon
	3.1 Start ipython notebook container in docker terminal run `docker run -d -p 443:8888 -v ~/Desktop/micro_rm_dev:/notebooks -e "PASSWORD=micro_rm_pw" varcall_notebook`
	3.2. open ipython notebook`localhost:8888` in browser, password "micro_rm_pw"
	3.3 merge bam files from individual datasets using samtools
		see notebook_2014_01_09_Pilon.pynb
	3.4 run pilon on merged bams - error running using Docker, memory issue, running from the commandline. Check out http://stackoverflow.com/questions/24422123/change-boot2docker-memory-assignment, blog post on changing memory usage for boot2docker
	