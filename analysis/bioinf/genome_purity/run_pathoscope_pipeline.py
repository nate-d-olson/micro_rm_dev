# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ## Processing micro RM datasets
# Extracts run information from input file (format defined below) and runs pathoqc, map, and id.  pathoscope_pipeline.py includes code for running the three python modules.
# 
# ### Parameter File Description
# Parameters and values seperated by '=', with each one individual lines in a file
# * qc_loc - path to pathoqc
# * path_loc - path to pathoscope.py
# * ref - location of reference sequence file
# * index_dir - directory with indexed reference files
# * analysis_out_dir - output directory for analysis
# * datasets - comma seperated list of SRR:plat for each dataset to run
# * n - number of threads to use
# * root_dir - project root directors - full system path
# 
# ### TODO
# * fix code to define the number of threads
# * Parallelize pathoid - possibly need second loop?


import sys
import re
import subprocess
from pathoscope_pipeline import * #need to place pathoscope_pipeline.py in directory

def define_run_params(run_info,parameters):
    # hash for storing run specific information
    from collections import defaultdict
    
    run_params = defaultdict(str)
    
    run_info_list = run_info.split(":")
    run_params['plat'] = run_info_list[1]
    run_params['run_id'] = run_info_list[0]
    
    fastq_root = parameters['root_dir'] + parameters['fastq_dir'] + \
                 run_info_list[1] + '/'+ 'fastq'+'/'

    if run_info_list[1] == "MiSeq":
        run_params['fastq1'] = fastq_root + run_info_list[0] + "_1.fastq"
        run_params['fastq2'] = fastq_root + run_info_list[0] + "_2.fastq"
    else:
        run_params['fastq1'] = fastq_root + run_info_list[0] + ".fastq"
        run_params['fastq2'] = None
    
    run_params['out_dir'] = parameters['root_dir'] + \
                            parameters['analysis_out_dir'] + \
                            '/' + run_params['run_id']
    return run_params

def run_pathoscope(run_params,patho_params):
    import glob
    # creating run directory
    subprocess.call(['mkdir',run_params['out_dir']])
    
    # creating file with run parameters
    run_log_file = open(run_params['out_dir']+"/run_parameters.txt", 'w')
    run_log_file.write("Parameter\tValue\n")
    for i in run_params.keys():
        run_log_file.write("%s\t%s\n" % (i, run_params[i]))
    for i in patho_params.keys():
        run_log_file.write("%s\t%s\n" % (i, patho_params[i]))
    run_log_file.close()
    
    ## running pathoqc
    pathoqc_command(plat=run_params['plat'],
                    fastq1=run_params['fastq1'],
                    fastq2=run_params['fastq2'],
                    out_dir=run_params['out_dir'],
                    path_pathoqc=patho_params['root_dir']+ patho_params['qc_loc'],
                    thread_num= 8)

    ## running pathomap
    if run_params['fastq2'] != None:
        trimmed_fastq1 = run_params['out_dir'] + '/' + run_params['run_id'] + '_1_tr.fq'
        trimmed_fastq2 = run_params['out_dir'] + '/' + run_params['run_id'] + '_2_tr.fq'
    else:
        trimmed_fastq1 = run_params['out_dir'] + '/' + run_params['run_id'] + '_tr.fq'
        trimmed_fastq2 = None
        
    pathomap_command(ref_path=patho_params['root_dir']+patho_params['ref'],
                     index_dir=patho_params['root_dir']+patho_params['index_dir'], 
                     exptag = run_params['run_id'], 
                     fastq1=trimmed_fastq1, 
                     fastq2=trimmed_fastq2, 
                     out_dir=run_params['out_dir'], 
                     path_pathoscope=patho_params['root_dir']+patho_params['patho_loc'])
    
    ## cleaning up pathomap files- removes find combined in appendAlign
    cleanup_log = open(run_params['out_dir']+"/"+"cleanup.log",'w')
    sam_files = glob.glob(run_params['out_dir']+"/"+run_params['run_id']+"-"+ \
                          re.sub('.fasta$|.fa$','',patho_params['ref'].split("/")[-1])+'*sam')
    cleanup_command = ['rm'] + sam_files
    subprocess.call(cleanup_command, stdout=cleanup_log)
    
    ## running pathoid
    if run_params['fastq2'] != None:
        sam_suffix = "_1_tr.sam"
    else:
        sam_suffix = "_tr.sam"
    pathomap_sam = run_params['out_dir']+"/" + run_params['run_id']+sam_suffix
    pathoid_command(path_pathoscope=patho_params['root_dir']+patho_params['patho_loc'], \
        input_sam=pathomap_sam, out_dir = run_params['out_dir'], exptag = run_params['run_id'])

def read_dat(filename):
    #process input file with configuration information
    from collections import defaultdict
    
    parameters = defaultdict(str)
    with open(filename,'r') as f:
        for line in f:
            param = line.strip().split("=")
            parameters[param[0]] = param[1]
    return parameters

def main(filename):
    #read run parameters from input file and process using pathoscope
    parameters = read_dat(filename)
    
    for i in parameters['datasets'].split(","):
        run_params = define_run_params(i,parameters)
        run_pathoscope(run_params, parameters)

if __name__ == '__main__':
    main(sys.argv[1])

