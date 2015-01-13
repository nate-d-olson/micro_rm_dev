
# coding: utf-8

# ## Mapping MiSeq Paired End Reads Using bwa
# To Do
# 1. remove tmp after run
# 2. add run_id to log files
# 3. add parameters/ help with verbose command

import sys
import re
import subprocess
from bwa_mem_pe import *

def define_run_params(run_info,parameters):
    # hash for storing run specific information
    from collections import defaultdict
    
    run_params = defaultdict(str)
    
    run_info_list = run_info.split(":")
    run_params['plat'] = run_info_list[1]
    run_params['run_id'] = run_info_list[0]
    run_params['vial'] = run_info_list[2]
    run_params['rep'] = run_info_list[3]
    
    fastq_root = parameters['root_dir'] + parameters['fastq_dir'] + run_info_list[1] + '/'+ 'fastq'+'/'
    if run_info_list[1] == "MiSeq":
        run_params['fastq1'] = fastq_root + run_info_list[0] + "_1.fastq"
        run_params['fastq2'] = fastq_root + run_info_list[0] + "_2.fastq"
    else:
        run_params['fastq1'] = fastq_root + run_info_list[0] + ".fastq"
        run_params['fastq2'] = None
    
    ref_name = parameters['ref'].split("/")[-1].replace(".fasta","_")
    run_params['out_dir'] = parameters['root_dir'] + parameters['analysis_out_dir']
    run_params['sam'] = run_params['out_dir'] + "/tmp/" + ref_name + run_params['run_id'] + ".sam"
    run_params['bam'] = run_params['out_dir'] + "/tmp/" + ref_name + run_params['run_id'] + ".bam"
    run_params['sorted_bam'] = run_params['out_dir'] + "/"+ ref_name + run_params['run_id'] + ".bam"
    run_params['read_group'] = "@RG\tID:%s\tCN:%s\tLB:%s\tPL:%s\tSM:%s" %(parameters['RM'],
                                                                          parameters['center'],
                                                                          run_params['vial']+"-"+run_params['rep'],
                                                                          run_params['plat'],
                                                                          parameters['run_id'])
    return run_params

def run_bwa_mem_pipeline(run_params,bwa_params):
    
    ## creating file with run parameters
    run_log_file = open(run_params['out_dir']+"/" + run_params['run_id'] +"-run_parameters.txt", 'w')
    run_log_file.write("Parameter\tValue\n")
    for i in run_params.keys():
        run_log_file.write("%s\t%s\n" % (i, run_params[i]))
    for i in bwa_params.keys():
        run_log_file.write("%s\t%s\n" % (i, bwa_params[i]))
    run_log_file.close()
    
    ## running bwa_mem
    bam_map_bwa_mem(ref = bwa_params['root_dir'] + bwa_params['ref'],
                        fq1 = run_params['fastq1'],
                        fq2 =run_params['fastq2'],
                        out_sam = run_params['sam'],
                        out_dir = run_params['out_dir'],
                        read_group = run_params['read_group'])
    sam_to_bam(in_sam = run_params['sam'], 
                     out_bam = run_params['bam'], 
                     out_dir = run_params['out_dir'])
    bam_sort(in_bam = run_params['bam'],
                 out_sort = run_params['sorted_bam'],
                 out_dir = run_params['out_dir'])
    bam_index(out_dir = run_params['out_dir'], bam = run_params['sorted_bam'])

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
    #read run parameters from input file and map reads to reference using bwa
    parameters = read_dat(filename)
    
    # creating temp directory
    subprocess.call(['mkdir','-p',parameters['root_dir']+ parameters['analysis_out_dir']+"/tmp/"])
    
    # indexing reference
    bwa_index_ref(parameters['root_dir'] + parameters['ref'],parameters['root_dir']+parameters['analysis_out_dir'])
    
    for i in parameters['datasets'].split(","):
        run_params = define_run_params(i,parameters)
        run_bwa_mem_pipeline(run_params, parameters)

if __name__ == '__main__':
    main(sys.argv[1])

