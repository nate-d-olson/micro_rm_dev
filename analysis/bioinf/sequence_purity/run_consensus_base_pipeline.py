
# coding: utf-8

# ###Whole genome vcf files
# 
# Input: Individually mapped bam files     
# Output: Single VCF file
# 
# TODO:
# * clean up input definitions
# * output directory 
# 
#     run_name - root directory
#     \tmp - intermediate files
#     \run_id
#         | final bam
#         \ logs - log files
#     \vcf - platform specific vcf files
# 
# * need reference dict file - create function, seperate python file with scripts for indexing and creating dict for the reference
# 
# Reference dict code
# java -jar ../../usr/local/bin/CreateSequenceDictionary.jar \
#     R=../data/RM8375/ref/CFSAN008157.HGAP.fasta \
#     O=../data/RM8375/ref/CFSAN008157.HGAP.dict

# In[3]:

import sys
import re
import subprocess
from bwa_mem_pe import *
from consensus_base_pipeline import *


# In[4]:

def define_run_params(run_info,pipeline_params):
    # hash for storing run specific information
    from collections import defaultdict
    
    run_params = defaultdict(str)
    run_info_list = run_info.split(":")
    run_params['plat'] = run_info_list[1]
    run_params['run_id'] = run_info_list[0] 
    run_params['vial'] = run_info_list[2]
    run_params['rep'] = run_info_list[3]       
    run_params['out_dir'] = pipeline_params['root_dir'] + pipeline_params['analysis_out_dir']

    # run specific log directory
    run_params['log_dir'] = run_params['out_dir'] + "logs" + "/" + run_params['run_id']
    subprocess.call(['mkdir','-p',run_params['log_dir']])

    bam_root = pipeline_params['root_dir'] + pipeline_params['analysis_out_dir'] + "tmp/" + \
                            pipeline_params['ref_root'] +"_"+ run_params['run_id']
    run_params['bam'] = pipeline_params['root_dir'] + pipeline_params['bam_dir'] + "/" + \
                            pipeline_params['ref_root'] +"_"+ run_params['run_id'] + ".bam"
    run_params['fix_file'] = bam_root + "_fix.bam"
    run_params['header_file'] = bam_root + "_header.bam"
    run_params['sort_file'] = bam_root + "_sort.bam"
    run_params['group_sort_file'] = bam_root + "_group_sort.bam"
    run_params['realign_file'] = bam_root + "_realign.bam"
    run_params['intervals_file'] = bam_root + ".intervals"
    run_params['markdup_file'] = bam_root + "_markdup.bam"
    run_params['metrics_file'] = bam_root + ".metrics"
    run_params['read_group'] = [(("RGID=%s") % pipeline_params['RM']),
                                (("RGLB=%s") % run_params['vial']+"-"+run_params['rep']),
                                (("RGPL=%s") % run_params['plat']),
                                (("RGPU=%s") % run_params['vial']+"-"+run_params['rep']),
                                (("RGSM=%s") % run_params['run_id']),
                                (("RGCN=%s") % pipeline_params['center'])]
    return run_params


# In[5]:

def dedup_realign_single_bam(run_params, pipeline_params):
    ''' Processing single bam file'''
    if run_params['plat'] == "MiSeq":
        bam_group_sort(in_bam = run_params['bam'], 
                       out_bam = run_params['group_sort_file'], 
                       log_dir = run_params['log_dir'])
        
        bam_fixmate(in_bam = run_params['group_sort_file'],
                    out_bam = run_params['fix_file'],
                    log_dir = run_params['log_dir'])

        bam_sort(in_bam = run_params['fix_file'], 
                out_sort = run_params['sort_file'], 
                out_dir = run_params['log_dir'])
    else:
        bam_add_header(in_bam=run_params['bam'], 
                       out_header=run_params['header_file'],
                       log_dir=run_params['log_dir'],
                       read_group=run_params['read_group'])

        bam_sort(in_bam = run_params['header_file'], 
                out_sort = run_params['sort_file'], 
                out_dir = run_params['log_dir'])

    bam_index(bam = run_params['sort_file'],
              out_dir = run_params['log_dir'])
    
    bam_realign(ref = pipeline_params['ref'],
                in_bam = run_params['sort_file'],
                out_bam = run_params['realign_file'], 
                intervals_file = run_params['intervals_file'],
                log_dir = run_params['log_dir'])
    
    bam_markdup(in_bam = run_params['realign_file'], 
                out_bam = run_params['markdup_file'], 
                metrics_file = run_params['metrics_file'],
                log_dir = run_params['log_dir'])
    
    bam_index(bam = run_params['markdup_file'],
              out_dir = run_params['log_dir'])


# In[6]:

def run_consensus_base_pipeline(run_params,pipeline_params):
    
    ## creating file with run parameters
    run_log_file = open(run_params['out_dir']+"/" + run_params['run_id'] +"-run_parameters.txt", 'w')
    run_log_file.write("Parameter\tValue\n")
    for i in run_params.keys():
        run_log_file.write("%s\t%s\n" % (i, run_params[i]))
    for i in pipeline_params.keys():
        run_log_file.write("%s\t%s\n" % (i, pipeline_params[i]))
    run_log_file.close()
    
    ## processing bam
    dedup_realign_single_bam(run_params, pipeline_params)


# In[28]:

def genome_pileups(parameters, markdup_files):
    
    for plat in ["MiSeq", "PGM"]:
        out_dir = parameters['root_dir'] + parameters['analysis_out_dir']
        vcf = out_dir + "/" + parameters['RM'] + "-" + plat + ".vcf"
        genome_calls_mpileup(bams=markdup_files[plat],ref=parameters['ref'],vcf_file=vcf,log_dir=out_dir)


# In[8]:

def read_dat(filename):
    #process input file with configuration information
    from collections import defaultdict
    
    parameters = defaultdict(str)
    with open(filename,'r') as f:
        for line in f:
            param = line.strip().split("=")
            parameters[param[0]] = param[1]
    return parameters


# In[29]:

def main(filename):
    #read run parameters from input file and process using pathoscope
    pipeline_params = read_dat(filename)
    
    # creating temp and log directories
    subprocess.call(['mkdir','-p',pipeline_params['root_dir']+ pipeline_params['analysis_out_dir']+"/tmp/"])
    subprocess.call(['mkdir','-p',pipeline_params['root_dir']+ pipeline_params['analysis_out_dir']+"/logs/"])
    
    
    # list of refined bams
    markdup_files = {'PGM':[], 'MiSeq':[]}
    for i in pipeline_params['datasets'].split(","):
        run_params = define_run_params(i,pipeline_params)
    
        markdup_files[run_params['plat']] += [run_params['markdup_file']]
            
        run_consensus_base_pipeline(run_params, pipeline_params)
       
    # pileups by platform
    genome_pileups(pipeline_params, markdup_files)

if __name__ == '__main__':
    main(sys.argv[1])

