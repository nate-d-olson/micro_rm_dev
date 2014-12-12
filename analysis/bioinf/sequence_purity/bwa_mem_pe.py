
# coding: utf-8

# ### Functions for mapping paired end reads to reference fasta using bwa-mem
# * map_bwa_mem - maps reads using bwa mem algorithm
# * bwa_index_ref - indexing reference genome
# * bam_sort - sorting bam file
# * bam_index - indexing bam
# * sam_to_bam - converting alignment format
# 

# In[13]:

import sys
import time
import subprocess


# In[14]:

def bwa_index_ref(ref,out_dir):
    ## log files for standard out and error
    log_file = open(out_dir + "/bwa_index_ref"+ time.strftime("-%Y-%m-%d-%H-%M-%S.log"),'w')
    stderr_file = open(out_dir + "/bwa_index_ref"+ time.strftime("-%Y-%m-%d-%H-%M-%S.stder"),'w')

    ## run command
    bwa_index_command = ["bwa","index",ref]
    subprocess.call(bwa_index_command, stdout=log_file, stderr=stderr_file)
    log_file.close(); stderr_file.close()


# In[15]:

def bam_map_bwa_mem(ref, fq1, fq2, out_sam, out_dir, read_group):
    ''' Mapping paired-end reads to reference'''
    
    ## log files for standard out and error
    sam_file = open(out_sam,'w')
    stderr_file = open(out_dir + "/bwa_mem"+ time.strftime("-%Y-%m-%d-%H-%M-%S.stder"),'w')
    
    ## run command
    bwa_mem_command = ["bwa","mem","-t","8","-R",read_group, ref,fq1,fq2]
    subprocess.call(bwa_mem_command, stdout=sam_file, stderr=stderr_file)  
    sam_file.close(); stderr_file.close()


# In[16]:

def sam_to_bam(in_sam, out_bam, out_dir):
    '''Convert sam to bam'''
    
    ## log files for standard out and error
    bam_file = open(out_bam,'w')
    stderr_file = open(out_dir + "/samtools_view"+ time.strftime("-%Y-%m-%d-%H-%M-%S.stder"),'w')
    
    ## run command
    sam_to_bam_command = ["samtools","view","-b",in_sam]
    subprocess.call(sam_to_bam_command, stdout=bam_file, stderr=stderr_file)  
    bam_file.close(); stderr_file.close()


# In[17]:

def bam_sort(in_bam, out_sort, out_dir):
    ''' Sorting bam'''
    
    ## log files for standard out and error
    out_file = open(out_sort,'w')
    stderr_file = open(out_dir + "/bwa_sort"+ time.strftime("-%Y-%m-%d-%H-%M-%S.stder"),'w')
    
    ## run command
    bam_sort_command = ["samtools","sort","-T","sort_temp","-O", "bam", in_bam]
    subprocess.call(bam_sort_command, stdout=out_file,stderr=stderr_file) 
    out_file.close(); stderr_file.close()


# In[18]:

def bam_index(out_dir, bam):
    ''' index bam file'''
    
    ## log files for standard out and error
    log_file = open(out_dir + "/bwa_index"+ time.strftime("-%Y-%m-%d-%H-%M-%S.log"),'w')
    stderr_file = open(out_dir + "/bwa_index"+ time.strftime("-%Y-%m-%d-%H-%M-%S.stder"),'w')
    
    ## run command
    subprocess.call(["samtools","index",bam],stdout=log_file,stderr=stderr_file)
    log_file.close(); stderr_file.close()

