#!/usr/bin/sh
# Bash script for running reapr in a docker container

# passing fastq and ref files from command line
REF=$1
FQ_1=$2
FQ_2=$3
# for single-end reads
# reapr facheck $REF ref.facheck
# reapr smaltmap ref.facheck.fa $FQ_1 $FQ_2 long_smalt_mapped.bam
reapr pipeline ref.facheck.fa long_smalt_mapped.bam reapr_out