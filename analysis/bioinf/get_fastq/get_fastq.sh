#!/usr/bin/sh
# downloads fastq files for accession numbers in an 
# input file, accession numbers on individual lines in file
while read p; do
  /sratoolkit.2.4.2-ubuntu64/bin/fastq-dump -A $p
done < $1