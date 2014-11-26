!#/usr/bin/bash
### Script for processing Seq data using pathoscope



### Input file
#paired end
FQ1=$1
FQ2=$2
SEQ_ROOT=sed($FQ1 's/_1.fq//')

#single end
FQ=$1
SEQ_ROOT=sed($FQ 's/.fq//')

### Input directory
DAT_DIR=

### Output directory
OUT_DIR_ROOT=analysis/bioinf/genome_purity
OUT_DIR=$OUT_DIR_ROOT/

### Pathoscope locations


## Pathoqc
python pathoqc.py -1 ../../../data/RM8375/MiSeq/fastq/SRR1555296_1.fastq -2 ../../../data/RM8375/MiSeq/fastq/SRR1555296_2.fastq -s Illumina -p 8 -o ../../../analysis/bioinf/genome_structure/SRR1555296

## Pathomap
python pathoscope.py MAP -U fastq -targetRefFiles nt.fasta  -outDir results -outAlign <input_name>.sam -expTag LT2_micro_rm
used -U for single end and -1, -2 for paired end

../../../utilities/pathoscope2/pathoscope/pathoscope.py MAP -1 ../../../data/RM8375/MiSeq/fastq/SRR1555296_1.fastq -2 ../../../data/RM8375/MiSeq/fastq/SRR1555296_2.fastq -targetRefFiles ../../../utilities/patho_utils/micro_rm_patho_db_ti.fa -indexDir ../../../utilities/patho_utils/ -outDir SRR1555296 -outAlign SRR155296.sam -expTag RM8375

## PathoID
nolson@Selenium:/media/nolson/second/current_projects/micro_rm_dev/analysis/bioinf/genome_purity/SRR1555296$ python ../../../../utilities/pathoscope2/pathoscope/pathoscope.py ID -alignFile SRR155296.sam -fileType sam -outDir ./ -expTag RM8375-ID
