## Loads vcf files into sqlite database
library(stringr) 
library(dplyr) 
library(tidyr) 
library(VariantAnnotation)

vcf_to_tbl <- function(ref, vcf_filename){
  vcf <- readVcf(vcf_filename, geno=ref)
  #converting to a datatable
}




# purity_stats()
# calculates general purity statistics for a dataframe

# loading vcfs
ref = "../../data/RM8375/ref/CFSAN008157.HGAP.fasta"
vcf_filename = "../../bioinf/sequence_purity/consensus_base_miseq/RM8375-MiSeq.vcf"