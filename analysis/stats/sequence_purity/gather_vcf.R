## Loads vcf files into sqlite database
#library(stringr) 
#library(tidyr)
library(data.table)
library(dplyr) 

# library(RSQLite)
# library(VariantAnnotation)

cbtsv_tosql <- function (tsv_file,db_conn, tbl_name) {
  pgm <- fread(tsv_file,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  pgm <- rename(pgm, PLATDP=DP) # added to fix issue with two DP columns
  copy_to(vcf_db, pgm, name = "full_pgm", temporary = FALSE, indexes = list("CHROM","POS","SAMPLE")) 
}


pgm_tsv <- "analysis/stats/sequence_purity/RM8375-PGM.tsv"
vcf_db <- src_sqlite("data/RM8375/RM8375.sqlite", create = TRUE)

cbtsv_tosql(pgm_tsv, vcf_db,"cb_pgm")
