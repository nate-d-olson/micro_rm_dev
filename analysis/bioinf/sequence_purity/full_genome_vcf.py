import sys
import time

def prep_bam_for_variant_calling(bam, ref, known_vcf):
    ''' Function for generating variant calls from bam files, 
        adapted from  
        http://www.htslib.org/workflow/#mapping_to_variant
    '''
    import subprocess
    bam_split = bam.split("/")[-1]
    ref_split = ref.split("/")[-1]
    bam_root = ref_split.replace(".fasta","") + "_" + bam_split.replace(".bam","")
    
    GATK_command = ["java","-Xmx4g","/micro_rm/utilities/GenomeAnalysisTK.jar"]
    
    ##%%## Need to work out logging code
    log_file = open(bam_root + ".log",'w')
    log_file.write("Fix mate pairs for: %s" % (bam))
    # fix mate pairs
    bam_fix = "tmp/" + bam_root + ".fix.bam"
    #not for PGM
    subprocess.call(["samtools","fixmate",bam,bam_fix], stdout="log/bam_fix" + timestamp() + log_file)
    
    # sort
    bam_sort = "tmp/" + bam_root + ".sort.bam"
    subprocess.call(["samtools", "sort", "-O", bam, "-o", bam_sort, bam_fix], \
                    stderr="log/bam_sort" + timestamp() + log_filelog_file)
    
    
    GATK_command = ["java","-Xmx4g","../../utilities/GenomeAnalysisTK.jar"]
    # realignment around indels
    intervals = "tmp/" + bam_root + ".intervals"
    bam_realign = "tmp/" + bam_root + "realign.bam"
    realigner_target_command = GATK_command + ["-T","RealignerTargetCreator", "-R",
                                             ref,"-I",bam_sort, "-o", intervals]
    subprocess.call(realigner_target_command, stderr="realign_target" + timestamp() + log_file)
    realigner_command = GATK_command + ["-T","IndelRealigner", "-R", ref,"-I",bam_sort,
                                "-targetIntervals", intervals, "-o", bam_realign]
    subprocess.call(realigner_command, stderr="log/bam_realign" + timestamp() + log_file)
    

def main(ref, known_vcf, vcf_filename, bams_to_process):
    # processing multiple bams
    processed_bams = []
    for i in bams_to_process:
        processed_bams.append(prep_bam_for_variant_calling(i, ref,known_vcf))
    genome_calls_mpileup(processed_bams,ref, vcf_filename)

# Need to check also incorporate command line parsing
def __name__ if name == '__main__':
    ref=sys.argv[1] #reference genome reads were mapped to
    known_vcf=sys.argv[2] # vcf of known variants
    vcf_filename=sys.argv[3] # name for output file
    bams = sys.argv[4:] #bams to process
    main(sys.argv)