## load all modules here
set -e
set -u
set -o

. /u/local/Modules/default/init/modules.sh
module load bwa
module load samtools
module load java/1.8.0_77
module load bcftools

#######################
## set up variables
scratch=/u/scratch/a/arun
working_dir=/u/home/a/arun/project-klohmuel/dominance/data
REF=/u/home/a/arun/project-klohmuel/dominance/data/reference/Arabidopsis_lyrata.v.1.0.dna.toplevel.fa
MEM=10G

#######################
#  1. call SNPs

java -jar ~/software/GenomeAnalysisTK.jar \
 -T UnifiedGenotyper \
 -R ${REF} \
 -I ${working_dir}/aligned/${1}.cleanedaligned.markdup.indel.bam \
 -o ${working_dir}/vcf/lyrata_single_sample/${1} \
 --output_mode EMIT_ALL_SITES 

bgzip -c ${working_dir}/vcf/lyrata_single_sample/${1} > ${working_dir}/vcf/lyrata_single_sample/${1}.vcf.b.gz