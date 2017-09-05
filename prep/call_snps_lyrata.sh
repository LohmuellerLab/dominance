## load all modules here
set -e
set -u
set -o

. /u/local/Modules/default/init/modules.sh
module load bwa
module load samtools
module load java/1.8.0_77


#######################
## set up variables
scratch=/u/scratch/a/arun
working_dir=/u/home/a/arun/project-klohmuel/dominance/data
data_dir=/u/home/a/arun/project-klohmuel/dominance/data/aligned
REF=/u/home/a/arun/project-klohmuel/dominance/data/reference/Arabidopsis_lyrata.v.1.0.dna.toplevel.fa
MEM=10G
mkdir -p $scratch/lyrata-temp

#######################
#  1. call SNPs

java -jar ~/software/GenomeAnalysisTK.jar \
 -T UnifiedGenotyper \
 -R ${REF} \
 -I ${data_dir}/SRR2040788.cleanedaligned.markdup.indel.bam \
 -I ${data_dir}/SRR2040789.cleanedaligned.markdup.indel.bam \
 -I ${data_dir}/SRR2040790.cleanedaligned.markdup.indel.bam \
 -I ${data_dir}/SRR2040791.cleanedaligned.markdup.indel.bam \
 -I ${data_dir}/SRR2040792.cleanedaligned.markdup.indel.bam \
 -I ${data_dir}/SRR2040793.cleanedaligned.markdup.indel.bam \
 -I ${data_dir}/SRR2040794.cleanedaligned.markdup.indel.bam \
 -I ${data_dir}/SRR2040795.cleanedaligned.markdup.indel.bam \
 -I ${data_dir}/SRR2040796.cleanedaligned.markdup.indel.bam \
 -I ${data_dir}/SRR2040797.cleanedaligned.markdup.indel.bam \
 -I ${data_dir}/SRR2040798.cleanedaligned.markdup.indel.bam \
 -I ${data_dir}/SRR2040826.cleanedaligned.markdup.indel.bam \
 -I ${data_dir}/SRR2040827.cleanedaligned.markdup.indel.bam \
 -I ${data_dir}/SRR2040828.cleanedaligned.markdup.indel.bam \
 -I ${data_dir}/SRR2040829.cleanedaligned.markdup.indel.bam \
 -I ${data_dir}/SRR2040830.cleanedaligned.markdup.indel.bam \
 -o ${working_dir}/vcf/lyrata.vcf \
 --output_mode EMIT_ALL_SITES