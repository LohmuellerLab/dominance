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
REF=/u/home/a/arun/project-klohmuel/dominance/data/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
MEM=10G
mkdir -p $scratch/thaliana-temp

#######################
#  1. call SNPs

java -jar ~/software/GenomeAnalysisTK.jar \
 -T UnifiedGenotyper \
 -R ${REF} \
 -I ${data_dir}/22004_AAGAGATC_C5E9DANXX_6_20140930B_20140930.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/22006_GAGTTAGC_C5E9DANXX_6_20140930B_20140930.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35513_GTGGCC_C8DE0ANXX_4_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35594_GCCAAT_C8DE0ANXX_5_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35599_GTTTCG_C8DE0ANXX_5_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35600_CGTACG_C8DE0ANXX_5_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35601_TGACCA_C8DE0ANXX_5_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35605_CAGATC_C8DE0ANXX_6_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35607_ACTGAT_C8DE0ANXX_6_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35608_CGATGT_C8DE0ANXX_6_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35610_GGCTAC_C8DE0ANXX_6_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35611_ATTCCT_C8DE0ANXX_6_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35615_TAGCTT_C8DE0ANXX_6_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35618_GCCAAT_C8DE0ANXX_7_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35619_AGTTCC_C8DE0ANXX_7_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35621_CGTACG_C8DE0ANXX_7_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35622_GTGAAA_C8DE0ANXX_7_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35623_CAGATC_C8DE0ANXX_7_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/35625_GTCCGC_C8DE0ANXX_7_20160129B_20160129.cleanedaligned.markdup.indel.bam.reheader.bam \
 -I ${data_dir}/37469_CAGATC_C89VWANXX_2_20160414B_20160414.cleanedaligned.markdup.indel.bam.reheader.bam \
 -o ${working_dir}/vcf/thaliana.vcf \
 -ploidy 1 \
 --output_mode EMIT_ALL_SITES

# -I ${data_dir}/18510_AAGACGGA_C41F2ACXX_6_20140423B_20140423.cleanedaligned.markdup.indel.bam \
# -I ${data_dir}/18512_CGAACTTA_C41F2ACXX_6_20140423B_20140423.cleanedaligned.markdup.indel.bam \
# -I ${data_dir}/18514_TGAAGAGA_C41F2ACXX_6_20140423B_20140423.cleanedaligned.markdup.indel.bam \
# -I ${data_dir}/22002_GAATCTGA_C5E9DANXX_6_20140930B_20140930.cleanedaligned.markdup.indel.bam \
