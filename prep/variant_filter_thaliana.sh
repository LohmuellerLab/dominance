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
data_dir=/u/home/a/arun/project-klohmuel/dominance/data/aligned
REF=/u/home/a/arun/project-klohmuel/dominance/data/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
MEM=10G
mkdir -p $scratch/thaliana-temp

#######################
#  1. filter

java -jar ~/software/GenomeAnalysisTK.jar \
 -T VariantFiltration \
 -R ${REF} \
 -o /dev/stdout/ \
 --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0 || HaplotypeScore > 13.0" \
 --filterName "Recommended" \
 --variant ${working_dir}/vcf/thaliana.vcf | gzip -c - > ${working_dir}/vcf/thaliana.filtered.vcf.gz 





### old lyrata stuff. need to change for thaliana
#gzip -cd ${working_dir}/vcf/lyrata_ss.vcf.gz | bgzip -c > ${working_dir}/vcf/lyrata_ss.vcf.b.gz
#tabix -p vcf ${working_dir}/vcf/lyrata_ss.vcf.b.gz

# filter single sample calls
#java -jar ~/software/GenomeAnalysisTK.jar \
# -T VariantFiltration \
# -R ${REF} \
# -o /dev/stdout/ \
# --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0 || HaplotypeScore > 13.0" \
# --filterName "Recommended" \
# --variant ${working_dir}/vcf/lyrata_ss.vcf.b.gz | bgzip -c > ${working_dir}/vcf/lyrata_ss.filtered.vcf.b.gz