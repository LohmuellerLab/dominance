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
sample=$1
REF=/u/home/a/arun/project-klohmuel/dominance/data/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
MEM=10G
mkdir -p $scratch/thaliana-temp

#####################
# 0. sort bam

#java -jar -Xmx${MEM} -jar ~/software/picard.jar RevertSam \
# I=${working_dir}/${sample}.bam \
# O=${scratch}/thaliana-temp/${sample}.ubam \
# SANITIZE=true \
# MAX_DISCARD_FRACTION=0.005 \
# ATTRIBUTE_TO_CLEAR=XT \
# ATTRIBUTE_TO_CLEAR=XN \
# ATTRIBUTE_TO_CLEAR=AS \
# ATTRIBUTE_TO_CLEAR=OC \
# ATTRIBUTE_TO_CLEAR=OP \
# TMP_DIR=${scratch}/thaliana-temp \
# 2> ./preprocessing_${sample}.err

######################
#  1. Mark adapters
#java -Xmx${MEM} -jar ~/software/picard.jar MarkIlluminaAdapters \
# I=${scratch}/thaliana-temp/${sample}.ubam \
# O=${scratch}/thaliana-temp/${sample}.ubam.adapters.bam \
# M=markadapters_${sample}_metric.txt \
# TMP_DIR=${scratch}/thaliana-temp \
# 2>> ./preprocessing_${sample}.err

#######################
# 2. align
# match parameters from Novikova et al 2016
#java -Xmx${MEM} -jar ~/software/picard.jar SamToFastq MAX_RECORDS_IN_RAM=750000 \
# I=${scratch}/thaliana-temp/${sample}.ubam.adapters.bam \
# FASTQ=/dev/stdout \
# CLIPPING_ATTRIBUTE=XT \
# CLIPPING_ACTION=2 \
# INTERLEAVE=true \
# NON_PF=true \
# TMP_DIR=${scratch}/thaliana-temp \
# 2>> ./preprocessing_${sample}.err | \
#bwa mem -U 15 -M -t 2 -p ${REF} /dev/stdin 2>> ./preprocessing_${sample}.err | \
#java -Xmx${MEM} -jar ~/software/picard.jar MergeBamAlignment MAX_RECORDS_IN_RAM=750000 \
# ALIGNED_BAM=/dev/stdin \
# UNMAPPED_BAM=${scratch}/thaliana-temp/${sample}.ubam \
# OUTPUT=${scratch}/thaliana-temp/${sample}.cleanedaligned.bam \
# R=${REF} CREATE_INDEX=true ADD_MATE_CIGAR=true \
# CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
# INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
# PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
# TMP_DIR=${scratch}/thaliana-temp 2>> ./preprocessing_${sample}.err

######################
# 3. Mark Duplicates

#java -Xmx${MEM} -jar ~/software/picard.jar MarkDuplicates INPUT=${scratch}/thaliana-temp/${sample}.cleanedaligned.bam OUTPUT=${scratch}/thaliana-temp/${sample}.cleanedaligned.markdup.bam M=${scratch}/thaliana-temp/${sample}.markdup.metrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=150000 CREATE_INDEX=true TMP_DIR=${scratch}/thaliana-temp 2>> ./preprocessing_${sample}.err

#mv ${scratch}/thaliana-temp/${sample}.cleanedaligned.markdup.bam ${working_dir}/aligned/

######################                                                                                                                   
# 4. Old school indel realignment                                                                                                        

samtools index ${working_dir}/aligned/${sample}.cleanedaligned.markdup.bam 

java -Xmx${MEM} -jar ~/software/GenomeAnalysisTK.jar \
 -T RealignerTargetCreator \
 -R ${REF} \
 -I ${working_dir}/aligned/${sample}.cleanedaligned.markdup.bam \
 -o ${scratch}/lyrata-temp/${sample}.intervals \
 2>> ./preprocessing_${sample}.err


java -Xmx${MEM} -jar ~/software/GenomeAnalysisTK.jar \
 -T IndelRealigner \
 -R ${REF} \
 -targetIntervals ${scratch}/lyrata-temp/${sample}.intervals \
 -I ${working_dir}/aligned/${sample}.cleanedaligned.markdup.bam \
 -o ${working_dir}/aligned/${sample}.cleanedaligned.markdup.indel.bam \
 2>> ./preprocessing_${sample}.err

