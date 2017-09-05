#!/bin/bash
#$ -cwd
#$ -V
#$ -N merge
#$ -l h_data=8G,time=4:00:00
#$ -M eplau
#$ -m bea

. /u/local/Modules/default/init/modules.sh
module load bcftools

bcftools merge data/vcf/lyrata_single_sample/*.gz | gzip  > data/vcf/lyrata_ss.vcf.gz