#!/bin/bash
#$ -cwd
#$ -V
#$ -N summarize
#$ -l h_data=8G,time=4:00:00
#$ -M eplau
#$ -m bea


# old multisample
#python scripts/summarize_filtered_vcf.py -v data/vcf/lyrata.filtered.annotated.CDS.vcf.gz -o results/lyrata_multisample.txt
#sort -k1,1 -k2,2 -n results/lyrata_multisample.txt | uniq > results/lyrata_multisample_sorted.txt

# try single sample
#python scripts/summarize_filtered_vcf.py -v data/vcf/lyrata_ss.filtered.annotated.CDS.vcf.gz -o results/lyrata_singlesample.txt
#sort -k1,1 -k2,2 -n results/lyrata_singlesample.txt | uniq > results/lyrata_singlesample_sorted.txt

. /u/local/Modules/default/init/modules.sh
module load bcftools

bcftools view -s ^SRR2040788,SRR2040795,SRR2040829 data/vcf/lyrata.filtered.annotated.CDS.vcf.gz | gzip  > data/vcf/lyrata.filtered.annotated.CDS.nostructure.vcf.gz

python scripts/summarize_filtered_vcf.py -v data/vcf/lyrata.filtered.annotated.CDS.nostructure.vcf.gz -o results/lyrata_multisample_nostructure.txt

sort -k1,1 -k2,2 -n results/lyrata_multisample_nostructure.txt | uniq > results/lyrata_multisample_nostructure_sorted.txt
