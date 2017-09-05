#!/bin/bash
#$ -cwd
#$ -V
#$ -N summarize
#$ -l h_data=8G,time=10:00:00
#$ -M eplau
#$ -m bea

# old multisample
#python scripts/summarize_filtered_vcf_haploid.py -v data/vcf/thaliana.filtered.annotated.CDS.vcf.gz -o results/thaliana_multisample.txt
#sort -k1,1 -k2,2 -n results/thaliana_multisample.txt | uniq > results/thaliana_multisample_sorted.txt

#. /u/local/Modules/default/init/modules.sh
#module load bcftools

#bcftools view -s ^35601,35513,35600,37469 data/vcf/thaliana.filtered.annotated.CDS.vcf.gz | gzip  > data/vcf/thaliana.filtered.annotated.CDS.nostructure.vcf.gz

python scripts/summarize_filtered_vcf_haploid.py -v data/vcf/thaliana.filtered.annotated.CDS.nostructure.vcf.gz -o results/thaliana_multisample_nostructure.txt 
sort -k1,1 -k2,2 -n results/thaliana_multisample_nostructure.txt | uniq > results/thaliana_multisample_nostructure_sorted.txt
