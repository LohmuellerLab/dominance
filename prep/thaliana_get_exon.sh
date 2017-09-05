#!/bin/bash
#$ -cwd
#$ -V
#$ -N exon
#$ -l h_data=8G,time=2:00:00
#$ -M eplau
#$ -m bea

## load all modules here
set -e
set -u
set -o

. /u/local/Modules/default/init/modules.sh
module load java/1.8.0_77
module load bcftools

zcat data/reference/Arabidopsis_thaliana.TAIR10.34.gff3.gz | awk '$3 == "CDS" {OFS="\t";print($1,$4,$5)}' > data/reference/thaliana_CDS.gff
java -Xmx4g -jar ~/software/snpEff/SnpSift.jar filter "(FILTER = 'PASS')" data/vcf/thaliana.filtered.annotated.vcf.gz| bgzip -c > data/vcf/thaliana.filtered.annotated.PASS.vcf.b.gz
tabix -p vcf data/vcf/thaliana.filtered.annotated.PASS.vcf.b.gz

bcftools view -R data/reference/thaliana_CDS.gff data/vcf/thaliana.filtered.annotated.PASS.vcf.b.gz | gzip -c - > data/vcf/thaliana.filtered.annotated.CDS.vcf.gz
