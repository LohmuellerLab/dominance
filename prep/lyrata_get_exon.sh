#!/bin/bash
#$ -cwd
#$ -V
#$ -N exon
#$ -l h_data=8G,time=4:00:00
#$ -M eplau
#$ -m bea

## load all modules here
set -e
set -u
set -o

. /u/local/Modules/default/init/modules.sh
module load java/1.8.0_77
module load bcftools

#zcat data/reference/Arabidopsis_lyrata.v.1.0.33.chr.gff3.gz | awk '$3 == "exon" {OFS="\t";print($1,$4,$5)}' > data/reference/lyrata_exon.gff
#zcat data/reference/Arabidopsis_lyrata.v.1.0.33.chr.gff3.gz | awk '$3 == "CDS" {OFS="\t";print($1,$4,$5)}' > data/reference/lyrata_CDS.gff
#java -Xmx4g -jar ~/software/snpEff/SnpSift.jar filter "(FILTER = 'PASS')" data/vcf/lyrata.filtered.annotated.vcf.gz| bgzip -c > data/vcf/lyrata.filtered.annotated.PASS.vcf.b.gz
#tabix -p vcf data/vcf/lyrata.filtered.annotated.PASS.vcf.b.gz

#bcftools view -R data/reference/lyrata_CDS.gff data/vcf/lyrata.filtered.annotated.PASS.vcf.b.gz | gzip -c - > data/vcf/lyrata.filtered.annotated.CDS.vcf.gz


## single sample
#java -Xmx4g -jar ~/software/snpEff/SnpSift.jar filter "(FILTER = 'PASS')" data/vcf/lyrata_ss.filtered.annotated.vcf.gz| bgzip -c > data/vcf/lyrata_ss.filtered.annotated.PASS.vcf.b.gz
#tabix -p vcf data/vcf/lyrata_ss.filtered.annotated.PASS.vcf.b.gz
bcftools view -R data/reference/lyrata_CDS.gff data/vcf/lyrata_ss.filtered.annotated.PASS.vcf.b.gz | gzip -c - > data/vcf/lyrata_ss.filtered.annotated.CDS.vcf.gz