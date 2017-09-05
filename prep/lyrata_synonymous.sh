## load all modules here
set -e
set -u
set -o

. /u/local/Modules/default/init/modules.sh
module load java/1.8.0_77

java -Xmx4g -jar ~/software/snpEff/SnpSift.jar filter " (FILTER = 'PASS') & (ANN[*].EFFECT has 'missense_variant')" data/vcf/lyrata.filtered.annotated.vcf.gz | gzip -c - > data/vcf/lyrata.filtered.annotated.nonsyn.vcf.gz
java -Xmx4g -jar ~/software/snpEff/SnpSift.jar filter " (FILTER = 'PASS') & (ANN[*].EFFECT has 'synonymous_variant')" data/vcf/lyrata.filtered.annotated.vcf.gz | gzip -c - > data/vcf/lyrata.filtered.annotated.syn.vcf.gz
