## load all modules here
set -e
set -u
set -o

. /u/local/Modules/default/init/modules.sh
module load java/1.8.0_77



#java -Xmx4g -jar ~/software/snpEff/snpEff.jar Arabidopsis_lyrata data/vcf/lyrata.filtered.vcf.gz | gzip -c - > data/vcf/lyrata.filtered.annotated.vcf.gz
#java -Xmx4g -jar /u/home/a/arun/software/snpEff-2/snpEff.jar Arabidopsis_lyrata data/vcf/lyrata.filtered.vcf.gz | gzip -c - > data/vcf/lyrata.filtered.annotated.3.6.vcf.gz
#java -Xmx4g -jar ~/software/snpEff/snpEff.jar Arabidopsis_lyrata data/vcf/lyrata.filtered.annotated.nonsyn.vcf.gz -s data/vcf/nonsyn_summary.html | /dev/null
#java -Xmx4g -jar ~/software/snpEff/snpEff.jar Arabidopsis_lyrata data/vcf/lyrata.filtered.annotated.syn.vcf.gz -s data/vcf/syn_summary.html | /dev/null

## single sample
#java -Xmx4g -jar ~/software/snpEff/snpEff.jar Arabidopsis_lyrata data/vcf/lyrata_ss.filtered.vcf.b.gz | gzip -c - > data/vcf/lyrata_ss.filtered.annotated.vcf.gz 


#thaliana
java -Xmx4g -jar ~/software/snpEff/snpEff.jar Arabidopsis_thaliana data/vcf/thaliana.filtered.vcf.gz | gzip -c - > data/vcf/thaliana.filtered.annotated.vcf.gz