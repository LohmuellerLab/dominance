#!/bin/bash
#$ -cwd
#$ -V
#$ -N GOterm
#$ -l highp,h_data=24G,time=300:00:00
#$ -M eplau
#$ -m bea


. /u/local/Modules/default/init/modules.sh
module load bedtools

#python scripts/fill_in_mono.py -s results/lyrata_multisample_nostructure_sorted.txt -g data/reference/Arabidopsis_lyrata.v.1.0.33.chr.gff3.gz -t data/reference/ATH_GO_GOSLIM.txt.gz > results/lyrata_multisample_nostructure_sorted_GO.txt

python scripts/fill_in_mono2.py -s results/thaliana_multisample_nostructure_sorted.txt -g data/reference/Arabidopsis_thaliana.TAIR10.34.gff3.gz -t data/reference/ATH_GO_GOSLIM.txt.gz > results/thaliana_multisample_nostructure_sorted_GO.txt