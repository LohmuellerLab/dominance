#!/bin/bash
#$ -cwd
#$ -V
#$ -N expr
#$ -l highp,h_data=4G,time=3:00:00
#$ -M eplau
#$ -m bea


. /u/local/Modules/default/init/modules.sh

python scripts/fill_in_expression.py -s results/lyrata_multisample_nostructure_sorted_GO.txt -e data/median_expression.tsv > results/lyrata_multisample_nostructure_sorted_GO_expr.txt
