echo "scripts/snp_eff_annotate.sh" | qsub -cwd -V -N annotate -l h_data=8G,time=33:00:00,highp -pe shared 2 -M eplau -m bea
