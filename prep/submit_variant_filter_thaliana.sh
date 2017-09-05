echo "scripts/variant_filter_thaliana.sh" | qsub -cwd -V -N variant_filter -l h_data=10G,time=33:00:00,highp -pe shared 2 -M eplau -m bea
