echo "scripts/call_snps_thaliana.sh" | qsub -cwd -V -N align${SAMPLE} -l h_data=10G,time=48:00:00,highp -pe shared 2 -M eplau -m bea
