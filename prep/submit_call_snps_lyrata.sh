echo "scripts/call_snps_lyrata.sh" | qsub -cwd -V -N align${SAMPLE} -l h_data=30G,time=330:00:00,highp -pe shared 2 -M eplau -m bea
