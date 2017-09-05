for SAMPLE in `ls data/aligned/SRR*.indel.bam | cut -f 3 -d "/" |cut -f 1 -d "."`
do
#    echo $SAMPLE
    echo "scripts/call_snps_lyrata_single.sh ${SAMPLE}" | qsub -cwd -V -N align${SAMPLE} -l h_data=30G,time=330:00:00,highp -pe shared 2 -M eplau -m bea

done