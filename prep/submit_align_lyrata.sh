for SAMPLE in `ls data/SRR*_1.fastq.gz | cut -f 2 -d "/" |cut -f 1 -d "_"`
do
    echo "scripts/align_lyrata.sh ${SAMPLE}" | qsub -cwd -V -N align${SAMPLE} -l h_data=30G,time=330:00:00,highp -pe shared 2 -M eplau -m bea

done