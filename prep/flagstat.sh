for sample in `ls data/aligned/*`
do
    out=`echo ${sample} | cut -f 3 -d '/'`
    samtools flagstat ${sample} > data/flagstats/${out}.flagstat
done