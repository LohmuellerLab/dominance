#need bcftools loaded


for file in `ls data/vcf/lyrata_single_sample/*`
do
    #echo ${file}
    bgzip -d < "$file" | head -n -4 | bgzip -c >"$file".tmp
    mv -- "$file".tmp "$file"
done