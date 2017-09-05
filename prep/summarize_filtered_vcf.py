# argument is a gzipped, filtered VCF of exons only (diploid)

import gzip
import sys
import csv
import argparse
import re
import numpy as np

parser = argparse.ArgumentParser(description="A script that summarizes a VCF.")
parser.add_argument("-v", "--vcf", action="store", required=True, help="Input VCF file.")
parser.add_argument("-o", "--out", action="store", required=True, help="Output filename")

args = parser.parse_args()

_file = args.vcf
out_name = args.out

re_syn = re.compile("syn*")
re_missense = re.compile("missense*")

with gzip.open(_file, 'rt') as vcf_in:
    vcf_in = csv.reader(vcf_in, delimiter="\t")
    tsv_out = csv.writer(open(out_name, 'w'), delimiter="\t", lineterminator="\n")
    for row in vcf_in:
        if any('##' in strings for strings in row):
            continue
        if any('#CHROM' in strings for strings in row):
            tsv_out.writerow(["chr", "pos", "frequency", "dac", "sample_size", "min_cov", "max_cov", "mean_cov", "function", "multiallelic", "gene_name"])
            continue

        # split up line into chunks
        chrom,pos=row[:2]
        alt = row[4]
        info = row[7]
        try:
            dp_field = row[8].split(":").index("DP")# need to get field corresponding to DP because it's not in the same place every time 
            gt_field = row[8].split(":").index("GT")
        except ValueError:
            print("Skipping:", chrom, " ", pos, ". DP field does not exist.")
            continue

        inds = row[9:]
        final_inds = []
        gts = []
        for field in inds:
            if "./." not in field:
                field_list = field.split(":")
                if int(field_list[dp_field]) > 10 and int(field_list[dp_field]) < 200:
                    final_inds.append(field)
                    gts.append(field_list[gt_field])
        if len(gts) == 0:
            continue # skip where we can't call anything

        g = [g.split("/") for g in gts]
        g2 = [item for sublist in g for item in sublist]
        dac = g2.count("1")

        cov = [int(c.split(":")[dp_field]) for c in final_inds]
        min_cov = np.min(cov)
        max_cov = np.max(cov)
        mean_cov = np.mean(cov)
        
        cov = np.array(cov)
        s = cov[cov > 10]
        sample_size = len(s[s < 200])*2
        
        freq = dac/sample_size
        if "," in alt:
            multiallelic = "True"
        else:
            multiallelic = "False"

        if alt == ".": # match to reference                                                                            
            tsv_out.writerow([chrom, pos, "0", "0", sample_size, min_cov, max_cov, mean_cov, "NA", multiallelic, "NA" ])
            continue

        ann_index = [idx for idx, s in enumerate(info.split(";")) if "ANN" in s][0]
        ann = info.split(";")[ann_index].split("|")[1] # 18 -> ann; 1 -> annotation

        if "syn" in ann:
            functional = "synonymous"
        elif "missense" in ann:
            functional = "nonsynonymous"
        else:
            functional = "NA"

        gene_name = info.split(";")[ann_index].split("|")[3] # 18 -> ann; 3 -> gene

        tsv_out.writerow([chrom, pos, freq, dac, sample_size, min_cov, max_cov, mean_cov, functional, multiallelic, gene_name])
