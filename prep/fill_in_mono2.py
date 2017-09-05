# script to fill in monomorphic sites in my frequency spectrum file
# use the GFF
# usage python fill_in_mono.py -f spectrum.txt -g gff -o  spectrum_monofilled.txt

import argparse
import csv
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("-s", action="store", required=True, help="Sites file")
parser.add_argument("-g", action="store", required=True, help="GFF")
parser.add_argument("-t", action="store", required=True, help="GO term list")

args = parser.parse_args()
sites = args.s
gff = args.g 
go_terms = args.t

g_dict = {}
g_list = []
# store GFF intervals
with gzip.open(gff, 'rt') as g:
    for line in g:
        if "#" in line:
            continue
        row = line.split("\t")
        if row[2] == "gene":
            #save interval
            try:
                chrom = int(row[0])
            except:
                continue # skip non integer chromosomes
            start = int(row[3])
            end = int(row[4])
            try:
                string = row[8].split(";")[0].split(":")[-1].split("_")[-1].split(".")[0] # attempts to grab the AT1 name, will grab something else if it fails
            except:
                string = "NA"
            pos_dict = {}
            for x in range(start, end+1):
                try:
                    g_dict[chrom][x] = string
                except KeyError:
                    g_dict[chrom] = {}
                    g_dict[chrom][x] = string


go_dict = {} # key = AT; value = GO
with gzip.open(go_terms, 'rt') as go:
    for line in go:
        row = line.split("\t")
        go_dict[row[0]] = row[5]

print("chr", "pos", "frequency", "dac", "sample_size", "min_cov", "max_cov", "mean_cov", "function", "multiallelic", "gene_name", "GO_term", sep="\t")
with open(sites, 'r') as s:
    for line in s:
        row = line.split("\t")
        if row[0] == "chr":
            continue
        chrom,pos = row[0],row[1]
 #       bed = pybedtools.BedTool([[chrom, pos, pos]])
 #       intersect = g_bed.intersect(bed)
 #       gene = [i[-1] for i in intersect][0]
        try:
            gene = g_dict[int(chrom)][int(pos)]
        except:
            gene = "NA"
#        if row[-1] == "NA\n":
        row[-1] = gene # replace last row with gene name regardless of if there used to be an NA
        try:
            row.append(go_dict[gene])
        except:
            row.append("NA")
        r2 = [s.rstrip() for s in row]
        print(*r2, sep="\t")
