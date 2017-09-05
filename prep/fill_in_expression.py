import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-s", action="store", required=True, help="Sites file")
parser.add_argument("-e", action="store", required=True, help="Expression file")

args = parser.parse_args()
sites = args.s
exp = args.e

exp_dict = {} # key = gene; value = media expression value
with open(exp, 'r') as e:
    for line in e:
        row = line.split(" ")
        exp_dict[row[0]] = row[1].rstrip()

print("chr", "pos", "frequency", "dac", "sample_size", "min_cov", "max_cov", "mean_cov", "function", "multiallelic", "gene_name", "GO_term", "med_expr", sep="\t")

with open(sites, 'r') as s:
    for line in s:
        row = line.split("\t")
        if row[0] == "chr":
            continue
        gene=row[10]
        try:
            med_exp = exp_dict[gene]
        except KeyError:
            med_exp = "NA"

        row.append(med_exp)
        r2 = [s.rstrip() for s in row]
        print(*r2, sep="\t")
