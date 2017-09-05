import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-s", action="store", required=True, help="Sites file")
parser.add_argument("-c", action="store", required=True, help="Connectivity file (3 columns)")

args = parser.parse_args()
sites = args.s
interact = args.c

i_dict = {} # key = gene; value = number of interactions
with open(interact, 'r') as e:
    for line in e:
        row = line.split(" ")
        if row[0] == "protein1":
            continue
        g1,g2,score = row[0].split(".")[1],row[1].split(".")[1],int(row[2].rstrip())
        if score > 700:
            try:
                i_dict[g1] += 1 # try to add 1 to the count
            except KeyError:
                i_dict[g1] = 1 # if the key doesn't exist, create the key
        else:
            continue

print("chr", "pos", "frequency", "dac", "sample_size", "min_cov", "max_cov", "mean_cov", "function", "multiallelic", "gene_name", "GO_term", "med_expr", "connections",sep="\t")

with open(sites, 'r') as s:
    for line in s:
        row = line.split("\t")
        if row[0] == "chr":
            continue
        gene=row[10]
        try:
            interactions = str(i_dict[gene])
        except KeyError:
            interactions = "NA"

        row.append(interactions)
        r2 = [s.rstrip() for s in row]
        print(*r2, sep="\t")
