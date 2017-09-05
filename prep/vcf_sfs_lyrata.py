#/bin/python2
# vcf_sfs.py
# checks for heterozygotes

import argparse
import gzip
import csv

parser = argparse.ArgumentParser(description="A script that makes a site frequency spectrum.")
parser.add_argument("-v","--vcf", action="store", required=True, help="Input VCF. Assumes that it is gzipped")
parser.add_argument("-i", "--interesting", action="store_true", required=False, help="Instead of SFS, output high frequency loci information. Argument is frequency cutoff for interesting SNPs in range: (0, 1)")
parser.add_argument("-b", "--bed", action="store_true", required=False, help="output bed format (should be used with -i)")

args = parser.parse_args()

vcf_in = args.vcf
interesting = args.interesting
bed = args.bed

SFS = []
high_freq_loci = []

THRESHOLD=0.7

with gzip.open(vcf_in, 'r') as tsvin:
    tsvin = csv.reader(tsvin, delimiter="\t")
    for row in tsvin:
        if any('##' in strings for strings in row):
            continue
        if any('#CHROM' in strings for strings in row):
            nsamples = len(row[9:])
            SFS = [0] * nsamples # initialize SFS with zeros
            continue

        chrom,pos,_id,ref,alt,qual,_filter,info,_format=row[0:9]
        
        if len(alt) > 1 or len(ref) > 1 or "," in alt: #skip non biallelic and indels
            continue
            
        gts1 = [geno.split(":")[0] for geno in row[9:]]
        
        gst2 = [geno.split("/")[:] for geno in gts1]
        
        print gts1
        print gts2


#        if "." in gts: # skip missing data
#            continue
#            
#        freq = gts.count("1|1")
#        if interesting:
#            if freq*1.0/len(gts) >= THRESHOLD:
#                high_freq_loci.append([chrom,pos,ref,alt, freq*1.0/len(gts)])
#        try:
#            SFS[freq] = SFS[freq] + 1 #add 1 to the frequency column
#        except IndexError:
#            continue

#if interesting:
#    if bed:
#        for locus in high_freq_loci:
#            print locus[0]+"\t"+str(int(locus[1])-10)+"\t"+str(int(locus[1])+10)
#    else:
#        for locus in high_freq_loci:
#            print locus[0]+":"+locus[1]+"\t"+str(locus[4])
#else:            
#    print SFS
