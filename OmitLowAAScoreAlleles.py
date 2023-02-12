#!/usr/bin/env python3

import sys
import warnings
import math
import argparse
from collections import OrderedDict
from pysam import VariantFile

# Masks alternate alleles falling below AAScore threhsold
# Assumes you've already removed all variants where there is one alt allele and it failed the threshold
# this saves up to hours of time by iterating through all samples only when absolutely necessary
# INFO wil lhave lots of multiple values now e.g. 'AC=' for that allele,
# pipe output using "-" as second arg to this script into bcftools to remove the freshly unused alt alleles and remove any sites that are all Ref:
# <this_script.py> $in_vcf stdout | bcftools view --trim-alt-alleles -Ou | bcftools view -m2 

if len(sys.argv) < 2 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
	sys.exit("Usage: OmitLowAAScoreAlleles.py <in_vcf> <out_vcf> <min_AA> \nIf no <out_vcf>, or stdout or - specified, outputs to stdout. Default minimum AAScore = 0.5" )
elif len(sys.argv) == 2 or len(sys.argv) == 3:
	AAthresh=0.5
	if len(sys.argv) == 2 or sys.argv[2].lower() == "stdout" or sys.argv[2] == "-":
		bcf_out_f = "-"
	else:
		bcf_out_f = sys.argv[2]
elif len(sys.argv) > 4:
	sys.exit("Too many arguments. Run script with '-h' or no arguments for help")
else:
	AAthresh=sys.argv[3]
vcf_f=sys.argv[1]
bcfIn=VariantFile(vcf_f)
bcf_out=VariantFile(bcf_out_f, 'w', header = bcfIn.header)


#for every variant,
#check if at least two aascores.
#If at least two, check every score[i]
#if score[i] < threshold
#for every gt, convert every matching allele to "."

for rec in bcfIn.fetch():
	if len(rec.info["AAScore"]) > 1:
		for i, aa in enumerate(rec.info["AAScore"]):
			if aa < AAthresh:
				bad_alt_allele = i+1
				#print("allele %d of %s failed: %s" % (i,rec.alts[i], aa))
				for sample in rec.samples:
					GT_tuple=()
					for j, allele in enumerate(rec.samples[sample]["GT"]):
						# Create a tuple with the good alleles, or
						# add the None value when hit a bad allele,
						# Then assign the appended genotype tuple.
						if allele == bad_alt_allele:
							GT_tuple += (None,)
							last_sample = sample
						else:
							GT_tuple += (allele,)
					rec.samples[sample]["GT"] = GT_tuple
			#else:
			#	print("allele %d passed with score: %d" % (i, aa))
	bcf_out.write(rec)

