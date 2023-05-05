#!/usr/bin/env python3

import sys
import warnings
import math
import argparse
from collections import OrderedDict
from pysam import VariantFile

# Convert any allele to missing data
# My use case is converting asterisks to missing, the default behavior
# Mostly copied from OmitLowAAScoreAlleles.py
# Usage:
# <this_script.py> $in_vcf stdout | bcftools view --trim-alt-alleles -Ou | bcftools view -m2

if len(sys.argv) < 2 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
        sys.exit("Usage: OmitLowAAScoreAlleles.py <in_vcf> <out_vcf> <bad_allele> \nIf no <out_vcf>, or stdout or - specified, outputs to stdout. Default bad allele = *" )
elif len(sys.argv) == 2 or len(sys.argv) == 3:
        bad_alt = "*"
        if len(sys.argv) == 2 or sys.argv[2].lower() == "stdout" or sys.argv[2] == "-":
                bcf_out_f = "-"
        else:
                bcf_out_f = sys.argv[2]
elif len(sys.argv) > 4:
        sys.exit("Too many arguments. Run script with '-h' or no arguments for help")
else:
        bad_alt=sys.argv[3]

vcf_f=sys.argv[1]
bcfIn=VariantFile(vcf_f)
bcf_out=VariantFile(bcf_out_f, 'w', header = bcfIn.header)


# For every record,
# check if at least one alt allele is asterisk.
# If it is, then check every genotype
# Convert every matching allele to "."
for rec in bcfIn.fetch():
	for i, alt in enumerate(rec.alts):
		if alt == bad_alt:
			#print("Bad ALT allele found")
			bad_alt_allele = i+1
			for sample in rec.samples:
				GT_tuple=()
				for j, allele in enumerate(rec.samples[sample]["GT"]):
					# Create tuple with good alelles, or
					# add None value when hit a bad allele
					# then assign appended genotype tuple
					if allele == bad_alt_allele:
						GT_tuple += (None,)
						last_sample = sample
						#print("Bad allele encountered in %s" % (sample))
					else:
						GT_tuple += (allele,)
				rec.samples[sample]["GT"] = GT_tuple
	bcf_out.write(rec)
