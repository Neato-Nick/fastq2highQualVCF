#!/usr/bin/env python3

# Script originally obtained from Upadhyay et al. 2021, modified to
# 1) to allow missing DP values
# 2) add parameter to optionally retain multiallelic variants
# Maulik Upadhyay sent Nick Cauldron this script 2022-02-23
# https://github.com/BioInf2305/trivial-python-scripts/blob/master/ReplaceGenoWithMissing.py

import sys
import warnings
import math
import argparse
from collections import OrderedDict
from pysam import VariantFile

def SetVarDepthMissing(bcfIn,sampleDepthFile,mapQ,nDepth,mode,multiall,bcfOut):
    nDepth=int(nDepth)
    mapQ=int(mapQ)
    if mode not in ["1","2"]:
        print("ERROR, mode should be set as either 1 or 2, refer to the help option")
        sys.exit(0)
    bcf_in=VariantFile(bcfIn)
    totalSamples=len(list((bcf_in.header.samples)))
    bcf_in.header.info.add("SD",str(totalSamples),"String","list of samples filter based on its respective depth")
    bcf_out = VariantFile(bcfOut, 'w', header=bcf_in.header)
    sampleDepth=OrderedDict()
    strNull="0"
    intNull=0
    with open(sampleDepthFile) as infile:
        for line in infile:
            line=line.split()
            sampleDepth[line[0]]=int(math.ceil(float(line[1])))
    for sample in list(bcf_in.header.samples):
        if sample not in sampleDepth:
            print("ERROR: "+sample+" is not present in file "+str(sampleDepthFile))
            sys.exit(0)
    for rec in bcf_in.fetch():
        variableDepth=0
        varCount=0
        sampleTuple=()
        if rec.info["MQ"]>mapQ:
            for sample in sampleDepth:
                # First condition allows sample dp file to have samples not in vcf
                if sample not in rec.samples.keys():
                    continue
	        # This if condition converts missing DP values (".") to 0
                if (type(rec.samples[sample]["DP"]) is not int):
                    rec.samples[sample]["DP"]=0
                    warnings.warn("Missing DP in VCF converted to 0")
                if (rec.samples[sample]["DP"]>(sampleDepth[sample]*nDepth)) or \
                (rec.samples[sample]["DP"]<(sampleDepth[sample]/nDepth)):
                    variableDepth+=1
                    rec.samples[sample]["GT"]=(None,None)
                    sampleTuple+=(sample,)
                    if mode!="1":
                        for key in rec.samples[sample]:
                            if key!="GT":
                                tupleEmpty=()
                                if type(rec.samples[sample][key])==type(intNull):
                                    rec.samples[sample][key]=intNull
                                else:
                                    for i in range(len(rec.samples[sample][key])):
                                        if type(rec.samples[sample][key][0])==type(strNull):
                                            addNull=strNull
                                        else:
                                            addNull=intNull
                                        tupleEmpty+=(addNull,)
                                    rec.samples[sample][key]=tupleEmpty
                elif rec.samples[sample]["GT"]!=(None,None) and rec.samples[sample]["GT"]!=(0,0):
                        varCount+=1
            if variableDepth!=len(sampleDepth) and varCount>0:
                if rec.alts!=None:
                    if len(rec.alts)>1 and multiall==1:
                        pass
                    else:
                        for i in range(totalSamples-len(sampleTuple)):
                            sampleTuple+=("",)
                        rec.info["SD"]=sampleTuple
                        bcf_out.write(rec)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="A small python script to set the geno field missing for \
    each individual sample based on user-defined threshold of depth. Positions with minor allele freq 0 as well as with \
    multi-allelic variants will be excluded by default. Additionally, the script will add the info tag SD listing the samples for \
    which the genotypes were set to missing at that particular position", epilog="author: Maulik Upadhyay (Upadhyay.maulik@gmail.com)")
    parser.add_argument("-v","--vcfF",metavar="File",help="compressed or uncompressed vcf or \
    bcf file",required=True)
    parser.add_argument("-s",'--sampleF',metavar="File",help="sample file with two columns:\
    sample_name average_depth; note that here sample names should match with those defined in \
    the bcf/vcf file (refer to the line starting with #CHROM in the file)",required=True)
    parser.add_argument("-q","--mapQ",metavar="Int",help="minimum mapping quality of SNPs to consider \
    (others will be filtered out,default=20). OPTIONAL FLAG",default=20,required=False)
    parser.add_argument("-n","--nDepth",metavar="Int",help="genotypes with depth below average_depth/n or \
    higher than average_depth*n for each resp. sample will be set to missing",required=True)
    parser.add_argument("-m","--mode",metavar="Str",help=" 1 or 2;whether to set only genotypes as missing or the entire format entries as missing \
    for example in mode 1, GT:PL:DP:SP:ADF:ADR:AD will be converted to ./.:0,120,227:40:0:27,0:13,0:40,0, while in mode 2, it will be set to \
    ./.:0,0,0:0:0:0,0:0,0:0,0. OPTIONAL FLAG, Default is 1",default="1",required=False)
    parser.add_argument("-a","--multiallelic",metavar="Str",help=" 1 or 2;whether to remove multi-allelic variants or not. \
    mode 1 removes these, mode 2 keeps them. OPTIONAL FLAG, Default is 1",default="1",required=False)
    parser.add_argument("-o","--outputF",metavar="Str",help="bcf or vcf file name for the output", required=True)
    args = parser.parse_args()
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        SetVarDepthMissing(args.vcfF,args.sampleF,args.mapQ,args.nDepth,args.mode,args.multiallelic,args.outputF)
