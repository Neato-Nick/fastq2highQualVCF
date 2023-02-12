#!/bin/bash

#$ -cwd
#$ -S /bin/bash
#$ -N fastq2cram_bwa
#$ -e align_callvars_err
#$ -o align_callvars_out
###$ -q *@!(nem*|megraw*|samwise*|fungi*|oryza*|cedro*)
###$ -q *@!(nem*|samwise*|fungi*|cedro*|symbiosis*|galls*)
###$ -q *@!(nem*|samwise*|cedro*|fungi*|amp*)
###$ -q *@!(nem*|samwise*|cedro*|amp*|symbiosis*)
###$ -q *@(nem|samwise)
#$ -q bpp
#$ -l mem_free=10G
#$ -V
#$ -t 1-669:1
#$ -tc 50
#$ -pe thread 4

# limit memory access, but vmem was causing all sorts of cryptic malloc errors in gatk
###$ -l h_vmem=10G
# 1 array task per sample

i=$(expr $SGE_TASK_ID - 1)

##### ##### ##### ##### #####
# Software

# http://bio-bwa.sourceforge.net/bwa.shtml
#BWA="~/bin/bwa-0.7.17/bwa"
# recompiled to fit my own file structure
BWA="$HOME/dfs_opt/bwa-0.7.17/bwa"

# SAMT="~/bin/samtools-1.9/samtools-1.9/samtools "
# upgrade to newest release on 5/20/2020 to utilize the new "coverage" sub-command
# this should help make sure we get the same number of total reads -> reads aligned+unaligned
# and can plot histogram across genome
# samtools stats also gets some nice upgrades
# lastly the misc/plot-bamstats script can be sued to check quality heatmaps for each bam
# there should be no conflicts with this - unless samtools changes allowable bam header chars
#SAMT="$HOME/dfs_opt/samtools-1.10/samtools"
SAMT="/local/cluster/samtools-1.16.1/bin/samtools"

##### ##### ##### ##### #####
# User provided materials

# Reference sequence
#REFDIR="/nfs5/BPP/Grunwald_Lab/lab_collabs/diagnostic_assays/genome_data"
REFDIR="."
BREF="${REFDIR}/PR-102_v4.fasta"

# The file samples.txt contains info about sample names and files.
# Each line is one sample and one job.
# The line is a semi colon delimited list.
# The first element is the sample name.
# The second element is the fastq file including any path info.
#
# t30-4;../fastqs/ATCGGC.fastq.gz
#
# cat ../illumina_mapping/samples_globaldiv_n667.txt ../illumina_mapping/samples_ubc_hybrids.txt | egrep -v "^#" | sed 's/$HOME/\/home\/bpp\/carleson/g' | sed 's/grunwald_lab\/raw/grunwald_lab_dfs\/raw/g' > samples_globaldiv_hybrids_n669.txt
#while read -r sample; do name=$(echo $sample | cut -d ';' -f 1); fwd=$(echo $sample | cut -d ';' -f 2); rev=$(echo $sample | cut -d ';' -f 3); echo $name; ls $fwd; ls $rev;done < samples_globaldiv_hybrids_n669.txt > find.out 2>&2
sm_file="samples_globaldiv_hybrids_n669.txt"

echo
echo $sm_file
echo
FILE=( `cat "$sm_file"`)
IFS=';' read -a arr <<< "${FILE[$i]}"
echo


##### ##### ##### ##### #####
# Report host info where we landed

echo -n "Running on: "
hostname
echo "SGE job id: ${JOB_ID}:$SGE_TASK_ID"
date
echo

myEpoch=(`date +%s`)
echo "Epoch start:" $myEpoch
startEpoch=$myEpoch

# Report files to be using, including refs and fastqs
echo "Reference for BWA:"
CMD="ls $BREF"
eval $CMD
echo
echo "Isolate name: ${arr[0]}"
echo "Sequence files:"
CMD="ls ${arr[1]} ${arr[2]}"
eval $CMD
echo

# http://bio-bwa.sourceforge.net/bwa.shtml
# Align reads with bwa.

# Report bwa version info.
echo "bwa info"
CMD="$BWA 2>&1"
echo
echo $CMD
eval $CMD
echo

# http://www.htslib.org/doc/
# Echo samtools version info.
echo "samtools info"
CMD="$SAMT --version"
echo
eval $CMD
echo

##### ##### ##### ##### #####
# Be a good custodian

# first we might need to clean up previous working dirs. Jobs failed at various states dependent on node
# so let's just erase everything from a particular isolate, before starting over
# affected dirs are align_callvars_tmp/ (bam-writing) bams/ (results) stats/ (dupmetrics + stats) gvcfs/

TEMPDIR="$PWD/align_callvars_tmp/"
# mkdir $TEMPDIR
echo "Files in temp dir related to this isolate:"
ls ${TEMPDIR}*${arr[0]}*

echo
echo "First, clearing temp files or any bams from the run with interrupted fq dump"
CMD="rm -f ${TEMPDIR}${arr[0]}*.sam* ${TEMPDIR}${arr[0]}*.cram*"
date
#echo $CMD
#eval $CMD
date


##### ##### ##### ##### #####
# Map reads

# The GATK needs read group info:
# https://software.broadinstitute.org/gatk/guide/article?id=6472
# SM: sample
# LB: library, may be sequenced multiple times
# ID: Read Group Identifier, a unique identifier
# PL: Platform/technology used
echo
echo "Aligning reads"
RG="@RG\tID:${arr[0]}\tLB:${arr[0]}\tPL:illumina\tSM:${arr[0]}\tPU:${arr[0]}"

CMD="$BWA mem -t $NSLOTS -M -R \"$RG\" $BREF ${arr[1]} ${arr[2]} > $TEMPDIR${arr[0]}.sam"
date
echo $CMD
eval $CMD
date

##### ##### ##### ##### #####
# Generate stats to validate the alignment went okay
# CMD="$SAMT stats $TEMPDIR${arr[0]}.sam | gzip -c > $TEMPDIR${arr[0]}_stats.txt.gz"
echo
echo "How did we do?"
CMD="$SAMT stats $TEMPDIR${arr[0]}.sam | gzip -c > ./stats/${arr[0]}_sam-raw_stats.txt.gz"
date
echo $CMD
eval $CMD
date

# convert sam to cram (to save most space) and move to crams dir
echo
CMD="$SAMT view -hC -@ $NSLOTS -T $BREF $TEMPDIR${arr[0]}.sam -o ./crams/${arr[0]}.cram"
date
echo $CMD
eval $CMD
date
echo

echo "Leaving SAM in $TEMPDIR until valid bam is created"
echo "This means, after gvcf script!"
echo "So, we are all done here!"
date
echo

# EOF
