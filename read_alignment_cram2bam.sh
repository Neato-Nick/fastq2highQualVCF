#!/bin/bash

#$ -cwd
#$ -S /bin/bash
#$ -N cram2bam
#$ -e align_callvars_err
#$ -o align_callvars_out
###$ -q *@!(nem*|megraw*|samwise*|fungi*|oryza*|cedro*)
###$ -q *@!(nem*|samwise*|fungi*|cedro*|symbiosis*|galls*|megraw*)
###$ -q *@!(nem*|samwise*|fungi*|cedro*|symbiosis*|galls*)
###$ -q *@!(nem*|samwise*|cedro*|fungi*)
###$ -q bpp@oryza0
###$ -q bpp@galls
###$ -q *@!(nem*|samwise*|cedro*)
###$ -q *@!(nem*|samwise*|amp*)
#$ -q *@!(amp*|symbiosis*|nem*)
#$ -l mem_free=15G
#$ -V
# #$ -h
#$ -t 501-669:1
#$ -tc 50
#$ -pe thread 1

# limit memory access, but vmem was causing all sorts of cryptic malloc errors in gatk
###$ -l h_vmem=10G
# 1 task per sample

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
# REFDIR="/nfs5/BPP/Grunwald_Lab/lab_collabs/diagnostic_assays/genome_data"
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
#sm_file="samples_cram2bam.20200828.txt"
#sm_file="samples_all.txt"
#sm_file="samples_fastq2cram.20201009.txt"
#sm_file="readsEU1-NA1_295.list" # 141+19=160
#sm_file="readsNA1_135.list" # 118+18-1=135
#sm_file="samples.SRA_558041_559872.no_overlap_samples-all.txt" # 82
#sm_file="samples_redo_20200503.txt"
#sm_file="KSondreli.samples.txt"
#sm_file="samples_re-fqdump.2022jan.txt"
#sm_file="samples_ubc_hybrids.txt"
sm_file="samples_globaldiv_hybrids_n669.txt"
echo
echo $sm_file
echo
FILE=( `cat "$sm_file"`)
IFS=';' read -a arr <<< "${FILE[$i]}"
echo

TEMPDIR="$PWD/align_callvars_tmp/"
# mkdir $TEMPDIR
ls $TEMPDIR*${arr[0]}*


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

# http://www.htslib.org/doc/
# Echo samtools version info.
echo "samtools info"
CMD="$SAMT --version"
echo
eval $CMD
echo

##### ##### ##### ##### #####
# Be a good custodian

# first we need to clean up previous working dirs. Jobs failed at various states dependent on node
# so let's just erase everything from a particular isolate, before starting over
# affected dirs are align_callvars_tmp/ (bam-writing) bams/ (results) stats/ (dupmetrics + stats) gvcfs/
echo
echo "First, clearing temp files or any bams from the previous run"
CMD="rm -f ${TEMPDIR}${arr[0]}*.bam*"
date
#echo $CMD
#eval $CMD
date


##### ##### ##### ##### #####
# Mark duplicates and sort

# `sort` writes to standard out, unless an output file is forced by -o
# `fixmate` and `markdup` write to a file (as a positional arg) unlesss forced to stdout by -
# `markdup` requires `fixmate` first, which requires ipnut sorted by name
# `markdup` requires coordinate-sorted alignments... annoying
# I think gatk requires coordinate-sorted bams so we'd need to sort again at the end
# gatk has iterated several times to markdup, but it says nothing about removing them,
# even though default behavior is only to mark... so for archival purposes, I will not remove duplicate reads
# https://qcb.ucla.edu/wp-content/uploads/sites/14/2016/03/GATKwr12-2-Marking_duplicates.pdf

#BAM="$TEMPDIR${arr[0]}_dupmrk.bam"
BAM="bams/${arr[0]}_dupmrk.bam"

echo
echo "Sort -name | fixmate | sort -coordinate | mark duplicates > file written for gatk"
ALN="crams/${arr[0]}.cram"
#ALN="$TEMPDIR${arr[0]}.sam"
CMD="samtools sort --reference $BREF -m 8g -T ${TEMPDIR}${arr[0]}srt1 --threads $NSLOTS -n $ALN | $SAMT \
        fixmate -m --threads $NSLOTS - - | $SAMT \
        sort -m 8g -T ${TEMPDIR}${arr[0]}srt2 --threads $NSLOTS - | $SAMT \
        markdup -T ${TEMPDIR}${arr[0]}mrkdup --threads $NSLOTS -f ./stats/${arr[0]}_dupmrk_metrics.txt - $BAM"
date
echo $CMD
eval $CMD
date
##### ##### ##### ##### #####
# Index

#CMD="$SAMT index bams/${arr[0]}_sorted.bam"
# BAM="bams/${arr[0]}_dupmrk.bam"
echo
echo "Creating index and calculating stats"
CMD="$SAMT index -@ $NSLOTS $BAM"
date
echo $CMD
eval $CMD
date


# Generate stats to validate alignment processing.
#CMD="$SAMT stats bams/${arr[0]}_sorted.bam | gzip -c > bams/${arr[0]}_sorted_stats.txt.gz"
# CMD="$SAMT stats $TEMPDIR${arr[0]}_dupmrk.bam | gzip -c > ./stats/${arr[0]}_bam-rmdup_stats.txt.gz"
echo
echo "Getting stats from the coord-sorted BAM with duplicates marked"
#CMD="$SAMT stats $TEMPDIR${arr[0]}_dupmrk.bam | gzip -c > $TEMPDIR${arr[0]}_dupmrk_stats.txt.gz"
CMD="$SAMT stats $BAM | gzip -c > ./stats/${arr[0]}_bam-dupmrk_stats.txt.gz"
date
echo $CMD
eval $CMD
date

echo "Finished getting BAM from CRAM, move onto variant calling"
# EOF.
