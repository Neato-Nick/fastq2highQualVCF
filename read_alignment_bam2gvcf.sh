#!/bin/bash

#$ -cwd
#$ -S /bin/bash
#$ -N bam2gvcf_ksondreli
#$ -e align_callvars_err
#$ -o align_callvars_out
###$ -q *@!(nem*|megraw*|samwise*|fungi*|oryza*|cedro*)
###$ -q *@!(nem*|samwise*|fungi*|cedro*|symbiosis*|galls*|megraw*)
###$ -q *@!(nem*|samwise*|fungi*|cedro*|symbiosis*|galls*)
###$ -q *@!(nem*|samwise*|cedro*|symbiosis*|oryza1*|galls*|fungi0*|cerebro*)
###$ -q bpp@galls
#$ -q *@!(nem*|samwise*|cedro*|fungi*|galls*|amp*)
#$ -l mem_free=18G
#$ -V
# #$ -h
###$ -t 121-160
###$ -t 61-70
#$ -t 1-14:1
#$ -tc 50
#$ -pe thread 1

i=$(expr $SGE_TASK_ID - 1)


##### ##### ##### ##### #####
# Software

# https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows
#GATK="~/bin/gatk4/gatk-4.1.4.1/gatk"
# upgrade to newest release on 5/20/2020 to solve "smith-waterman alignment error" bugs (4.1.7.0)
# also bug fixes to GenotypeGVCFs (4.1.5.0) and to use port of DepthOfCoverage (4.1.6.0)
GATK="$HOME/dfs_opt/gatk-4.1.7.0/gatk"

JAVA="/home/bpp/knausb/bin/javadir/jre1.8.0_25/bin/java"
# when I run $gatk -version I get HTSJDK v2.21.2


##### ##### ##### ##### #####
# User provided materials

# Reference sequence
#REF="/home/bpp/knausb/Grunwald_Lab/home/knausb/pinf_bwa/bwaref/pinf_super_contigs.fa"
#
#REFDIR="/dfs/Grunwald_Lab/final_data/ramorum_pacbio_assemblies"
#BREF="${REFDIR}/PR-102.fasta.gz"
REFDIR="/dfs/Grunwald_Lab/lab_collabs/diagnostic_assays/genome_data"
BREF="${REFDIR}/PR-102_v3.1.fasta"

# GATK reference
GREF="$BREF"


# The file samples.txt contains info about sample names and files.
# Each line is one sample and one job.
# The line is a semi colon delimited list.
# The first element is the sample name.
# The second element is the fastq file including any path info.
#
# t30-4;../fastqs/ATCGGC.fastq.gz
#

#sm_file="samples_all.txt"
#sm_file="samples_bam2gvcf.20200830.txt"
#sm_file="samples_fastq2cram.20201009.txt"
#sm_file="readsEU1_160.list" # 141+19=160
#sm_file="readsNA1_135.list" # 118+18-1=135
#sm_file="readsEU1-NA1_295.list" # 160+135=295
#sm_file="samples.SRA_558041_559872.no_overlap_samples-all.txt" # 82
#sm_file="samples_redo_20200503.txt"
#sm_file="samples_re-fqdump.2022jan.txt"
sm_file="KSondreli.samples.txt"
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
# Report what we ended up with

echo -n "Running on: "
hostname
echo "SGE job id: $JOB_ID"
date
echo

myEpoch=(`date +%s`)
echo "Epoch start:" $myEpoch
startEpoch=$myEpoch

# Report files to be using, including refs and fastqs
echo "Reference for BWA:"
CMD="ls $BREF"
eval $CMD
echo "Reference for GATK:"
CMD="ls $GREF"
eval $CMD
echo
echo "Isolate name: ${arr[0]}"
echo "Sequence files:"
CMD="ls ${arr[1]} ${arr[2]}"
eval $CMD
echo
echo "Alignment file: "
BAM="bams/${arr[0]}_dupmrk.bam"
ls $BAM

##### ##### ##### ##### #####
# Create gvcf

# the java option limits how much memory is used,
# ensure "-Xmx${n}g" matches "#$ -l mem_free=${n}G" in the qsub params
# UPDATE 5/28: consider leaving unzipped, and compressing later
# UPDATE 6/2: moved parallel threads from NSLOTS to 1 b/c of how many sms I have
#CMD="$GATK --java-options \"-Djava.io.tmpdir=/data/ -Xmx4g\" HaplotypeCaller \
echo
echo "Calling variants!"
gvcf="$TEMPDIR${arr[0]}.g.vcf"
# CMD="$GATK HaplotypeCaller \
CMD="$GATK --java-options \"-Xmx12g\" HaplotypeCaller \
   -R $GREF \
   -O $gvcf \
   -I $BAM \
   --tmp-dir $TEMPDIR \
   --native-pair-hmm-threads 1 \
   -ERC GVCF"
date
echo $CMD
eval $CMD
date

### Validate gvcf
echo
echo "Validating new gvcf"
CMD="$GATK --java-options \"-Xmx8g\" ValidateVariants \
	-gvcf true -R $GREF -V $gvcf"
date
echo $CMD
eval $CMD
date

### Count variants in there
echo
echo "Counting variants stats, for statistics"
CMD="$GATK --java-options \"-Xmx8g\" CountVariants \
	-R $GREF -V $gvcf"
date
echo $CMD
eval $CMD
date

myEpoch=(`date +%s`)
echo "Epoch start:" $myEpoch

echo "Finished calling variants! If it worked, the following gvcf and its tabix index exist"
CMD="ls ${gvcf}*"
echo $CMD
eval $CMD
date

##### ##### ##### ##### #####
# Copy files and clean up.
mv_src_dest="rsync -avz --remove-source-files"

echo
echo "Cleaning temp files, moving to ./bams and ./gvcfs"
date


# CMD="$mv_src_dest $TEMPDIR${arr[0]}.g.vcf.gz ./gvcfs/"
# gvcf and its index
echo
CMD="$mv_src_dest ${gvcf}* ./gvcfs/"
date
echo $CMD
eval $CMD
date

echo
CMD="rm -f ${TEMPDIR}${arr[0]}*.sam*"
date
#echo $CMD
#eval $CMD
date

echo
echo "Finished cleaning up files"
echo
date
echo

# EOF.
