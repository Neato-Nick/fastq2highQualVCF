#!/bin/bash
#$ -V
#$ -N graphtyper_n669_allscafs
#$ -e graphtyper_err
#$ -o graphtyper_out
#$ -q *@!(nem*|cedro*|amp*|symbiosis*)
#$ -v TMPDIR=genotype_tmp
#$ -l mem_free=20G
#$ -S /bin/bash
#$ -cwd
#$ -t 1-24:1
#$ -tc 10
#$ -pe thread 12

# 1 array task per contig in ref genome
i=$(expr $SGE_TASK_ID - 1)

echo -n "Running on: "
hostname
echo "SGE job id: ${JOB_ID}:$SGE_TASK_ID"
date
echo

REFDIR="."
REF="${REFDIR}/PR-102_v4.fasta"
graphtyper --version

# cut -d ';' -f 1 readsEU1-NA1_295.list | sed 's/^/bams\//' | sed 's/$/_dupmrk.bam/' > bamsEU1-NA1_295.list
# bamlist="samples_globaldiv_n669.bams.txt"
bamlist="samples_globaldiv_hybrids_n669.bams.list"
outdir="vcfs/graphtyper/globaldiv_n669"

echo "Will use reference:$REF"

FAI=( `cut -f 1 "${REF}.fai" `)
read scaffold <<< "${FAI[$i]}"

echo
echo "Genotyping all bams in $bamlist"
CMD="graphtyper genotype $REF --sams=$bamlist --region=$scaffold --threads=$NSLOTS --output=$outdir"
date
echo $CMD
eval $CMD
date

echo
echo "Concatenating segmented results"
cd $outdir/$scaffold
CMD="graphtyper vcf_concatenate ./*.vcf.gz | bgzip -c > ${scaffold}.vcf.gz"
date
echo $CMD
eval $CMD
date

echo
echo "Indexing concatenated VCF"
CMD="$GATK IndexFeatureFile -I ${scaffold}.vcf.gz"
date
echo $CMD
eval $CMD
echo "exit $?"
echo
date


