#!/bin/bash

#$ -V
#$ -N upad_vcftools_multiallelic_lowAA
#$ -e filter_err
#$ -o filter_out
###$ -q *@!(nem*|samwise*|amp*|debary*|galls*|anduin*|symbiosis*)
#$ -q samwise
#$ -l mem_free=20G
#$ -S /bin/bash
#$ -cwd
#$ -t 1-24:1
#$ -tc 15
#$ -pe thread 8

# Filter each scaffold and phase

i=$(expr $SGE_TASK_ID - 1)

echo -n "Running on: "
hostname
echo "SGE job id: ${JOB_ID}:$SGE_TASK_ID"
date
echo

REFDIR="/nfs5/BPP/Grunwald_Lab/home/carleson/ramorum/illumina_mapping_pr102v4"
REF="${REFDIR}/PR-102_v4.fasta"


FAI=( `cut -f 1 "${REF}.fai" `)
read scaffold <<< "${FAI[$i]}"

GATK="$HOME/dfs_opt/gatk-4.1.7.0/gatk"

# Files set up
workdir=$scaffold
big_vcf_name=$scaffold
big_vcf=$workdir/${big_vcf_name}.vcf.gz

# Ideally, before phasing we would remove any alternate allele with an AAScore < 0.5
# As previously implemented, these would be included as long as at least one alt allele *did* meet that threshold
# Sites with mult alleles were removed because the Upadhyay et al. script removes any multi-allelic site!
# Unfavorable - ideally remove AAInfo alts BEFOREhand, which would rescue more SNPs.
# 500 SNPs were removed from pchr7 that shouldn't have been.
# A minor difference but means I'm throwing out real diversity
# Fixed by making the script not discard multiallelic sites,
# and omitting low AAScores so I can filter multi-allelics out later if I want
# bcftools view --threads 4 --trim-alt-alleles -Ou pchr7.vcf.gz | \
# bcftools view --threads 4 -f .,PASS \
#--include '(INFO/AAScore[0] >= 0.5 & INFO/AAScore[1] < 0.5) | (INFO/AAScore[1] >= 0.5 & INFO/AAScore[0] < 0.5)' \
# -o pchr7.splitAA.pass.vcf.gz


# Type of average used and Minimum mapping qual
avg="median"
mapQ=0
avg_depths="mean-dp_globaldiv_n669.mosdepth_mapQ${mapQ}_${avg}"
avg_depths_f="${REFDIR}/stats/${avg_depths}.txt"
# initial awk-ing of min_dp_i only neds to be performed once
min_dp_i=15
# Include List:
CMD="awk -v min_dp=$min_dp_i '\$2 >= min_dp' $avg_depths_f | cut -f 1 > ${avg_depths}.minDP${min_dp_i}.include.list"
CMD="awk -v min_dp=$min_dp_i '\$2 <= min_dp' $avg_depths_f | cut -f 1 > ${avg_depths}.minDP${min_dp_i}.exclude.list"
echo $CMD
eval $CMD
# Exclude list:
CMD="awk -v min_dp=$min_dp_i '\$2 <= min_dp' $avg_depths_f | cut -f 1 > ${avg_depths}.minDP${min_dp_i}.exclude.list"
echo $CMD
eval $CMD

initial_filt_vcf="$workdir/""$(basename --suffix ".vcf.gz" ${big_vcf})"".bcft_Idp${min_dp_i}_aa05_PASS.vcf.gz"
CMD="bcftools view --threads 4 -Ou \
	-S ${avg_depths}.minDP${min_dp_i}.include.list \
	$big_vcf | \
	bcftools view --threads 4 --trim-alt-alleles \
	-f .,PASS --include '(INFO/AAScore[*] >= 0.5)' \
	-o $initial_filt_vcf"
TABIX="tabix $initial_filt_vcf"
echo $CMD
eval $CMD
echo $TABIX
eval $TABIX

# Censor based on missing data quantile
# x/3 or 3x
gt_dp=3
vcf_upad_out="$workdir/""$(basename --suffix ".vcf.gz" ${initial_filt_vcf})"".upad_MQ${mapQ}_${avg}_${gt_dp}x.vcf.gz"
CMD="python3 ~/dfs_opt/scripts/ReplaceGenoWithMissing_nc.py -v $initial_filt_vcf \
	-s $avg_depths_f -q $mapQ -n $gt_dp -o $vcf_upad_out -m 1 -a 2"
echo $CMD
eval $CMD
date
echo

# Censor based on hard cut-off of minimum cov, and remove samples averaging below X coverage
# target is to remove 7881-E12 but there may be others we want to get rid of
# Assign that ($min_dp_i) based on awk of whole-genome stats so every chr has the same indivs
min_dp_v=4
# Filter on minQ, minDP per genotype, remove samples
remove_sms_out="$workdir/""$(basename --suffix ".vcf.gz" ${vcf_upad_out})"".vcft_Vdp${min_dp_v}_PASS"
CMD="vcftools --gzvcf $vcf_upad_out --out $remove_sms_out --recode --recode-INFO-all \
	--minDP $min_dp --minQ $min_Q --remove-filtered-all"
echo $CMD
eval $CMD
CMD="bgzip -f ${remove_sms_out}.recode.vcf"
echo $CMD
eval $CMD
CMD="tabix -f ${remove_sms_out}.recode.vcf.gz"
echo $CMD
eval $CMD

# Remove any variants where number of Homref+MissingGTs=number of samples,
# possible for any samples we removed containing singletons
filter_aa_out="$workdir/""$(basename --suffix ".vcf.gz" ${remove_sms_out})"".bcft_RRmis_lowAA"
#CMD="bcftools view --trim-alt-alleles --threads $NSLOTS  -o ${filter_aa_out}.vcf \
CMD="bcftools view --trim-alt-alleles --threads $NSLOTS  -Ou \
	${remove_sms_out}.recode.vcf.gz | \
	bcftools view --threads $NSLOTS --trim-alt-alleles -o ${filter_aa_out}.vcf
	--include '(COUNT(GT=\"RR\")+(COUNT(GT=\"mis\"))) < N_SAMPLES & INFO/AAScore[*] >= 0.5'"
echo $CMD
eval $CMD
date
echo
CMD="bgzip -f ${filter_aa_out}.vcf"
echo $CMD
eval $CMD
CMD="tabix -f ${filter_aa_out}.vcf.gz"
echo $CMD
eval $CMD

# Remove alt alleles that don't meet AAScore, trim unused alts, and remove variants with no ALT
# Rarely, this duplicated some variant IDs. Of those, take only the first record.
omit_lowAAScore_out="${filter_aa_out}.omit_lowAAS"
CMD="$HOME/dfs_opt/scripts/OmitLowAAScoreAlleles.py ${filter_aa_out}.vcf.gz - | \
	bcftools view --threads $NSLOTS --trim-alt-alleles -Ou | \
	bcftools view --threads $NSLOTS -m2 -Ou |
	bcftools norm -d exact -o ${omit_lowAAScore_out}.vcf"
echo $CMD
eval $CMD
date
echo

CMD="bgzip -f ${omit_lowAAScore_out}.vcf"
echo $CMD
eval $CMD
CMD="tabix -f ${omit_lowAAScore_out}.vcf.gz"
echo $CMD
eval $CMD

# What did we end up with?
VCF="${omit_lowAAScore_out}.vcf.gz"
VCF_name=$(basename --suffix ".vcf.gz" $VCF)
CMD="bcftools stats $VCF | egrep \"^SN\"" 
echo
echo $CMD
eval $CMD

echo "Phasing VCF"
beagle="/nfs5/BPP/Grunwald_Lab/home/carleson/opt/bin/beagle.jar"
out_pre="phased/${VCF_name}.phased"
echo
CMD="java -Xmx100g -jar $beagle gt=$VCF out=$out_pre nthreads=$NSLOTS impute=FALSE"
date
echo $CMD
eval $CMD


echo "exit $?"
date
