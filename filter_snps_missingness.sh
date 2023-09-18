#!/bin/bash

#$ -hold_jid_ad phase_unimpute
#$ -V
#$ -N bcffilter_to-prams
#$ -e filter_err
#$ -o filter_out
###$ -q *@!(nem*|samwise*|amp*|debary*|galls*|anduin*|symbiosis*)
###$ -q nem
###$ -q nem
#$ -q nem
#$ -l mem_free=20G
#$ -S /bin/bash
#$ -cwd
#$ -t 1-24:1
#$ -tc 8
#$ -pe thread 6

i=$(expr $SGE_TASK_ID - 1)

echo -n "Running on: "
hostname
echo "SGE job id: ${JOB_ID}.${SGE_TASK_ID}"
date
echo

REFDIR="/nfs5/BPP/Grunwald_Lab/home/carleson/ramorum/illumina_mapping_pr102v4"
REF="${REFDIR}/PR-102_v4.fasta"

FAI=( `cut -f 1 "${REF}.fai" `)
read scaffold <<< "${FAI[$i]}"

# Files set up
workdir=$scaffold
group_filt_vcf=$workdir/${scaffold}.bcft_Idp15_aa05_PASS.upad_MQ0_median_3x.vcft_Vdp4_PASS.bcft_RRmis_lowAA.omit_lowAAS.vcf.gz
group_filt_name=$(basename --suffix ".vcf.gz" $group_filt_vcf)

#egrep -v "Plat|Phib|Pfoli" mean-dp_globaldiv_n669.mosdepth_mapQ0_median.minDP15.include.list > mean-dp_globaldiv_n669.mosdepth_mapQ0_median.minDP15.include.Prams.list
#grep "Plat" sample_of_ramorum_div.samples.median_dp_mapQ0.txt | cut -f 1 > mean-dp_globaldiv_n669.mosdepth_mapQ0_median.minDP15.include.Prams.Plat_rep.list
#cat mean-dp_globaldiv_n669.mosdepth_mapQ0_median.minDP15.include.Prams.list >> mean-dp_globaldiv_n669.mosdepth_mapQ0_median.minDP15.include.Prams.Plat_rep.list
# egrep -v "Plat|Phib|Pfoli" mean-dp_curryco_freshNA1_n144.mosdepth_mapQ0_median.minDP15.include.list > mean-dp_curryco_freshNA1_n144.mosdepth_mapQ0_median.minDP15.include.Prams.list

# First include only Prams and snps
# Trim unused alternate alleles, to ensure accurate biallelic filtering in next step
# Get only biallelics, exclude vars where all sms are either hom ref or hom alt
# Remove variants below a minimum mapping quality or above a max % missing data (baseline)
# (doesn't need to be in a sep step but I did it anyway for readability)
min_mapQ=20
max_missing_vars_1=0.25 # 0 allows sites all missing, 1 allows no missing data
pram_basic="$workdir/${group_filt_name}.prams_poly-bi-snps.mq20_vmiss0.25"
pram_basic1="$workdir/${group_filt_name}.prams"
pram_basic2="$workdir/${group_filt_name}.prams_poly-bi-snps"
pram_basic3="$workdir/${group_filt_name}.prams_poly-bi-snps.mq20"
pram_basic4="$workdir/${group_filt_name}.prams_poly-bi-snps.vmiss0.25"
pram_basic5="$workdir/${group_filt_name}.prams_poly-bi-snps.mq20_vmiss0.25"
pram_basic_vcf="${pram_basic}.vcf.gz"
bcftools view --threads $NSLOTS -Ou -v snps \
	-S mean-dp_curryco_freshNA1_n144.mosdepth_mapQ0_median.minDP15.include.Prams.list \
	--trim-alt-alleles $group_filt_vcf | \
	bcftools view --threads $NSLOTS -Ou -m2 -M2 \
	-e 'COUNT(GT="AA")+COUNT(GT="mis")=N_SAMPLES || COUNT(GT="RR")+COUNT(GT="mis")=N_SAMPLES' | \
	bcftools view --threads $NSLOTS \
	-e 'INFO/MQ<20 || F_MISSING>0.75' \
	-o $pram_basic_vcf
CMD="bcftools stats $pram_basic_vcf | egrep \"^SN\""
echo $CMD
eval $CMD
#bcftools view --threads $NSLOTS -o ${pram_basic1}.vcf.gz -v snps \
#	-S mean-dp_globaldiv_n669.mosdepth_mapQ0_median.minDP15.include.Prams.list \
#	--trim-alt-alleles $group_filt_vcf
#CMD="bcftools stats ${pram_basic1}.vcf.gz | egrep \"^SN\""
#echo $CMD
#eval $CMD
#bcftools view --threads $NSLOTS -o ${pram_basic2}.vcf.gz -m2 -M2 \
#	-e 'COUNT(GT="AA")+COUNT(GT="mis")=N_SAMPLES || COUNT(GT="RR")+COUNT(GT="mis")=N_SAMPLES' \
#	${pram_basic1}.vcf.gz
#CMD="bcftools stats ${pram_basic2}.vcf.gz | egrep \"^SN\""
#echo $CMD
#eval $CMD
#bcftools view --threads $NSLOTS \
#	-e 'INFO/MQ>=20' \
#	-o ${pram_basic3}.vcf.gz ${pram_basic2}.vcf.gz
#CMD="bcftools stats ${pram_basic3}.vcf.gz | egrep \"^SN\""
#echo $CMD
#eval $CMD
#bcftools view --threads $NSLOTS \
#	-e 'F_MISSING>0.75' \
#	-o ${pram_basic4}.vcf.gz ${pram_basic2}.vcf.gz
#CMD="bcftools stats ${pram_basic4}.vcf.gz | egrep \"^SN\""
#echo $CMD
#eval $CMD
#bcftools view --threads $NSLOTS \
#	-e 'INFO/MQ>=20 || F_MISSING>0.75' \
#	-o ${pram_basic5}.vcf.gz ${pram_basic2}.vcf.gz
#CMD="bcftools stats ${pram_basic5}.vcf.gz | egrep \"^SN\""
#echo $CMD
#eval $CMD



# separate job that's held until array finishes:
# concatenate all chroms,
# remove samples >50% missing data
# remove variants >25% missing data, recode info.
# vcfA: remove variants minor allele count < 2
# vcfB: remove variants mac < 2 or maf < 0.05
# Concatenate all chroms and select only these samples and loci in vcfA |
# trim unused alt alleles |
# Remove multi-allelic SNPs again - some will newly have a 3rd allele from an outgroup
# bcftools stats -s - > file.vchk
# plot-vcfstats -p outdir file.vchk
