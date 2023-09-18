#!/bin/bash

#$ -hold_jid bcffilter_to-prams
#$ -V
#$ -N concat_pramselect_macfilter_freshNA1 
#$ -e filter_err
#$ -o filter_out
#$ -q bpp@anduin
#$ -l mem_free=20G
#$ -S /bin/bash
#$ -cwd
#$ -pe thread 8

# Filter out ramorum samples -> vars based on missingness %
# then re-include outgroups at the end
# Two VCFs generated, one just no singletons and the also with higher MAF

echo "Landed on:"
hostname

# Remove samples on missingness % can't be done in bcftools alone
GATK="~/grunwald_lab_dfs_me/opt/bin/gatk"
vcfdir="/nfs5/BPP/Grunwald_Lab/home/carleson/ramorum/illumina_mapping_pr102v4/vcfs/graphtyper/freshNA1_currycoNA1"
group_filt_name="bcft_Idp15_aa05_PASS.upad_MQ0_median_3x.vcft_Vdp4_PASS.bcft_RRmis_lowAA.omit_lowAAS"
pram_basic_concat="concat/${group_filt_name}.prams_poly-bi-snps.mq20_vmiss0.25"
vcf_list=poly-bi_prams.list
CMD="bcftools concat --file-list $vcf_list --threads $NSLOTS -o ${pram_basic_concat}.vcf.gz"
if [ ! -f ${pram_basic_concat}.vcf.gz ]
then
	echo $CMD
	eval $CMD
	CMD="$GATK IndexFeatureFile -I ${pram_basic_concat}.vcf.gz"
	echo $CMD
	eval $CMD
fi

tmpdir=/data/carleson/bcftools_concat/${group_filt_name}.prams_poly-bi-snps.mq20_vmiss0.25
echo $tmpdir
vcf_list_base=concat_vcfs_filt-lowAAS.list
mkdir -p $tmpdir && cd $tmpdir
mkdir concat
pwd -P
ln -s $vcfdir/$vcf_list_base ./
ln -s $vcfdir/${pram_basic_concat}.vcf.gz* concat/
for vcf in $(cat $vcf_list_base); do
	chr=$(dirname $vcf)
	vcf_f=$(basename $vcf)
	mkdir $chr
	#ln -sf $vcfdir/$vcf $chr/$vcf_f
	ln -sf $vcfdir/${vcf}* $chr/
done

high_miss_sms="${pram_basic_concat}.rm-sms-50perc"
CMD="vcftools --gzvcf ${pram_basic_concat}.vcf.gz --missing-indv --out $high_miss_sms"
if [ ! -f ${high_miss_sms}.imiss ]
then
	echo $CMD
	eval $CMD
fi
# awk messes up if I cmd this in quotes
# while here get include list with outgroup P lat
cat ${high_miss_sms}.imiss | awk '{if($5>0.5)print $1}' | grep -v INDV> ${high_miss_sms}.exclude.args
cat ${high_miss_sms}.imiss | awk '{if($5<=0.5)print $1}' | grep -v INDV> ${high_miss_sms}.include.args
# Add a few lateralis samples and outgroups
#printf "Plat_RH_5\nPlat_SMST21\nPlat_SMSTG\nPhibe\nPfoli\n" >> ${high_miss_sms}.include.args
#grep "Plat" sample_of_ramorum_div.samples.median_dp_mapQ0.txt | cut -f 1 >> ${high_miss_sms}.include.args
CMD="$GATK SelectVariants --variant ${pram_basic_concat}.vcf.gz --output ${high_miss_sms}.vcf.gz \
        --exclude-non-variants true --remove-unused-alternates true \
        --exclude-sample-name ${high_miss_sms}.exclude.args"
if [ ! -f ${high_miss_sms}.vcf.gz ]
then
	echo $CMD
	eval $CMD
	date
	echo
fi

# remove variants >25% missing data, recode info.
# vcfA: remove variants >25% missing data, minor allele count < 2
no_singles="${high_miss_sms}.vmiss0.75_mac2"
if [ ! -f ${no_singles}.vcf.gz ]
then
	bcftools view --threads $NSLOTS \
		-e 'F_MISSING>0.25 || MAC<2' -o ${no_singles}.vcf.gz ${high_miss_sms}.vcf.gz
	# vcfB: remove variants >25% missing data, mac < 2 or maf < 0.05
	no_singles_lowmaf="${high_miss_sms}.vmiss0.75_mac2_maf0.05"
	bcftools view --threads $NSLOTS \
		-e 'F_MISSING>0.25 || MAC<2 || MAF<0.05' -o ${no_singles_lowmaf}.vcf.gz ${high_miss_sms}.vcf.gz
fi

# Concatenate all chroms from base and select only these samples+Plat and vcfA loci --trim-alt-alleles 
# Remove multi-allelic SNPs again - some will newly have a 3rd allele from an outgroup 
no_singles_nm=$(basename --suffix ".vcf.gz" ${no_singles}.vcf.gz)
no_singles_lowmaf_nm=$(basename --suffix ".vcf.gz" $no_singles_lowmaf)
# link regions file
CMD="ln -s $vcfdir/${no_singles}.vcf.gz* concat/"
echo $CMD
eval $CMD
CMD="ln -s $vcfdir/${high_miss_sms}.include.args concat/"
echo $CMD
eval $CMD
no_singles_base="${no_singles}.re-PlatPhibePfoli"
no_singles_lowmaf_base="${no_singles_lowmaf}.re-PlatPhibePfoli"
echo
CMD="bcftools concat --threads $NSLOTS --file-list $vcf_list_base --threads $NSLOTS -O z -o ${no_singles_base}.allconcat.vcf.gz -a && \
	tabix ${no_singles_base}.allconcat.vcf.gz"
date
echo $CMD
eval $CMD
date
echo
CMD="bcftools view --threads $NSLOTS -Ou -S ${high_miss_sms}.include.args -R ${no_singles}.vcf.gz ${no_singles_base}.allconcat.vcf.gz | \
	bcftools view --threads $NSLOTS -m2 -M2 --trim-alt-alleles -o ${no_singles_base}.vcf.gz"
echo $CMD
eval $CMD
date

echo "Finished concatenating full SNPs, running stats then will do MAF snps."
no_singles_base_nm=$(basename --suffix ".vcf.gz" ${no_singles_base}.vcf.gz)
mkdir stats
echo
CMD="bcftools stats ${no_singles_base}.vcf.gz > stats/${no_singles_base_nm}.vchk"
echo $CMD
eval $CMD
egrep "^SN" stats/${no_singles_base_nm}.vchk
echo
CMD="plot-vcfstats -p stats/${no_singles_base_nm}_stats stats/${no_singles_base_nm}.vchk"
echo $CMD
eval $CMD
date
echo

# link regions file
# Run again for low maf
echo
CMD="ln -s $vcfdir/${no_singles_lowmaf}.vcf.gz* concat/"
echo $CMD
eval $CMD
CMD="bcftools view --threads $NSLOTS -Ou -S ${high_miss_sms}.include.args -R ${no_singles_lowmaf}.vcf.gz ${no_singles_base}.allconcat.vcf.gz | \
	bcftools view --threads $NSLOTS -m2 -M2 --trim-alt-alleles -o ${no_singles_lowmaf_base}.vcf.gz"
date
echo $CMD
eval $CMD
date

echo
no_singles_lowmaf_base_nm=$(basename --suffix ".vcf.gz" ${no_singles_lowmaf_base}.vcf.gz)
mkdir stats
CMD="bcftools stats ${no_singles_lowmaf_base}.vcf.gz > stats/${no_singles_lomwaf_base_nm}.vchk"
echo $CMD
eval $CMD
egrep "^SN" stats/${no_singles_lowmaf_base_nm}.vchk
CMD="plot-vcfstats -p stats/${no_singles_lowmaf_base_nm}_stats stats/${no_singles_base_nm}.vchk"
echo $CMD
eval $CMD
date
echo

# copy back to /nfs5
rsync -av stats/ $vcfdir/stats
rsync ${no_singles_base}.vcf.gz $vcfdir/concat
rsync ${no_singles_lowmaf_base}.vcf.gz $vcfdir/concat

# EOF
