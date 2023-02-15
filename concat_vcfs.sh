#!/bin/bash

#vcf_list=concat_vcfs_phased_lightfilt.list
#vcf_list=concat_vcfs_phased_filt-lowAAS.list
#vcf_list=concat_vcfs_filt-lowAAS.list
vcf_list=concat_vcfs_filt-lowAAS_re-na.list
#ls -1 phased/*.upad_MQ0_median_3x.vcft_Vdp4_Idp15_PASS.bcft_aa05.phased.vcf.gz > $vcf_list
#ls -1 phased/*bcft_Idp15_aa05_PASS.upad_MQ0_median_3x*lowAAS.phased.vcf.gz > $vcf_list
#ls -1 */*bcft_Idp15_aa05_PASS.upad_MQ0_median_3x.vcft_Vdp4_PASS.bcft_RRmis_lowAA.omit_lowAAS.vcf.gz > $vcf_list
#ls -1 phased/*bcft_Idp15_aa05_PASS.upad_MQ0_median_3x.vcft_Vdp4_PASS.bcft_RRmis_lowAA.omit_lowAAS.impute.phased.re-NA.vcf.gz > $vcf_list

prefix="ramorum_div_n49"
filt_strat="bcft_Idp15_aa05_PASS.upad_MQ0_median_3x.vcft_Vdp4_PASS.bcft_RRmis_lowAA.omit_lowAAS.impute.phased.re-NA"
#concat_vcf="phased/ramorum_div_n49.upad_MQ0_median_3x.vcft_Vdp4_Idp15_PASS.bcft_aa05.phased.vcf.gz"
#concat_vcf="phased/ramorum_div_n49.bcft_Idp15_aa05_PASS.upad_MQ0_median_3x.vcft_Vdp4_PASS.bcft_RRmis_lowAA.omit_lowAAS.phased.vcf.gz"
concat_vcf="concat/${prefix}.${filt_strat}.vcf.gz"
samples="sample_of_ramorum_div.samples.list"

# bcftools needs vcfs to be indexed
#for vcf in $(cat $vcf_list); do tabix $vcf; done

# Combine all chromosomes / phased for the samples we want
# Because we're subsetting samples, any singletons they harbored will no longer be informative.
# Remove these by excluding any variant where all genotypes are missing + homo-ref
bcftools concat -Ou --file-list $vcf_list --threads 4 | \
	bcftools view -Ou --samples-file $samples | \
	bcftools view --trim-alt-alleles -o $concat_vcf \
	--exclude '(COUNT(GT="RR")+(COUNT(GT="mis")))=N_SAMPLES'

bcftools stats -F ../../../PR-102_v4.fasta $concat_vcf > "stats/${prefix}.${filt_strat}.vchk"
plot-vcfstats -p "stats/${prefix}.${filt_strat}" -s "stats/${prefix}.${filt_strat}.vchk"

# EOF
