#!/bin/bash

vcf_list=concat_vcfs_phased_lightfilt.list
#ls -1 phased/*.upad_MQ0_median_3x.vcft_Vdp4_Idp15_PASS.bcft_aa05.phased.vcf.gz > $vcf_list
# Manually move pchr10:pchr13 to be after pchr9
#vim $vcf_list
concat_vcf="phased/ramorum_div_n49.upad_MQ0_median_3x.vcft_Vdp4_Idp15_PASS.bcft_aa05.phased.vcf.gz"
samples="sample_of_ramorum_div.samples.list"

# bcftools needs vcfs to be indexed
for vcf in $(cat $vcf_list); do tabix $vcf; done

# Combine all chromosomes / phased for the samples we want
# Because we're subsetting samples, any singletons they harbored will no longer be informative.
# Remove these by excluding any variant where all genotypes are missing + homo-ref
bcftools concat -Ou --file-list $vcf_list --threads 4 | \
	bcftools view -Ou --samples-file $samples | \
	bcftools view --trim-alt-alleles -o $concat_vcf \
	--exclude '(COUNT(GT="RR")+(COUNT(GT="mis")))=N_SAMPLES'

# EOF