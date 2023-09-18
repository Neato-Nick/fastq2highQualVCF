#!/bin/bash

# Usage: SGE_Batch -c 'roh_plink.sh' -q samwise -P 3 -r roh-plink

GATK="$HOME/dfs_opt/gatk-4.1.7.0/gatk"

# Subset original VCF to NA1 or EU1
# 2023-09-11: Not necessary for this dataset
VCF_dir="/nfs5/BPP/Grunwald_Lab/home/carleson/ramorum/illumina_mapping_pr102v4/vcfs/graphtyper/freshNA1_currycoNA1/concat"
VCF_orig="$VCF_dir/bcft_Idp15_aa05_PASS.upad_MQ0_median_3x.vcft_Vdp4_PASS.bcft_RRmis_lowAA.omit_lowAAS.prams_poly-bi-snps.mq20_vmiss0.25.rm-sms-50perc.vmiss0.75_mac2.re-PlatPhibePfoli.vcf.gz"
#VCF_orig="bcft_Idp15_aa05_PASS.upad_MQ0_median_3x.vcft_Vdp4_PASS.bcft_RRmis_lowAA.omit_lowAAS.prams_poly-bi-snps.mq20_vmiss0.25.rm-sms-50perc.vmiss0.75_mac2.re-PlatPhibePfoli.vcf.gz"
VCF_orig_pre=$(basename --suffix ".vcf.gz" $VCF_orig)
#curryco_NA1="${VCF_orig_pre}.curryco_NA1"
#curryco_EU1="${VCF_orig_pre}.curryco_EU1"
#asia_IC1="${VCF_orig_pre}.asia_IC1"
curryco_NA1="${VCF_orig_pre}"
#curryco_NA1_samples=samplesNA1_135.args
#curryco_EU1_samples=samplesEU1_160.args
# awk '$2 ~ /IC1/ {print $1}' /nfs5/BPP/Grunwald_Lab/home/carleson/ramorum/illumina_mapping_pr102v4/vcfs/graphtyper/globaldiv_n669/pram_lookuptable_clean_added_info_lineages.tsv > samplesIC1_55.args
#asia_IC1_samples=samplesIC1_55.args

echo
CMD="$GATK SelectVariants --variant $VCF_orig --output ${curryco_NA1}.vcf.gz --sample-name=$curryco_NA1_samples \
	--exclude-non-variants true --restrict-alleles-to BIALLELIC --remove-unused-alternates true \
	--select-type-to-include SNP --allow-nonoverlapping-command-line-samples"
date
#echo $CMD
#eval $CMD &
echo
CMD="$GATK SelectVariants --variant $VCF_orig --output ${curryco_EU1}.vcf.gz --sample-name=$curryco_EU1_samples \
	--exclude-non-variants true --restrict-alleles-to BIALLELIC --remove-unused-alternates true \
	--select-type-to-include SNP --allow-nonoverlapping-command-line-samples"
#echo $CMD
#eval $CMD &
echo
CMD="$GATK SelectVariants --variant $VCF_orig --output ${asia_IC1}.vcf.gz --sample-name=$asia_IC1_samples \
	--exclude-non-variants true --restrict-alleles-to BIALLELIC --remove-unused-alternates true \
	--select-type-to-include SNP --allow-nonoverlapping-command-line-samples"
#echo $CMD
#eval $CMD &
echo

#echo "Waiting for all selectvariants to finish"
date
jobs
wait
jobs
date

hws=20 # homozyg-window-snp
hk=5 # homozyg-kb
hd=250 # homozyg-density
CMD_1="plink --vcf $VCF_dir/${curryco_NA1}.vcf.gz --vcf-half-call m --homozyg --homozyg-window-snp $hws --homozyg-kb $hk --homozyg-density $hd --out ${curryco_NA1}.roh_plink_homozyg-window-snp${hws}_homozyg-kb${hk}_homozyg-density${hd} --allow-extra-chr 0"
CMD_2="plink --vcf ${curryco_EU1}.vcf.gz --vcf-half-call m --homozyg --homozyg-window-snp $hws --homozyg-kb $hk --homozyg-density $hd --out ${curryco_EU1}.roh_plink_homozyg-window-snp${hws}_homozyg-kb${hk}_homozyg-density${hd} --allow-extra-chr 0"
CMD_3="plink --vcf ${asia_IC1}.vcf.gz --vcf-half-call m --homozyg --homozyg-window-snp $hws --homozyg-kb $hk --homozyg-density $hd --out ${asia_IC1}.roh_plink_homozyg-window-snp${hws}_homozyg-kb${hk}_homozyg-density${hd} --allow-extra-chr 0"

echo
date
echo $CMD_1
eval $CMD_1
echo
#echo $CMD_2
#eval $CMD_2 &
#echo
#echo $CMD_3
#eval $CMD_3 &
#echo

echo "Waiting for all PLINK ROH inference to finish"
jobs
wait
echo "All jobs done"
date
jobs

# PLINK run setup before I edited

#VCF_1="data/vcfs/curry_county/genotyped_CurryCo_EU2_NA2reps_n322_wg.gatk_select_curryco.upad_MQ20_3x.vcft_dp4_Q20.recode.gatk_imiss50perc.vcft_rm-snps25perc_mac2.recode.vcf.gz"
#VCF_2="data/vcfs/curry_county/genotyped_CurryCo_NA1_n135_wg.gatk_select.upad_MQ20_3x.vcft_dp4_Q20.recode.gatk_imiss50perc.vcft_rm-snps25perc_mac2.recode.vcf.gz"
#VCF_3="vcfs/genotyped_refEU1_CurryCo_EU2_NA2reps_n322_wg.gatk_select.upad_MQ20_3x.vcft_dp4_Q20.recode.gatk_imiss50perc.vcft_rm-snps25perc_mac2.recode.polarized.vcf.gz"
#vcf_pre_1=$(basename --suffix ".vcf.gz" $VCF_1)
#vcf_pre_2=$(basename --suffix ".vcf.gz" $VCF_2)
#vcf_pre_3=$(basename --suffix ".vcf.gz" $VCF_3)
#VCF_EU1="data/vcfs/curry_county/genotyped_CurryCo_EU2_NA2reps_n322_wg.gatk_select_curryco.upad_MQ20_3x.vcft_dp4_Q20.recode.gatk_imiss50perc.vcft_rm-snps25perc_mac2.recode.curryco_EU1.vcf.gz"
#vcf_pre_EU1=$(basename --suffix ".vcf.gz" $VCF_EU1)
#VCF_NA1="data/vcfs/curry_county/genotyped_CurryCo_EU2_NA2reps_n322_wg.gatk_select_curryco.upad_MQ20_3x.vcft_dp4_Q20.recode.gatk_imiss50perc.vcft_rm-snps25perc_mac2.recode.curryco_NA1.vcf.gz"
#vcf_pre_NA1=$(basename --suffix ".vcf.gz" $VCF_NA1)

hws=20 # homozyg-window-snp
hk=10 # homozyg-kb
hd=250 # homozyg-density
CMD_1="plink --vcf $VCF_1 --homozyg --homozyg-window-snp $hws --homozyg-kb $hk --homozyg-density $hd --out data/roh/${vcf_pre_1}.roh_plink_homozyg-window-snp${hws}_homozyg-kb${hk}_homozyg-density${hd} --allow-extra-chr"
CMD_2="plink --vcf $VCF_2 --homozyg --homozyg-window-snp $hws --homozyg-kb $hk --homozyg-density $hd --out data/roh/${vcf_pre_2}.roh_plink_homozyg-window-snp${hws}_homozyg-kb${hk}_homozyg-density${hd} --allow-extra-chr"
CMD_3="plink --vcf $VCF_3 --homozyg --homozyg-window-snp $hws --homozyg-kb $hk --homozyg-density $hd --out roh/${vcf_pre_3}.roh_plink_homozyg-window-snp${hws}_homozyg-kb${hk}_homozyg-density${hd} --allow-extra-chr"
#plink --vcf $VCF_NA1 --homozyg --homozyg-window-snp 20 --homozyg-kb 10 --homozyg-density 250 --out data/roh/${vcf_pre_NA1}.roh_plink_homozyg-window-snp20_homozyg-kb10_homozyg-density250 --allow-extra-chr
#plink --vcf $VCF_EU1 --homozyg --homozyg-window-snp 20 --homozyg-kb 10 --homozyg-density 250 --out data/roh/${vcf_pre_EU1}.roh_plink_homozyg-window-snp20_homozyg-kb10_homozyg-density250 --allow-extra-chr
#echo $CMD_1
#eval $CMD_1
#echo
#echo $CMD_2
#eval $CMD_2
#echo
#echo $CMD_3
#eval $CMD_3
#echo

hws=20 # homozyg-window-snp
hk=5 # homozyg-kb
hd=250 # homozyg-density
CMD_1="plink --vcf $VCF_1 --homozyg --homozyg-window-snp $hws --homozyg-kb $hk --homozyg-density $hd --out data/roh/${vcf_pre_1}.roh_plink_homozyg-window-snp${hws}_homozyg-kb${hk}_homozyg-density${hd} --allow-extra-chr"
CMD_2="plink --vcf $VCF_2 --homozyg --homozyg-window-snp $hws --homozyg-kb $hk --homozyg-density $hd --out data/roh/${vcf_pre_2}.roh_plink_homozyg-window-snp${hws}_homozyg-kb${hk}_homozyg-density${hd} --allow-extra-chr"
CMD_3="plink --vcf $VCF_3 --homozyg --homozyg-window-snp $hws --homozyg-kb $hk --homozyg-density $hd --out roh/${vcf_pre_3}.roh_plink_homozyg-window-snp${hws}_homozyg-kb${hk}_homozyg-density${hd} --allow-extra-chr"
#echo $CMD_1
#eval $CMD_1
#echo
#echo $CMD_2
#eval $CMD_2
#echo
#echo $CMD_3
#eval $CMD_3
#echo
