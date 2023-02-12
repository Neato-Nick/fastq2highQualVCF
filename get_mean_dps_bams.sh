#!/bin/bash
# Author: Nick Carleson
# Usage: SGE_Batch -c './get_mean_dps_bams.sh' -q bpp -r median_dps_n669 -P 2

# Mean depths required for filtering based on average depth.
# Calculate for a list of bams.

mapQ=0
avg="median"

#grep -v -f bamsEU1_160.list samples_globaldiv_n667.bams.txt | grep -v -f bamsNA1_135.list > samples_globaldiv_nocurryco_n372.list
bam_list=samples_globaldiv_hybrids_n669.bams.list
mean_dps_out=stats/mean-dp_globaldiv_n669.mosdepth_mapQ${mapQ}_${avg}.txt

touch $mean_dps_out
while read bam; do
	#echo "Calculating mean cov in:"
	#ls $bam
	sm_name=$(basename --suffix "_dupmrk.bam" $bam)
	#mean_cov=$(samtools coverage $bam | tail -n +2 | awk '{w = w + $3; e = e + $7 * $3;} END {print e/w}')
	#echo $mean_cov
	#printf "${sm_name}\t${mean_cov}\n" >> $mean_dps_out
    CMD="~/dfs_opt/mosdepth --no-per-base -t $NSLOTS --use-median -x --mapq $mapQ stats/${sm_name}_mapQ${mapQ}_${avg} $bam"
    echo $CMD
    eval $CMD
    #mean_cov_mosdp=$(tail -n +2 ${sm_name}_mapQ${mapQ}.mosdepth.summary.txt | awk '{w = w + $2; e = e + $4 * $2;} END {print e/w}')
    ls stats/${sm_name}_mapQ${mapQ}.mosdepth.summary.txt
    # get weighted average across all chromosomes
	mean_cov_mosdp=$(tail -n +2 stats/${sm_name}_mapQ${mapQ}_${avg}.mosdepth.summary.txt | awk '{w = w + $2; e = e + $4 * $2;} END {print e/w}')
    printf "${sm_name}\t${mean_cov_mosdp}\n"
    printf "${sm_name}\t${mean_cov_mosdp}\n" >> $mean_dps_out
done < $bam_list
wc -l $mean_dps_out
