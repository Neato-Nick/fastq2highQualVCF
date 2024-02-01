#!/bin/bash
# Author: Nick Carleson
# Usage: SGE_Batch -c './get_mean_dps_bams.sh' -q bpp -r median_dps_n669 -P 2

# Input: List of bams
# Outputs:
#   * One stats file per sample with the median read depth of every contig
#   * One two-columnm, tab-separated summary file: sample_name	average_depth
# Average depths are required for filtering based on average depth.
# Calculate for a list of bams.

mapQ=0
avg="median"

# bamlist should be a list of all the BAMs, one file per line.
# Absolute paths preferred but relative paths relative to execution dir works
#grep -v -f bamsEU1_160.list samples_globaldiv_n667.bams.txt | grep -v -f bamsNA1_135.list > samples_globaldiv_nocurryco_n372.list
bam_list="/oscar/data/ccuomo1/Cauris_Res/allterra_allsra.n_1606.bams.list"
mean_dps_out="stats/mean-dp_global-Cauris_n1606.B11205.mosdepth_mapQ${mapQ}_${avg}.txt"


which mosdepth
mosdepth --version

touch $mean_dps_out
while read bam; do
        #echo "Calculating mean cov in:"
        #ls $bam
        sm_name=$(basename --suffix ".reordered.bam" $bam)
	mosdp_out_stem="stats/${sm_name}_mapQ${mapQ}_${avg}.B11205"
        mosdp_summ="${mosdp_out_stem}.mosdepth.summary.txt"
        CMD="mosdepth --no-per-base -t $SLURM_CPUS_PER_TASK --use-median -x --mapq $mapQ ${mosdp_out_stem} $bam"
        if [ ! -f $mosdp_summ ]
        then
                echo $CMD
                eval $CMD
        else
                echo "Mosdepth already ran on ${bam}. Simply outputting weighted mean"
        fi
        
        ls -lht $mosdp_summ
        # get weighted average across all chromosomes
        mean_cov_mosdp=$(tail -n +2 $mosdp_summ | awk '{w = w + $2; e = e + $4 * $2;} END {print e/w}')
        printf "${sm_name}\t${mean_cov_mosdp}\n"
        printf "${sm_name}\t${mean_cov_mosdp}\n" >> $mean_dps_out

        #date
        #echo
done < $bam_list
wc -l $mean_dps_out
