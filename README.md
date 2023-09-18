# Raw reads to filtered VCF
This collection of bash and python scripts are how I obtain VCF data from FASTQ files.

## FASTQ to VCF

They are submitted to SGE in this order,
after jobs from a previous step are completed.

1. read_alignment_fastq2crams.sh
2. read_alignment_cram2bams.sh

##### GraphTyper

3. read_alignment_bam2vcf_graphtyper.sh

##### GATK v4

3. read_alignment_bam2gvcf.sh
4. combine_gvcfs_by_chrom.sh

Steps 2-4 can be queued up after submitting Step 1.
1-2 (and Step 3 of GATK) are parallelized per sample, and can submit while waiting for earlier steps to finish e.g. `#$ -hold_jid_ad fastq2cram`
Step 3 Graphtyper (4 of GATK) is parallelized differently and shouldn't be performed until all samples are done.
Therefore, use e.g. `-hold_jid cram2bam` instead.
For details and correct usage, see: https://stackoverflow.com/a/28884945

## VCF Filtering

There are lots of ways to filter, based on the objective of the study

### For markers

Other times I just wanted a rigorous, high quality set of SNPs.
Since I was interested in one specific population but also wanted outgroups in the analysis,
I selected only variants polymorphic within a species,
filtered and censored genotypes,
extracted those loci from the full VCF,
and re-censored genotypes.


### For haplotypes 

To run PAML and phylogeography tools I wanted haplotype sequences from the VCF,
So I was less aggressive on filtering and did not discard based on minor allele count.
There will be additional scripts here to go from the phased VCF to these haplotypes.

1. filter_bcf_upadhyay.sh
     i. ReplaceGenoWithMissing_nc.py
     ii. OmitLowAAScoreAlleles.py
2. concat_vcfs.sh 

## LOH / ROH

From fastq to LOH calls, my workflow is:

1. read_alignment_fastq2cram.sh
2. read_alignment_cram2bam.sh
3. read_alignment_bam2vcf_graphtyper.sh
4. get_mean_dps_bams.sh
5. Filter VCF lightly
6. roh_plink.sh

### LOH inference tool

There are model-based and rules-based tools.
For a review, see Ceballos et al. 2018, *Nat. Rev. Genet.*
roh_plink.sh wraps PLINK v1.9, a rules-based tool.
I have also used bcftools roh, a model-based tool.
It produced similar calls to PLINK.
Clonal populations violate the assumptions made by the model-based bcftools roh.
That being said, the results were very similar for P. ramorum so either would probably would be fine.

There are a lot of parameters to PLINK roh. Important considerations are in Meyermans et al. 2020, *BMC Genomics*

### Filtering VCFs for LOH/ROH

False positive SNPs likely reduce LOHs inferred, while false negatives could both inflate and reduce the number inferred.
I tend to err on the cautious side, and lightly filter the VCF.
My filtering scripts focus on the in-population, then apply the filters to the VCF containing outgroups.
This results in high-quality SNPs for inferences on populations of inference, just using outgroups as anchor points on SNPs informative within-species.

n.b. These scripts are fairly project-specific.
To re-use them takes me ~1 hour to refamiliarize myself with inputs, outputs, and software dependencies.
The lightly filtered and phased SNPs could be used for other purposes, including haplotype reconstruction.
Phasing information **is not used** for LOH inference, but I needed the filtered VCF to be phased for other tools.

1. filter_bcf_upadhyay.sh
2. filter_snps_missingness.sh
3. Generate VCF lists
	1. ls -1 \*/\*.bcft_Idp15_aa05_PASS.upad_MQ0_median_3x.vcft_Vdp4_PASS.bcft_RRmis_lowAA.omit_lowAAS.vcf.gz > concat_vcfs_filt-lowAAS.list
	2. vim concat_vcfs_filt-lowAAS.list # manually natural sort chrom/scaffold order
	3. ls -1 \*/\*.bcft_Idp15_aa05_PASS.upad_MQ0_median_3x.vcft_Vdp4_PASS.bcft_RRmis_lowAA.omit_lowAAS.prams_poly-bi-snps.mq20_vmiss0.25.vcf.gz  > poly-bi_prams.list
	4. vim poly-bi_prams.list # manually natural sort chrom/scaffold order
4. filter_snps_missingness_concat_scratchdir.sh


