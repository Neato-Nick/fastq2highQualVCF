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
