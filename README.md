## Raw reads to filtered VCF
This collection of bash and python scripts are how I obtain VCF data from FASTQ files.

1. read_alignment_fastq2crams.sh
2. read_alignment_cram2bams.sh
*GraphTyper*
3. read_alignment_bam2vcf_graphtyper.sh
*GATK v4*
3. read_alignment_bam2gvcf.sh
4. combine_gvcfs_by_chrom.sh

There are lots of ways to filter, based on the objective of the study

To run PAML and phylogeography tools I wanted haplotype sequences from the VCF,
So I was less aggressive on filtering and did not discard based on minor allele count.
There will be additional scripts here to go from the phased VCF to these haplotypes.
1. filter_bcf_upadhyay.sh
2. concat_vcfs.sh 

Other times I just wanted a rigorous, high quality set of markers.
Since I was interested in one specific population,
I selected only variants polymorphic within a species,
filtered and censored genotypes,
extracted those loci from the full VCF,
and re-censored genotypes.
