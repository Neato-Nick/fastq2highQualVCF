#! /bin/bash

#$ -N combine_gvcfs_allpops_20210507_n653
#$ -e genotype_err
#$ -o genotype_out
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -l mem_free=450g
###$ -q bpp
###$ -q *@!(cedro*|megraw*|cerebro*|samwise*|nem*|galls*|oryza2*)
#$ -q *@!(cedro*|samwise*|nem*|galls*|fungi*)
#$ -t 3,8
#$ -tc 3
#$ -m ea
#$ -M carleson@oregonstate.edu

# easiest to combine a gvcf first before genotyping, with gatk4+
# combine gvcfs here in one thread
# then genotype the combined gvcf per chromosome!

i=$(expr $SGE_TASK_ID - 1)

GATK="$HOME/dfs_opt/gatk-4.1.7.0/gatk"
#REFDIR="/dfs/Grunwald_Lab/final_data/ramorum_pacbio_assemblies"
REFDIR="/dfs/Grunwald_Lab/lab_collabs/diagnostic_assays/genome_data"
#REF="${REFDIR}/PR-102.fasta.gz"
REF="${REFDIR}/PR-102_v3.1.fasta"
FAI=( `cut -f 1 "${REF}.fai" `)
read scaffold <<< "${FAI[$i]}"

# Set file parameters - where files come from and go from
echo "Using $GATK"
TEMP_DIR="genotype_tmp"
out_folder=$(echo "$PWD/vcfs")
out_name="all_pops_n653.scaf${SGE_TASK_ID}"
combined_gvcf_out="${out_folder}/combined_${out_name}.g.vcf.gz"
genotyped_vcf_out="${out_folder}/genotyped_${out_name}.vcf.gz"
polymorph_vcf_out="${out_folder}/poly_bi.genotyped_${out_name}.vcf.gz"

# list of gvcfs to use!
#gvcfs_list_file="vcfs/tmp/all_pops_n276.files.list"
gvcfs_list_file="vcfs/tmp/global-div-allpops-n653_gvcfs.list"

echo "Landed as $JOB_ID on:"
hostname
echo "Executing $SGE_TASK_ID"

echo
echo "Now combining gvcfs"
over_intervals="-L $scaffold"
CMD="$GATK --java-options \"-Xmx400g\" CombineGVCFs \
	-R $REF -V $gvcfs_list_file  -O ${combined_gvcf_out} \
	$over_intervals \
	--annotation TandemRepeat \
	--tmp-dir $TEMP_DIR"
date
#echo $CMD
#eval $CMD
date

echo
echo "Genotyping ith contig $SGE_TASK_ID $scaffold"
CMD="$GATK --java-options '-Xmx250g -Djava.io.tmpdir=genotype_tmp -XX:ParallelGCThreads=1' GenotypeGVCFs \
	-R $REF -V $combined_gvcf_out -O $genotyped_vcf_out"
date
echo $CMD
eval $CMD
date

echo
# Remove monomorphic and multi-allelic variants
CMD="$GATK SelectVariants --variant $genotyped_vcf_out --output $polymorph_vcf_out \
	--exclude-non-variants true --restrict-alleles-to BIALLELIC"
date
echo $CMD
eval $CMD
date

echo
# how many samples and variants did we end up up with?
zgrep -m 1 "CHROM" $polymorph_vcf_out | cut -f 10- | sed 's/\t/\n/g' | wc -l
CMD="$GATK CountVariants --variant $polymorph_vcf_out"
date
echo $CMD
eval $CMD
date


# EOF
