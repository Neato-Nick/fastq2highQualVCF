#!/usr/bin/env Rscript

# Undo imputation - in the beagle output,
# convert any genotypes to missing that were missing in the unphased data

args = commandArgs(trailingOnly = TRUE)

usage_text <- "Usage: CensorImputedGTs.R <contig_name> <filtering strategy>\nExpects unphased data are in dir <contig_name>/<contig_name>.<filt_strategy>.vcf.gz\nPhased data in phased/<contig_name>.<filt_strat>.impute.phased.vcf.gz"
if ((length(args) == 0) || args[1] == "-h" || args[1] == "--help") {
        stop(paste0("Must provide input sequence ID.\n", usage_text))
} else if (length(args) == 1) {
        filter_strat <- ".bcft_Idp15_aa05_PASS.upad_MQ0_median_3x.vcft_Vdp4_PASS.bcft_RRmis_lowAA.omit_lowAAS"
        print(paste0("Filtering ", args[1], " with last filtering strategy tried:\n", filter_strat))
} else if (length(args) == 2) {
        filter_strat <- args[2]
} else {
        stop(paste0("Too many arguments supplied.\n", usage_text))
}

suppressMessages(library(vcfR))

ref_seq <- args[1]

message(paste0("Filtering contig: ", args[1]))

vcf_f_unphased <- paste0(ref_seq, "/", ref_seq, filter_strat, ".vcf.gz")
vcf_f_phased <- paste0("phased/", ref_seq, filter_strat, ".impute.phased.vcf.gz")
vcf_f_phased_out <- paste0("phased/", ref_seq, filter_strat, ".impute.phased.re-NA.vcf.gz")

vcf_unphased <- read.vcfR(vcf_f_unphased, verbose = FALSE)
vcf_phased <- read.vcfR(vcf_f_phased, verbose = FALSE)

gt_unphased <- extract.gt(vcf_unphased, element = "GT")
gt_phased <- extract.gt(vcf_phased, element = "GT")

# Convert missing genotypes in unphased data to missing in phased
vcf_phased@gt[,-1][is.na(gt_unphased) == TRUE ] <- NA
# This would get us missingness stats but not necesssary
# vcf_unphased@gt[,-1][is.na(gt_unphased) == TRUE ] <- NA

# Testing one specific locus I know has one known and one missing allele
#scaf_6  62167   scaf_6:62167:SG T       A,G     33043   PASS    AAScore=0.2845,0.8172;
#gt_unphased[rownames(gt_unphased) == "scaf_6:62167:SG",]
#            Plat_MPF4             Plat_MPF6             Plat_RH_5
#                "0/."                    NA                 "0/."
# BEAGLE5 imputed missing allele and overwrote what I had.
#gt_phased[rownames(gt_phased) == "scaf_6:62167:SG",]
#            Plat_MPF4             Plat_MPF6             Plat_RH_5
#                "1|1"                    NA                 "1|1"
#> grep("(./\\.)|(\\./.)", gt_unphased[rownames(gt_phased) == "scaf_6:62167:SG",], perl = TRUE)
# [1] 409 411
#> sum(grepl("(./\\.)|(\\./.)", gt_unphased[rownames(gt_phased) == "scaf_6:62167:SG",], perl = TRUE))
#[1] 2
# Convert to missing b/c there were actually a lot of similar cases
# Can't phase one missing allele, so censor the whole GT

vcf_phased@gt[,-1][grepl("(./\\.)|(\\./.)", gt_unphased[,], perl = TRUE)] <- NA

write.vcf(vcf_phased, vcf_f_phased_out)
