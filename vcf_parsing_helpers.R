# Parses snpeff annotations from VCF into a dataframe
# optionally, parse HGVS amino acide code to more familiar variant codes.
# solving the snpeff data cleanup took 1.5hrs. Fun stuff.
# remove_GT_AAs passed to parse_HGVS.p
tidy_snpeff <- function(vcf, AA_abbr = TRUE, remove_GT_AAs = TRUE) {
  # INFO has Key but not chrom/pos/ref/alt
  # and I need alt to ensure only correct annotations are shown
  tidy_fix <- vcfR2tidy(vcf, info_only = TRUE)$fix %>%
    select(CHROM, POS, REF, ALT) %>%
    mutate(Key = row_number())
  
  tidy_eff <- extract_info_tidy(vcf) %>%
    left_join(tidy_fix, by = c("Key")) %>%
    separate_longer_delim(c(AC, AF, MLEAC, MLEAF, ALT), delim = stringr::regex(","))
  
  tidy_eff <- parse_ANN(tidy_eff, filter_alt_allele = TRUE,
                        AA_abbr = AA_abbr, remove_GT_AAs = remove_GT_AAs)
  
  return(tidy_eff)
}

# Parses SnpEff annotations from a string, assuming ANN format
# Two main use-cases: in VCF, and in output of vcf2matrix python script
parse_ANN <- function(df, ann_col = "ANN", too_few = "error", filter_alt_allele = FALSE, AA_abbr = TRUE, remove_GT_AAs = TRUE) {
  # can be multiple annotations per mutation, need to separate those before clever separating
  # echo 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' | sed 's/ |/","/g' | sed 's/" /"/g'
  snpeff_cols <- c("Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID","Feature_Type","Feature_ID","Transcript_BioType","Rank","HGVS.c","HGVS.p","cDNA.pos / cDNA.length","CDS.pos / CDS.length","AA.pos / AA.length","Distance","ERRORS / WARNINGS / INFO")
  # echo "Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected" | sed 's/ |/","/g' | sed 's/" /"LOF_/g'
  LOF_cols <- c("LOF_Gene_Name","LOF_Gene_ID","LOF_Number_of_transcripts_in_gene","LOF_Percent_of_transcripts_affected")
  NMD_cols <- c("NMD_Gene_Name","NMD_Gene_ID","NMD_Number_of_transcripts_in_gene","NMD_Percent_of_transcripts_affected")
  
  ann_df <- df %>%
    separate_longer_delim(ANN, delim = stringr::regex("(?<=\\|),")) %>%
    separate_longer_delim(ANN, delim = stringr::regex(",(?=[ACTGN]{1,1000})")) %>%
    separate_wider_delim(cols = ANN, names = c(snpeff_cols), delim = "|", too_few = too_few)
  if(filter_alt_allele) {
    ann_df <- filter(ann_df, Allele == ALT)
  }
  if("LOF" %in% colnames(ann_df)) {
    ann_df <- ann_df %>%
      mutate(LOF = str_remove(LOF, "^\\(")) %>%
      mutate(LOF = str_remove(LOF, "\\)$")) %>%
      separate_wider_delim(LOF, names = LOF_cols, delim = "|")
  }
  if("NMD" %in% colnames(ann_df)) {
    ann_df <- ann_df %>%
      mutate(NMD = str_remove(NMD, "^\\(")) %>%
      mutate(NMD = str_remove(NMD, "\\)$")) %>%
      separate_wider_delim(NMD, names = NMD_cols, delim = "|")
  }
  ann_df <- ann_df %>%
    mutate(across(where(is.character), ~ na_if(., ""))) %>% # blank cells to NA
    type_convert()
  
  if(AA_abbr) {
    ann_df <- parse_HGVS.p(ann_df, "HGVS.p", remove_GT_AAs)
  }
}

# Parse HGVS amino acid code to the more familiar format of
# <AA_ref><AA_pos><AA_alt>
# Not all HGVS codes are supported for easy conversion,
# this was originally designed for simple substitutions.
parse_HGVS.p <- function(df, col = "HGVS.p", remove_GT_AAs = TRUE) {
  tidy_eff <- df %>% # now convert to old nomenclature
    separate_wider_delim(all_of(col), delim = ".", names = c(NA, "AA"),
                         cols_remove = FALSE, too_few = "align_start") %>%
    mutate(AA = str_replace(AA, "\\*$", "Xaa")) %>%
    mutate(AA = str_replace(AA, "(?<=\\d)\\?", "Del")) %>%
    separate_wider_delim(AA, delim = stringr::regex("(?=\\w{3})(?=\\d)"),
                         names = c("AA_ref", "AA_posalt"),
                         cols_remove = TRUE, too_many = "merge") %>%
    separate_wider_delim(AA_posalt, delim = stringr::regex("(?=[A-Z,a-z]{3,100})|(?=fs)|_"), 
                         names = c("AA_pos", "AA_alt"), cols_remove = TRUE,
                         too_many = "merge", too_few = "align_end")
  
  tidy_eff$AA_ref_ab <- names(Biostrings::AMINO_ACID_CODE)[match(tidy_eff$AA_ref, Biostrings::AMINO_ACID_CODE)]
  tidy_eff$AA_alt_ab <- names(Biostrings::AMINO_ACID_CODE)[match(tidy_eff$AA_alt, Biostrings::AMINO_ACID_CODE)]
  # Get abbreviation of amino acid if possible, otherwise use from HGVS e.g. "del"
  # Then join the reference and alternate abbreviations
  tidy_eff <- mutate(tidy_eff, AA_alt_ab = coalesce(AA_alt_ab, AA_alt)) %>%
    unite("AA_var",
          AA_ref_ab, AA_pos, AA_alt_ab,
          sep = "", na.rm = TRUE, remove = FALSE) %>%
    mutate(AA_var = na_if(AA_var, ""))
  if(remove_GT_AAs) {
    tidy_eff <- tidy_eff %>%
      select(-c(AA_alt_ab, AA_ref, AA_alt,
                AA_ref_ab, AA_pos))
  }
  return(tidy_eff)
}

# Parse HGVS codon code to a simpler format.
# Mostly intended for interpreting variants with "modifier" effect,
# e.g. promoter regions
parse_HGVS.c <- function(df, col = "HGVS.c", coalesce_HGVS.a = TRUE) {
  tidy_eff <- df %>%
    separate_wider_delim(all_of("HGVS.c"), delim = ".", names = c(NA, "codon"),
                         cols_remove = FALSE, too_few = "align_start") %>%
    mutate(codon = str_replace(codon, "ins", ".>")) %>%
    mutate(codon = str_remove(codon, "\\*")) %>%
    separate_wider_delim(codon, names = c("nt_pos_ref", "nt_alt"), delim = regex(">"),
                         too_few = "align_start", cols_remove = FALSE)
  # TODO: Support deletions
  
  if(coalesce_HGVS.a) {
    tidy_eff <- mutate(tidy_eff, AA_var = coalesce(AA_var, codon))
  }
  
  return(tidy_eff)
}

# Prepare genotype matrix using VCF or tidy-extracted data
# To reduce memory overhead when bug-testing, use prefiltered info and GT
# Final output is matrix with rownames according to individual IDs
# TODO: need to check if any loci within a gene have Alt at a locus in goodmaf
vcf_to_gt_matr <- function(vcf, tidy_filt_info = NULL, tidy_filt_gts = NULL, max_indel_size = 3, max_prop_amb_samples = 0.1, plot_prop_amb = FALSE,
                           coding_only = TRUE, ignore_low = FALSE, collapse_freq = 0.01, collapse_singletons = FALSE, separate_singletons = FALSE) {
  
  # Prep info
  if(is.null(tidy_filt_info)) {
    tidy_filt_info <- vcfr_to_gwas_info(vcf, max_indel_size,
                                        coding_only, ignore_low, collapse_freq)
  } 
  
  # Prep GTs
  if(is.null(tidy_filt_gts)) {
    tidy_filt_gts <- vcfr_to_gwas_gts(vcf, tidy_filt_info, max_prop_amb_samples, plot_prop_amb)
  }
  
  # Measure AF. Assumes output order is same as input variant order
  maf_tb <- as_tibble(maf(vcf_raw), rownames = "chr_pos") %>%
    separate_wider_delim(chr_pos, delim = "_", names = c("chr", "pos")) %>%
    rowid_to_column("Key") %>%
    filter(Key %in% tidy_filt_info$Key) %>%
    rename("MAC" = Count, "MAF" = Frequency, "missingGTcount" = `NA`)
  
  key_chr_pos <- vcfR2tidy(vcf_raw, info_only = TRUE)$fix %>%
    select(CHROM, POS) %>%
    rowid_to_column("Key")
  
  # Join GT matrix and chr / pos info
  tidy_filt_gts <- tidy_filt_gts %>%
    left_join(maf_tb, by = "Key") %>%
    type_convert()
  # Somehow collapse variants based on collapse_freq
  # Overall strategy: Make separate matrices of regular and collapsed, then combine
  # 1. Filter down to variants that DO pass MAF threshold
  # These are the full variants we can look at
  # Allow all alt alleles to be combined per locus, it's easier for now
  goodmaf_vars <- filter(tidy_filt_gts, MAF > collapse_freq) %>%
    mutate(bin_allele = case_when(
      gt_GT == 0 ~ 0,
      gt_GT >= 1 ~ 1)) %>%
    unite(locus, c(chr, pos)) %>%
    select(Indiv, bin_allele, locus)
  # 2. Collapse genes for rare variants.
  # OK to filter on gene name even if only alt alleles have that since thats all we care about
  # should I exclude indivs that have common variants? nah leave in for now
  # TODO: Need to verify on a larger dataset. Not enough 0s for me to be satisfied.
  lowmaf_vars_multiples <- filter(tidy_filt_gts, MAF < collapse_freq) %>%
    left_join(tidy_filt_info, by = c("Key")) %>%
    filter(!is.na(Gene_ID)) %>%
    distinct(Indiv, gt_GT, .keep_all = TRUE) %>%
    mutate("locus" = paste(Gene_ID, "rare", sep = "_"))# %>%
    # mutate(bin_allele = 1)
  
  # Collapse singletons separately
  # TODO: This data structure wont work
  # It's not compatible with vars_multiples for the individualizing-func
  # to work properly
  # XXX Do not use these options
  if(collapse_singletons) {
    lowmaf_vars_singletons_multiples <- filter(tidy_filt_gts, MAC == 1) %>%
      left_join(tidy_filt_info, by = c("Key")) %>%
      filter(!is.na(Gene_ID)) %>%
      distinct(Indiv, gt_GT, .keep_all = TRUE) %>%
      mutate("locus" = paste(Gene_ID, "rare_single", sep = "_")) #%>%
      # mutate(bin_allele = 1)
    
    # If desired, singleton genotypes would appear as Ref in the regular _rare locus
    if(separate_singletons) {
      filter(lowmaf_vars_multiples, gt_GT >= 1 & MAC == 1) %>%
        mutate(gt_GT = 0)
    }
    
    lowmaf_vars_multiples <- bind_rows(lowmaf_vars_multiples, lowmaf_vars_singletons_multiples)
  }
  
  lowmaf_vars <- lowmaf_vars_multiples %>%
    group_by(Gene_ID, Indiv) %>%
    # summarize(rare.all  = max(gt_GT, na.rm = TRUE),
    #           rare.maybe= min(gt_GT, na.rm = TRUE),
    #           .groups = "keep") %>%
    summarize(rare = max(gt_GT, na.rm = TRUE),
              .groups = "keep") %>%
    # TODO: need to check if any loci within a gene have Alt at a locus in goodmaf
    # use case_when with the rare column and check if the Indiv is in goodmaf
    # problem is that goodmaf doesnt have gene_ID, so it's not that simple
    # mutate(rare.maybe = case_when(goodmaf_vars[goodmaf_vars$Indiv == Indiv, bin_allele > 0] )) %>%
    # mutate(rare.maybe = case_when())
    mutate(across(contains("rare"),
                  ~ str_replace(.x, stringr::regex("-?Inf"), "NA")
    )) %>%
    type_convert() %>%
    mutate(across(contains("rare"),
                  ~ case_when(rare >= 1 ~ 1,
                              rare == 0 ~ 0,
                              .default = NA))) %>%
    pivot_longer(contains("rare"),
                 names_to = "rare_kind", values_to = "bin_allele") %>%
    unite("locus", c(Gene_ID, rare_kind))
  
  # bind, munge into weird matrix format we need
  # Example: see scripts/aina_treewas.ppt, or from treeWAS: data("snps)
  final_long_tb <- bind_rows(goodmaf_vars, lowmaf_vars) %>%
    pivot_wider(names_from = locus, values_from = bin_allele)
  
  final_long_matr <- as.matrix(final_long_tb[,-1])
  rownames(final_long_matr) <- final_long_tb$Indiv
  
  return(final_long_matr)
}

# Filter variants down to target loci for GWAS
vcfr_to_gwas_info <- function(vcf, max_indel_size = 3, coding_only = TRUE, ignore_low = FALSE,
                              collapse_freq = 0.01) {
  tidy_filt_info <- tidy_snpeff(vcf)
  
  tidy_filt_info <- tidy_filt_info %>%
    mutate(Alt_len = stringr::str_length(ALT)) %>%
    mutate(Ref_len = stringr::str_length(REF)) %>%
    filter(abs(Alt_len - Ref_len) <= max_indel_size)
  
  if(coding_only) {
    tidy_filt_info <- tidy_filt_info %>%
      filter(Annotation_Impact != "MODIFIER")
  }
  if(ignore_low) {
    tidy_filt_info <- tidy_filt_info %>%
      filter(Annotation_Impact != "LOW")
  }
  
  return(tidy_filt_info)
}

# Filter GTs to only the locations in info
vcfr_to_gwas_gts <- function(vcf, tidy_filt_info, max_prop_amb_samples = 0.2, plot_prop_amb = FALSE) {
  tidy_filt_gts <- extract_gt_tidy(vcf, verbose = FALSE) %>%
    filter(Key %in% tidy_filt_info$Key)
  
  tidy_amb_samp_distr <- tidy_filt_gts %>%
    group_by(Key) %>%
    summarize(prop_missing = sum(is.na(gt_GT))/n())
  
  tidy_keep_vars <- tidy_amb_samp_distr %>%
    filter(prop_missing < max_prop_amb_samples )
    
  tidy_filt_gts <- tidy_filt_gts %>%
    filter(Key %in% tidy_keep_vars$Key)
  
  if(plot_prop_amb) {
    tidy_amb_samp_distr_gg <- ggplot(tidy_amb_samp_distr) +
      geom_histogram(aes(prop_missing), binwidth = 0.005) +
      labs(x = "Proportion of samples with missing GT", y = "Variant count")
    
    print(tidy_amb_samp_distr_gg)
  }
  
  return(tidy_filt_gts)
}

# Reads and parses interproscan TSV file
read_iprscan <- function(filename, trim_cols = TRUE, run_GO = TRUE, run_pathways = FALSE,
                         split_GO = FALSE, convert_hyphens = FALSE) {
  # column descriptions from: https://interproscan-docs.readthedocs.io/en/latest/OutputFormats.html
  ipr_colnames <- c("prot_name", "prot_md5", "prot_length",
                    "database", "domain_acc", "domain_desc",
                    "domain_start", "domain_stop", "domain_score",
                    "match_status", "run_date",
                    "ipr_acc", "ipr_desc", "GO_terms", "pathways")

  # for some reason, read_tsv doesnt like the last column being go terms,
  # even though excel picks it up just fine.
  # Oh this is probably because most lines have 13 columns,
  # and only lines with interpro domains have 14 cols
  # The above debugging was only true for ipr 5.61. In 5.65 the tabs are correct
  # And theres even a last col that's empty, I suspect for Pathways data if it was present
  # if(run_GO) {ipr_colnames <- c(ipr_colnames, "GO_terms")}
  # XXX Idk how my workaround will impact Pathways parsing since I dont use that interproscan output
  # it may get munged together with GO terms
  # if(run_pathways) {ipr_colnames <- c(ipr_colnames, "pathways")}
  ipr <- readr::read_tsv(filename, col_names = ipr_colnames,
                         show_col_types = FALSE, na = c("-", "", "NA"))
  if(!run_GO) {ipr <- select(ipr, -GO_terms)}
  if(!run_pathways) {ipr <- select(ipr, -pathways)}
  # if(run_GO) {
  #   ipr <- tidyr::separate_wider_delim(ipr, ipr_desc,
  #                               delim = stringr::regex("\\s(?=(GO)|-)"),
  #                               names = c("ipr_desc", "GO_terms"),
  #                               too_few = "align_start", too_many = "merge")
  # }
  
  # column names that are not that useful and clutter the object
  # if removing match_status, for accuracy of results first filter to matches == TRUE
  if(trim_cols) {
    ipr <- ipr %>%
      dplyr::filter(match_status) %>%
      dplyr::select(-c(prot_md5, match_status, run_date))
  }
  
  # Interproscan outputs hyphens as missing data, convert to NA
  # Not necessary if columns are correctly parsed, as when using iprscan >=5.66
  if(convert_hyphens) {
    ipr <- dplyr::mutate(ipr, dplyr::across(tidyselect::where(is.character), ~ dplyr::na_if(.x, "-"))) %>%
      readr::type_convert()
  }
  
  # unique line for every GO term
  if(split_GO) {
    if(!run_GO) {
      print("Cannot split GO terms that were not parsed.
            If you did run interproscan with GO terms, set function parameter run_GO = TRUE")
    } else { ipr <- split_ipr_goterms(ipr) }
  }
  
  return(ipr)
}

# Input interproscan tsv parsed with at least the beginning of read_iprscan()
# Outputs the same dataframe but separated by GO terms
split_ipr_goterms <- function(x, GO_col = "GO_terms") {
  x_split <- tidyr::separate_longer_delim(x, tidyselect::all_of(GO_col), delim = "|") %>%
    tidyr::separate_wider_regex(GO_terms, patterns = c(GO_terms = "GO:\\d+", "\\(",
                                                       GO_database = "\\w+", "\\)"))
  
  return(x_split)
}

# Parses output from bcftools stats runs
# Designed around runs of exactly 2 input vcf files
# TODO: Implement for a run of 1 vcfs
# TODO: Implement for runs with any n vcfs
parse_bcf_stats_compare_1iso <- function(filename, return_DP_distr = FALSE) {
  ID_cols <- c("file1", "file2")
  SN_cols <- c("key", "value")
  TSTV_cols <- c("ts", "tv", "ts_tv", "ts_alt1", "tv_alt1", "ts_tv_alt1")
  DP_cols <- c("bin", "n_genotypes", "frac_genotypes", "n_sites", "frac_sites")
  
  full_df <- read_tsv(filename, col_names = c("category", "ID", "data"),
                      id = "bcf_filename", comment = "#", show_col_types = FALSE)
  
  ID_df <- filter(full_df, category == "ID") %>%
    separate_wider_delim(data, "\t", names = ID_cols, too_few = "align_start") %>%
    unite("files", file1, file2, sep = "-", remove = FALSE) %>%
    mutate(records_in = case_when(grepl("^vcf_links", files) ~ "Unique_truth",
                                  grepl("^split_vcfs", files) & grepl("-NA$", files) ~ "Unique_tool",
                                  .default = "Shared"))
  ID_dict <- select(ID_df, ID, records_in)
  
  # Any section I parse including type_convert() produces a ton of messages.
  # Wrapping it in the middle of the pipes caused an error, so I split up any pipe as needed.
  SN_df <- filter(full_df, category == "SN") %>%
    left_join(ID_dict, by = "ID") %>%
    separate_wider_delim(data, "\t", names = SN_cols)
  SN_df <- suppressMessages(type_convert(SN_df)) %>%
    filter(grepl("SNPs|indels", key)) %>%
    mutate(value = na_if(value, 0)) %>%
    separate_wider_delim(key, delim = " ", names = c(NA, NA, "variant_type")) %>%
    mutate(variant_type = str_remove(variant_type, ":")) %>%
    pivot_wider(names_from = variant_type, values_from = value)
  
  TSTV_df <- filter(full_df, category == "TSTV") %>%
    left_join(ID_dict, by = "ID") %>%
    separate_wider_delim(data, "\t", names = TSTV_cols)
  TSTV_df <- suppressMessages(type_convert(TSTV_df)) %>%
    select(!contains("_alt1"))
  
  # Assumes default parameters, which is that the highest DP category is ">500"
  # All of those will be shrunk to simply "501", potentially dangerous if ceiling used is higher
  # I *wanted* to use Inf but that's actually not statistically valid, better to shrink to max+1
  DP_df <- filter(full_df, category == "DP") %>%
    left_join(ID_dict, by = "ID") %>%
    separate_wider_delim(data, "\t", names = DP_cols) %>%
    filter(records_in == "Unique_tool") %>%
    mutate("bin" = str_replace(bin, ">.*", "501"))
  DP_df <- suppressMessages(type_convert(DP_df))
  DP_summ <- DP_df %>%
    group_by(bcf_filename, ID, records_in) %>%
    summarize(mean_dp = weighted.mean(bin, n_genotypes), .groups = "drop")
  
  bcfstats_df <- left_join(SN_df, TSTV_df, by = c("ID", "bcf_filename", "records_in")) %>%
    left_join(DP_summ, by = c("ID", "bcf_filename", "records_in")) %>%
    mutate("isolate" = str_extract(bcf_filename, "(?<=\\.)B[[:digit:]]+(?=\\.)")) %>%
    select(-starts_with("category"))
  
  if(return_DP_distr) {
    bcfstats <- list("summary" = bcfstats_df, "DP" = DP_df)
  } else {
    bcfstats <- bcfstats_df
  }
  
  return(bcfstats)
}


# Pair of functions to map 2-fold dilution schemes to observed MIC values
# Generated in cooperation with chatGPT
# Dynamically generate the 2-fold step scheme
# Usage:
# dilution_scheme <- generate_dilution_scheme(min_value, max_value)
generate_dilution_scheme <- function(min_value, max_value) {
  scheme <- numeric()
  value <- min_value
  while (value <= max_value) {
    scheme <- c(scheme, value)
    value <- value * 2
  }
  return(scheme)
}

# Assign step indices with proper handling of `>` values
# values below min_value are assigned a proportional fraction between 0 and 1
# Usage:
# mutate(phen, FLU_MIC_step = assign_step_index(FLU_MIC_val.ineq, dilution_scheme))
assign_step_index <- function(values, scheme, min_value) {
  sapply(values, function(x) {
    if (grepl("^>", x)) {
      # Handle fake-numeric values prefixed with `>`
      numeric_part <- as.numeric(sub("^>", "", x))
      if (!is.na(numeric_part)) {
        return(which(scheme == numeric_part) + 1)  # Increment by 1
      } else {
        stop("Invalid format for >value.")
      }
    } else {
      x <- as.numeric(x)
      if (!is.na(x) && x %in% scheme) {
        return(which(scheme == x))  # Exact match
      } else if (!is.na(x)) {
        if (x < min_value) {
          # Assign fractional step index for values < min_value
          return(round(x / min_value, 2))
        } else {
          # Interpolation for numeric values not in the scheme
          lower <- max(scheme[scheme < x])
          upper <- min(scheme[scheme > x])
          lower_index <- which(scheme == lower)
          upper_index <- which(scheme == upper)
          interpolated <- lower_index + (x - lower) / (upper - lower)
          return(round(interpolated, 2))  # Rounded result
        }
      } else if (is.na(x)) {
        return(NA)
      } else {
        stop("Unrecognized value in observed values.")
      }
    }
  })
}

# Get variant position and SnpEff annotations from Genotype Matrix
# Expected inputs:
# gt_file is the unmodified output of variantMatrixToPlink.pl.
# genes_bed_file is necessary to get the BP-coordinate
# at which collapsed variants of each gene are shown.
read_gt_matr_annots <- function(gt_file, genes_bed_file = NULL) {
  # Read matrix from vcf2matrix collapsing and convert to long format
  gt_col_df <- read_tsv(gt_file, n_max = 1,
                        show_col_types = FALSE)
  gt_colnames <- colnames(gt_col_df)
  gt_annots <- read_tsv(gt_file, skip = 1, n_max = 1, col_names = gt_colnames,
                        show_col_types = FALSE) %>%
    select(-1) %>%
    pivot_longer(everything(),
                 names_to = "variant", values_to = "ANN") %>%
    mutate(plink_index = row_number()) %>%
    separate_wider_delim("variant", delim = "_", names = c("chr", "pos_variable"),
                         too_few = "align_end", too_many = "merge",
                         cols_remove = FALSE) %>%
    mutate(pos_variable_num = as.numeric(pos_variable))
  
  # Parse SnpEff annotations to get variant info
  gt_annots_ann <- parse_ANN(gt_annots, too_few = "align_start") %>%
    group_by(variant) %>%
    mutate(across(c(Allele,HGVS.c, AA_var, HGVS.p,
                    `cDNA.pos / cDNA.length`, `CDS.pos / CDS.length`,
                    `AA.pos / AA.length`, Distance),
                  ~paste(.x, collapse = ','))) %>%
    arrange(plink_index) %>%
    distinct(variant, .keep_all = TRUE) %>%
    ungroup()
  
  # Read in all gene coordinates to add more positional data,
  # join with snpeff annotations
  if(!is.null(genes_bed_file)) {
    ref_genes_bed <- read_tsv(genes_bed_file,
                              col_names = c("chr", "begin", "gene_stop", "Gene_ID")) %>%
      mutate("gene_start" = begin+1) %>%
      select(chr, gene_start, gene_stop, Gene_ID)
    
    gt_annots_ann <- gt_annots_ann %>%
      left_join(ref_genes_bed, by = c("chr", "Gene_ID")) %>%
      mutate(sort_pos = coalesce(pos_variable_num, gene_start))
  } else {
    gt_annots_ann <- mutate(gt_annots_ann, sort_pos = pos_variable_num)
  }
  
  return(gt_annots_ann)
}

# Read GEMMA file
# Select p-value column on which to calculate expected p-value (for QQ plot)
# Extra arguments for read_tsv go to ...
read_gemma <- function(filename, n_indivs, use_p_val = "p_lrt", ...) {
  gemma_in <- read_tsv(filename, show_col_types = FALSE)
  
  # proceed only if Requested P-value is *actually* in the data
  if(! use_p_val %in% colnames(gemma_in)) {
    stop("Requested p-value not found in GEMMA output. Select different value")
  }
  
  # Calculate expected p-value
  gemma_out <- gemma_in %>%
    arrange(!!sym(use_p_val)) %>%
    mutate("exp_pval" = 1:n_indivs/n_indivs) %>%
    # mutate("exp_pval" = 1:37599/37599) %>%
    mutate("negLog10_expected" = -log10(exp_pval)) %>%
    # mutate("negLog10_{{use_p_val}}" := -log10({{ use_p_val }})) %>%
    mutate(across(starts_with("p_", ignore.case = FALSE), ~ -log10(.x), .names = "negLog10_{.col}")) %>%
    arrange(chr, ps)
  return(gemma_out)
}

# Manhattan plot from gemma dataframe (parsed using function above)
# Allows use of any p value in dataset
# Options: "p_wald", "p_lrt", "p_score"
# To get the manhattan plot sorted, associate model results with snpeff annotations.
# XXX expects gt_annots_ann to be present in environment
# ... are extra arguments to geom_text_repel
gemma2manhattan <- function(df, use_p_val = "p_lrt", label_quants = 0.9999, label_percentile = NULL,
                            y_trans = NULL, show_QQ = FALSE, coords = TRUE, return_annots = FALSE,
                            quant_probs = c(0, 0.9, 0.95, 0.99, 0.999, 0.9999, 1.0),...) {
  # proceed only if Requested P-value is *actually* in the data
  if(! use_p_val %in% colnames(df)) {
    stop("Requested p-value not found in GEMMA output. Select different value")
  }
  
  # If no labelling required and just want to show distribution of P-values, return a simple plot
  if(!coords) {
    gg_upsidedown_manhattan <- ggplot(df, aes(x = ps, y = !!sym(use_p_val))) +
      geom_point() +
      scale_y_continuous(trans= "log10") +
      labs(x = "position") +
      theme_bw() +
      theme(panel.grid = element_blank())
    
    return(gg_upsidedown_manhattan)
  }
  # Assign each p value to a quantile of the observed distribution
  # df_quants <- quantile(df$negLog10_p_wald, probs = c(0, 0.9, 0.95, 0.99, 0.999, 0.9999, 1.0))
  df_quants <- quantile(pull(df, !!sym(paste0("negLog10_", use_p_val))), probs = quant_probs)
  df_quants_df <- as_tibble(df_quants) %>%
    mutate(!!glue::glue("negLog10_{use_p_val}_min") := lag(value)) %>%
    mutate(!!glue::glue("negLog10_{use_p_val}_quant") := names(df_quants)) %>%
    mutate(!!glue::glue("negLog10_{use_p_val}_prob") := as.numeric(str_remove(names(df_quants), "%"))*.01) %>%
    rename(!!glue::glue("negLog10_{use_p_val}_max") := value)
  
  negLog10_p_min = paste0("negLog10_", use_p_val, "_min")
  negLog10_p_quant = paste0("negLog10_", use_p_val, "_quant")
  negLog10_p_prob = paste0("negLog10_", use_p_val, "_prob")
  negLog10_p_max = paste0("negLog10_", use_p_val, "_max")
  
  # by <- join_by("negLog10_p_wald" > "negLog10_p_wald_min", "negLog10_p_wald" <= "negLog10_p_wald_max")
  # by <- join_by("negLog10_{{use_p_val}}" > "negLog10_{{use_p_val}}_min", "negLog10_{{use_p_val}}" <= "negLog10_{{use_p_val}}_max")
  by <- join_by(!!sym(paste0("negLog10_", use_p_val)) > !!sym(negLog10_p_min),
                !!sym(paste0("negLog10_", use_p_val)) <= !!sym(negLog10_p_max))
  df_pos <- df %>%
    select(-c(chr, rs)) %>%
    left_join(gt_annots_ann, by = c("ps" = "plink_index")) %>%
    distinct(ps, .keep_all = TRUE) %>%
    arrange(chr, sort_pos, ps) %>%
    mutate(manhattan_pos = row_number()) %>%
    left_join(df_quants_df, by)
  if(return_annots) {
    return(df_pos)
  }
  
  # filter data labels to show - 
  # play around to get quantile threshold desired
  # df_pos_filt <- filter(df_pos, negLog10_p_wald_prob > label_quants)
  df_pos_filt <- filter(df_pos, !!sym(negLog10_p_prob) > label_quants)
  
  # Show horizontal line at chosen threshold
  negLog10_p_sigVal <- unlist(df_quants_df[
    df_quants_df[,negLog10_p_prob] == label_quants, negLog10_p_max])
  # doesn't work with quantile 0.999 so use percentile instead
  if(purrr::is_empty(negLog10_p_sigVal)) {
    if(is.null(label_percentile)) {
      stop("Quantile selection for plotting above significance threshold failed, and no percentile to fall back on.")
    }
    else {
      negLog10_p_sigVal <- unlist(df_quants_df[
        df_quants_df[,negLog10_p_quant] == label_percentile, negLog10_p_max])
    }
  }
  # print(negLog10_p_wald_sigVal)
  
  df_gg <- ggplot(df_pos, aes(x = manhattan_pos, y = !!sym(paste0("negLog10_", use_p_val)))) +
    geom_hline(aes(yintercept = negLog10_p_sigVal), linetype = "dashed", color = "blue4") +
    geom_point(aes(color = chr)) +
    geom_text_repel(data = df_pos_filt,
                    mapping = aes(label = paste(Gene_ID, AA_var, sep = ":")),
                    min.segment.length = 0.1, ...) +
    labs(x = "Relative position") +
    scale_x_continuous(expand = c(0.01,0.01)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none")
  
  if(!is.null(y_trans)) {
    df_gg <- df_gg + scale_y_continuous(expand = c(0.01,0.01), transform = y_trans)
  } else {
    df_gg <- df_gg + scale_y_continuous(expand = c(0.01,0.01))
  }
  
  if(show_QQ) {
    df_qq_gg <- gemma2QQ(df)
    
    df_gg <- list(df_gg, df_qq_gg)
  }
  return(df_gg)
}

# Gemma dataframe to qq plot - mostly used internally by gemma2manhattan
gemma2QQ <- function(df, use_p_val = "p_lrt") {
  # proceed only if Requested P-value is *actually* in the data
  if(! use_p_val %in% colnames(df)) {
    stop("Requested p-value not found in GEMMA output. Select different value")
  }
  ggplot(df, aes(x = negLog10_expected, y = !!sym(paste0("negLog10_", use_p_val)))) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_classic()
}

# Fisher's exact test on any variant within a summarized count matrix
# input format example:
# Variant           count FLU_32_res_named gt_num_named
# <chr>             <int> <chr>            <chr>       
# CP060341.1_215724   135 Flu32.S          Ref.allele  
# CP060341.1_215724   491 Flu32.R          Ref.allele  
# CP060341.1_215724     8 Flu32.S          Alt.allele  
# CP060341.1_215724   228 Flu32.R          Alt.allele  
# Example usage: fisher_outs <- map(check_specific_vars, contingency_fisher) %>% bind_rows() 
contingency_fisher <- function(check_var, x = gt_polymorphic_interest_phen_counts) {
  # filter to only known alleles in the variant being tested
  var_subset <- filter(x, Variant == check_var) %>%
    filter(!grepl("Mis.allele", gt_num_named))
  # contingency tibble
  contingency_tbl <- var_subset %>%
    pivot_wider(names_from = gt_num_named, values_from = count) %>%
    select(-Variant) %>%
    mutate(across(-1, ~ replace_na(.x, 0)))
  # Make contingency df - rownames needed
  tmp_df <- data.frame(contingency_tbl)
  rownames(tmp_df) <- tmp_df[,1]
  contingency_df <- tmp_df[,-1]
  # print(contingency_df)
  fisher_res <- fisher.test(contingency_df, alternative = "two.sided")
  fisher_broom <- broom::tidy(fisher_res) %>%
    mutate(Variant = check_var)
  # Make contingency matrix for fisher test - fails when a value is missing
  # contingency_matr <- matrix(pull(var_subset, count),
  #                            nrow = 2, dimnames = list(
  #                              Pheno = unique(var_subset$FLU_32_res_named),
  #                              Geno = unique(var_subset$gt_num_named)))
  # print(contingency_matr)
  # fisher_out <- fisher.test(contingency_matr, alternative = "two.sided")
  return(fisher_broom)
}

# Convert true/false into simplified gemma input;
# binary phenotypes (0/1) or forced continuous (1/2)
# Expected input: Named boolean vector of phenotypes
phen_logistic2numeric <- function(x, binary = TRUE) {
  # Default is binary output (0/1)
  if(binary) {
    x_converted <- str_replace_all(as.character(x), "FALSE", '0')
    x_converted <- str_replace_all(as.character(x_converted), "TRUE", '1')
    x_converted <- as.numeric(x_converted)
    names(x_converted) <- names(x)
  }
  # Alternative is continous distribution (1/2)
  if(!binary) {
    x_converted <- str_replace_all(as.character(x), "FALSE", '1')
    x_converted <- str_replace_all(as.character(x_converted), "TRUE", '2')
    x_converted <- as.numeric(x_converted)
    names(x_converted) <- names(x)
  }
  
  return(x_converted)
}

# Parse relatedness matrix from GEMMA
# inds: a character vector of isolate names
read_rel_matrix <- function(file, inds, long = TRUE) {
  x <- read_tsv(file, col_names = FALSE, show_col_types = FALSE)
  colnames(x) <- inds
  if(long) {
    x <- x %>% 
      mutate(IID_1 = inds) %>%
      pivot_longer(-IID_1, names_to = "IID_2", values_to = "relatedness") 
    # Filter to non-redundant positions
    # essentially uses a decorate-sort-undecorate pattern
    x <- x %>%
      # filter(IID_1 != IID_2) %>%
      mutate(DECORATE_comp_lg = case_when(IID_1 > IID_2 ~ IID_1,.default = IID_2),
             DECORATE_comp_sm = case_when(IID_1 > IID_2 ~ IID_2,.default = IID_1)) %>%
      distinct(DECORATE_comp_lg, DECORATE_comp_sm, .keep_all = TRUE) %>%
      select(-c(DECORATE_comp_lg, DECORATE_comp_sm))
  }
  
  return(x)
}