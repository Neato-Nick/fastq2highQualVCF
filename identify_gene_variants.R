#!/usr/bin/env Rscript
# Author: Nicholas Cauldron (modified by GitHub Copilot / Neato-Nick)
# Date: 2025-10-29 (refactor 2025-10-31)
# identify_gene_variants.R
#
# FUNCTIONAL refactor of the original script:
# - Breaks large steps into functions so they can be moved into an R package (e.g. gemmar).
# - Preserves original behavior and CLI options.
# - Exposes building blocks:
#     * read_and_tidy_vcf()
#     * parse_required_and_flip_keys()
#     * build_flip_df()
#     * compute_cc_variants()
#     * apply_flips_and_join_cc()
#     * assign_hotspots()
#     * summarize_per_isolate()
#     * enforce_required_positions()
# - Minimal behavior changes: output formatting, NA handling, verbose/debug preserved.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(vcfR)
  if (!requireNamespace("gemmar", quietly = TRUE)) {
    remotes::install_github("Neato-Nick/gemmar", quiet = TRUE)
  }
  library(gemmar)
})

# ---- CLI options ----
option_list <- list(
  make_option(c("-g", "--gene"), type = "character", default = NULL, help = "Gene ID (e.g. 'CJI82_03103')"),
  make_option(c("-v", "--vcf"), type = "character", default = NULL, help = "Path to VCF file (snpEff annotated)."),
  make_option(c("-c", "--clade"), type = "character", default = NULL, help = "Optional: path to clade file (Analysis_ID,Clade)."),
  make_option(c("-s", "--hotspots"), type = "character", default = NULL, help = "Optional: hotspots file with Start_aa,End_aa,Segment."),
  make_option(c("-o", "--out"), type = "character", default = NULL, help = "Output TSV path. Default: <gene>_variants.tsv"),
  make_option(c("-l", "--flip_variants"), type = "character", default = NULL, help = "comma-separated AA positions to flip with respect to reference."),
  make_option(c("-r", "--required_positions"), type = "character", default = NULL, help = "'all' or comma-separated AA positions to require."),
  make_option(c("--verbose"), action = "store_true", default = FALSE, help = "Print verbose debug messages."),
  make_option(c("--debug_isolates"), type = "character", default = NULL, help = "Comma-separated isolates to dump debug info for when --verbose is set.")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$gene) || is.null(opt$vcf)) {
  print_help(opt_parser)
  stop("You must provide --gene and --vcf", call. = FALSE)
}

# ---- convenience ----
`%||%` <- function(a, b) if (is.null(a)) b else a

gene_id <- opt$gene
vcf_path <- opt$vcf
clade_path <- opt$clade
hotspots_path <- opt$hotspots
out_path <- opt$out %||% paste0(gene_id, "_variants.tsv")
req_pos_arg <- opt$required_positions
flip_variants_arg <- opt$flip_variants
verbose <- opt$verbose
debug_isolates_arg <- opt$debug_isolates

vmessage <- function(...) { if (verbose) message(..., appendLF = TRUE) }

vmessage("Running identify_gene_variants.R")
vmessage("  gene: ", gene_id)
vmessage("  vcf : ", vcf_path)
if (!is.null(clade_path)) vmessage("  clade: ", clade_path)
if (!is.null(hotspots_path)) vmessage("  hotspots: ", hotspots_path)
vmessage("  required_positions: ", ifelse(is.null(req_pos_arg), "<none>", req_pos_arg))
vmessage("  out: ", out_path)
if (!is.null(debug_isolates_arg)) vmessage("  debug_isolates: ", debug_isolates_arg)

# ---- Function definitions ----

# Read VCF and tidy snpEff annotations, filter for gene of interest
read_and_tidy_vcf <- function(vcf_file, gene_id, verbose = FALSE) {
  if (verbose) message("Reading VCF...")
  vcf_raw <- read.vcfR(vcf_file, verbose = FALSE)
  vcf_samples <- colnames(vcf_raw@gt)[-1]
  if (verbose) message("  samples in VCF: ", length(vcf_samples), "  variant rows: ", nrow(vcf_raw@fix))
  if (verbose) message("Tidying snpEff INFO...")
  vcf_info_locus_raw <- suppressMessages(tidy_snpeff(vcf_raw))
  vcf_info_locus <- vcf_info_locus_raw %>%
    filter(Annotation_Impact != "MODIFIER") %>%
    filter(Gene_ID == gene_id) %>%
    separate_wider_delim(`AA.pos / AA.length`, delim = "/", names = c("AA_pos", "AA_length"))
  list(vcf_raw = vcf_raw, vcf_info_locus = vcf_info_locus, vcf_samples = vcf_samples)
}

# Parse required_positions and flip_variants args into Keys for the gene based on vcf_info_locus
parse_required_and_flip_keys <- function(vcf_info_locus, req_pos_arg, flip_variants_arg, verbose = FALSE) {
  keys_for_gene <- unique(vcf_info_locus$Key)
  required_keys_for_positions <- integer(0)
  require_all_keys_flag <- FALSE
  if (!is.null(req_pos_arg)) {
    req_trim <- trimws(tolower(req_pos_arg))
    if (req_trim == "all") {
      require_all_keys_flag <- TRUE
      required_keys_for_positions <- keys_for_gene
    } else {
      req_positions <- str_split(req_pos_arg, ",", simplify = TRUE) %>% as.character() %>% str_trim()
      req_positions_num <- suppressWarnings(as.integer(req_positions))
      if (any(is.na(req_positions_num))) stop("required_positions must be 'all' or comma-separated integers like '143,640'")
      mapping <- vcf_info_locus %>% filter(!is.na(AA_pos) & AA_pos %in% req_positions_num) %>% select(Key, AA_pos)
      missing_pos <- setdiff(req_positions_num, mapping$AA_pos)
      if (length(missing_pos) > 0) stop("Requested required_positions not found in annotations for this gene: ", paste(missing_pos, collapse = ","))
      required_keys_for_positions <- unique(mapping$Key)
    }
  }
  flip_keys_for_positions <- integer(0)
  if (!is.null(flip_variants_arg)) {
    flip_variants <- str_split(flip_variants_arg, ",", simplify = TRUE) %>% as.character() %>% str_trim()
    flip_variants_num <- suppressWarnings(as.integer(flip_variants))
    if (any(is.na(flip_variants_num))) stop("flip_variants must not be specified or only contain comma-separated integers like '143,640'")
    mapping <- vcf_info_locus %>% filter(!is.na(AA_pos) & AA_pos %in% flip_variants_num) %>% select(Key, AA_pos)
    missing_pos <- setdiff(flip_variants_num, mapping$AA_pos)
    if (length(missing_pos) > 0) stop("Requested flip_variants not found in annotations for this gene: ", paste(missing_pos, collapse = ","))
    flip_keys_for_positions <- unique(mapping$Key)
  }
  if (require_all_keys_flag && verbose) message("required_positions = all -> requiring genotype calls at ALL Keys for this gene.")
  if (length(required_keys_for_positions) > 0 && verbose) message("Requiring genotype calls at Keys: ", paste(required_keys_for_positions, collapse = ", "))
  if (length(flip_keys_for_positions) > 0 && verbose) message("Flipping genotype calls at Keys: ", paste(flip_keys_for_positions, collapse = ", "))
  list(keys_for_gene = keys_for_gene,
       required_keys_for_positions = required_keys_for_positions,
       require_all_keys_flag = require_all_keys_flag,
       flip_keys_for_positions = flip_keys_for_positions)
}

# Build flip mapping (auto + manual fallback). Returns tibble(Key, AA_pos, original_label, flipped_label)
build_flip_df <- function(vcf_info_tbl, flip_keys_vec, verbose = FALSE) {
  if (length(flip_keys_vec) == 0) return(tibble())
  if (verbose) message("Searching for genotype positions to flip")
  reverse_AAvar <- function(lbl) {
    if (is.na(lbl)) return(NA_character_)
    m <- str_match(lbl, "^([A-Za-z]+)(\\d+)([A-Za-z]+)$")
    if (is.na(m[1,1])) return(NA_character_)
    ref <- m[1,2]; pos <- m[1,3]; alt <- m[1,4]
    paste0(alt, pos, ref)
  }
  auto_map_tbl <- vcf_info_tbl %>%
    filter(Key %in% flip_keys_vec) %>%
    group_by(Key) %>%
    summarize(original_label = first(na.omit(AA_var)), AA_pos = first(na.omit(AA_pos)), .groups = "drop") %>%
    filter(!is.na(original_label)) %>%
    mutate(flipped_label = map_chr(original_label, ~ reverse_AAvar(.x))) %>%
    mutate(flipped_label = if_else(flipped_label == "NA" | is.na(flipped_label), NA_character_, flipped_label))
  # Manual mapping placeholder (keeps behavior of original script)
  target_switches <- tibble("AA_pos" = c(143, 640),
                            "original_label" = c("R143K", "V640A"),
                            "flipped_label" = c("K143R", "A640V"))
  manual_flip_map_df <- left_join(target_switches, vcf_info_tbl %>% select(Key, AA_pos, AA_var) %>% { suppressMessages(readr::type_convert(.)) }, by = "AA_pos") %>% filter(!is.na(Key))
  if (nrow(manual_flip_map_df) > 0) {
    manual_norm <- manual_flip_map_df
    if (!"Key" %in% names(manual_norm) && "AA_pos" %in% names(manual_norm)) {
      manual_norm <- manual_norm %>%
        left_join(select(vcf_info_tbl, Key, AA_pos), by = "AA_pos") %>%
        select(Key, flipped_label)
    } else if (!"flipped_label" %in% names(manual_norm) && "final_variants" %in% names(manual_norm)) {
      manual_norm <- manual_norm %>% rename(flipped_label = final_variants)
    }
    manual_norm <- manual_norm %>% filter(!is.na(Key) & !is.na(flipped_label)) %>%
      mutate(Key = as.integer(Key), flipped_label = as.character(flipped_label))
    auto_map_tbl <- auto_map_tbl %>%
      left_join(manual_norm %>% select(Key, flipped_label) %>% rename(manual_flipped = flipped_label), by = "Key") %>%
      mutate(flipped_label = coalesce(manual_flipped, flipped_label)) %>%
      select(-manual_flipped)
  }
  flip_df <- auto_map_tbl %>%
    select(Key, AA_pos, original_label, flipped_label) %>%
    mutate(Key = as.integer(Key))
  # drop incomplete mappings
  if (nrow(flip_df) > 0) {
    flip_df <- flip_df %>% filter(!is.na(flipped_label) & !is.na(original_label))
    if (verbose) {
      if (nrow(flip_df) == 0) message("build_flip_df: no usable flip mappings detected (manual or auto).")
      else message("build_flip_df: flip mapping produced for Keys: ", paste(flip_df$Key, collapse = ", "))
    }
  }
  flip_df
}

# Extract tidy genotypes from vcf_raw and filter to gene Keys and join annotation info
extract_and_prepare_gts <- function(vcf_raw, vcf_info_locus, keys_for_gene) {
  vcf_gts_raw <- suppressMessages(extract_gt_tidy(vcf_raw))
  if (!all(c("Key","Indiv","gt_GT","gt_GT_alleles") %in% names(vcf_gts_raw))) stop("extract_gt_tidy() must return Key, Indiv, gt_GT, gt_GT_alleles")
  vcf_gts <- vcf_gts_raw %>%
    filter(Key %in% keys_for_gene) %>%
    left_join(select(vcf_info_locus, Key, Annotation, Annotation_Impact, AC, HGVS.c, HGVS.p, Allele, AA_pos, AA_var), by = c("Key", "gt_GT_alleles" = "Allele")) %>%
    mutate(Variant = case_when(
      is.na(gt_GT) ~ NA_character_,
      gt_GT == 0 ~ "WT",
      gt_GT >= 1 ~ AA_var
    ))
  vcf_gts
}

# Compute clade-conserved (CC) variants using per-clade/per-Key unanimity
compute_cc_variants <- function(vcf_gts, clade_df, flip_df = NULL, vcf_info_locus = NULL, verbose = FALSE) {
  if (is.null(clade_df)) return(list(cc_variants = NULL, cc_variants_summ = tibble(Clade = character(), CC_variants = character())))
  if (verbose) message("Computing clade-conserved variants (CC) using per-clade/per-Key unanimity...")
  cc_variants <- select(clade_df, Analysis_ID, Clade) %>%
    right_join(vcf_gts, by = c("Analysis_ID" = "Indiv"))
  if (!is.null(flip_df) && nrow(flip_df) > 0) {
    cc_variants <- cc_variants %>%
      left_join(flip_df %>% select(Key, original_label, flipped_label), by = "Key") %>%
      mutate(
        Variant_oref = case_when(
          !is.na(flipped_label) & Variant == "WT" ~ flipped_label,
          !is.na(original_label) & Variant == original_label ~ "WT",
          TRUE ~ Variant
        )
      ) %>%
      select(-Variant) %>% rename(Variant = Variant_oref)
  }
  cc_variants <- filter(cc_variants, !is.na(gt_GT), !is.na(Clade)) %>%
    group_by(Clade, Key, Variant) %>%
    mutate(n_Clade_Key_Variant_isos = n()) %>%
    group_by(Clade, Key) %>%
    mutate(N_clade_Key_isos = n()) %>%
    filter(n_Clade_Key_Variant_isos == N_clade_Key_isos) %>%
    ungroup() %>%
    distinct(Clade, Variant, .keep_all = FALSE) %>%
    filter(Variant != "WT") %>%
    mutate("Var_CC" = "CC")
  cc_variants_summ <- tibble(Clade = character(), CC_variants = character())
  if (nrow(cc_variants) > 0) {
    cc_variants_summ <- cc_variants %>%
      group_by(Clade) %>%
      summarize(CC_variants = paste(sort(unique(Variant)), collapse = ","), .groups = "drop")
  }
  if (verbose) message("  Clade-conserved entries: ", nrow(cc_variants), "; Across clades: ", nrow(cc_variants_summ))
  list(cc_variants = cc_variants, cc_variants_summ = cc_variants_summ)
}

# Apply flips and join precomputed CC variants; produce vcf_rows ready for hotspot assignment and summarization
apply_flips_and_join_cc <- function(vcf_gts, clade_df = NULL, flip_df = NULL, cc_variants = NULL, verbose = FALSE) {
  vcf_rows <- vcf_gts %>%
    rename(Analysis_ID = Indiv) %>%
    left_join(clade_df %||% tibble(Analysis_ID = character(), Clade = character()), by = "Analysis_ID")
  if (!is.null(flip_df) && nrow(flip_df) > 0) {
    if (verbose) message("Joining gt data to flip_df (found).")
    # bring in original/flipped labels AND the AA_pos from the flip map so we can
    # treat WT genotypes as effective flipped variants for hotspot assignment.
    vcf_rows <- left_join(vcf_rows, flip_df %>% select(Key, original_label, flipped_label, AA_pos) %>% rename(flip_AA_pos = AA_pos), by = "Key")
  } else {
    if (verbose) message("No flip_df found; adding placeholder original_label/flipped_label columns.")
    vcf_rows <- mutate(vcf_rows, original_label = NA_character_, flipped_label = NA_character_)
  }
  if (!is.null(flip_df) && nrow(flip_df) > 0) {
    vcf_rows <- vcf_rows %>%
      mutate(
        Variant_oref = case_when(
          !is.na(flipped_label) & Variant == "WT" ~ flipped_label,
          !is.na(original_label) & Variant == original_label ~ "WT",
          TRUE ~ Variant
        )
      )
  } else {
    vcf_rows <- vcf_rows %>% mutate(Variant_oref = Variant)
  }
  # --- NEW: make flipped-aware AA_pos + Variant for downstream hotspot detection ---
  # If the VCF shows WT at a position but a flipped_label exists for the Key,
  # then treat this row as having the flipped amino-acid (and the AA_pos from the flip map)
  # for the purposes of hotspot assignment. We keep Variant (reported) as Variant_oref
  # so other logic (Variants, Variants_CC) continues to behave as before.
  vcf_rows <- vcf_rows %>%
    mutate(
      Variant_for_hotspot = if_else(!is.na(flipped_label) & !is.na(gt_GT) & gt_GT == 0, flipped_label, Variant_oref),
      AA_pos_for_hotspot = if_else(!is.na(flipped_label) & !is.na(gt_GT) & gt_GT == 0, flip_AA_pos, AA_pos)
    )
  if (!is.null(cc_variants) && nrow(cc_variants) > 0) {
    if (verbose) message("Joining gt data to pre-computed clade variants.")
    vcf_rows <- vcf_rows %>% left_join(., cc_variants %>% select(Clade, Variant, Var_CC), by = c("Clade", "Variant_oref" = "Variant"))
  } else {
    if (verbose) message("No pre-computed clade variants found; skipping join")
    vcf_rows <- vcf_rows %>% mutate(Var_CC = NA_character_)
  }
  vcf_rows <- vcf_rows %>% mutate(
    is_cc = !is.na(Var_CC),
    Variant = Variant_oref,
    Variant_CC = if_else(is_cc, "CC", Variant)
  ) %>%
    select(Key, Analysis_ID, Clade, gt_GT, Variant, Variant_CC, is_cc, AA_pos, original_label, flipped_label, Variant_for_hotspot, AA_pos_for_hotspot) %>%
    { suppressMessages(readr::type_convert(.)) }
  vcf_rows
}

# Assign hotspots per row based on AA_pos and hotspots table
assign_hotspots <- function(vcf_rows, hotspots = NULL, verbose = FALSE) {
  if (is.null(hotspots)) return(vcf_rows)
  if (verbose) message("Identifying mutational hot-spot regions.")
  # join hotspots using the flipped-aware AA_pos_for_hotspot column (falls back to AA_pos)
  vcf_rows <- left_join(
    vcf_rows %>% mutate(.AA_pos_for_hs = coalesce(AA_pos_for_hotspot, AA_pos)),
    hotspots %>% mutate(Start_aa = as.integer(Start_aa), End_aa = as.integer(End_aa)),
    by = join_by(.AA_pos_for_hs >= Start_aa, .AA_pos_for_hs <= End_aa)
  ) %>%
    # decide whether this row is effectively non-WT for hotspot purposes using the flipped-aware Variant_for_hotspot
    mutate(
      has_effective_nonWT = case_when(
        !is.na(Variant_for_hotspot) & Variant_for_hotspot != "WT" ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    mutate(
      # use the flipped-aware AA/Variant when deciding hotspot membership:
      # - if no segment and no AA_pos info -> WT (called but outside annotated AA positions)
      # - if no segment but AA_pos present and effective nonWT -> NonHS
      # - if Variant is NA -> NA (missing genotype/annotation)
      # - otherwise the hotspot Segment (HS label)
      Hotspot = case_when(
        is.na(Segment) & is.na(.AA_pos_for_hs) & !is.na(gt_GT) ~ "WT",
        is.na(Segment) & !is.na(.AA_pos_for_hs) & !is.na(gt_GT) & has_effective_nonWT ~ "NonHS",
        is.na(Variant) ~ NA_character_,
        .default = Segment
      )
    ) %>%
    # Hotspot_CC should reflect hotspot for non-CC rows but preserve NA for rows that are CC
    mutate(Hotspot_CC = case_when(
      Variant_CC == "CC" ~ NA_character_,
      TRUE ~ Hotspot
    )) %>%
    select(-.AA_pos_for_hs, -has_effective_nonWT)
  vcf_rows
}

# Summarize per-isolate values with the CC/Hotspot logic preserved
summarize_per_isolate <- function(vcf_rows, vcf_samples, clade_df = NULL, cc_variants_summ = NULL, verbose = FALSE) {
  isolate_list <- sort(unique(c(vcf_samples, vcf_rows$Analysis_ID, if (!is.null(clade_df)) clade_df$Analysis_ID else character())))
  summaries <- map_dfr(isolate_list, function(iso) {
    rows_iso <- filter(vcf_rows, Analysis_ID == iso)
    n_called <- sum(!is.na(rows_iso$gt_GT))
    called_nonwt_variants <- unique(rows_iso$Variant[!is.na(rows_iso$Variant) & rows_iso$Variant != "WT" & !is.na(rows_iso$gt_GT)])
    Variants_val <- if (n_called == 0) {
      NA_character_
    } else if (length(called_nonwt_variants) > 0) {
      paste(called_nonwt_variants, collapse = ",")
    } else {
      "WT"
    }
    # Variants_CC
    if (n_called == 0) {
      Variants_CC_val <- NA_character_
    } else {
      rows_called <- rows_iso %>% filter(!is.na(gt_GT))
      if (nrow(rows_called) == 0) {
        Variants_CC_val <- NA_character_
      } else {
        noncc_nonwt <- unique(rows_called$Variant[!is.na(rows_called$Variant) & rows_called$Variant != "WT" & !rows_called$is_cc])
        has_cc <- any(rows_called$is_cc, na.rm = TRUE)
        if (length(noncc_nonwt) > 0) {
          outvals <- unique(noncc_nonwt)
          if (has_cc) outvals <- c(outvals, "CC")
          Variants_CC_val <- paste(outvals, collapse = ",")
        } else if (has_cc) {
          Variants_CC_val <- "CC"
        } else {
          Variants_CC_val <- "WT"
        }
      }
    }
    # Hotspots
    hs_rows <- rows_iso %>% filter(!is.na(gt_GT) & !is.na(Variant) & Variant != "WT")
    hs_vals <- unique(hs_rows$Hotspot[!is.na(hs_rows$Hotspot)])
    Hotspots_val <- if (n_called == 0) {
      NA_character_
    } else if (length(hs_vals) > 0) {
      paste(hs_vals, collapse = ",")
    } else if (n_called > 0 && all((rows_iso$Variant[!is.na(rows_iso$gt_GT)]) == "WT", na.rm = TRUE)) {
      "WT"
    } else if (any(!is.na(rows_iso$Variant) & rows_iso$Variant != "WT" & !is.na(rows_iso$gt_GT))) {
      "NonHS"
    } else {
      "WT"
    }
    # Hotspots_CC
    # - Use Variant_CC (coalesced) to identify non-CC contributors so we avoid mismatches
    #   between the per-row hotspot assignment and the is_cc flag arising from joins.
    # - Treat rows with NA Hotspot but a present non-CC non-WT call as "NonHS".
    hs_cc_rows <- rows_iso %>%
      filter(!is.na(gt_GT) & !is.na(Variant) & Variant != "WT" &
               coalesce(Variant_CC, Variant) != "CC")
    hs_cc_vals <- unique(coalesce(hs_cc_rows$Hotspot, "NonHS")[!is.na(coalesce(hs_cc_rows$Hotspot, "NonHS"))])
    has_cc <- any(coalesce(rows_iso$is_cc, FALSE), na.rm = TRUE)
    Hotspots_CC_val <- if (n_called == 0) {
      NA_character_
    } else if (length(hs_cc_vals) > 0) {
      vals <- unique(hs_cc_vals)
      if (has_cc) vals <- c(vals, "CC")
      paste(vals, collapse = ",")
    } else if (has_cc) {
      "CC"
    } else {
      if (n_called > 0 && all((rows_iso$Variant[!is.na(rows_iso$gt_GT)]) == "WT", na.rm = TRUE)) {
        "WT"
      } else if (any(!is.na(rows_iso$Variant) & rows_iso$Variant != "WT" & !is.na(rows_iso$gt_GT))) {
        "NonHS"
      } else {
        "WT"
      }
    }
    clade_val <- if (!is.null(clade_df)) {
      clade_df %>% filter(Analysis_ID == iso) %>% pull(Clade) %>% .[1] %||% NA_character_
    } else NA_character_
    CC_variants_val <- NA_character_
    if (!is.null(clade_df) && !is.na(clade_val)) {
      cc_val <- cc_variants_summ %>% filter(Clade == clade_val) %>% pull(CC_variants)
      if (length(cc_val) == 1 && !is.na(cc_val) && cc_val != "") CC_variants_val <- cc_val
    }
    tibble(
      Analysis_ID = iso,
      Clade = clade_val,
      Variants = Variants_val,
      Variants_CC = Variants_CC_val,
      Hotspots = Hotspots_val,
      Hotspots_CC = Hotspots_CC_val,
      CC_variants = CC_variants_val,
      n_called = n_called
    ) %>%
      mutate(across(contains("Hotspot"), ~str_remove(., "WT,"))) %>%
      mutate(across(contains("Hotspot"), ~str_remove(., ",WT")))
  })
  summaries
}

# Enforce required positions -> set outputs to NA for isolates missing required Keys
enforce_required_positions <- function(summaries, vcf_rows, required_keys_for_positions, require_all_keys_flag, keys_for_gene, verbose = FALSE) {
  if (!(require_all_keys_flag || length(required_keys_for_positions) > 0)) return(summaries)
  if (verbose) message("Enforcing required positions constraint -> isolates missing required genotype positions will get NA for variants/hotspots.")
  if (require_all_keys_flag) {
    isolates_missing <- summaries %>% filter(n_called < length(keys_for_gene)) %>% pull(Analysis_ID)
  } else {
    req_df <- vcf_rows %>% filter(Key %in% required_keys_for_positions) %>%
      group_by(Analysis_ID) %>% summarize(n_req_called = sum(!is.na(gt_GT)), n_req = n_distinct(Key), .groups = "drop")
    isolates_missing <- req_df %>% filter(n_req_called < n_req) %>% pull(Analysis_ID)
    isolates_missing <- unique(c(isolates_missing, setdiff(summaries$Analysis_ID, req_df$Analysis_ID)))
  }
  if (length(isolates_missing) > 0) {
    if (verbose) message("Isolates missing required positions: ", paste(head(isolates_missing, 10), collapse = ", "), ifelse(length(isolates_missing) > 10, paste0(" ...(", length(isolates_missing) - 10, " more)"), ""))
    summaries <- summaries %>% mutate(
      Variants = if_else(Analysis_ID %in% isolates_missing, NA_character_, Variants),
      Variants_CC = if_else(Analysis_ID %in% isolates_missing, NA_character_, Variants_CC),
      Hotspots = if_else(Analysis_ID %in% isolates_missing, NA_character_, Hotspots),
      Hotspots_CC = if_else(Analysis_ID %in% isolates_missing, NA_character_, Hotspots_CC)
    )
  } else {
    if (verbose) message("No isolates missing the required positions.")
  }
  summaries
}

# ---- Main execution flow using functions ----

# Read & tidy VCF and grab annotations for the gene
rv <- read_and_tidy_vcf(vcf_path, gene_id, verbose = verbose)
vcf_raw <- rv$vcf_raw
vcf_info_locus <- rv$vcf_info_locus
vcf_samples <- rv$vcf_samples

if (nrow(vcf_info_locus) == 0) stop("No annotations found for gene ", gene_id, call. = FALSE)
vmessage("Found ", nrow(vcf_info_locus), " annotated Keys for gene ", gene_id)

# Read clade/hotspots if provided
clade_df <- NULL
if (!is.null(clade_path)) {
  vmessage("Reading clade file...")
  clade_df <- readr::read_tsv(clade_path, col_types = cols(.default = "c")) %>%
    rename_with(~ make.names(.x)) %>%
    filter(Analysis_ID %in% vcf_samples)
  if (!("Analysis_ID" %in% names(clade_df)) || !("Clade" %in% names(clade_df))) stop("Clade file must contain Analysis_ID and Clade")
  vmessage("  clade rows: ", nrow(clade_df))
}

hotspots <- NULL
if (!is.null(hotspots_path)) {
  vmessage("Reading hotspots file...")
  hotspots <- readr::read_tsv(hotspots_path, col_types = cols(.default = "c")) %>%
    mutate(Start_aa = as.integer(Start_aa), End_aa = as.integer(End_aa))
  if (!all(c("Start_aa","End_aa","Segment") %in% names(hotspots))) stop("Hotspots file must contain Start_aa, End_aa, Segment")
  vmessage("  hotspots segments: ", nrow(hotspots))
}

# Parse required & flip keys
pf <- parse_required_and_flip_keys(vcf_info_locus, req_pos_arg, flip_variants_arg, verbose = verbose)
keys_for_gene <- pf$keys_for_gene
required_keys_for_positions <- pf$required_keys_for_positions
require_all_keys_flag <- pf$require_all_keys_flag
flip_keys_for_positions <- pf$flip_keys_for_positions

# Extract genotypes and prepare gts table
vcf_gts <- extract_and_prepare_gts(vcf_raw, vcf_info_locus, keys_for_gene)
vmessage("vcf_gts rows for this gene: ", nrow(vcf_gts), " distinct samples: ", length(unique(vcf_gts$Indiv)))

# Build flip mapping and perform multi-allelic safety check
flip_df <- build_flip_df(vcf_info_locus, flip_keys_for_positions, verbose = verbose)
if (nrow(flip_df) > 0) {
  cc_unique_variants_check <- distinct(vcf_gts, Key, AA_pos, Variant) %>% filter(!Variant %in% c("WT")) %>% { suppressMessages(readr::type_convert(.)) }
  multi_allelic <- cc_unique_variants_check %>%
    filter(Key %in% flip_df$Key) %>%
    filter(!is.na(AA_pos)) %>%
    group_by(Key) %>%
    summarize(n_alleles = n(), alleles = paste(sort(unique(Variant)), collapse = ","), .groups = "drop") %>%
    filter(n_alleles > 1)
  if (nrow(multi_allelic) > 0) {
    stop("Cannot safely flip Keys that are multi-allelic. Problematic Keys:\n",
         paste0(multi_allelic$Key, " (", multi_allelic$alleles, ")\n"))
  }
}

# Compute CC variants if clade provided
cc_res <- compute_cc_variants(vcf_gts, clade_df, flip_df = flip_df, vcf_info_locus = vcf_info_locus, verbose = verbose)
cc_variants <- cc_res$cc_variants
cc_variants_summ <- cc_res$cc_variants_summ

# Apply flips, join CC info, and produce vcf_rows
vcf_rows <- apply_flips_and_join_cc(vcf_gts, clade_df = clade_df, flip_df = flip_df, cc_variants = cc_variants, verbose = verbose)

# Assign hotspots if provided
vcf_rows <- assign_hotspots(vcf_rows, hotspots = hotspots, verbose = verbose)

# Debug dump for requested isolates if verbose
if (!is.null(debug_isolates_arg)) {
  debug_isolates <- str_split(debug_isolates_arg, ",", simplify = TRUE) %>% as.character() %>% str_trim()
  if (verbose) {
    for (iso in debug_isolates) {
      vmessage("==== DEBUG: vcf_rows for isolate ", iso, " ====")
      print(filter(vcf_rows, Analysis_ID == iso))
      vmessage("==== DEBUG: vcf_info_locus rows for AA positions present in that isolate ====")
      poslist <- filter(vcf_rows, Analysis_ID == iso) %>% pull(AA_pos) %>% unique() %>% sort(na.last = TRUE)
      if (length(poslist) > 0) {
        print(filter(vcf_info_locus, AA_pos %in% poslist))
      } else {
        vmessage("  no AA positions found for this isolate in gene Keys")
      }
      vmessage("==== DEBUG: clade info for isolate ", iso, " ====")
      if (!is.null(clade_df)) print(filter(clade_df, Analysis_ID == iso)) else vmessage("  no clade file provided")
      vmessage("==== END DEBUG for ", iso, " ====")
    }
  }
}

# Summarize per isolate
summaries <- summarize_per_isolate(vcf_rows, vcf_samples, clade_df = clade_df, cc_variants_summ = cc_variants_summ, verbose = verbose)

# Enforce required positions
summaries <- enforce_required_positions(summaries, vcf_rows, required_keys_for_positions, require_all_keys_flag, keys_for_gene, verbose = verbose)

# Final output write: use literal "NA" for missing values
out_tbl <- summaries %>% select(Analysis_ID, any_of("Clade"), Variants, Variants_CC, Hotspots, Hotspots_CC, CC_variants) %>% arrange(Analysis_ID)
if (is.null(clade_df)) out_tbl <- out_tbl %>% select(-any_of(c("Clade","Variants_CC","Hotspots_CC","CC_variants")))
out_tbl <- out_tbl %>% select(-any_of("n_called"))
if (is.null(hotspots_path)) {
  out_tbl <- select(out_tbl, -contains("hotspot"))
}

vmessage("Writing output to: ", out_path, " (NULL/NA values will be written as literal 'NA')")
readr::write_tsv(out_tbl, out_path, na = "NA")
vmessage("Done. Wrote ", nrow(out_tbl), " isolates. Columns: ", paste(names(out_tbl), collapse = ", "))
