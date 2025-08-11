#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import sys
from collections import defaultdict
import argparse # Added for command-line argument parsing

def process_mosdepth_summary(file_path):
    """
    Parses a single mosdepth summary file and returns contig mean depths.
    Excludes the 'total' row, as it represents an aggregate value, not a specific contig.
    """
    try:
        df = pd.read_csv(file_path, sep='\t')
        df = df[df['chrom'] != 'total']
        # Return a Series with 'chrom' as index and 'mean' as values
        return df.set_index('chrom')['mean']
    except Exception as e:
        print(f"Error processing {file_path}: {e}", file=sys.stderr)
        return None

def generate_contig_ploidy_priors(mosdepth_files, output_file_path):
    """
    Generates the contig copy number priors table based on mosdepth summary files.

    Args:
        mosdepth_files (list): A list of paths to the _mosdepth.summary.txt files.
        output_file_path (str): The path for the output contig-ploidy-priors.tsv file.

    Returns:
        tuple: A tuple containing:
            - list: aneuploidy_observations (list of (isolate_name, contig, observed_ploidy))
            - dict: max_ploidy_tracker (dict of {contig: {'max_val': int, 'isolate': str}})
    """

    # GATK's default max-copy-number is 5 [1].
    MAX_PLOIDY = 4

    # A small probability floor (pseudo-count strength) for any ploidy state to ensure non-zero prior [Conversation, 100].
    # This value ensures that even unobserved ploidy states have a non-zero, but very small, probability,
    # consistent with GATK's example prior tables (e.g., 0.000001) [Conversation, 100].
    # Example table here has a min of 0.01: https://gatk.broadinstitute.org/hc/en-us/articles/13832680852635-DetermineGermlineContigPloidy
    #MIN_PROB_FLOOR = 1e-6
    MIN_PROB_FLOOR = 1e-3

    # Dictionary to store observed ploidy counts for each contig
    # Format: {contig_name: {ploidy_0: count, ploidy_1: count, ...}}
    contig_ploidy_counts = defaultdict(lambda: defaultdict(int))
    all_contigs = set()

    # New data structures for additional requests
    aneuploidy_observations = [] # Stores (isolate_name, contig, observed_ploidy) for non-1 ploidy
    max_ploidy_tracker = defaultdict(lambda: {'max_val': -1, 'isolate': None}) # Stores {contig: {'max_val': int, 'isolate': str}}

    print(f"Processing {len(mosdepth_files)} mosdepth summary files...")

    for i, file_path in enumerate(mosdepth_files):
        # Extract isolate name from filename (e.g., "isolate1" from "isolate1_mapQ10_median.B11205.mosdepth.summary.txt")
        #isolate_name = os.path.basename(file_path).split('_mapQ10_median') if '_mapQ10_median' in os.path.basename(file_path) else os.path.basename(file_path).split('.')
        isolate_name_parts = os.path.basename(file_path).split('_mapQ10_median')
        isolate_name = isolate_name_parts[0]
        print(f" Processing isolate {isolate_name} ({i+1}/{len(mosdepth_files)})...")

        contig_depths = process_mosdepth_summary(file_path)

        if contig_depths is None:
            continue

        # For Candida auris, all 7 chromosomes are considered core and haploid [Conversation, 102].
        # Calculate the median depth across all contigs for this isolate to establish baseline coverage.
        # This median coverage defines the expected haploid read depth for that specific isolate [85, 103, Conversation].
        if not contig_depths.empty:
            baseline_coverage = contig_depths.median()
        else:
            print(f" Warning: No contig data found for {isolate_name}. Skipping.", file=sys.stderr)
            continue

        # Ensure baseline coverage is positive to prevent division by zero or invalid calculations.
        if baseline_coverage <= 0:
            print(f" Warning: Baseline coverage for {isolate_name} is zero or negative ({baseline_coverage:.2f}). Skipping to avoid division by zero or invalid relative depths.", file=sys.stderr)
            continue

        for chrom, mean_depth in contig_depths.items():
            all_contigs.add(chrom)

            # Infer observed relative ploidy for each contig in the isolate [Conversation]
            relative_depth = mean_depth / baseline_coverage
            observed_ploidy = int(round(relative_depth))

            # Cap the observed_ploidy at MAX_PLOIDY (GATK's default is 5) and ensure it's not negative [49, 104, Conversation].
            observed_ploidy = max(0, min(observed_ploidy, MAX_PLOIDY))

            # Aggregate ploidy observations across the entire cohort [Conversation]
            contig_ploidy_counts[chrom][observed_ploidy] += 1

            # --- Store aneuploidy observations ---
            # For Candida auris, expected ploidy for core chromosomes is 1 [Conversation, 105].
            # An observed_ploidy different from 1 is considered an aneuploidy.
            # Aneuploidy, such as a whole chromosome duplication of chromosome 5, has been observed
            # in C. auris in vitro evolution experiments, conferring azole resistance [2, 3].
            # However, this mechanism can be unstable after drug pressure relief, which might explain
            # why it's not commonly observed in clinical isolates [4].
            # It's important to note that while aneuploidy can confer fitness benefits under stress,
            # it might also act as an "evolutionary diversion" rather than a "stepping stone" to adaptation,
            # potentially delaying adaptation at the lineage level [5-11].
            # The rate of aneuploidy (mis-segregation/non-disjunction) can be higher under heat stress [12].
            if observed_ploidy != 1:
                aneuploidy_observations.append((isolate_name, chrom, observed_ploidy))

            # --- Track maximum observed ploidy per contig ---
            if observed_ploidy > max_ploidy_tracker[chrom]['max_val']:
                max_ploidy_tracker[chrom]['max_val'] = observed_ploidy
                max_ploidy_tracker[chrom]['isolate'] = isolate_name

    print("\nAggregating observations and calculating prior probabilities...")

    # Prepare data for the output table
    output_data = []

    # Ensure all contigs are sorted for consistent output, which is generally good practice for GATK inputs.
    sorted_contigs = sorted(list(all_contigs))

    # Define the header row for the TSV file.
    header = ['CONTIG_NAME'] + [f'PLOIDY_PRIOR_{p}' for p in range(MAX_PLOIDY + 1)]
    output_data.append(header)

    # Calculate probabilities for each contig.
    for chrom in sorted_contigs:
        row = [chrom]
        # Get the total number of valid samples that contributed data for this specific contig.
        total_observations_for_contig = sum(contig_ploidy_counts[chrom].values())

        # Calculate a 'prior strength' for smoothing. This scales the MIN_PROB_FLOOR
        # based on the number of actual observations for the contig.
        # This approach helps ensure that observed frequencies dominate, while still giving
        # a small non-zero probability to unobserved states [Conversation].
        prior_strength = MIN_PROB_FLOOR * total_observations_for_contig

        # The denominator for normalization includes the total observed counts plus
        # the sum of prior strengths across all possible ploidy states.
        denominator_sum = total_observations_for_contig + (MAX_PLOIDY + 1) * prior_strength

        # Handle the edge case where a contig had no valid observations across the entire cohort.
        if total_observations_for_contig == 0 or denominator_sum == 0:
            # If no data, distribute probability uniformly across all possible ploidy states (0 to MAX_PLOIDY).
            # This is a reasonable default when no empirical data is available.
            uniform_prob = 1.0 / (MAX_PLOIDY + 1)
            for _ in range(MAX_PLOIDY + 1):
                #row.append(f"{uniform_prob:.6f}")
                row.append(f"{uniform_prob:.3f}")
        else:
            # Calculate probabilities for each ploidy state (0 through MAX_PLOIDY).
            for p in range(MAX_PLOIDY + 1):
                count = contig_ploidy_counts[chrom][p]
                # The probability is (observed_count + prior_strength) / (total_observed_count + total_prior_strength_for_all_ploidies).
                probability = (count + prior_strength) / denominator_sum
                # row.append(f"{probability:.6f}") # Format to 6 decimal places as seen in GATK examples.
                row.append(f"{probability:.3f}") # Format to 6 decimal places as seen in GATK examples.

        output_data.append(row)

    # Write the assembled table to the specified output file path.
    print(f"Writing contig ploidy priors to {output_file_path}...")
    with open(output_file_path, 'w') as f:
        for row in output_data:
            f.write('\t'.join(row) + '\n')

    print(f"Contig ploidy priors table generated successfully at {output_file_path}.")

    return aneuploidy_observations, max_ploidy_tracker

# Main execution block: This part runs when the script is executed from the command line.
if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description=(
            "Generates a contig copy number priors table for GATK's DetermineGermlineContigPloidy tool.\n"
            "It processes mosdepth summary files, calculates a baseline coverage per isolate, "
            "infers contig ploidy, and aggregates observations across the cohort to determine "
            "prior probabilities for each contig. A small pseudo-count is applied to ensure "
            "all plausible ploidy states have a non-zero prior probability. "
            "It also reports observed aneuploidies and maximum ploidy per contig for sanity checks."
        ),
        formatter_class=argparse.RawTextHelpFormatter # For multiline description
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path for the output contig-ploidy-priors.tsv file."
    )
    parser.add_argument(
        "-i", "--input-list",
        required=True,
        help="Path to a file containing a list of mosdepth summary file paths, one filename per line."
    )

    # Parse arguments
    args = parser.parse_args()

    output_tsv = args.output
    input_file_list_path = args.input_list

    # Read mosdepth summary file paths from the input list file
    mosdepth_files_to_process = []
    try:
        with open(input_file_list_path, 'r') as f:
            for line in f:
                # Remove leading/trailing whitespace and add to list
                file_path = line.strip()
                if file_path: # Ensure line is not empty
                    mosdepth_files_to_process.append(file_path)
    except FileNotFoundError:
        print(f"Error: Input file list '{input_file_list_path}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading input file list '{input_file_list_path}': {e}", file=sys.stderr)
        sys.exit(1)

    if not mosdepth_files_to_process:
        print("Error: No mosdepth summary files found in the input list file. Exiting.", file=sys.stderr)
        sys.exit(1)

    # Call the main function to generate the priors table and get the additional data.
    aneuploidies, max_ploidies = generate_contig_ploidy_priors(mosdepth_files_to_process, output_tsv)

    # --- Print Aneuploidy Table ---
    if aneuploidies:
        print("\n--- Isolates with Aneuploidies (Observed Ploidy != 1) ---")
        print(f"{'Isolate':<20} {'Contig':<10} {'Observed Ploidy':<15}")
        print("-" * 45)
        for iso, contig, ploidy in sorted(aneuploidies): # Sort for consistent output
            print(f"{iso:<20} {contig:<10} {ploidy:<15}")
    else:
        print("\n--- No Aneuploidies (Observed Ploidy != 1) Detected ---")
        print("All contigs in all isolates were observed as ploidy 1.")

    # --- Print Max Ploidy Observed per Contig ---
    if max_ploidies:
        print("\n--- Maximum Ploidy Observed per Contig ---")
        print(f"{'Contig':<10} {'Max Ploidy':<12} {'Isolate (with Max)':<30}")
        print("-" * 52)
        # Sort by contig name for consistent output
        for contig in sorted(max_ploidies.keys()):
            max_val = max_ploidies[contig]['max_val']
            isolate_with_max = max_ploidies[contig]['isolate']
            print(f"{contig:<10} {max_val:<12} {isolate_with_max:<30}")
    else:
        print("\n--- No Contigs Processed to Determine Max Ploidy ---")
