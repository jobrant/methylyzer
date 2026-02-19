#!/usr/bin/env python3

"""
bsDiagnostic.py - Diagnostic tool for identifying problematic reads

Identifies "dark phase" reads that show no methylation signal at either
endogenous (HCG/CG) or accessibility (GCH/GC) sites. These reads may indicate:
- Bisulfite over-conversion
- Failed enzyme treatment
- Contaminating DNA
- Alignment artifacts

Outputs:
- List of flagged read IDs for extraction
- Per-read metrics TSV for further analysis
- Summary statistics

Author: Jason Orr Brant, 2026
"""

import argparse
import os
import sys
import glob
import re
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


def get_cytosine_offset(site):
    """Determine the cytosine position (0-indexed) in a methylation site pattern."""
    site_upper = site.upper()
    if 'C' not in site_upper:
        raise ValueError(f"No cytosine found in site pattern: {site}")
    return site_upper.index('C')


def find_site_positions(sequence, site):
    """Find all positions where a methylation site occurs in the sequence."""
    positions = []
    seq_str = str(sequence).upper()
    
    # Handle IUPAC ambiguity codes
    site_pattern = site.upper()
    site_pattern = site_pattern.replace('H', '[ACT]')
    site_pattern = site_pattern.replace('W', '[AT]')
    site_pattern = site_pattern.replace('S', '[CG]')
    site_pattern = site_pattern.replace('M', '[AC]')
    site_pattern = site_pattern.replace('K', '[GT]')
    site_pattern = site_pattern.replace('R', '[AG]')
    site_pattern = site_pattern.replace('Y', '[CT]')
    site_pattern = site_pattern.replace('B', '[CGT]')
    site_pattern = site_pattern.replace('D', '[AGT]')
    site_pattern = site_pattern.replace('V', '[ACG]')
    site_pattern = site_pattern.replace('N', '[ACGT]')
    
    site_regex = re.compile(site_pattern)
    
    for match in site_regex.finditer(seq_str):
        positions.append(match.start())
    
    return positions


def analyze_read(read_seq, ref_seq, site_positions, c_offset):
    """
    Analyze methylation status at site positions for a single read.
    
    Returns dict with:
        - methylated: count of C (methylated)
        - unmethylated: count of T (unmethylated)
        - wrong_base: count of A/G/other (potential shift)
        - no_data: count of gaps/N
        - total: total sites analyzed
    """
    stats = {
        'methylated': 0,
        'unmethylated': 0,
        'wrong_base': 0,
        'no_data': 0,
        'total': 0
    }
    
    read_str = str(read_seq).upper()
    ref_str = str(ref_seq).upper()
    
    for start_pos in site_positions:
        c_pos = start_pos + c_offset
        
        if c_pos >= len(read_str) or c_pos >= len(ref_str):
            stats['no_data'] += 1
            stats['total'] += 1
            continue
        
        ref_base = ref_str[c_pos]
        read_base = read_str[c_pos]
        
        if ref_base != 'C':
            continue  # Not actually a cytosine position
        
        stats['total'] += 1
        
        if read_base in {'-', 'N'}:
            stats['no_data'] += 1
        elif read_base == 'C':
            stats['methylated'] += 1
        elif read_base == 'T':
            stats['unmethylated'] += 1
        else:
            stats['wrong_base'] += 1
    
    return stats


def calculate_conversion_at_other_c(read_seq, ref_seq, exclude_positions):
    """
    Calculate bisulfite conversion rate at cytosines NOT in monitored sites.
    High conversion (>95%) is expected; very high (>99%) might indicate over-conversion.
    
    exclude_positions: set of reference positions to skip (the monitored sites)
    """
    read_str = str(read_seq).upper()
    ref_str = str(ref_seq).upper()
    
    c_count = 0
    t_count = 0
    
    for i, ref_base in enumerate(ref_str):
        if ref_base != 'C':
            continue
        if i in exclude_positions:
            continue
        if i >= len(read_str):
            continue
        
        read_base = read_str[i]
        if read_base == 'C':
            c_count += 1
        elif read_base == 'T':
            t_count += 1
    
    total = c_count + t_count
    if total == 0:
        return None
    
    return t_count / total * 100


def process_fasta(fasta_file, left_site='HCG', right_site='GCH', 
                  min_fold_length=50, fold_identity=80):
    """
    Process a single FASTA file and return per-read metrics.
    
    Returns list of dicts, one per read.
    """
    records = list(SeqIO.parse(fasta_file, 'fasta'))
    
    if len(records) < 2:
        print(f"  Warning: {fasta_file} has fewer than 2 sequences, skipping")
        return []
    
    ref_seq = records[0].seq
    reads = records[1:]
    
    # Find site positions
    left_positions = find_site_positions(ref_seq, left_site)
    right_positions = find_site_positions(ref_seq, right_site)
    
    left_offset = get_cytosine_offset(left_site)
    right_offset = get_cytosine_offset(right_site)
    
    # Build set of all monitored cytosine positions (for conversion rate calc)
    monitored_positions = set()
    for pos in left_positions:
        monitored_positions.add(pos + left_offset)
    for pos in right_positions:
        monitored_positions.add(pos + right_offset)
    
    print(f"  Reference: {len(ref_seq)} bp")
    print(f"  {left_site} sites: {len(left_positions)}")
    print(f"  {right_site} sites: {len(right_positions)}")
    print(f"  Reads: {len(reads)}")
    
    results = []
    
    for read in reads:
        left_stats = analyze_read(read.seq, ref_seq, left_positions, left_offset)
        right_stats = analyze_read(read.seq, ref_seq, right_positions, right_offset)
        
        # Analyze truncation
        trunc_stats = analyze_truncation(read.seq, ref_seq)
        
        # Detect folded reads
        fold_stats = detect_fold(read.seq, ref_seq, min_fold_length, fold_identity)
        
        # Calculate methylation percentages (of informative sites only)
        left_informative = left_stats['methylated'] + left_stats['unmethylated']
        right_informative = right_stats['methylated'] + right_stats['unmethylated']
        
        left_meth_pct = (left_stats['methylated'] / left_informative * 100) if left_informative > 0 else None
        right_meth_pct = (right_stats['methylated'] / right_informative * 100) if right_informative > 0 else None
        
        # Calculate conversion rate at non-site cytosines
        other_c_conversion = calculate_conversion_at_other_c(read.seq, ref_seq, monitored_positions)
        
        # Calculate read coverage (ungapped length / ref length)
        ungapped_len = len(str(read.seq).replace('-', ''))
        coverage_pct = ungapped_len / len(ref_seq) * 100
        
        results.append({
            'read_id': read.id,
            'fasta_file': os.path.basename(fasta_file),
            
            # Left site (endogenous) stats
            'left_methylated': left_stats['methylated'],
            'left_unmethylated': left_stats['unmethylated'],
            'left_wrong_base': left_stats['wrong_base'],
            'left_no_data': left_stats['no_data'],
            'left_total': left_stats['total'],
            'left_meth_pct': left_meth_pct,
            
            # Right site (accessibility) stats
            'right_methylated': right_stats['methylated'],
            'right_unmethylated': right_stats['unmethylated'],
            'right_wrong_base': right_stats['wrong_base'],
            'right_no_data': right_stats['no_data'],
            'right_total': right_stats['total'],
            'right_meth_pct': right_meth_pct,
            
            # Truncation metrics
            'first_aligned': trunc_stats['first_aligned'],
            'last_aligned': trunc_stats['last_aligned'],
            'five_prime_missing': trunc_stats['five_prime_missing'],
            'three_prime_missing': trunc_stats['three_prime_missing'],
            'truncation_type': trunc_stats['truncation_type'],
            
            # Fold metrics
            'is_folded': fold_stats['is_folded'],
            'fold_point_ref': fold_stats['fold_point_ref'],
            'fold_point_read': fold_stats['fold_point_read'],
            'fold_identity': fold_stats['fold_identity'],
            'fold_length': fold_stats['fold_length'],
            
            # Other metrics
            'other_c_conversion': other_c_conversion,
            'coverage_pct': coverage_pct,
        })
    
    return results


def classify_read(row, dark_threshold=10, shift_threshold=20, truncation_threshold=90):
    """
    Classify a read based on its metrics.
    
    dark_threshold: max % methylation at both sites to be called "dark phase"
    shift_threshold: min % wrong bases to be called "shifted"
    truncation_threshold: min % coverage to NOT be called "truncated"
    
    Returns classification string.
    """
    left_meth = row['left_meth_pct']
    right_meth = row['right_meth_pct']
    left_wrong = row['left_wrong_base']
    right_wrong = row['right_wrong_base']
    left_total = row['left_total']
    right_total = row['right_total']
    coverage = row['coverage_pct']
    is_folded = row.get('is_folded', False)
    
    # Calculate wrong base percentages
    left_wrong_pct = (left_wrong / left_total * 100) if left_total > 0 else 0
    right_wrong_pct = (right_wrong / right_total * 100) if right_total > 0 else 0
    
    # Check for truncated reads first (can combine with other classifications)
    is_truncated = coverage < truncation_threshold
    
    # Check for folded reads first - these are a distinct artifact
    if is_folded:
        return 'FOLDED'
    
    # Check for shifted reads
    if left_wrong_pct >= shift_threshold or right_wrong_pct >= shift_threshold:
        return 'TRUNCATED_SHIFTED' if is_truncated else 'SHIFTED'
    
    # Check for dark phase (low methylation at BOTH sites)
    if left_meth is not None and right_meth is not None:
        if left_meth <= dark_threshold and right_meth <= dark_threshold:
            return 'TRUNCATED_DARK' if is_truncated else 'DARK_PHASE'
    
    # Check for partial dark (only one site type affected)
    if left_meth is not None and left_meth <= dark_threshold:
        return 'TRUNCATED' if is_truncated else 'LOW_ENDOGENOUS'
    if right_meth is not None and right_meth <= dark_threshold:
        return 'TRUNCATED' if is_truncated else 'LOW_ACCESSIBILITY'
    
    # Truncated but otherwise normal
    if is_truncated:
        return 'TRUNCATED'
    
    return 'NORMAL'


def analyze_truncation(read_seq, ref_seq):
    """
    Analyze where truncation occurs in a read.
    
    Returns dict with:
        - first_aligned: first position with non-gap base (0-indexed)
        - last_aligned: last position with non-gap base (0-indexed)
        - five_prime_missing: bases missing from 5' end
        - three_prime_missing: bases missing from 3' end
        - truncation_type: '5prime', '3prime', 'both', or 'internal'
    """
    read_str = str(read_seq).upper()
    ref_len = len(ref_seq)
    
    # Find first and last non-gap positions
    first_aligned = None
    last_aligned = None
    
    for i, base in enumerate(read_str):
        if base not in {'-', 'N'}:
            if first_aligned is None:
                first_aligned = i
            last_aligned = i
    
    if first_aligned is None:
        # All gaps
        return {
            'first_aligned': None,
            'last_aligned': None,
            'five_prime_missing': ref_len,
            'three_prime_missing': ref_len,
            'truncation_type': 'empty'
        }
    
    five_prime_missing = first_aligned
    three_prime_missing = ref_len - last_aligned - 1 if last_aligned < ref_len else 0
    
    # Classify truncation type
    if five_prime_missing > 10 and three_prime_missing > 10:
        truncation_type = 'both'
    elif five_prime_missing > 10:
        truncation_type = '5prime'
    elif three_prime_missing > 10:
        truncation_type = '3prime'
    else:
        truncation_type = 'none'
    
    return {
        'first_aligned': first_aligned,
        'last_aligned': last_aligned,
        'five_prime_missing': five_prime_missing,
        'three_prime_missing': three_prime_missing,
        'truncation_type': truncation_type
    }


def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                  'N': 'N', '-': '-', 'R': 'Y', 'Y': 'R',
                  'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq.upper()))


def calculate_identity(seq1, seq2):
    """Calculate percent identity between two sequences of equal length."""
    if len(seq1) != len(seq2) or len(seq1) == 0:
        return 0
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a not in {'-', 'N'})
    total = sum(1 for a, b in zip(seq1, seq2) if a not in {'-', 'N'} and b not in {'-', 'N'})
    
    return (matches / total * 100) if total > 0 else 0


def detect_fold(read_seq, ref_seq, min_fold_length=50, identity_threshold=80):
    """
    Detect if a read appears to be "folded" - where the 3' portion is the 
    reverse complement of the 5' portion.
    
    This can happen due to:
    - Incomplete hairpin adapter removal in PacBio libraries
    - Intramolecular template switching
    - Strong secondary structure causing polymerase to loop back
    
    Args:
        read_seq: The aligned read sequence (may contain gaps)
        ref_seq: The reference sequence
        min_fold_length: Minimum length of fold region to consider (bp)
        identity_threshold: Minimum % identity to call a fold
    
    Returns dict with:
        - is_folded: True if fold detected
        - fold_point_ref: Reference position where fold occurs (0-indexed)
        - fold_point_read: Position in ungapped read where fold occurs
        - fold_identity: % identity between 5' and revcomp(3')
        - fold_length: Length of the folded region
        - fold_type: 'hairpin' (folds back on itself) or 'none'
    """
    read_str = str(read_seq).upper()
    ref_len = len(ref_seq)
    
    # Extract ungapped sequence and track reference positions
    ungapped_seq = ''
    ref_positions = []  # Maps ungapped position -> reference position
    
    for i, base in enumerate(read_str):
        if base not in {'-'}:
            ungapped_seq += base
            ref_positions.append(i)
    
    ungapped_len = len(ungapped_seq)
    
    # Need sufficient length to detect fold
    if ungapped_len < min_fold_length * 2:
        return {
            'is_folded': False,
            'fold_point_ref': None,
            'fold_point_read': None,
            'fold_identity': 0,
            'fold_length': 0,
            'fold_type': 'none'
        }
    
    best_fold = {
        'is_folded': False,
        'fold_point_ref': None,
        'fold_point_read': None,
        'fold_identity': 0,
        'fold_length': 0,
        'fold_type': 'none'
    }
    
    # Scan for fold points - check if 3' portion reverse complements to match 5'
    # Start from middle and work outward to find best fold point
    mid_point = ungapped_len // 2
    
    # Try different fold points around the middle
    for offset in range(0, min(mid_point - min_fold_length, 100), 5):
        for fold_point in [mid_point - offset, mid_point + offset]:
            if fold_point < min_fold_length or fold_point > ungapped_len - min_fold_length:
                continue
            
            # Get 5' portion and 3' portion
            five_prime = ungapped_seq[:fold_point]
            three_prime = ungapped_seq[fold_point:]
            
            # Reverse complement the 3' portion
            three_prime_rc = reverse_complement(three_prime)
            
            # Compare: align from the fold point going outward
            # The 3' revcomp should match the 5' portion going backward from fold point
            compare_len = min(len(five_prime), len(three_prime_rc))
            
            if compare_len < min_fold_length:
                continue
            
            # Compare the end of 5' with the start of 3' revcomp
            five_prime_end = five_prime[-compare_len:]
            three_prime_rc_start = three_prime_rc[:compare_len]
            
            identity = calculate_identity(five_prime_end, three_prime_rc_start)
            
            if identity >= identity_threshold and identity > best_fold['fold_identity']:
                best_fold = {
                    'is_folded': True,
                    'fold_point_ref': ref_positions[fold_point] if fold_point < len(ref_positions) else None,
                    'fold_point_read': fold_point,
                    'fold_identity': identity,
                    'fold_length': compare_len,
                    'fold_type': 'hairpin'
                }
    
    return best_fold


def main():
    parser = argparse.ArgumentParser(
        description="Diagnostic tool for identifying problematic methylation reads",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Classifications:
  NORMAL            - Typical methylation patterns, full-length read
  DARK_PHASE        - Low methylation at BOTH endogenous AND accessibility sites
  SHIFTED           - High proportion of wrong bases (A/G) at cytosine positions
  LOW_ENDOGENOUS    - Low endogenous methylation only (accessibility normal)
  LOW_ACCESSIBILITY - Low accessibility only (endogenous normal)
  TRUNCATED         - Read covers <90% of reference (otherwise normal)
  TRUNCATED_DARK    - Truncated AND dark phase
  TRUNCATED_SHIFTED - Truncated AND shifted
  FOLDED            - Read appears folded (3' is revcomp of 5')

Output files:
  <o>_metrics.tsv     - Per-read metrics for all reads
  <o>_flagged.txt     - Read IDs flagged as problematic (for extraction)
  <o>_truncated.txt   - Read IDs of truncated reads only (with truncation details)
  <o>_folded.txt      - Read IDs of folded reads (with fold point details)
  <o>_summary.txt     - Summary statistics

Examples:
  # Analyze extracted FASTA files
  python bsDiagnostic.py extracted/ -o diagnostic_results
  
  # Custom sites and thresholds
  python bsDiagnostic.py extracted/ -o results --left-site CG --right-site GC --dark-threshold 5
  
  # Stricter truncation threshold (95% coverage required)
  python bsDiagnostic.py extracted/ -o results --truncation-threshold 95
  
  # Adjust fold detection sensitivity
  python bsDiagnostic.py extracted/ -o results --fold-identity 75 --min-fold-length 40
  
  # Extract flagged reads for further analysis
  seqtk subseq original.fa diagnostic_results_flagged.txt > flagged_reads.fa
        """
    )
    
    parser.add_argument('input_path',
                       help='Directory containing extracted FASTA files or single FASTA file')
    parser.add_argument('-o', '--output', default='diagnostic',
                       help='Output prefix for result files (default: diagnostic)')
    parser.add_argument('--left-site', default='HCG',
                       help='Endogenous methylation site pattern (default: HCG)')
    parser.add_argument('--right-site', default='GCH',
                       help='Accessibility site pattern (default: GCH)')
    parser.add_argument('--dark-threshold', type=float, default=10,
                       help='Max %% methylation to classify as dark phase (default: 10)')
    parser.add_argument('--shift-threshold', type=float, default=20,
                       help='Min %% wrong bases to classify as shifted (default: 20)')
    parser.add_argument('--truncation-threshold', type=float, default=90,
                       help='Min %% coverage to NOT be classified as truncated (default: 90)')
    parser.add_argument('--fold-identity', type=float, default=80,
                       help='Min %% identity to detect a folded read (default: 80)')
    parser.add_argument('--min-fold-length', type=int, default=50,
                       help='Min length (bp) of fold region to detect (default: 50)')
    parser.add_argument('--strand', default='ab',
                       help='Which strands to analyze: a, b, or ab (default: ab)')
    
    args = parser.parse_args()
    
    # Find input files
    if os.path.isfile(args.input_path):
        fasta_files = [args.input_path]
    elif os.path.isdir(args.input_path):
        patterns = []
        if 'a' in args.strand.lower():
            patterns.append(os.path.join(args.input_path, "A_*.fa"))
        if 'b' in args.strand.lower():
            patterns.append(os.path.join(args.input_path, "B_*.fa"))
        
        fasta_files = []
        for pattern in patterns:
            fasta_files.extend(glob.glob(pattern))
        fasta_files = sorted(fasta_files)
    else:
        print(f"Error: {args.input_path} is not a valid file or directory")
        sys.exit(1)
    
    if not fasta_files:
        print(f"No FASTA files found in {args.input_path}")
        print("Expected pattern: A_*.fa or B_*.fa")
        sys.exit(1)
    
    print("=" * 60)
    print("Methylyzer Diagnostic: Dark Phase, Truncation & Fold Detection")
    print("=" * 60)
    print(f"Input: {args.input_path}")
    print(f"Left site (endogenous): {args.left_site}")
    print(f"Right site (accessibility): {args.right_site}")
    print(f"Dark phase threshold: ≤{args.dark_threshold}% methylation at both sites")
    print(f"Shift threshold: ≥{args.shift_threshold}% wrong bases")
    print(f"Truncation threshold: <{args.truncation_threshold}% coverage")
    print(f"Fold detection: ≥{args.fold_identity}% identity, ≥{args.min_fold_length} bp")
    print(f"Files to process: {len(fasta_files)}")
    print()
    
    # Process all files
    all_results = []
    
    for fasta_file in fasta_files:
        print(f"Processing: {os.path.basename(fasta_file)}")
        results = process_fasta(fasta_file, args.left_site, args.right_site,
                                args.min_fold_length, args.fold_identity)
        all_results.extend(results)
        print()
    
    if not all_results:
        print("No reads found to analyze")
        sys.exit(1)
    
    # Classify all reads
    for row in all_results:
        row['classification'] = classify_read(row, args.dark_threshold, args.shift_threshold, args.truncation_threshold)
    
    # Count classifications
    class_counts = defaultdict(int)
    for row in all_results:
        class_counts[row['classification']] += 1
    
    # Write metrics TSV
    metrics_file = f"{args.output}_metrics.tsv"
    with open(metrics_file, 'w') as f:
        # Header
        columns = [
            'read_id', 'fasta_file', 'classification',
            'coverage_pct', 'truncation_type', 'five_prime_missing', 'three_prime_missing',
            'is_folded', 'fold_point_ref', 'fold_identity', 'fold_length',
            'left_meth_pct', 'left_methylated', 'left_unmethylated', 'left_wrong_base', 'left_no_data', 'left_total',
            'right_meth_pct', 'right_methylated', 'right_unmethylated', 'right_wrong_base', 'right_no_data', 'right_total',
            'other_c_conversion'
        ]
        f.write('\t'.join(columns) + '\n')
        
        for row in all_results:
            values = []
            for col in columns:
                val = row.get(col)
                if val is None:
                    values.append('NA')
                elif isinstance(val, float):
                    values.append(f'{val:.2f}')
                else:
                    values.append(str(val))
            f.write('\t'.join(values) + '\n')
    
    print(f"Wrote per-read metrics: {metrics_file}")
    
    # Write flagged read IDs
    flagged_file = f"{args.output}_flagged.txt"
    flagged_count = 0
    with open(flagged_file, 'w') as f:
        for row in all_results:
            if row['classification'] != 'NORMAL':
                f.write(f"{row['read_id']}\t{row['classification']}\n")
                flagged_count += 1
    
    print(f"Wrote flagged read IDs: {flagged_file} ({flagged_count} reads)")
    
    # Write truncated read IDs with details
    truncated_file = f"{args.output}_truncated.txt"
    truncated_count = 0
    with open(truncated_file, 'w') as f:
        f.write("# read_id\tclassification\tcoverage_pct\ttruncation_type\t5prime_missing\t3prime_missing\n")
        for row in all_results:
            if 'TRUNCATED' in row['classification'] or row['coverage_pct'] < args.truncation_threshold:
                f.write(f"{row['read_id']}\t{row['classification']}\t{row['coverage_pct']:.1f}\t{row['truncation_type']}\t{row['five_prime_missing']}\t{row['three_prime_missing']}\n")
                truncated_count += 1
    
    print(f"Wrote truncated read IDs: {truncated_file} ({truncated_count} reads)")
    
    # Write folded read IDs with details
    folded_file = f"{args.output}_folded.txt"
    folded_count = 0
    with open(folded_file, 'w') as f:
        f.write("# read_id\tfold_point_ref\tfold_identity\tfold_length\tcoverage_pct\n")
        for row in all_results:
            if row.get('is_folded', False):
                f.write(f"{row['read_id']}\t{row['fold_point_ref']}\t{row['fold_identity']:.1f}\t{row['fold_length']}\t{row['coverage_pct']:.1f}\n")
                folded_count += 1
    
    print(f"Wrote folded read IDs: {folded_file} ({folded_count} reads)")
    
    # Write summary
    summary_file = f"{args.output}_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("METHYLYZER DIAGNOSTIC SUMMARY\n")
        f.write("=" * 60 + "\n\n")
        
        f.write(f"Input: {args.input_path}\n")
        f.write(f"Left site: {args.left_site}\n")
        f.write(f"Right site: {args.right_site}\n")
        f.write(f"Dark threshold: ≤{args.dark_threshold}%\n")
        f.write(f"Shift threshold: ≥{args.shift_threshold}%\n")
        f.write(f"Truncation threshold: <{args.truncation_threshold}%\n")
        f.write(f"Fold detection: ≥{args.fold_identity}% identity, ≥{args.min_fold_length} bp\n\n")
        
        f.write(f"Total reads analyzed: {len(all_results)}\n\n")
        
        f.write("CLASSIFICATION COUNTS:\n")
        f.write("-" * 40 + "\n")
        for cls in ['NORMAL', 'DARK_PHASE', 'SHIFTED', 'LOW_ENDOGENOUS', 'LOW_ACCESSIBILITY', 
                    'TRUNCATED', 'TRUNCATED_DARK', 'TRUNCATED_SHIFTED', 'FOLDED']:
            count = class_counts.get(cls, 0)
            pct = count / len(all_results) * 100
            f.write(f"  {cls:20s}: {count:6d} ({pct:5.1f}%)\n")
        
        f.write("\n")
        
        # Fold statistics
        folded_reads = [r for r in all_results if r.get('is_folded', False)]
        if folded_reads:
            f.write("FOLD STATISTICS:\n")
            f.write("-" * 40 + "\n")
            
            # Average fold point
            fold_points = [r['fold_point_ref'] for r in folded_reads if r['fold_point_ref'] is not None]
            if fold_points:
                avg_fold = sum(fold_points) / len(fold_points)
                min_fold = min(fold_points)
                max_fold = max(fold_points)
                f.write(f"  Fold point (ref position):\n")
                f.write(f"    Min: {min_fold}\n")
                f.write(f"    Max: {max_fold}\n")
                f.write(f"    Avg: {avg_fold:.1f}\n")
            
            # Average identity
            identities = [r['fold_identity'] for r in folded_reads]
            avg_identity = sum(identities) / len(identities)
            f.write(f"  Avg fold identity: {avg_identity:.1f}%\n")
            
            # Average fold length
            lengths = [r['fold_length'] for r in folded_reads]
            avg_length = sum(lengths) / len(lengths)
            f.write(f"  Avg fold length: {avg_length:.1f} bp\n")
            
            f.write("\n")
        
        # Truncation statistics
        truncated_reads = [r for r in all_results if 'TRUNCATED' in r['classification'] or r['coverage_pct'] < args.truncation_threshold]
        if truncated_reads:
            f.write("TRUNCATION STATISTICS:\n")
            f.write("-" * 40 + "\n")
            
            # Count by truncation type
            trunc_types = defaultdict(int)
            for r in truncated_reads:
                trunc_types[r['truncation_type']] += 1
            
            f.write("  By truncation location:\n")
            for ttype in ['5prime', '3prime', 'both', 'none']:
                count = trunc_types.get(ttype, 0)
                f.write(f"    {ttype:10s}: {count:6d}\n")
            
            # Average missing bases
            five_prime = [r['five_prime_missing'] for r in truncated_reads if r['five_prime_missing'] is not None]
            three_prime = [r['three_prime_missing'] for r in truncated_reads if r['three_prime_missing'] is not None]
            
            if five_prime:
                f.write(f"  Avg 5' missing: {sum(five_prime)/len(five_prime):.1f} bp\n")
            if three_prime:
                f.write(f"  Avg 3' missing: {sum(three_prime)/len(three_prime):.1f} bp\n")
            
            f.write("\n")
        
        # Stats for dark phase reads specifically
        dark_reads = [r for r in all_results if r['classification'] in ('DARK_PHASE', 'TRUNCATED_DARK')]
        if dark_reads:
            f.write("DARK PHASE READ STATISTICS:\n")
            f.write("-" * 40 + "\n")
            
            # Average conversion at other C's
            conversions = [r['other_c_conversion'] for r in dark_reads if r['other_c_conversion'] is not None]
            if conversions:
                avg_conv = sum(conversions) / len(conversions)
                f.write(f"  Avg conversion at non-site C's: {avg_conv:.1f}%\n")
            
            # Average coverage
            coverages = [r['coverage_pct'] for r in dark_reads]
            avg_cov = sum(coverages) / len(coverages)
            f.write(f"  Avg read coverage: {avg_cov:.1f}%\n")
            
            f.write("\n")
        
        # Compare dark phase vs normal
        normal_reads = [r for r in all_results if r['classification'] == 'NORMAL']
        if dark_reads and normal_reads:
            f.write("COMPARISON: DARK PHASE vs NORMAL\n")
            f.write("-" * 40 + "\n")
            
            # Conversion rates
            dark_conv = [r['other_c_conversion'] for r in dark_reads if r['other_c_conversion'] is not None]
            normal_conv = [r['other_c_conversion'] for r in normal_reads if r['other_c_conversion'] is not None]
            
            if dark_conv and normal_conv:
                f.write(f"  Non-site C conversion (dark):   {sum(dark_conv)/len(dark_conv):.1f}%\n")
                f.write(f"  Non-site C conversion (normal): {sum(normal_conv)/len(normal_conv):.1f}%\n")
            
            f.write("\n")
    
    print(f"Wrote summary: {summary_file}")
    
    # Print summary to console
    print()
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total reads: {len(all_results)}")
    print()
    for cls in ['NORMAL', 'DARK_PHASE', 'SHIFTED', 'LOW_ENDOGENOUS', 'LOW_ACCESSIBILITY',
                'TRUNCATED', 'TRUNCATED_DARK', 'TRUNCATED_SHIFTED', 'FOLDED']:
        count = class_counts.get(cls, 0)
        pct = count / len(all_results) * 100
        print(f"  {cls:20s}: {count:6d} ({pct:5.1f}%)")
    print()
    
    if class_counts.get('DARK_PHASE', 0) > 0 or class_counts.get('TRUNCATED_DARK', 0) > 0:
        print(f"Dark phase read IDs written to: {flagged_file}")
    if truncated_count > 0:
        print(f"Truncated read IDs written to: {truncated_file}")
    if folded_count > 0:
        print(f"Folded read IDs written to: {folded_file}")
    print("Use with seqtk to extract: seqtk subseq original.fa <file>.txt > extracted.fa")


if __name__ == "__main__":
    main()
