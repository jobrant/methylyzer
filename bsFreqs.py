#!/usr/bin/env python3

"""
bsFreqs.py - Methylation frequency analysis from bsExtract FASTA files

Analyzes methylation patterns at specified sequence contexts (e.g., HCG, GCH, CG, GC).
Generates frequency tables and optional CSV matrices for visualization.

Supported site patterns:
  - Dinucleotides: CG, GC
  - Trinucleotides with IUPAC codes: HCG, GCH, WCG, etc.
  
IUPAC ambiguity codes supported:
  H = A, C, or T (not G)
  W = A or T
  S = C or G
  etc.

The cytosine position is automatically detected from the pattern.

Author: Jason Orr Brant, 2026
"""

import argparse
import os
import sys
import glob
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq
import re


def get_cytosine_offset(site):
    """
    Determine the cytosine position (0-indexed) in a methylation site pattern.
    
    Args:
        site: Methylation context pattern (e.g., 'CG', 'GC', 'HCG', 'GCH')
        
    Returns:
        Integer offset of the cytosine within the pattern
        
    Raises:
        ValueError: If no cytosine found in pattern
        
    Examples:
        get_cytosine_offset('CG')  -> 0
        get_cytosine_offset('GC')  -> 1  
        get_cytosine_offset('HCG') -> 1
        get_cytosine_offset('GCH') -> 1
        get_cytosine_offset('WCG') -> 1
    """
    site_upper = site.upper()
    
    if 'C' not in site_upper:
        raise ValueError(f"No cytosine found in site pattern: {site}")
    
    # Find first occurrence of C
    c_index = site_upper.index('C')
    
    # Check for multiple C's and warn
    c_count = site_upper.count('C')
    if c_count > 1:
        print(f"WARNING: Site pattern '{site}' contains {c_count} cytosines. "
              f"Using first C at position {c_index} (0-indexed). "
              f"Multi-cytosine patterns like CCG may need special handling.")
    
    return c_index


def validate_sites(sites):
    """
    Validate a list of methylation site patterns.
    
    Args:
        sites: List of site pattern strings
        
    Returns:
        Dict mapping site -> cytosine_offset
        
    Prints warnings for any problematic patterns.
    """
    offsets = {}
    
    for site in sites:
        try:
            offset = get_cytosine_offset(site)
            offsets[site] = offset
            print(f"  Site '{site}': cytosine at position {offset} (0-indexed)")
        except ValueError as e:
            print(f"ERROR: {e}")
            sys.exit(1)
    
    return offsets


class MethylationAnalyzer:
    def __init__(self):
        self.sites = ['HCG', 'GCH']  # Default methylation contexts
        self.ref_seq = None
        self.sequences = []
        self.site_positions = {}  # site -> list of positions
        self.cytosine_offsets = {}  # site -> cytosine offset (0-indexed)
        
    def set_sites(self, sites):
        """
        Set methylation sites to analyze and validate/compute cytosine offsets.
        
        Args:
            sites: List of site patterns (e.g., ['HCG', 'GCH'] or ['CG', 'GC'])
        """
        self.sites = sites
        print("Validating methylation site patterns:")
        self.cytosine_offsets = validate_sites(sites)
        print()
        
    def find_site_positions(self, sequence, site):
        """Find all positions where a methylation site occurs in the sequence."""
        positions = []
        site_len = len(site)
        seq_str = str(sequence).upper()
        
        # Handle ambiguous bases in site patterns (IUPAC codes)
        site_pattern = site.upper()
        site_pattern = site_pattern.replace('H', '[ACT]')  # H = not G
        site_pattern = site_pattern.replace('W', '[AT]')   # W = weak (A or T)
        site_pattern = site_pattern.replace('S', '[CG]')   # S = strong (C or G)
        site_pattern = site_pattern.replace('M', '[AC]')   # M = amino (A or C)
        site_pattern = site_pattern.replace('K', '[GT]')   # K = keto (G or T)
        site_pattern = site_pattern.replace('R', '[AG]')   # R = purine (A or G)
        site_pattern = site_pattern.replace('Y', '[CT]')   # Y = pyrimidine (C or T)
        site_pattern = site_pattern.replace('B', '[CGT]')  # B = not A
        site_pattern = site_pattern.replace('D', '[AGT]')  # D = not C
        site_pattern = site_pattern.replace('V', '[ACG]')  # V = not T
        site_pattern = site_pattern.replace('N', '[ACGT]') # N = any
        
        site_regex = re.compile(site_pattern)
        
        for match in site_regex.finditer(seq_str):
            positions.append(match.start())
        
        return positions
    
    def classify_methylation(self, ref_base, seq_base, site, position_in_site):
        """
        Classify methylation state based on reference and sequence bases.
        
        Uses dynamic cytosine offset detection rather than hardcoded positions.
        
        Returns: 'M' (methylated), 'U' (unmethylated), 'N' (ambiguous/gap), 'H' (context)
        """
        if ref_base == '-' or seq_base == '-' or seq_base == 'N':
            return 'N'
        
        ref_base = ref_base.upper()
        seq_base = seq_base.upper()
        
        # Get cytosine offset for this site
        c_offset = self.cytosine_offsets.get(site, get_cytosine_offset(site))
        
        if position_in_site == c_offset:
            # This is the cytosine position
            if ref_base == 'C':
                return 'M' if seq_base == 'C' else ('U' if seq_base == 'T' else 'N')
        
        # For non-cytosine positions, check if it's a context base or informative
        # Context bases (like H in HCG) are not directly informative
        return 'H'  # Context/flanking base
    
    def analyze_sequence(self, seq_record, ref_seq):
        """Analyze methylation pattern for one sequence."""
        results = {}
        
        for site in self.sites:
            site_data = []
            positions = self.site_positions.get(site, [])
            
            for pos in positions:
                site_len = len(site)
                site_result = []
                
                for i in range(site_len):
                    ref_pos = pos + i
                    if ref_pos < len(ref_seq) and ref_pos < len(seq_record.seq):
                        ref_base = str(ref_seq)[ref_pos]
                        seq_base = str(seq_record.seq)[ref_pos]
                        meth_state = self.classify_methylation(ref_base, seq_base, site, i)
                        site_result.append(meth_state)
                    else:
                        site_result.append('N')
                
                site_data.append({
                    'position': pos,
                    'states': site_result,
                    'pattern': ''.join(site_result)
                })
            
            results[site] = site_data
        
        return results
    
    def calculate_frequencies(self, output_dir, file_info):
        """Calculate methylation frequencies at methylatable cytosines."""
        strand = file_info['strand']
        locus = file_info['locus']
        
        print(f"Calculating methylation frequencies for {len(self.sequences)} sequences...")
        
        for site_type in self.sites:
            print(f"Processing {site_type} sites")
            start_positions = self.site_positions.get(site_type, [])
            
            if not start_positions:
                print(f"No {site_type} sites found")
                continue
            
            # Dynamic cytosine offset detection
            cytosine_offset = self.cytosine_offsets.get(site_type, get_cytosine_offset(site_type))
            
            print(f"DEBUG: Using cytosine_offset = {cytosine_offset} for {site_type}")
            print(f"DEBUG: Site start positions (0-based): {start_positions[:5]}")  # show first 5
            
            # Calculate actual cytosine positions for verification
            c_positions = [start + cytosine_offset for start in start_positions]
            print(f"DEBUG: Cytosine positions (0-based): {c_positions[:5]}")
            print(f"DEBUG: Cytosine positions (1-based): {[p + 1 for p in c_positions[:5]]}")
            
            # Verify reference bases at cytosine positions
            ref_str = str(self.ref_seq).upper()
            for pos in c_positions[:5]:
                if pos < len(ref_str):
                    print(f"DEBUG: Pos {pos+1} ref base = {ref_str[pos]} (should be 'C')")
            
            position_base_counts = defaultdict(Counter)
            total_reads_at_pos = defaultdict(int)

            for seq_record in self.sequences:
                seq_str = str(seq_record.seq).upper()
                
                for start in start_positions:
                    c_pos = start + cytosine_offset  # 0-based position of the C
                    
                    if c_pos >= len(self.ref_seq) or c_pos >= len(seq_str):
                        continue
                    
                    ref_base = self.ref_seq[c_pos].upper()
                    read_base = seq_str[c_pos]
                    
                    if ref_base != 'C':
                        print(f"WARNING: Pos {c_pos+1} ref is {ref_base} — skipping")
                        continue
                    
                    if read_base in {'-', 'N'}:
                        continue
                    
                    # Always use forward logic: reads are already oriented correctly by alignment
                    # Methylated → C, Unmethylated → T
                    if read_base == 'C':
                        methylated_base = 'C'
                    elif read_base == 'T':
                        methylated_base = 'T'
                    else:
                        methylated_base = read_base  # rare errors
                    
                    # Optional: Debug first few positions
                    # if c_pos in c_positions[:2]:  # first two sites
                    #     print(f"DEBUG: Pos {c_pos+1} | read_base = {read_base} → counted as {methylated_base}")
                    
                    position_base_counts[c_pos][methylated_base] += 1
                    total_reads_at_pos[c_pos] += 1
            
            # Prepare output file
            out_file = os.path.join(output_dir, f"{site_type}-cytosines-{strand}-{locus}.tsv")
            print(f"  → Writing to {out_file}")
            
            with open(out_file, 'w') as f:
                f.write(f"# Methylation at {site_type} cytosines\n")
                f.write(f"# Locus: {locus}\n")
                f.write(f"# Strand: {strand} ({'top' if strand == 'A' else 'bottom'} strand capture)\n")
                f.write(f"# Reference length: {len(self.ref_seq)} bp\n")
                f.write(f"# Sequences analyzed: {len(self.sequences)}\n")
                f.write(f"# {site_type} sites found: {len(start_positions)}\n")
                f.write(f"# Cytosine offset in pattern: {cytosine_offset}\n\n")
                
                f.write("Pos\tA\tC\tG\tT\n")
                
                for pos in sorted(position_base_counts.keys()):  # pos is 0-based cytosine position
                    counts = position_base_counts[pos]
                    total = total_reads_at_pos[pos]
                    
                    if total == 0:
                        continue
                    
                    fa = counts['A'] / total
                    fc = counts['C'] / total
                    fg = counts['G'] / total
                    ft = counts['T'] / total
                    
                    # Print pos as 1-based genomic position
                    f.write(f"{pos + 1}\t{fa:.6f}\t{fc:.6f}\t{fg:.6f}\t{ft:.6f}\n")

    def write_csv_data(self, output_dir, file_info):
        """Write coded CSV files for plotting (optional)."""
        strand = file_info['strand']
        locus = file_info['locus']
        
        print("Writing CSV matrix files...")
        
        for site in self.sites:
            positions = self.site_positions.get(site, [])
            if not positions:
                continue
            
            csv_file = os.path.join(output_dir, f"{site}-{strand}-{locus}.matrix.csv")
            print(f"Writing {site} matrix to {csv_file}")
            
            with open(csv_file, 'w') as f:
                # Header
                f.write("Sequence_ID")
                for pos in sorted(positions):
                    f.write(f",{site}_{pos + 1}")
                f.write("\n")
                
                # Data rows
                for seq_record in self.sequences:
                    seq_results = self.analyze_sequence(seq_record, self.ref_seq)
                    site_data = seq_results[site]
                    
                    f.write(seq_record.id)
                    
                    # Create position lookup
                    pos_lookup = {info['position']: info['pattern'] for info in site_data}
                    
                    for pos in sorted(positions):
                        pattern = pos_lookup.get(pos, 'N' * len(site))
                        # Convert pattern to numeric code (M=1, U=0, N/H=-1)
                        numeric_pattern = []
                        for state in pattern:
                            if state == 'M':
                                numeric_pattern.append('1')
                            elif state == 'U':
                                numeric_pattern.append('0')
                            else:
                                numeric_pattern.append('-1')
                        f.write(f",{'_'.join(numeric_pattern)}")
                    f.write("\n")
    
    def write_coded_map_csv(self, output_dir, file_info):
        """
        Write coded CSV files for plotting methylation maps.
        Each row is a molecule, each column is a reference position.
        Creates full-length vectors for visualization.
        """
        strand = file_info['strand']
        locus = file_info['locus']
        
        print("Writing coded map CSV files...")
        
        # Coding scheme matching methylscaper expectations
        codes = {
            'HCG_methylated': 2.0,      # Red patch
            'HCG_unmethylated': -2.0,   # Black patch  
            'GCH_methylated': 1.0,      # Yellow patch
            'GCH_unmethylated': 1.0,    # Yellow patch (same for GCH)
            'boundary': 0.0,            # Grey/background
            'no_data': -1.0,            # Missing data
            'context': -1.0             # Context bases within sites
        }
        
        ref_length = len(self.ref_seq)
        
        # Create position maps for each site type with dynamic cytosine offset
        position_to_site = {}  # pos -> (site_type, offset_in_site, is_cytosine)
        
        for site_type in self.sites:
            positions = self.site_positions.get(site_type, [])
            c_offset = self.cytosine_offsets.get(site_type, get_cytosine_offset(site_type))
            
            for start_pos in positions:
                for offset in range(len(site_type)):
                    abs_pos = start_pos + offset
                    if abs_pos < ref_length:
                        is_cytosine = (offset == c_offset)
                        position_to_site[abs_pos] = (site_type, offset, is_cytosine)
        
        # Use tab-separated format and .csv extension
        csv_file = os.path.join(output_dir, f"{strand}-{locus}_coded_map.csv")
        print(f"  Writing coded map to {csv_file}")
        
        with open(csv_file, 'w') as f:
            # Write header with all reference positions 
            f.write("#Seq")
            for pos in range(1, ref_length + 1):  # 1-based numbering
                f.write(f"\tC{pos}")
            f.write("\n")
            
            # Process each sequence/molecule
            for seq_record in self.sequences:
                seq_str = str(seq_record.seq).upper()
                f.write(seq_record.id)
                
                # Analyze each reference position
                for ref_pos in range(ref_length):  # 0-based for indexing
                    if ref_pos >= len(seq_str):
                        # Beyond read length
                        f.write(f"\t{codes['no_data']}")
                        continue
                    
                    read_base = seq_str[ref_pos]
                    ref_base = str(self.ref_seq)[ref_pos].upper()
                    
                    if read_base in {'-', 'N'}:
                        f.write(f"\t{codes['no_data']}")
                        continue
                    
                    # Check if this position is part of a methylation site
                    if ref_pos in position_to_site:
                        site_type, offset, is_cytosine = position_to_site[ref_pos]
                        
                        # Only analyze methylatable cytosines (using dynamic offset)
                        if is_cytosine and ref_base == 'C':
                            if site_type == 'HCG':
                                if read_base == 'C':
                                    code_value = codes['HCG_methylated']    # 2.0 (Red)
                                elif read_base == 'T':
                                    code_value = codes['HCG_unmethylated']  # -2.0 (Black)
                                else:
                                    code_value = codes['no_data']           # -1.0 (ambiguous)
                            
                            elif site_type == 'GCH':
                                # GCH sites get yellow (1.0) regardless of methylation state
                                code_value = codes['GCH_methylated']       # 1.0 (Yellow)
                            
                            else:
                                # For other site types (CG, GC, etc.), use HCG-style coding
                                if read_base == 'C':
                                    code_value = codes['HCG_methylated']    # 2.0
                                elif read_base == 'T':
                                    code_value = codes['HCG_unmethylated']  # -2.0
                                else:
                                    code_value = codes['no_data']
                        
                        else:
                            # Non-methylatable position within a site (context bases)
                            code_value = codes['context']                  # -1.0
                    
                    else:
                        # Position not in any methylation site
                        code_value = codes['boundary']                     # 0.0 (Grey)
                    
                    f.write(f"\t{code_value}")
                
                f.write("\n")

    def write_separate_site_maps(self, output_dir, file_info):
        """
        Write separate tab-separated CSV files for each site type,
        WITH patch-filling between consecutive sites.
        
        Format matches what the R plotting script (batch_methylscaper_plots.R) expects:
        - Tab-separated with #Seq header and C1, C2, ... column names
        - Values:
            2  = methylated site (C retained)
           -2  = unmethylated site (C→T conversion)
            1  = fill between two methylated sites
           -1  = fill between two unmethylated sites
            0  = fill between mixed sites OR non-site background
            .  = no data (gap, N, beyond read, ambiguous)
        - Same coding for ALL site types; the R recode handles color mapping
          (HCG file → red/black, GCH file → yellow/black)
        """
        strand = file_info['strand']
        locus = file_info['locus']
        
        print("Writing separate site map CSV files (with patch filling)...")
        
        ref_length = len(self.ref_seq)
        
        for site_type in self.sites:
            positions = self.site_positions.get(site_type, [])
            if not positions:
                continue
            
            # Get dynamic cytosine offset for this site type
            c_offset = self.cytosine_offsets.get(site_type, get_cytosine_offset(site_type))
            
            # Build SORTED list of cytosine positions for this site type
            cytosine_positions_list = sorted(set(
                start_pos + c_offset for start_pos in positions
                if start_pos + c_offset < ref_length
            ))
            cytosine_positions_set = set(cytosine_positions_list)
            
            csv_file = os.path.join(output_dir, f"{site_type}-{strand}-{locus}_map.csv")
            print(f"  Writing {site_type} map to {csv_file}")
            print(f"    {len(cytosine_positions_list)} cytosine sites to fill between")
            
            with open(csv_file, 'w') as f:
                # Tab-separated header matching coded_map format
                f.write("#Seq")
                for pos in range(1, ref_length + 1):
                    f.write(f"\tC{pos}")
                f.write("\n")
                
                # Process each sequence
                for seq_record in self.sequences:
                    seq_str = str(seq_record.seq).upper()
                    
                    # STEP 1: Determine methylation state at each cytosine site
                    site_states = {}  # pos -> 2 (meth), -2 (unmeth), 3 (wrong base), or None (no data)
                    for c_pos in cytosine_positions_list:
                        if c_pos >= len(seq_str):
                            site_states[c_pos] = None
                            continue
                        read_base = seq_str[c_pos]
                        ref_base = str(self.ref_seq)[c_pos].upper()
                        if read_base in {'-', 'N'} or ref_base != 'C':
                            site_states[c_pos] = None
                        elif read_base == 'C':
                            site_states[c_pos] = 2   # Methylated
                        elif read_base == 'T':
                            site_states[c_pos] = -2  # Unmethylated
                        else:
                            site_states[c_pos] = 3
                    
                    # STEP 2: Compute fill values between consecutive sites
                    # fill_regions[pos] = fill value for positions between site i and site i+1
                    fill_regions = {}  # ref_pos -> fill value
                    
                    for idx in range(len(cytosine_positions_list) - 1):
                        left_pos = cytosine_positions_list[idx]
                        right_pos = cytosine_positions_list[idx + 1]
                        left_state = site_states.get(left_pos)
                        right_state = site_states.get(right_pos)
                        
                        # Determine fill value
                        if left_state is None or right_state is None:
                            fill_val = None  # No data → "."
                        elif left_state == 2 and right_state == 2:
                            fill_val = 1     # Both methylated → fill
                        elif left_state == -2 and right_state == -2:
                            fill_val = -1    # Both unmethylated → fill
                        elif left_state == 3 and right_state == 3:
                            fill_val = 4
                        else:
                            fill_val = 0     # Mixed → boundary/gray
                        
                        # Apply fill to all positions between left and right
                        for fill_pos in range(left_pos + 1, right_pos):
                            fill_regions[fill_pos] = fill_val
                    
                    # STEP 3: Write the row
                    f.write(seq_record.id)
                    
                    for ref_pos in range(ref_length):
                        if ref_pos >= len(seq_str):
                            f.write("\t.")  # Beyond read
                            continue
                        
                        read_base = seq_str[ref_pos]
                        
                        if read_base in {'-', 'N'}:
                            f.write("\t.")  # Gap or ambiguous
                            continue
                        
                        # Is this a cytosine site?
                        if ref_pos in cytosine_positions_set:
                            state = site_states.get(ref_pos)
                            if state is not None:
                                f.write(f"\t{state}")  # 2 or -2
                            else:
                                f.write("\t.")
                        
                        # Is this in a fill region between two sites?
                        elif ref_pos in fill_regions:
                            fill_val = fill_regions[ref_pos]
                            if fill_val is not None:
                                f.write(f"\t{fill_val}")  # 1, -1, or 0
                            else:
                                f.write("\t.")
                        
                        # Outside any site region (before first site or after last)
                        else:
                            f.write("\t0")  # Background
                    
                    f.write("\n")

    def write_all_csv_formats(self, output_dir, file_info):
        """Write both coded and pattern CSV formats."""
        self.write_coded_map_csv(output_dir, file_info)
        self.write_separate_site_maps(output_dir, file_info)


def parse_filename(filename):
    """Parse bsExtract filename to get strand and locus info."""
    basename = os.path.splitext(os.path.basename(filename))[0]
    
    # Expected format: A_locus or B_locus
    parts = basename.split('_', 1)
    if len(parts) >= 2:
        strand = parts[0]  # A or B
        locus = parts[1]   # Everything after first underscore
        return {'strand': strand, 'locus': locus}
    else:
        # Fallback
        return {'strand': 'unknown', 'locus': basename}


def find_fasta_files(input_path):
    """Find all A_*.fa and B_*.fa files in the input path."""
    if os.path.isfile(input_path):
        # Single file provided
        return [input_path]
    elif os.path.isdir(input_path):
        # Directory provided - find all A_*.fa and B_*.fa files
        pattern_a = os.path.join(input_path, "A_*.fa")
        pattern_b = os.path.join(input_path, "B_*.fa")
        files = glob.glob(pattern_a) + glob.glob(pattern_b)
        return sorted(files)
    else:
        print(f"Error: {input_path} is not a valid file or directory")
        return []


def main():
    parser = argparse.ArgumentParser(
        description="Generate methylation frequency tables from bsExtract FASTA files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Supported methylation site patterns:
  Dinucleotides:  CG, GC
  Trinucleotides: HCG, GCH, WCG, etc.

The cytosine position is automatically detected from the pattern:
  CG  -> cytosine at position 0
  GC  -> cytosine at position 1
  HCG -> cytosine at position 1
  GCH -> cytosine at position 1

IUPAC ambiguity codes:
  H = A, C, or T (not G)
  W = A or T (weak)
  S = C or G (strong)
  M = A or C (amino)
  K = G or T (keto)
  R = A or G (purine)
  Y = C or T (pyrimidine)

Examples:
  # Standard MAPit analysis (HCG + GCH, excludes GCG)
  python bsFreqs.py /path/to/extracted/ -o frequencies
  
  # Include GCG sites (when no endogenous CG methylation)
  python bsFreqs.py /path/to/extracted/ -s GC CG -o frequencies
  
  # CG-only analysis
  python bsFreqs.py /path/to/extracted/ -s CG -o frequencies
  
  # Generate CSV matrices for methylscaper plotting
  python bsFreqs.py /path/to/extracted/ -o frequencies --csv
"""
    )
    parser.add_argument('input_path', 
                       help="Directory containing bsExtract output files (A_*.fa, B_*.fa) or single FASTA file")
    parser.add_argument('-o', '--output', default='methylation_analysis',
                       help="Output directory (default: methylation_analysis)")
    parser.add_argument('-s', '--sites', nargs='+', default=['HCG', 'GCH'],
                       help="Methylation sites to analyze (default: HCG GCH). "
                            "Cytosine position is auto-detected. "
                            "Examples: CG, GC, HCG, GCH, WCG")
    parser.add_argument('--csv', action='store_true',
                       help="Also generate CSV matrices for methylscaper plotting")
    parser.add_argument('--pattern', default="[AB]_*.fa",
                       help="File pattern to match (default: [AB]_*.fa)")
    
    args = parser.parse_args()
    
    # Find input files
    fasta_files = find_fasta_files(args.input_path)
    
    if not fasta_files:
        print("No FASTA files found!")
        print(f"Looking for pattern: {args.pattern}")
        if os.path.isdir(args.input_path):
            all_files = os.listdir(args.input_path)
            fa_files = [f for f in all_files if f.endswith('.fa')]
            print(f"Available .fa files: {fa_files}")
        sys.exit(1)
    
    print(f"Found {len(fasta_files)} FASTA files to process:")
    for f in fasta_files:
        print(f"  {f}")
    print()
    
    # Create main output directory (all files go here)
    os.makedirs(args.output, exist_ok=True)
    
    # Process each input file
    for fasta_file in fasta_files:
        file_info = parse_filename(fasta_file)
        print(f"=== Processing {os.path.basename(fasta_file)} ===")
        print(f"  Strand: {file_info['strand']}, Locus: {file_info['locus']}")
        
        analyzer = MethylationAnalyzer()
        analyzer.set_sites(args.sites)  # Validates and computes cytosine offsets
        
        # Read sequences
        sequences = list(SeqIO.parse(fasta_file, 'fasta'))
        if not sequences:
            print(f"Warning: No sequences found in {fasta_file}")
            continue
            
        # First sequence should be reference
        analyzer.ref_seq = sequences[0].seq
        analyzer.sequences = sequences[1:]  # Rest are reads
        
        print(f"  Reference: {len(analyzer.ref_seq)} bp")
        print(f"  Sequences: {len(analyzer.sequences)}")
        
        if len(analyzer.sequences) == 0:
            print("  No read sequences found, skipping...")
            continue
        
        # Find methylation sites in reference
        for site in analyzer.sites:
            positions = analyzer.find_site_positions(analyzer.ref_seq, site)
            analyzer.site_positions[site] = positions
            print(f"  Found {len(positions)} {site} sites")
        
        # Generate frequency tables directly in main output directory
        analyzer.calculate_frequencies(args.output, file_info)
        
        # Generate CSV if requested
        if args.csv:
            analyzer.write_coded_map_csv(args.output, file_info)
            analyzer.write_separate_site_maps(args.output, file_info)
        
        print()
    
    print(f"Analysis complete! Results in: {args.output}")
    print("\nOutput files:")
    freq_files = glob.glob(os.path.join(args.output, "*.tsv"))
    matrix_files = glob.glob(os.path.join(args.output, "*.csv"))
    
    for f in sorted(freq_files + matrix_files):
        print(f"  {os.path.basename(f)}")


if __name__ == "__main__":
    main()