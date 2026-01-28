#!/usr/bin/env python3
"""
bsFreqs.py - Generate methylation frequency tables and coded maps from bsExtract FASTA files

PATCHED VERSION: Fixed coded CSV generation with proper patch-filling logic
to match the original meTools.py / MethMap.py behavior for methylscaper visualization.

Value encoding (matching methylscaper's preprocessSingleMolecule output):
   2  = methylated site (dark tick at actual C)
  -2  = unmethylated site (dark tick at actual C)  
   1  = fill between consecutive methylated sites (colored patch)
  -1  = fill between consecutive unmethylated sites (colored patch)
   0  = boundary (transition between meth/unmeth) or background (grey)
  '.'  = no data (gap, N, or beyond read)
"""

import argparse
import os
import sys
import glob
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq
import re

class MethylationAnalyzer:
    def __init__(self):
        self.sites = ['HCG', 'GCH']  # Default methylation contexts
        self.ref_seq = None
        self.sequences = []
        self.site_positions = {}  # site -> list of positions
        
    def find_site_positions(self, sequence, site):
        """Find all positions where a methylation site occurs in the sequence."""
        positions = []
        site_len = len(site)
        seq_str = str(sequence).upper()
        
        # Handle ambiguous bases in site patterns
        site_pattern = site.replace('H', '[ACT]')  # H = not G
        site_regex = re.compile(site_pattern)
        
        for match in site_regex.finditer(seq_str):
            positions.append(match.start())
        
        return positions
    
    def classify_methylation(self, ref_base, seq_base, site, position_in_site):
        """
        Classify methylation state based on reference and sequence bases.
        
        For HCG sites: H(0) C(1) G(2) - only positions 1,2 are methylatable
        For GCH sites: G(0) C(1) H(2) - only positions 0,1 are methylatable
        
        Returns: 'M' (methylated), 'U' (unmethylated), 'N' (ambiguous/gap), 'H' (context)
        """
        if ref_base == '-' or seq_base == '-' or seq_base == 'N':
            return 'N'
        
        ref_base = ref_base.upper()
        seq_base = seq_base.upper()
        site_upper = site.upper()
        
        if site_upper == 'HCG':
            if position_in_site == 0:  # H position
                return 'H'  # Context base, not methylatable
            elif position_in_site == 1:  # C position
                if ref_base == 'C':
                    return 'M' if seq_base == 'C' else ('U' if seq_base == 'T' else 'N')
            elif position_in_site == 2:  # G position
                if ref_base == 'G':
                    return 'M' if seq_base == 'G' else ('U' if seq_base == 'A' else 'N')
                    
        elif site_upper == 'GCH':
            if position_in_site == 0:  # G position
                if ref_base == 'G':
                    return 'M' if seq_base == 'G' else ('U' if seq_base == 'A' else 'N')
            elif position_in_site == 1:  # C position
                if ref_base == 'C':
                    return 'M' if seq_base == 'C' else ('U' if seq_base == 'T' else 'N')
            elif position_in_site == 2:  # H position
                return 'H'  # Context base, not methylatable
                
        elif site_upper == 'CG':
            if position_in_site == 0:  # C position
                if ref_base == 'C':
                    return 'M' if seq_base == 'C' else ('U' if seq_base == 'T' else 'N')
            elif position_in_site == 1:  # G position
                if ref_base == 'G':
                    return 'M' if seq_base == 'G' else ('U' if seq_base == 'A' else 'N')
        
        return 'N'  # Default
    
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
            
            # Cytosine is always at offset 1 in both HCG (H C G) and GCH (G C H)
            cytosine_offset = 1
            
            # Calculate actual cytosine positions
            c_positions = [start + cytosine_offset for start in start_positions]
            
            # Verify reference bases at cytosine positions (first 5 for debug)
            ref_str = str(self.ref_seq).upper()
            
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
                f.write(f"# {site_type} sites found: {len(start_positions)}\n\n")
                
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

    def write_coded_map_csv(self, output_dir, file_info):
        """
        Write SEPARATE coded CSV files for each site type with proper patch-filling.
        
        This replicates the original meTools.py __score__() and MethMap.py logic:
        - Sites marked with 2 (methylated) or -2 (unmethylated)
        - Positions between sites filled based on flanking site states
        - Boundaries (mixed states) marked as 0 (grey)
        
        Value encoding (matching methylscaper's preprocessSingleMolecule output):
          2  = methylated site (dark tick)
         -2  = unmethylated site (dark tick)
          1  = fill between methylated sites (patch color)
         -1  = fill between unmethylated sites (patch color)
          0  = boundary or background (grey)
          .  = missing data
        """
        strand = file_info['strand']
        locus = file_info['locus']
        
        print("Writing coded map CSV files with patch-filling...")
        
        # Value encoding matching methylscaper's output (INTEGERS, not floats!)
        METH_SITE = 2         # methylated cytosine site (dark tick)
        UNMETH_SITE = -2      # unmethylated cytosine site (dark tick)
        METH_FILL = 1         # fill between methylated sites (patch color)
        UNMETH_FILL = -1      # fill between unmethylated sites (patch color)
        BOUNDARY = 0          # boundary or background (grey)
        NO_DATA = '.'         # Missing data (gap, N, beyond read)
        
        ref_length = len(self.ref_seq)
        ref_str = str(self.ref_seq).upper()
        
        # Process each site type separately (creates separate files)
        for site_type in self.sites:
            site_start_positions = self.site_positions.get(site_type, [])
            if not site_start_positions:
                print(f"  No {site_type} sites found, skipping...")
                continue
            
            # Get cytosine positions (offset 1 for both HCG and GCH)
            # These are the actual methylatable positions (0-based)
            cytosine_positions = sorted([start + 1 for start in site_start_positions 
                                         if start + 1 < ref_length])
            
            # Output file per site type
            csv_file = os.path.join(output_dir, f"{site_type}-{strand}-{locus}_map.csv")
            print(f"  Writing {site_type} map to {csv_file}")
            print(f"    {len(cytosine_positions)} cytosine positions")
            
            with open(csv_file, 'w') as f:
                # Header
                f.write("#Seq")
                for pos in range(1, ref_length + 1):
                    f.write(f"\tC{pos}")
                f.write("\n")
                
                # Process each sequence
                for seq_record in self.sequences:
                    seq_str = str(seq_record.seq).upper()
                    f.write(seq_record.id)
                    
                    # Step 1: Determine methylation state at each cytosine position
                    # state: 1 = methylated (C), -1 = unmethylated (T), 0 = ambiguous/no data
                    site_states = {}
                    for c_pos in cytosine_positions:
                        if c_pos >= len(seq_str):
                            site_states[c_pos] = 0  # Beyond read
                        else:
                            read_base = seq_str[c_pos]
                            ref_base = ref_str[c_pos]
                            
                            if read_base in {'-', 'N'} or ref_base != 'C':
                                site_states[c_pos] = 0  # No data or not a C in reference
                            elif read_base == 'C':
                                site_states[c_pos] = 1   # Methylated
                            elif read_base == 'T':
                                site_states[c_pos] = -1  # Unmethylated
                            else:
                                site_states[c_pos] = 0   # Ambiguous
                    
                    # Step 2: Find informative sites (non-zero states)
                    informative_sites = [(pos, site_states[pos]) for pos in cytosine_positions 
                                         if site_states[pos] != 0]
                    
                    # Step 3: Build the coded vector for this sequence
                    coded_vector = [BOUNDARY] * ref_length  # Start with all background
                    
                    if informative_sites:
                        # Mark site positions and fill between them
                        for idx, (site_pos, state) in enumerate(informative_sites):
                            # Mark the site itself
                            if state == 1:
                                coded_vector[site_pos] = METH_SITE      # 2.0
                            else:  # state == -1
                                coded_vector[site_pos] = UNMETH_SITE    # -2.0
                            
                            # Fill positions between this site and the next
                            if idx < len(informative_sites) - 1:
                                next_pos, next_state = informative_sites[idx + 1]
                                
                                # Determine fill value based on flanking states
                                # This matches meTools.py: ' +-'[(J[i]+J[j])//2]
                                if state == 1 and next_state == 1:
                                    fill_value = METH_FILL      # 1 - both methylated
                                elif state == -1 and next_state == -1:
                                    fill_value = UNMETH_FILL    # -1 - both unmethylated
                                else:
                                    fill_value = BOUNDARY       # 0 - mixed states
                                
                                # Fill positions between sites (exclusive of site positions)
                                for fill_pos in range(site_pos + 1, next_pos):
                                    coded_vector[fill_pos] = fill_value
                    
                    # Step 4: Write the vector, handling missing data
                    for ref_pos in range(ref_length):
                        if ref_pos >= len(seq_str):
                            f.write(f"\t{NO_DATA}")
                        elif seq_str[ref_pos] in {'-', 'N'}:
                            f.write(f"\t{NO_DATA}")
                        else:
                            f.write(f"\t{coded_vector[ref_pos]}")
                    
                    f.write("\n")
            
            print(f"    Written {len(self.sequences)} sequences")

    def write_csv_data(self, output_dir, file_info):
        """Write coded CSV files for plotting (optional) - DEPRECATED, use write_coded_map_csv instead."""
        # Redirect to the new method
        self.write_coded_map_csv(output_dir, file_info)


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
        epilog="Example: python bsFreqs.py /path/to/extracted_files/ -o methylation_analysis"
    )
    parser.add_argument('input_path', 
                       help="Directory containing bsExtract output files (A_*.fa, B_*.fa) or single FASTA file")
    parser.add_argument('-o', '--output', default='methylation_analysis',
                       help="Output directory (default: methylation_analysis)")
    parser.add_argument('-s', '--sites', nargs='+', default=['HCG', 'GCH'],
                       help="Methylation sites to analyze (default: HCG GCH)")
    parser.add_argument('--csv', action='store_true',
                       help="Also generate CSV matrices for plotting")
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
        analyzer.sites = args.sites
        
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
        
        print()
    
    print(f"Analysis complete! Results in: {args.output}")
    print("\nOutput files:")
    freq_files = glob.glob(os.path.join(args.output, "*.tsv"))
    matrix_files = glob.glob(os.path.join(args.output, "*_map.csv"))
    
    for f in sorted(freq_files + matrix_files):
        print(f"  {os.path.basename(f)}")

if __name__ == "__main__":
    main()