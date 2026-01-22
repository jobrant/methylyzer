#!/usr/bin/env python3

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
        """Calculate methylation frequencies for each site type."""
        strand = file_info['strand']
        locus = file_info['locus']
        
        print(f"Calculating frequencies for {len(self.sequences)} sequences...")
        
        for site in self.sites:
            print(f"Processing site: {site}")
            positions = self.site_positions.get(site, [])
            
            if not positions:
                print(f"No {site} sites found in reference sequence")
                continue
            
            # Initialize counters for each position
            position_counts = {}
            for pos in positions:
                position_counts[pos] = Counter()
            
            # Analyze each sequence
            valid_sequences = 0
            for seq_record in self.sequences:
                seq_results = self.analyze_sequence(seq_record, self.ref_seq)
                site_data = seq_results[site]
                
                has_data = False
                for site_info in site_data:
                    pos = site_info['position']
                    pattern = site_info['pattern']
                    if 'M' in pattern or 'U' in pattern:  # Has some readable data
                        has_data = True
                    position_counts[pos][pattern] += 1
                
                if has_data:
                    valid_sequences += 1
            
            # Write frequency file with new naming convention
            freq_file = os.path.join(output_dir, f"{site}-{strand}-{locus}.freqs.csv")
            print(f"Writing {site} frequencies to {freq_file}")
            
            with open(freq_file, 'w') as f:
                f.write(f"# Methylation frequencies for {site} sites\n")
                f.write(f"# Locus: {locus}\n")
                f.write(f"# Strand: {strand}\n")
                f.write(f"# Reference sequence: {len(self.ref_seq)} bp\n")
                f.write(f"# Total sequences analyzed: {len(self.sequences)}\n")
                f.write(f"# Sequences with data: {valid_sequences}\n")
                f.write(f"# {site} sites found: {len(positions)}\n\n")
                
                # Header
                f.write("Position,Site_Pattern,Count,Frequency\n")
                
                # Data for each position
                for pos in sorted(positions):
                    total_reads = sum(position_counts[pos].values())
                    if total_reads == 0:
                        continue
                        
                    for pattern, count in position_counts[pos].most_common():
                        freq = count / total_reads if total_reads > 0 else 0
                        f.write(f"{pos + 1},{pattern},{count},{freq:.4f}\n")
    
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
        epilog="Example: python bsFreq.py /path/to/extracted_files/ -o methylation_analysis"
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
            analyzer.write_csv_data(args.output, file_info)
        
        print()
    
    print(f"Analysis complete! Results in: {args.output}")
    print("\nOutput files:")
    freq_files = glob.glob(os.path.join(args.output, "*.freqs.csv"))
    matrix_files = glob.glob(os.path.join(args.output, "*.matrix.csv"))
    
    for f in sorted(freq_files + matrix_files):
        print(f"  {os.path.basename(f)}")

if __name__ == "__main__":
    main()