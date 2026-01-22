#!/usr/bin/env python3
# 2026, Jason Orr Brant, University of Florida Health Cancer Institute

"""
bsExtract.py - Modern minimal extractor for bsAlign database
Extracts filtered, deduplicated single-molecule FASTA files per locus
and generates a useful QC report.tsv
For downstream use with methylscaper (R package)
"""

import argparse
import os
import sqlite3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from datetime import datetime
from os import mkdir, path

##############################################################################

class meTable:
    def __init__(self, contig, BS=95, unambig=True, **kwargs):
        self.BS = BS / 100.0
        self.__head__ = []
        self.__offsets__ = {}
        self.__sites__ = {}        
        self.ref = contig[0]
        self.seqs = {seq.id: seq for seq in contig[1:]}
        
        # DICT of LAMBDA functions
        self.__rules__ = {
            'locus': lambda x: self.ref.id,
            'seqID': lambda x: x.id
        }
        
        # default sites
        if not kwargs: 
            kwargs = dict(CG=1, GC=2)
            
        for k, v in kwargs.items():
            try: 
                assert type(v) is int and v != 0
            except: 
                raise IOError('Site parameters must be in format GC=2, etc.')
            # add site to Header
            self.__head__ += [k]
            # get list of instances of site
            self.__sites__[k] = self.__match__(k, v)
            self.__offsets__[k] = v
            # make new rule for methylation site k
            self.__rules__[k] = lambda x, y=k: self.__score__(x, y)
            
        self.__head__ = ['locus', 'seqID'] + sorted(self.__head__)
        self.__checkBS__()
        
        # remove overlapping sites from each list of positions
        if unambig:
            self.__sites__ = {k: self.__pare__(k, self.__sites__[k])
                            for k in kwargs.keys()}
    
    def __checkBS__(self):
        """This is the key bisulfite conversion filter that was more permissive"""
        B = self.BS / (1 - self.BS)
        C = self.__pare__('C', self.__match__('C', 1))
        D = lambda x, y: sum([x.seq[c].upper() == y for c in C if c < len(x.seq)])        
        self.seqs = {i: j for i, j in self.seqs.items() if D(j, 'T') >= B * D(j, 'C')}
    
    def __getitem__(self, i):
        if type(i) == str:
            return [self.__rules__[j](self.seqs[i]) for j in self.__head__]
        else:
            seq_values = list(self.seqs.values())  # Python 3 fix
            return [self.__rules__[j](seq_values[i]) for j in self.__head__]
    
    def __len__(self):
        return len(self.seqs)
    
    def __match__(self, M, N):
        """Get bp positions of all instances of site M (methylated at bp N)"""
        ref = str(self.ref.seq).upper()
        return [n + N - 1 for n in range(len(ref)) if ref.find(M, n) == n]
    
    def __pare__(self, K, I):
        """List all bp position from list I that do not overlap with a methylation site other than K"""
        return sorted(set(I) - set(
            [j for k in self.__head__[2:] if k != K for j in self.__sites__.get(k, [])]))

##############################################################################

def snowflake(seqs, sites):
    """Deduplicate reads by unique endogenous methylation pattern at monitored sites."""
    U = {''.join([str(seq.seq[c].upper()) for c in sites if c < len(seq.seq)]): seq for seq in seqs}
    K = []
    for u in sorted(U.keys()):
        for n, k in enumerate(K):
            # k in u: replace k 
            if all([i in ('-', j) for i, j in zip(k, u)]): 
                K[n] = u
                break
            # mismatch: try next
            elif not all([j in ('-', i) for i, j in zip(k, u)]): 
                continue
            break
        else:
            K.append(u)
    return [U[k] for k in K]
            
##############################################################################

def unpack(db_path, dest_dir=None,
           loci=None, strands=['ab'],
           methyl=None, min_len='100', min_bs=95, 
           uniques=False, chrom=None, sample_name_override=None):
    loci = loci or []

    # Extract sample name
    if sample_name_override:
        sample_name = sample_name_override
    else:
        # Extract from directory path
        db_directory = os.path.dirname(os.path.abspath(db_path))
        sample_name = os.path.basename(db_directory)
    
    print(f"Sample name: {sample_name}")
    
    # Handle strands argument (default 'ab')
    if isinstance(strands, str):
        strands = strands.lower().strip()
        if strands == 'ab':
            strands = ['a', 'b']
        elif strands in ('a', 'b'):
            strands = [strands]
        else:
            # Accept comma-separated or other separators
            strands = [s.strip() for s in strands.replace(',', ' ').split() if s.strip() in ('a', 'b')]
    else:
        # If list, normalize
        strands = [s.lower() for s in strands if s.lower() in ('a', 'b')]

    # Remove duplicates and ensure at least something
    strands = list(set(strands))
    if not strands:
        strands = ['a', 'b']  # fallback to both

    print(f"Processing strands: {strands}")
    
    methyl = methyl or {'CG': 1, 'GC': 2}
    
    # prepare target directory
    if not dest_dir: 
        dest_dir = os.path.splitext(db_path)[0] + "_extracted"
    os.makedirs(dest_dir, exist_ok=True)    
    
    conn = sqlite3.connect(db_path)
    curs = conn.cursor()
    
    # Auto-populate loci if not specified
    if not loci:
        curs.execute("SELECT DISTINCT locus FROM records")
        loci = [row[0] for row in curs.fetchall()]
        
    if chrom:
        loci = [l for l in loci if l.startswith(f"{chrom}-")]
    
    # Debug: Print what we found
    print(f"Found loci: {loci[:5]}{'...' if len(loci) > 5 else ''} (total: {len(loci)})")
    print(f"Processing strands: {strands}")
    
    # For report: track counts
    report = {}  # key: (locus, strand) → [aligned, failed_len, failed_bs, removed_dupe, final]
    
    # Positions to monitor for snowflake dedup
    monitored_positions = []
    if uniques:
        curs.execute("SELECT sequence FROM loci ORDER BY LENGTH(sequence) DESC LIMIT 1")
        ref_row = curs.fetchone()
        if ref_row:
            longest_ref = ref_row[0].upper()
            for context, offset in methyl.items():
                clen = len(context)
                for i in range(len(longest_ref) - clen + 1):
                    if longest_ref[i:i+clen] == context:
                        monitored_positions.append(i + offset - 1)
    
    total_files = 0
    total_reads_final = 0
    
    for locus in loci:
        for strand in strands:
            key = (locus, strand)
            print(f"\n--- Processing: locus={locus}, strand={strand} ---")
            print(f"  Strands being processed: {strands}")
            
            # Fetch reference - try exact match first, then without strand suffix
            ref_seq = None
            ref_candidates = [locus, locus.rstrip('.ab')]
            for ref_id in ref_candidates:
                curs.execute("SELECT sequence FROM loci WHERE locus = ?", (ref_id,))
                ref_row = curs.fetchone()
                if ref_row:
                    ref_seq = Seq(ref_row[0].upper())
                    break
                    
            if not ref_seq:
                print(f"  No reference found for locus '{locus}' (tried: {ref_candidates})")
                continue
                
            print(f"  Reference length: {len(ref_seq)} bp")
            ref_rec = SeqRecord(ref_seq, id=locus, description="")
            
            # Fetch all aligned reads for this locus and strand
            curs.execute("""
                SELECT read, sequence 
                FROM records 
                WHERE locus = ? AND strand = ?
            """, (locus, strand))
            
            raw_reads = curs.fetchall()
            aligned_count = len(raw_reads)
            print(f"  Aligned reads found: {aligned_count}")
            
            if aligned_count == 0:
                continue
            
            # Length filter
            length_passed = []
            for read_id, seq_str in raw_reads:
                ungapped_len = len(seq_str.replace('-', ''))
                if '%' in str(min_len):
                    required = (int(min_len.strip('%')) / 100.0) * len(ref_seq)
                else:
                    required = int(min_len)
                    
                if ungapped_len >= required:
                    length_passed.append((read_id, seq_str))
                    
            print(f"  Passed length filter ({min_len}): {len(length_passed)}")
            
            # BS conversion filter 
            if len(length_passed) > 0:

                contig = [ref_rec] + [SeqRecord(Seq(seq_str), id=read_id) 
                            for read_id, seq_str in length_passed]
    
                # Apply meTable filtering with the same parameters as old pipeline
                try:
                    data = meTable(contig, BS=min_bs, **methyl)
                    bs_passed = [(seq.id, str(seq.seq)) for seq in data.seqs.values()]
                except Exception as e:
                    print(f"  Warning: meTable filtering failed: {e}")
                    # Fallback to passing all reads if meTable fails
                    bs_passed = length_passed
            else:
                bs_passed = []

            print(f"  Passed BS filter (min {min_bs}%): {len(bs_passed)} (meTable filtering)")
            
            # Create SeqRecords
            seq_records = []
            for read_id, seq_str in bs_passed:
                full_seq = Seq(seq_str.upper())
                seq_records.append(SeqRecord(full_seq, id=read_id, description=""))
            
            final_before_dedup = len(seq_records)
            
            # Deduplicate if requested
            if uniques and monitored_positions:
                before_dedup = len(seq_records)
                seq_records = snowflake(seq_records, monitored_positions)
                print(f"  Deduplicated: {before_dedup} → {len(seq_records)} unique patterns")
            
            final_count = len(seq_records)
            print(f"  Final reads written: {final_count}")
            
            # Store report counts
            report[key] = [aligned_count, aligned_count - len(length_passed),
                           len(length_passed) - len(bs_passed), final_before_dedup - final_count, final_count]
            
            if not seq_records:
                continue
            
            # Write FASTA - sample name and clean locus name
            clean_locus = locus.split(':')[0] if ':' in locus else locus
            filename = f"{strand.upper()}_{sample_name}_{clean_locus}.fa"
            out_path = path.join(dest_dir, filename)
            
            with open(out_path, 'w') as handle:
                SeqIO.write([ref_rec] + seq_records, handle, 'fasta')
            
            total_reads_final += final_count
            total_files += 1
            print(f"Wrote {final_count} reads → {out_path}")
    
    # Write report
    report_path = path.join(dest_dir, 'report.tsv')
    with open(report_path, 'w') as rpt:
        rpt.write(f"# Generated: {datetime.now()}\n")
        rpt.write(f"# Source DB: {db_path}\n")
        rpt.write(f"# Filters: min_len={min_len}, min_bs={min_bs}%, uniques={uniques}\n")
        rpt.write("#LOCUS\tSTRAND\tALIGNED\tFAILED_LEN\tFAILED_BS\tREMOVED_DUPE\tFINAL\n")
        
        for (locus, strand) in sorted(report.keys(), key=lambda x: (x[0], x[1])):
            counts = report[(locus, strand)]
            rpt.write('\t'.join([locus, strand.upper()] + [str(c) for c in counts]) + '\n')
    
    print(f"\nExtraction complete:")
    print(f"  {total_reads_final} final reads in {total_files} FASTA files")
    print(f"  Report written to: {report_path}")
        
##############################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extract filtered/deduplicated single-molecule FASTA from bsAlign DB + QC report"
    )
    parser.add_argument('database', help="bsAlign output SQLite database (.db)")
    parser.add_argument('--dest', default=None, help="Output directory")
    parser.add_argument('--loci', nargs='*', default=[], help="Locus IDs")
    parser.add_argument('--strand', default='ab', help="Strands to extract (a/b/ab)")
    parser.add_argument('--min-len', default='100', help="Min length (bp or e.g. '80%%')")
    parser.add_argument('--min-bs', type=int, default=95, help="Min bisulfite conversion rate (0-100)")
    parser.add_argument('--uniques', action='store_true', help="Deduplicate by methylation pattern")
    parser.add_argument('--sites', type=str, default=None,
                        help="Methylation contexts, e.g. 'CG=1,GC=2,CC=1'")
    parser.add_argument('--chrom', default=None, help="Only loci from this chromosome")
    parser.add_argument('--sample-name', default=None, 
                        help="Sample name for output files (default: use directory name)")
    
    args = parser.parse_args()
    
    sites = {'CG': 1, 'GC': 2}
    if args.sites:
        sites = dict(item.split('=') for item in args.Sites.split(','))
        sites = {k: int(v) for k, v in sites.items()}
    
    unpack(
        db_path=args.database,
        dest_dir=args.dest,
        loci=args.loci or None,
        strands=args.strand,
        methyl=sites,
        min_len=args.min_len,
        min_bs=args.min_bs,
        uniques=args.uniques,
        chrom=args.chrom,
        sample_name_override=args.sample_name
    )