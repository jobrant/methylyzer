#!/usr/bin/env python3
# Rewritten 2025, Jason Orr Brant, University of Florida Health Cancer Institute

"""
* aligns sequences to references using BLAST
* masks cytosines for unbiased alignment
* returns SQL table
* user can choose to align on either or both strands
* works for cytosine methylation in any context (CG, GC, etc.)
* uses makeblastdb, blastn
* RPD4 2013
"""

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
import os
from os import environ, path, remove
from subprocess import call, run, PIPE, CalledProcessError
from time import time
import argparse
import sqlite3
import warnings
import tempfile
from io import StringIO
import glob
from shutil import which

##############################################################################

def find_executable(name, env_var=None):
    """
    Find executable in PATH, with optional environment variable override.
    """
    # Check if user provided custom path via environment variable
    if env_var and os.environ.get(env_var):
        custom_path = os.environ.get(env_var)
        if os.path.isfile(custom_path) and os.access(custom_path, os.X_OK):
            return custom_path
        else:
            raise RuntimeError(f"Custom path {custom_path} from {env_var} is not executable")
    
    # Use shutil.which to find in PATH
    exe_path = which(name)
    if exe_path:
        return exe_path
    else:
        raise RuntimeError(f"{name} not found in PATH. Please load the appropriate module (e.g., 'module load blast') or set {env_var} environment variable.")

# Find required executables
try:
    bn = find_executable('blastn', 'BLASTN_PATH')
    c2bm = find_executable('convert2blastmask', 'CONVERT2BLASTMASK_PATH') 
    mbdb = find_executable('makeblastdb', 'MAKEBLASTDB_PATH')
    
    print(f"Using blastn: {bn}")
    print(f"Using convert2blastmask: {c2bm}")
    print(f"Using makeblastdb: {mbdb}")
    
except RuntimeError as e:
    print(f"Error: {e}")
    print("\nTo fix this:")
    print("1. Load the BLAST module: module load blast")
    print("2. Or set custom paths: export BLASTN_PATH=/path/to/blastn")
    exit(1)

##############################################################################

# to swap extensions
suffix = lambda x, y: path.splitext(x)[0] + '.' + y

# to suppress warnings from use of tempnam
warnings.filterwarnings('ignore')

# for ambiguous base calls
ambigIUPAC = {'A': 'MRW', 'C': 'MSY', 'G': 'KRS', 'T': 'KWY',
              'K': 'GT', 'M': 'AC', 'R': 'AG', 'S': 'CG', 'W': 'AT', 'Y': 'CT'}
ambigIUPAC = {i + j: (list(set(list(ambigIUPAC[i])) & set(list(ambigIUPAC[j]))) + ['N'])[0]
              for i in 'ACGKMRSTWY' for j in 'ACGKMRSTWY' if i != j}

##############################################################################

def check_IDs(records):
    # Extract the base name and create the indexed file name
    indexed = os.path.splitext(records)[0] + '_indexed.fa'
    
    # Extract sequence IDs
    names = [record.id for record in SeqIO.parse(records, 'fasta')]
    
    # Check for duplicates
    if len(set(names)) < len(names):
        N = len(str(len(names)))
        
        # Open the source file and the new indexed file
        with open(records, 'r') as source, open(indexed, 'w') as handle:
            for n, rec in enumerate(SeqIO.parse(source, 'fasta')):
                rec.id = str(n).zfill(N) + '|' + rec.id
                rec.description = rec.id  # Ensure the description is updated as well
                SeqIO.write(rec, handle, 'fasta')
        return indexed
    else:
        return records
    
##############################################################################

def deaminate(source, strand='ab'):
    """
    Perform in silico bisulfite conversion using the modern faconvert tool.
    """
    seqs = []

    # Process each requested strand ('a', 'b', or both)
    for i in strand:
        # Determine conversion flag for faconvert
        # Assuming: empty string or nothing = C->T (Watson/'a' strand)
        # "GA" = G->A (Crick/'b' strand)

        # Run faconvert and capture its stdout directly
        try:
            if i == 'b':
                # G to A conversion for b strand
                result = run(
                    ['faconvert', source, 'GA'],
                    stdout=PIPE,
                    stderr=PIPE,
                    check=True,
                    text=True
                )
            else:
                # C to T conversion for a strand (default)
                result = run(
                    ['faconvert', source, 'CT'],
                    stdout=PIPE,
                    stderr=PIPE,
                    check=True,
                    text=True
                )
    
        except CalledProcessError as e:
            raise RuntimeError(f"faconvert failed for strand '{i}': {e.stderr.strip()}")

        # Parse the converted sequences directly from stdout
        for seq in SeqIO.parse(StringIO(result.stdout), 'fasta'):
            # For 'b' strand: reverse complement after conversion
            if i == 'b':
                seq.seq = seq.seq.reverse_complement()

            # Append strand label to ID
            seq.id += f".{i}"
            seq.description = ''

            seqs.append(seq)

    # Write final combined output 
    output_file = suffix(source, strand) 

    with open(output_file, 'w') as handle:
        SeqIO.write(seqs, handle, 'fasta')

    print(f"Deaminated sequences written to: {output_file}")

##############################################################################

def parse(whichFile):
    total_alignments = 0
    for record in NCBIXML.parse(open(whichFile, 'r')):
        # Get best alignments
        EI = [(alignment.hsps[0].expect, i)
              for i, alignment in enumerate(record.alignments)]
        EI.sort()
        for e, i in EI:
            if e / 10 <= EI[0][0]:
                total_alignments += 1
                alignment = record.alignments[i]
                hsp = alignment.hsps[0]
                
                Q = record.query.split()[0].split('.')
                S = alignment.accession.split('.')
                qStr, sStr = Q[-1], S[-1]
                qID, sID = '.'.join(Q[:-1]), '.'.join(S[:-1])

                yield sID, qID, qStr, sStr, hsp.expect, hsp.sbjct_start, hsp.sbjct_end, hsp.query_start, hsp.query_end, hsp.sbjct, hsp.query
            else:  
                break
    print(f"Debug: parse() yielded {total_alignments} total alignments")

##############################################################################

def put_data(seqs, refs, dest, mask=False, strand='ab', update=False):
    try:
        # Ensure all identifiers are unique & count input sequences        
        seqs, refs = check_IDs(seqs), check_IDs(refs)
        
        # IMPORTANT: Index original references BEFORE deamination
        original_refs_index = SeqIO.to_dict(SeqIO.parse(refs, 'fasta'))
        
        # files: (D)atabase, (X)ML blastn output, (Q)uery deaminated, (S)ubject deaminated
        D, X = suffix(dest, 'db'), suffix(seqs, 'xml')
        Q, S = suffix(seqs, strand), suffix(refs, 'ab')  # Q and S are the .ab files
        
        # Initial setup (skip for second phase of paired-end)
        if not update:
            if path.exists(D):
                remove(D)
                
        if not path.exists(X):
            start = time()
            print('Converting bases ...')
            deaminate(seqs, strand)
            deaminate(refs, 'ab')
            print(f"DEBUG: Query file created: {Q}")
            print(f"DEBUG: Subject file created: {S}")
            
            # Check if the files actually exist and have content
            if os.path.exists(Q):
                q_count = sum(1 for _ in SeqIO.parse(Q, 'fasta'))
                print(f"DEBUG: Query sequences in {Q}: {q_count}")
            else:
                print(f"DEBUG: Query file {Q} does not exist!")
                
            if os.path.exists(S):
                s_count = sum(1 for _ in SeqIO.parse(S, 'fasta'))
                print(f"DEBUG: Subject sequences in {S}: {s_count}")
            else:
                print(f"DEBUG: Subject file {S} does not exist!")
                
            print('in', time() - start)
            print('Aligning ...')
            
            start = time()
            run_BLAST(Q, S, X, mask)
            print("BLAST finished — starting filing data")
            
            # Use in-memory indexing instead of SQLite-based indexing
            print("Creating sequence index...")
            seqs_index = {}

            # Index original sequences (needed for reamination)
            for record in SeqIO.parse(seqs, 'fasta'):
                seqs_index[record.id] = record
                # Also index with strand suffixes for lookup
                seqs_index[record.id + '.a'] = record
                seqs_index[record.id + '.b'] = record

            for record in SeqIO.parse(refs, 'fasta'):
                seqs_index[record.id] = record
                seqs_index[record.id + '.a'] = record
                seqs_index[record.id + '.b'] = record

            # Additionally, index deaminated sequences for strand detection
            for record in SeqIO.parse(Q, 'fasta'):
                seqs_index[record.id + '_deaminated'] = record

            for record in SeqIO.parse(S, 'fasta'):
                seqs_index[record.id + '_deaminated'] = record
            
            # Index Q file
            if os.path.exists(Q):
                for record in SeqIO.parse(seqs, 'fasta'):
                    seqs_index[record.id] = record
            
            # Index S file  
            if os.path.exists(S):
                for record in SeqIO.parse(refs, 'fasta'):
                    seqs_index[record.id] = record
                    
            print(f"Total sequences in index: {len(seqs_index)}")
            print('in', time() - start)
            
        else:
            print('Warning: XML file already exists')
            print('(delete & run again for fresh alignment).')
            # Still need the index for filing data
            print("Creating sequence index...")
            seqs_index = {}
            
            # Index Q file
            if os.path.exists(Q):
                for record in SeqIO.parse(Q, 'fasta'):
                    seqs_index[record.id] = record
            
            # Index S file  
            if os.path.exists(S):
                for record in SeqIO.parse(S, 'fasta'):
                    seqs_index[record.id] = record
            
        print('Filing data ...')
        start = time()
        
        # Set up SQL database
        conn = sqlite3.connect(D)
        curs = conn.cursor()
        curs.execute('CREATE TABLE IF NOT EXISTS '
                     'records (id INTEGER PRIMARY KEY NOT NULL, '
                     'read, expt, locus, expect REAL, strand, sequence)')
        curs.execute('CREATE TABLE IF NOT EXISTS '
                     'loci (id INTEGER PRIMARY KEY NOT NULL, locus, sequence)')
        curs.execute('CREATE INDEX IF NOT EXISTS loc_idx ON loci (locus)')
        curs.execute('CREATE INDEX IF NOT EXISTS els_idx ON records '
                     '(expt, locus, strand)')
        curs.execute('CREATE INDEX IF NOT EXISTS rls_idx ON records '
                     '(read, locus, strand)')
        
        failed_lookups = 0
        for S_hit, Q_hit, qStr, sStr, e, s5, s3, q5, q3, sBS, qBS in parse(X):
            # Use the clean IDs from parse()
            qID = Q_hit.split(' ', 1)[0].strip()
            sID = S_hit.split(' ', 1)[0].strip()
            
            # Get original sequences for reamination (strand info already extracted)
            qSeq = seqs_index.get(qID)
            sSeq = seqs_index.get(sID)
                        
            if qSeq is None or sSeq is None:
                failed_lookups += 1
                if failed_lookups <= 10:  # Show first 10 failures
                    print(f"Failed lookup: qID='{qID}', sID='{sID}'")
                    print(f"  Available query IDs: {[k for k in seqs_index.keys() if 'query_pattern' in k][:5]}")
                # Only show warning for truly missing sequences
                continue
            
            if s5 > s3:
                continue
            
            # Reverse complement if it was the 'b' strand
            if qStr == 'b':
                qSeq = qSeq.reverse_complement()
            if sStr == 'b':
                sSeq = sSeq.reverse_complement()
            
            # Reaminate (convert back from bisulfite space)
            seq = reaminate(*[list(x.upper()) for x in [
                sSeq[s5 - 1:s3], qSeq[q5 - 1:q3], sBS, qBS]])
            # Pad with dashes to full reference length
            seq = '-' * (s5 - 1) + seq + '-' * (len(sSeq) - s3)
            # Return to original reference orientation
            if sStr == 'b':
                seq = seq.reverse_complement()
            # Extract read ID and experiment (if present)
            qID_parts = qID.split('|')
            read_id, expt = (qID_parts + [''])[:2]
            if update:
                # Update existing record (for paired-end merging)
                done = False
                for n, old in curs.execute(
                    'SELECT id, sequence FROM records WHERE read=? '
                    'AND locus=? AND strand=?', (read_id, sID, sStr)):
                    new = list(old)
                    for i, j in enumerate(str(seq).upper()):
                        if j in '-N' or j == new[i]:
                            continue
                        elif new[i] in '-N':
                            new[i] = j
                        else:
                            new[i] = ambigIUPAC.get(''.join(sorted(new[i] + j)), 'N')
                    curs.execute('UPDATE records SET sequence = ? WHERE id = ?',
                                 (''.join(new), n))
                    done = True
                if done:
                    continue
            # Insert new record
            curs.execute('INSERT INTO records VALUES (NULL,?,?,?,?,?,?)',
                         (read_id, expt, sID, e, sStr, str(seq).upper()))
        print('in', time() - start)
        print(f"Total failed sequence lookups: {failed_lookups}")
        print('Adding refs ...')
        
        start = time()
        existing_loci = {row[0] for row in curs.execute('SELECT locus FROM loci')}
        used_loci = {row[0] for row in curs.execute('SELECT locus FROM records GROUP BY locus')}
        
        for locus_id in set(used_loci) - set(existing_loci):
            # Use the original reference sequences (before deamination)
            if locus_id in original_refs_index:
                ref_seq = str(original_refs_index[locus_id].seq)
                curs.execute('INSERT INTO loci VALUES (NULL, ?, ?)', (locus_id, ref_seq))
                print(f"  Added original reference for locus: {locus_id}")
            else:
                print(f"  Warning: Could not find original reference for locus '{locus_id}'")
                
        print('in', time() - start)
        
        # Final summary
        read_count = curs.execute('SELECT COUNT(DISTINCT read) FROM records').fetchone()[0]
        locus_count = curs.execute('SELECT COUNT(DISTINCT locus) FROM records').fetchone()[0]
        print(f'Aligned {read_count} sequences to {locus_count} loci')
        print()
        conn.commit()
        conn.close()
        
    except Exception as e:
        print(f"Error in put_data: {e}")
        raise

##############################################################################
            
def reaminate(seq1, seq2, conv1, conv2):
    I = range(len(conv1))

    # align original sequences by BLAST of fully-converted sequences
    for i in I:
        if conv1[i] == '-':
            seq1[i:i] = ['-']
        if conv2[i] == '-':
            seq2[i:i] = ['-']    
    assert len(conv1) == len(seq1) == len(conv2) == len(seq2)

    # clone seq2 as base for final output
    final = list(seq2)

    # resolve indels
    for i in I:
        if seq1[i] == 'C' != seq2[i]:

            # deletion
            if seq2[i] == '-':
                if 'C' in seq2[:i + 2][-3:]:
                    final[i] = 'C'
                elif 'T' in seq2[:i + 2][-3:]:
                    final[i] = 'T'
                else:
                    final[i] = 'N'
            elif seq2[i] == 'T':

                # insertion
                if 'C' in seq2[:i + 2][-3:] and '-' in seq1[:i + 2][-3:]:
                    final[i] = 'C'

                # neither
                else:
                    final[i] = 't'

    # remove insertions
    return Seq(''.join([final[i] for i in I if seq1[i] != '-']))

##############################################################################

import glob
import subprocess

def run_BLAST(query, subject, XML, mask):
    # incorporate lowercase masking
    if mask:
        # Clean old DB files (keep your existing code)
        old_db_files = glob.glob(subject + '.*')
        for old_file in old_db_files:
            if old_file != subject:
                if os.path.exists(old_file):
                    os.remove(old_file)

        call([c2bm, '-in', subject, '-out', suffix(subject, 'asnb'),
              '-masking_algorithm', 'repeat',
              '-masking_options', 'repeatmasker,default',
              '-outfmt', 'maskinfo_asn1_bin', '-parse_seqids'])
        call([mbdb, '-dbtype', 'nucl', '-parse_seqids', '-hash_index', '-in', subject,
              '-mask_data', suffix(subject, 'asnb')])

    # makeblastdb (non-masked branch)
    else:
        # Clean old DB files
        old_db_files = glob.glob(subject + '.*')
        for old_file in old_db_files:
            if old_file != subject:
                if os.path.exists(old_file):
                    os.remove(old_file)

        # Add -hash_index here
        call([mbdb, '-dbtype', 'nucl', '-parse_seqids', '-hash_index', '-in', subject])

    cmd = [
        bn, 
        '-task', 'megablast', 
        '-dust', '50 64 1',
        '-evalue', '0.001', 
        '-max_target_seqs', '10',
        '-outfmt', '5',
        '-query', query,
        '-db', subject,
        '-out', XML,
        '-num_threads', '8'
    ]

    if mask:
        cmd += ['-db_soft_mask', '40']
    
    print("Running BLAST command:")
    print(' '.join(cmd))

    result = call(cmd)
    if result != 0:
        raise RuntimeError(f"blastn failed with return code {result}")
            
##############################################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('seqs', help='FASTA-format file of read sequences')
    parser.add_argument('refs', help='FASTA-format file of reference sequences')
    parser.add_argument('-Save', help='output file name')
    parser.add_argument('-mask', action='store_true', help='enable lowercase masking in references')
    parser.add_argument('-pair', type=str, help='2nd FASTA-format file for paired-end sequencing')
    parser.add_argument('-update', action='store_true')
    parser.add_argument('-1', '--a1_or_b1', action='append_const', const='a',
                        dest='strand', help='convert BS-read G to A, i.e. seq. primer=a1 or b1')
    parser.add_argument('-2', '--a2_or_b2', action='append_const', const='b',
                        dest='strand', help='convert BS-read C to T, i.e. seq. primer=a2 or b2')
    args = parser.parse_args()

    # if strand unspecified, try both
    if not args.strand:
        args.strand = 'ab'
    else:
        args.strand = ''.join(sorted(args.strand)).lower()

    # label output
    if not args.Save:
        args.Save = path.join(path.dirname(args.seqs), '.'.join([
            path.split(args.seqs)[-1].split('.')[0][:10],
            path.split(args.refs)[-1].split('.')[0][:10], '']))
    for arg, val in sorted(vars(args).items()):
        print(arg, val)
    print()

    put_data(args.seqs, args.refs, args.Save, args.mask, args.strand, args.update)
    
    if args.pair:
        args.strand = {'a': 'b', 'b': 'a', 'ab': 'ab'}[args.strand]
        put_data(args.pair, args.refs, args.Save, args.mask, args.strand, True)

