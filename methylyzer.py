#!/usr/bin/env python3
"""
Methylyzer - Bisulfite Sequencing Analysis Workflow

Successor to reAminator, this tool orchestrates:
1. bsAlign_updated.py  - Alignment of bisulfite reads
2. bsExtract.py - Extraction of quality reads  
3. bsFreqs.py   - Frequency analysis of methylation sites

Output structure:
  sample_dir/
  ├── *.db                    # bsAlign database
  ├── extracted/              # Full filtered dataset
  │   ├── A_locus_sample.fa
  │   └── report.tsv
  ├── subsampled/             # Random subsample for plotting
  │   ├── A_locus_sample.fa
  │   └── report.tsv
  ├── frequencies/            # Frequency tables from full data
  │   ├── HCG-cytosines-A-locus.tsv
  │   └── GCH-cytosines-A-locus.tsv
  └── frequencies_subsampled/ # Freq tables + CSV maps from subsampled
      ├── HCG-cytosines-A-locus.tsv
      ├── GCH-cytosines-A-locus.tsv
      ├── HCG-A-locus_map.csv
      └── GCH-A-locus_map.csv

Author: Jason Orr Brant, 2026
"""

import argparse
import os
import sys
import subprocess
import time
import tempfile
from pathlib import Path

class MethylyzerSLURM:
    def __init__(self):
        self.script_dir = Path(__file__).parent
        self.steps = [1, 2, 3]
        # Default SLURM settings for Marie in Kladde lab (can be overridden)
        self.slurm_settings = {
            'account': 'kladde',
            'qos': 'kladde-b',
            'partition': None,
            'time': '12:00:00',
            'mem': '32G',
            'cpus': 8,
            'email': 'm.gauthier@ufl.edu'
        }
        # bsExtract settings
        self.extract_settings = {
            'min_len': '90%',
            'min_bs': 95,
            'uniques': False,
            'strand': 'ab',
            'subsample': 1000, # Default: subsample 1000 reads
            'subsample_seed': 42 # The answer!
        }
        #bsFreqs settings
        self.freqs_settings = {
            'freq_source': 'both',
            'sites': ['HCG', 'GCH']
        }
        
    def create_sbatch_script(self, job_name, commands, dependencies=None):
        """Create a temporary sbatch script."""
        script_content = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --time={self.slurm_settings['time']}
#SBATCH --mem={self.slurm_settings['mem']}
#SBATCH --cpus-per-task={self.slurm_settings['cpus']}
#SBATCH -o LOGS/{job_name}.%A_%a.out
#SBATCH -e LOGS/{job_name}.%A_%a.err
"""
        
        # Add optional SLURM settings
        if self.slurm_settings['account']:
            script_content += f"#SBATCH --account={self.slurm_settings['account']}\n"
        if self.slurm_settings['qos']:
            script_content += f"#SBATCH --qos={self.slurm_settings['qos']}\n"
        if self.slurm_settings['partition']:
            script_content += f"#SBATCH --partition={self.slurm_settings['partition']}\n"
        if self.slurm_settings['email']:
            script_content += f"#SBATCH --mail-user={self.slurm_settings['email']}\n"
            script_content += f"#SBATCH --mail-type=END,FAIL\n"
        if dependencies:
            script_content += f"#SBATCH --dependency=afterok:{':'.join(map(str, dependencies))}\n"
            
        script_content += """
# Load required modules
module purge
module load python/3.8
module load ncbi_blast

# Set up environment
set -e  # Exit on any error

echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"

"""
        
        # Add the actual commands
        for cmd in commands:
            script_content += f"{cmd}\n"
            
        script_content += '\necho "Job completed at: $(date)"\n'
        
        return script_content
    
    def submit_job(self, script_content, job_name):
        """Submit a job to SLURM and return job ID."""
        # Create temporary script file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as f:
            f.write(script_content)
            script_file = f.name
        
        try:
            # Submit job
            cmd = ['sbatch', script_file]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Extract job ID from output 
            job_id = result.stdout.strip().split()[-1]
            print(f"✓ Submitted {job_name} (Job ID: {job_id})")
            return int(job_id)
            
        except subprocess.CalledProcessError as e:
            print(f"✗ Failed to submit {job_name}: {e.stderr}")
            return None
        finally:
            # Clean up temporary script
            os.unlink(script_file)
    
    def wait_for_jobs(self, job_ids, check_interval=30):
        """Wait for all jobs to complete."""
        if not job_ids:
            return True
            
        print(f"\nWaiting for {len(job_ids)} job(s) to complete...")
        print(f"Job IDs: {job_ids}")
        
        while True:
            # Check job status
            try:
                cmd = ['squeue', '-j', ','.join(map(str, job_ids)), '--noheader', '--format=%i,%T']
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode != 0:
                    # Jobs are done (not in queue anymore)
                    break
                    
                running_jobs = []
                for line in result.stdout.strip().split('\n'):
                    if line:
                        job_id, status = line.split(',')
                        running_jobs.append((job_id, status))
                
                if not running_jobs:
                    break
                    
                print(f"Still running: {len(running_jobs)} jobs", end='')
                for job_id, status in running_jobs:
                    print(f" {job_id}({status})", end='')
                print()
                
                time.sleep(check_interval)
                
            except KeyboardInterrupt:
                print("\nCancelling remaining jobs...")
                for job_id in job_ids:
                    subprocess.run(['scancel', str(job_id)], capture_output=True)
                return False
        
        # Check if jobs completed successfully
        failed_jobs = []
        for job_id in job_ids:
            cmd = ['sacct', '-j', str(job_id), '--format=JobID,State', '--noheader', '--parsable2']
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if 'FAILED' in result.stdout or 'CANCELLED' in result.stdout or 'TIMEOUT' in result.stdout:
                failed_jobs.append(job_id)
        
        if failed_jobs:
            print(f"✗ Jobs failed: {failed_jobs}")
            return False
        else:
            print("✓ All jobs completed successfully!")
            return True
    
    def process_directory(self, input_dir, refs_file):
        """Process one directory through the workflow."""
        input_path = Path(input_dir).resolve()
        refs_file_path = Path(refs_file).resolve() 
        
        print(f"\n{'='*60}")
        print(f"Processing directory: {input_path}")
        print(f"Using reference: {refs_file_path}")
        print(f"{'='*60}")
        
        if not input_path.exists():
            print(f"Error: Directory does not exist: {input_path}")
            return False
        
        if not refs_file_path.exists():
            print(f"Error: Reference file does not exist: {refs_file_path}")
            return False
        
        job_ids = []
        
        # Step 1: bsAlign
        if 1 in self.steps:
            # Find input FASTA (CCS reads only)
            ccs_files = list(input_path.glob("ccs.*.fasta")) + list(input_path.glob("ccs.*.fa"))

            if not ccs_files:
                print(f"Error: No CCS files found in {input_path}")
                print("Expected naming: ccs.samplename.fasta or ccs.samplename.fa")
                # Show what files ARE there for debugging
                all_fasta = list(input_path.glob("*.fa")) + list(input_path.glob("*.fasta"))
                if all_fasta:
                    print(f"Found FASTA files: {[f.name for f in all_fasta]}")
                return False
            
            if len(ccs_files) > 1:
                print(f"Multiple CCS files found: {[f.name for f in ccs_files]}")
                print(f"Using: {ccs_files[0].name}")
            
            input_fasta = str(ccs_files[0])
            print(f"  CCS reads: {ccs_files[0].name}")
            print(f"  Reference: {refs_file_path.name}")
            
            commands = [
                f"cd {input_path}",
                f"python /orange/kladde/KLADDE_SCRIPTS/Methylyzer/bsAlign_updated.py {input_fasta} {refs_file_path }"
            ]
            
            script = self.create_sbatch_script(f"bsAlign_{input_path.name}", commands)
            job_id = self.submit_job(script, f"bsAlign_{input_path.name}")
            if job_id:
                job_ids.append(job_id)
            else:
                return False
        
        # Step 2: bsExtract (depends on step 1)
        if 2 in self.steps:
            extract_dir = input_path / "extracted"
            
            # Build bsExtract command with custom arguments
            extract_cmd = f"python /orange/kladde/KLADDE_SCRIPTS/Methylyzer/bsExtract.py *.db"
            
            # Add optional arguments
            extract_cmd += f" --min-len {self.extract_settings['min_len']}"
            extract_cmd += f" --min-bs {self.extract_settings['min_bs']}"
            extract_cmd += f" --strand {self.extract_settings['strand']}"
            
            if self.extract_settings['uniques']:
                extract_cmd += " --uniques"

            # Add subsampling if enabled
            if self.extract_settings['subsample'] and self.extract_settings['subsample'] > 0:
                extract_cmd += f" --subsample {self.extract_settings['subsample']}"
                extract_cmd += f" --subsample-seed {self.extract_settings['subsample_seed']}"
            
            commands = [
                f"cd {input_path}",
                extract_cmd
            ]

            print(f"  Debug - self.extract_settings: {self.extract_settings}")
            print(f"  Debug - uniques setting: {self.extract_settings['uniques']}")
            print(f"  Debug - extract_cmd: {extract_cmd}")
            
            dependencies = job_ids if 1 in self.steps else None
            script = self.create_sbatch_script(f"bsExtract_{input_path.name}", commands, dependencies)
            job_id = self.submit_job(script, f"bsExtract_{input_path.name}")
            if job_id:
                job_ids.append(job_id)
            else:
                return False
        
        # Step 3: bsFreq (depends on step 2)
        if 3 in self.steps:
            freq_source = self.freqs_settings['freq_source']
            freq_commands = []
            
            # Run on full data (extracted/)
            if freq_source in ('full', 'both'):
                extracted_dir = input_path / "extracted"
                freq_dir = input_path / "frequencies"
                
                # Frequency tables only (no CSV maps for full data)
                freq_cmd = f"python /orange/kladde/KLADDE_SCRIPTS/Methylyzer/bsFreqs.py {extracted_dir} -o {freq_dir}"
                freq_cmd += f" -s {' '.join(self.freqs_settings['sites'])}"
                freq_commands.append(f"echo 'Running bsFreqs on full data...'")
                freq_commands.append(freq_cmd)
                print(f"  bsFreqs (full): {extracted_dir} → {freq_dir}")
            
            # Run on subsampled data (subsampled/)
            if freq_source in ('subsampled', 'both'):
                subsampled_dir = input_path / "subsampled"
                freq_subsampled_dir = input_path / "frequencies_subsampled"
                
                # Frequency tables + CSV maps for subsampled data
                freq_cmd_sub = f"python /orange/kladde/KLADDE_SCRIPTS/Methylyzer/bsFreqs.py {subsampled_dir} -o {freq_subsampled_dir} --csv"
                freq_cmd_sub += f" -s {' '.join(self.freqs_settings['sites'])}"
                
                freq_commands.append(f"echo 'Running bsFreqs on subsampled data (with CSV maps)...'")
                freq_commands.append(freq_cmd_sub)
                print(f"  bsFreqs (subsampled + CSV): {subsampled_dir} → {freq_subsampled_dir}")
            
            if freq_commands:
                commands = [f"cd {input_path}"] + freq_commands
                
                dependencies = [job_ids[-1]] if job_ids else None
                script = self.create_sbatch_script(f"bsFreqs_{input_path.name}", commands, dependencies)
                job_id = self.submit_job(script, f"bsFreqs_{input_path.name}")
                if job_id:
                    job_ids.append(job_id)
                else:
                    return False
        
        return job_ids
    
    def run(self, input_dirs, refs_file):
        """Run the workflow on all directories."""
        print("Methylyzer: Bisulfite Sequencing Analysis Workflow")
        print(f"Steps to run: {self.steps}")
        print(f"Reference file: {refs_file}")
        print(f"Input directories: {len(input_dirs)}")
        print(f"Subsample: {self.extract_settings['subsample']} reads")
        print(f"Freq source: {self.freqs_settings['freq_source']}")
        
        if not Path(refs_file).exists():
            print(f"Error: Reference file does not exist: {refs_file}")
            return False
        
        all_job_ids = []
        
        for input_dir in input_dirs:
            job_ids = self.process_directory(input_dir, refs_file)
            if job_ids:
                all_job_ids.extend(job_ids)
            else:
                print(f"Failed to submit jobs for {input_dir}")
                return False
        
        if all_job_ids:
            print(f"\nSubmitted {len(all_job_ids)} total jobs")
            success = self.wait_for_jobs(all_job_ids)
            return success
        else:
            print("No jobs were submitted")
            return False

def main():
    parser = argparse.ArgumentParser(
        description="Methylyzer: Bisulfite sequencing analysis workflow with SLURM job submission",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Output structure:
  sample_dir/
  ├── extracted/              # Full filtered dataset
  ├── subsampled/             # Random subsample for plotting  
  ├── frequencies/            # Frequency tables from full data
  └── frequencies_subsampled/ # Freq tables + CSV maps from subsampled

Examples:
  # Process single directory (default: subsample=1000, freq-source=both)
  python methylyzer.py -r references.fa sample_001/
  
  # Process multiple directories
  python methylyzer.py -r references.fa sample_*/
  
  # Only run extraction and frequencies (skip alignment)  
  python methylyzer.py -r references.fa -s 23 sample_001/
  
  # Disable subsampling (run full pipeline only)
  python methylyzer.py -r references.fa --subsample 0 --freq-source full sample_001/
  
  # Custom subsample size
  python methylyzer.py -r references.fa --subsample 500 sample_001/
  
  # Custom SLURM settings
  python methylyzer.py -r references.fa --account mylab --time 8:00:00 --mem 32G sample_001/
        """
    )
    
    parser.add_argument('directories', nargs='+',
                       help='Input directories containing FASTA files')
    parser.add_argument('-r', '--references', required=True,
                       help='Reference sequences FASTA file')
    parser.add_argument('-s', '--steps', default='123',
                       help='Steps to run (1=bsAlign, 2=bsExtract, 3=bsFreq)')
    
    # bsExtract arguments
    parser.add_argument('--min-len', default='90%',
                       help='Minimum read length for bsExtract (bp or percentage, e.g. "80%%")')
    parser.add_argument('--min-bs', type=int, default=95,
                       help='Minimum bisulfite conversion rate for bsExtract (0-100, default: 95)')
    parser.add_argument('--uniques', action='store_true',
                       help='Deduplicate reads by methylation pattern in bsExtract')
    parser.add_argument('--strand', default='ab',
                       help='Strands to extract (a/b/ab, default: ab)')
    parser.add_argument('--subsample', type=int, default=1000,
                       help='Randomly subsample reads for plotting (default: 1000). '
                            'Set to 0 to disable subsampling.')
    parser.add_argument('--subsample-seed', type=int, default=42,
                       help='Random seed for reproducible subsampling (default: 42)')
    
    # bsFreqs arguments
    parser.add_argument('--freq-source', choices=['full', 'subsampled', 'both'], default='both',
                       help='Run bsFreqs on full, subsampled, or both (default: both). '
                            'CSV map files are only generated for subsampled data.')
    parser.add_argument('--sites', nargs='+', default=['HCG', 'GCH'],
                        help='Methylation sites for bsFreqs (default: HCG GCH)')

    
    # SLURM options
    parser.add_argument('--account', help='SLURM account')
    parser.add_argument('--qos', help='SLURM QOS')
    parser.add_argument('--partition', help='SLURM partition')
    parser.add_argument('--time', default='12:00:00', help='Job time limit (default: 12:00:00)')
    parser.add_argument('--mem', default='32G', help='Memory per job (default: 32G)')
    parser.add_argument('--cpus', type=int, default=8, help='CPUs per job (default: 8)')
    parser.add_argument('--email', help='Email for job notifications')
    
    args = parser.parse_args()
    
    # Parse steps
    valid_steps = set('123')
    if not all(s in valid_steps for s in args.steps):
        print("Error: Steps must contain only digits 1, 2, and/or 3")
        sys.exit(1)
    
    workflow = MethylyzerSLURM()
    workflow.steps = [int(s) for s in args.steps]
    
    # Update SLURM settings
    if args.account:
        workflow.slurm_settings['account'] = args.account
    if args.qos:
        workflow.slurm_settings['qos'] = args.qos
    if args.partition:
        workflow.slurm_settings['partition'] = args.partition
    if args.time:
        workflow.slurm_settings['time'] = args.time
    if args.mem:
        workflow.slurm_settings['mem'] = args.mem
    if args.cpus:
        workflow.slurm_settings['cpus'] = args.cpus
    if args.email:
        workflow.slurm_settings['email'] = args.email

    # Update bsExtract settings
    workflow.extract_settings['min_len'] = args.min_len
    workflow.extract_settings['min_bs'] = args.min_bs
    workflow.extract_settings['uniques'] = args.uniques
    workflow.extract_settings['strand'] = args.strand

    # Update bsFreqs settings
    workflow.freqs_settings['freq_source'] = args.freq_source
    workflow.freqs_settings['sites'] = args.sites

    # subsample=0 means no subsampling, subsample>0 means subsample to that number
    workflow.extract_settings['subsample'] = args.subsample if args.subsample > 0 else None
    workflow.extract_settings['subsample_seed'] = args.subsample_seed
    
    # Run workflow
    success = workflow.run(args.directories, args.references)
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
