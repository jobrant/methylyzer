# Methylyzer

**Bisulfite sequencing analysis pipeline for single-molecule methylation mapping**

Methylyzer is a Python-based pipeline for processing, aligning, and analyzing bisulfite-converted sequencing reads from single-molecule methylation assays such as MAPit (Methyltransferase Accessibility Protocol for individual Templates). It produces per-molecule methylation maps and frequency tables suitable for visualization with [methylscaper](https://bioconductor.org/packages/methylscaper/) or custom plotting tools.

Methylyzer is the successor to the reAminator pipeline (Darst & Riva, University of Florida).

![Example methylation map](mapit_image.png)
*Example output: Side-by-side heatmaps of endogenous CG methylation (HCG, left, red=methylated) and accessibility (GCH, right, yellow=accessible) at a single locus. Each row is one molecule; molecules are ordered by PCA.*

---

## Overview

The pipeline consists of three core steps, orchestrated by `methylyzer.py` for SLURM-based HPC environments:

| Step | Script | Description |
|------|--------|-------------|
| 1 | `bsAlign_updated.py` | Align bisulfite reads to references via BLAST |
| 2 | `bsExtract.py` | Filter, QC, and extract single-molecule FASTA files |
| 3 | `bsFreqs.py` | Calculate methylation frequencies and generate visualization matrices |

Optional downstream plotting is handled by `batch_methylscaper_plots.R` and `plot_methylscaper.sh`.

---

## Installation

### Dependencies

**Python 3.8+**
- [Biopython](https://biopython.org/)

**External tools (Step 1)**
- NCBI BLAST+ (`blastn`, `makeblastdb`, `convert2blastmask`)
- `faconvert` (for in silico bisulfite conversion)

**R 4.0+ (for plotting, optional)**
- `seriation`
- `data.table`
- `Rfast`

On HPC systems with module environments:
```bash
module load python/3.8
module load ncbi_blast/2.15.0
module load R/4.2  # for plotting
```

### Setup

Clone the repository and ensure the scripts are accessible:
```bash
git clone https://github.com/your-username/Methylyzer.git
cd Methylyzer
```

No installation step is needed beyond having the dependencies available. The pipeline scripts can be run directly or called through `methylyzer.py`.

---

## Quick Start

### Full pipeline via SLURM wrapper

```bash
# Standard MAPit analysis (HCG + GCH)
python methylyzer.py -r references.fa sample_001/

# Process multiple samples
python methylyzer.py -r references.fa sample_001/ sample_002/ sample_003/

# Wildcard
python methylyzer.py -r references.fa sample_*/
```

### Running steps individually

```bash
# Step 1: Align
python bsAlign_updated.py ccs.sample.fasta references.fa

# Step 2: Extract
python bsExtract.py sample.db --dest sample_dir/ --min-len 90% --min-bs 95

# Step 3: Frequencies
python bsFreqs.py extracted/ -o frequencies/
```

---

## Pipeline Details

### Step 1: Alignment (`bsAlign_updated.py`)

Aligns bisulfite-converted reads to reference sequences using BLAST in a strand-aware manner.

**What it does:**
1. Performs in silico bisulfite conversion of both reads and references (C→T for top strand, G→A for bottom strand)
2. Aligns converted reads to converted references via megaBLAST
3. "Reaminates" — reconstructs the original methylation information by comparing aligned sequences back to unconverted originals
4. Stores alignments in a SQLite database (`.db` file)

**Input:** FASTA file of CCS reads (expected naming: `ccs.samplename.fasta`) and a reference FASTA.

**Output:** SQLite database containing aligned reads with methylation state preserved.

```bash
python bsAlign_updated.py reads.fasta references.fa

# Options
python bsAlign_updated.py reads.fasta references.fa -mask      # Enable lowercase masking
python bsAlign_updated.py reads.fasta references.fa -pair R2.fa # Paired-end
python bsAlign_updated.py reads.fasta references.fa -1          # Top strand only
python bsAlign_updated.py reads.fasta references.fa -2          # Bottom strand only
```

### Step 2: Extraction (`bsExtract.py`)

Extracts quality-filtered, single-molecule FASTA files from the alignment database.

**Filters applied:**
- **Read length** — Minimum ungapped length (bp or percentage of reference, default: `90%`)
- **Bisulfite conversion** — Minimum conversion rate at non-CG/GC cytosines (default: `95%`)
- **Deduplication** — Optional removal of duplicate methylation patterns (`--uniques`)
- **Subsampling** — Optional random downsampling for visualization (default: `1000` when called via methylyzer)

**Output structure:**
```
sample_dir/
├── extracted/                  # Full filtered dataset
│   ├── A_locus_sample.fa       # Top strand reads (reference + reads)
│   ├── B_locus_sample.fa       # Bottom strand reads
│   └── report.tsv              # QC report
└── subsampled/                 # Random subsample (if --subsample > 0)
    ├── A_locus_sample.fa
    ├── B_locus_sample.fa
    └── report.tsv
```

Each FASTA file contains the reference sequence as the first record, followed by aligned reads. This format is directly compatible with methylscaper and bsFreqs.

```bash
python bsExtract.py sample.db --min-len 90% --min-bs 95 --strand ab

# With deduplication
python bsExtract.py sample.db --uniques

# With subsampling for plotting
python bsExtract.py sample.db --subsample 1000 --subsample-seed 42
```

### Step 3: Frequency Analysis (`bsFreqs.py`)

Generates methylation frequency tables and optional coded CSV matrices for visualization.

**Frequency tables (.tsv):** Per-position base frequencies (A, C, G, T) at each methylatable cytosine. Methylated cytosines retain C; unmethylated cytosines are converted to T by bisulfite treatment.

**CSV map files (.csv):** Full-length molecule × position matrices with numeric codes for methylscaper-style visualization.

```bash
# Standard analysis
python bsFreqs.py extracted/ -o frequencies/

# With CSV maps for plotting
python bsFreqs.py subsampled/ -o frequencies_subsampled/ --csv

# Custom methylation contexts (e.g., no endogenous CG methylation)
python bsFreqs.py extracted/ -o frequencies/ --sites CG GC
```

---

## Methylation Contexts

By default, Methylyzer analyzes **HCG** (endogenous CG methylation) and **GCH** (enzyme accessibility) contexts, which exclude ambiguous GCG sites.

### Supported site patterns

The cytosine position is automatically detected from the pattern:

| Pattern | Regex | Cytosine Position | Description |
|---------|-------|:-----------------:|-------------|
| `HCG` | `[ACT]CG` | 1 | CG methylation, excluding GCG |
| `GCH` | `GC[ACT]` | 1 | GC accessibility, excluding GCG |
| `CG` | `CG` | 0 | All CG dinucleotides (including GCG) |
| `GC` | `GC` | 1 | All GC dinucleotides (including GCG) |
| `WCG` | `[AT]CG` | 1 | CG preceded by A or T only |

### IUPAC ambiguity codes

| Code | Bases | Meaning |
|------|-------|---------|
| H | A, C, T | Not G |
| W | A, T | Weak |
| S | C, G | Strong |
| M | A, C | Amino |
| K | G, T | Keto |
| R | A, G | Purine |
| Y | C, T | Pyrimidine |
| B | C, G, T | Not A |
| D | A, G, T | Not C |
| V | A, C, G | Not T |
| N | A, C, G, T | Any |

### When to use dinucleotide vs trinucleotide contexts

For standard MAPit experiments, use the defaults (`HCG` + `GCH`). These exclude GCG sites, which are ambiguous because both endogenous methylation and enzyme accessibility could produce a signal.

If your experiment has **no endogenous CG methylation** (e.g., certain assay designs), you can safely include GCG sites:
```bash
python bsFreqs.py extracted/ --sites CG GC
# or via methylyzer:
python methylyzer.py -r refs.fa sample_001/ --sites CG GC
```

---

## SLURM Workflow (`methylyzer.py`)

`methylyzer.py` orchestrates all three steps as SLURM jobs with dependency chaining.

### Default behavior

```bash
python methylyzer.py -r references.fa sample_001/
```

This will:
1. **Align** CCS reads to references (Step 1)
2. **Extract** filtered reads to `extracted/` and subsample 1000 reads to `subsampled/` (Step 2)
3. **Calculate frequencies** on both full and subsampled data (Step 3)

### Output directory structure

```
sample_001/
├── ccs.sample.fasta            # Input reads
├── sample.db                   # Alignment database (Step 1)
├── extracted/                  # Full filtered reads (Step 2)
│   ├── A_locus_sample.fa
│   ├── B_locus_sample.fa
│   └── report.tsv
├── subsampled/                 # Subsampled reads (Step 2)
│   ├── A_locus_sample.fa
│   ├── B_locus_sample.fa
│   └── report.tsv
├── frequencies/                # From full data (Step 3)
│   ├── HCG-cytosines-A-locus.tsv
│   └── GCH-cytosines-A-locus.tsv
└── frequencies_subsampled/     # From subsampled data (Step 3)
    ├── HCG-cytosines-A-locus.tsv
    ├── GCH-cytosines-A-locus.tsv
    └── A-locus_coded_map.csv   # CSV maps for visualization
```

### Full argument reference

```
usage: methylyzer.py [-h] -r REFERENCES [-s STEPS] [--min-len MIN_LEN]
                     [--min-bs MIN_BS] [--uniques] [--strand STRAND]
                     [--subsample SUBSAMPLE] [--subsample-seed SUBSAMPLE_SEED]
                     [--freq-source {full,subsampled,both}]
                     [--sites SITES [SITES ...]] [--csv]
                     [--account ACCOUNT] [--qos QOS] [--partition PARTITION]
                     [--time TIME] [--mem MEM] [--cpus CPUS] [--email EMAIL]
                     directories [directories ...]
```

**Pipeline arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `-r`, `--references` | *required* | Reference sequences FASTA file |
| `-s`, `--steps` | `123` | Steps to run (1=align, 2=extract, 3=freqs) |

**bsExtract arguments (Step 2):**

| Argument | Default | Description |
|----------|---------|-------------|
| `--min-len` | `90%` | Minimum read length (bp or percentage) |
| `--min-bs` | `95` | Minimum bisulfite conversion rate (0-100) |
| `--uniques` | off | Deduplicate by methylation pattern |
| `--strand` | `ab` | Strands to extract (a/b/ab) |
| `--subsample` | `1000` | Random subsample size (0 to disable) |
| `--subsample-seed` | `42` | Random seed for reproducibility |

**bsFreqs arguments (Step 3):**

| Argument | Default | Description |
|----------|---------|-------------|
| `--sites` | `HCG GCH` | Methylation contexts to analyze |
| `--csv` | off | Generate CSV matrices for plotting |
| `--freq-source` | `both` | Run on `full`, `subsampled`, or `both` |

**SLURM arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `--account` | `kladde` | SLURM account |
| `--qos` | `kladde-b` | SLURM QOS |
| `--partition` | none | SLURM partition |
| `--time` | `12:00:00` | Job time limit |
| `--mem` | `32G` | Memory per job |
| `--cpus` | `8` | CPUs per job |
| `--email` | none | Email for notifications |

### Common usage patterns

```bash
# Skip alignment, just re-extract and recompute frequencies
python methylyzer.py -r refs.fa -s 23 sample_001/

# Only run frequencies on existing extraction
python methylyzer.py -r refs.fa -s 3 sample_001/

# Disable subsampling (full data only)
python methylyzer.py -r refs.fa --subsample 0 --freq-source full sample_001/

# Custom subsample size with deduplication
python methylyzer.py -r refs.fa --subsample 500 --uniques sample_001/

# Different methylation contexts
python methylyzer.py -r refs.fa --sites CG GC sample_001/

# Custom SLURM settings
python methylyzer.py -r refs.fa --account mylab --time 8:00:00 --mem 64G sample_*/
```

> **Note:** When using `--sites` with `nargs='+'`, place the positional `directories` argument first or use `--` to separate:
> ```bash
> python methylyzer.py -r refs.fa sample_001/ --sites CG GC
> # or
> python methylyzer.py -r refs.fa --sites CG GC -- sample_001/
> ```

---

## Visualization

### Batch plotting with R

`batch_methylscaper_plots.R` generates methylscaper-style heatmaps from the coded CSV map files. It is self-contained and does not require the methylscaper R package.

```bash
# Basic usage
Rscript batch_methylscaper_plots.R /path/to/frequencies_subsampled/

# With options
Rscript batch_methylscaper_plots.R /path/to/frequencies_subsampled/ /path/to/plots/ \
    --format both --dpi 300 --width 10 --height 12
```

| Option | Default | Description |
|--------|---------|-------------|
| `--width` | 8 | Plot width (inches) |
| `--height` | 10 | Plot height (inches) |
| `--format` | `png` | Output format: `png`, `pdf`, or `both` |
| `--dpi` | 150 | Resolution for PNG |

### SLURM batch plotting

`plot_methylscaper.sh` automates plotting across samples via SLURM:

```bash
# From project directory
sbatch plot_methylscaper.sh sample_001
sbatch plot_methylscaper.sh sample_001 sample_002 sample_003
sbatch plot_methylscaper.sh sample_*
```

Output is placed in `frequencies_subsampled/plots/` within each sample directory.

### Using methylscaper directly

The extracted FASTA files and CSV maps are also compatible with the [methylscaper](https://bioconductor.org/packages/methylscaper/) R/Bioconductor package and its Shiny app for interactive exploration.

---

## File Descriptions

### Core pipeline

| File | Description |
|------|-------------|
| `methylyzer.py` | SLURM workflow orchestrator |
| `bsAlign_updated.py` | Bisulfite-aware BLAST alignment (Python 3) |
| `bsExtract.py` | Read extraction, filtering, QC, and subsampling |
| `bsFreqs.py` | Methylation frequency and visualization matrix generation |

### Plotting

| File | Description |
|------|-------------|
| `batch_methylscaper_plots.R` | Command-line R script for batch heatmap generation |
| `plot_methylscaper.sh` | SLURM wrapper for batch plotting |
| `seqPlot.R` | `plotSequence()` function for methylscaper-style heatmaps |

### Legacy / Reference

These files are from the original reAminator pipeline (Python 2.7, Darst & Riva 2013) and the MethylMapper tool (Riva 2017). They are retained for reference but are not used by the current pipeline.

| File | Description |
|------|-------------|
| `bsAlign_fixed_recent.py` | Original bsAlign (Python 2, uses fastools + ConfigParser) |
| `bsDraw.py` | Original extraction + omelet visualization |
| `meTools.py` | Methylation analysis library (Python 2) |
| `meMapper.py` | MethylMapper wrapper |
| `MethMap.py` | MethylMapper map/frequency engine |
| `methylmapper.py` | MethylMapper main application |
| `preprocessSingleMolecule.R` | R-based alignment for methylscaper |
| `initialOrder.R` | Molecule ordering (PCA/seriation) for methylscaper |

---

## Citation

If you use Methylyzer in your research, please cite the underlying tools and methods:

- MAPit: Darst RP, Pardo CE, Ai L, Brown KD, Kladde MP. *Bisulfite sequencing of DNA from fixed tissues.* Methods Mol Biol. 2010.
- methylscaper: Parker MA, et al. *methylscaper: an R/Bioconductor package for joint visualization of DNA methylation and nucleosome occupancy at the single-molecule level.* Bioinformatics. 2022.

---

## Authors

- **Jason Orr Brant** — Pipeline modernization and Methylyzer development (2025–2026), University of Florida Health Cancer Institute
- **Russell Darst** — Original reAminator pipeline (2013), University of Florida
- **Alberto Riva** — MethylMapper (2017), ICBR Bioinformatics, University of Florida

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
