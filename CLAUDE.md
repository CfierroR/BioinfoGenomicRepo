# CLAUDE.md — BioinfoGenomicRepo

This file provides context for AI assistants (Claude Code and similar tools) working in this repository. It describes the project structure, development workflows, and conventions to follow.

---

## Project Overview

**BioinfoGenomicRepo** is a bioinformatics pipeline repository providing end-to-end analysis workflows for three major types of genomic sequencing data:

| Module | Purpose | Primary Language |
|--------|---------|-----------------|
| `RNA/` | RNA-Seq differential expression analysis | Bash, R, Python |
| `ChIPSeq/` | ChIP-Seq peak calling | Snakemake, YAML |
| `Exome/` | Whole Exome Sequencing (WES) visualization | R |
| `test/` | Demo/test scripts and synthetic data | R |

This repository is intended for researchers running genomic analyses. Scripts are modular and often need to be run manually in sequence rather than as a fully automated pipeline (except ChIP-Seq, which uses Snakemake).

---

## Repository Structure

```
BioinfoGenomicRepo/
├── CLAUDE.md                        # This file
├── README.md                        # Minimal project overview
├── .gitignore                       # Ignores test/data and test/output
│
├── RNA/                             # RNA-Seq analysis pipeline
│   ├── README.md
│   ├── bash_pipeline/               # Sequential processing scripts (run in order)
│   │   ├── 00_fastqscreen.sh        # QC screening (FastQ Screen)
│   │   ├── 01_fastqc.sh             # Quality assessment (FastQC)
│   │   ├── 02_trimgalore.sh         # Adapter/quality trimming (Trim Galore)
│   │   ├── 02_trimgalore_parallel.sh # Parallel wrapper for trimming
│   │   ├── 03_star.sh               # Alignment to reference genome (STAR)
│   │   ├── 04_featurecounts.sh      # Read counting (featureCounts)
│   │   ├── 05_SarTools.r            # Differential expression (DESeq2/SARTools)
│   │   ├── 05_TPM-normalization.sh  # TPM normalization
│   │   └── 06_annotateGene.R        # Annotate Ensembl IDs with gene symbols
│   ├── R/                           # Visualization scripts
│   │   ├── MAPlot.R                 # MA plot for DEG results
│   │   └── VolcanoPlot.R            # Volcano plot for DEG results
│   └── python/
│       └── joinTables.py            # Merge DEG tables from multiple studies
│
├── ChIPSeq/                         # ChIP-Seq pipeline (Snakemake automated)
│   ├── README.md                    # Comprehensive setup and usage docs
│   ├── Snakefile                    # Snakemake workflow definition
│   └── config.yaml                  # Pipeline configuration
│
├── Exome/                           # WES visualization scripts
│   ├── README.md
│   └── R/
│       ├── plotAnnotation.R         # Annotation bar chart
│       ├── plot_chrom.r             # Chromosome coverage ideogram
│       └── plot_chrom_diff.r        # Differential chromosome coverage
│
└── test/                            # Demo and testing
    ├── codes/
    │   ├── create_dummy.R           # Generates synthetic test data
    │   └── dotplot.R                # Scatter plot with regression
    ├── data/
    │   └── genes_data.csv           # Synthetic test dataset (gitignored)
    └── outputs/                     # Generated plots (gitignored)
        ├── output.png
        ├── output_2.png
        ├── output_3.png
        └── output_4.png
```

---

## Analysis Pipelines

### RNA-Seq Pipeline

Scripts run **sequentially** in numbered order. Each script expects specific input directory structures.

```
Raw FASTQ files
  → 00_fastqscreen.sh   (outputs to outputs/00_fastqScreen/)
  → 01_fastqc.sh        (outputs to outputs/01_fastqc/)
  → 02_trimgalore.sh    (outputs to outputs/02_trimgalore/)
  → 03_star.sh          (outputs to outputs/03_star/, requires hg38 index)
  → 04_featurecounts.sh (outputs to outputs/04_featurecounts/)
  → 05_SarTools.r       (DESeq2 differential expression analysis)
  → 05_TPM-normalization.sh (produces .tpm.txt files)
  → MAPlot.R / VolcanoPlot.R (visualization of DEG results)
```

**Key parameters in bash scripts:**
- Threads: 30 (STAR, FastQC), 8 (Trim Galore), 10 (featureCounts)
- Reference genome: hg38 (hardcoded paths — must be updated per environment)
- Paired-end sequencing assumed throughout

**Running SARTools DESeq2 (`05_SarTools.r`):**
```bash
Rscript RNA/bash_pipeline/05_SarTools.r \
  --targetFile targets.txt \
  --rawDir counts/ \
  --projectName MyProject \
  --author "Researcher Name"
```
Supports 24 command-line options via `optparse`. Key defaults: `alpha=0.05`, `pAdjustMethod="BH"`, `typeTrans="VST"`.

**Joining DEG tables across studies (`joinTables.py`):**
```bash
python RNA/python/joinTables.py file1.txt file2.txt output.txt
```
Performs an outer join on the `Gene` column from two TSV files.

### ChIP-Seq Pipeline (Snakemake)

Fully automated via Snakemake. Edit `ChIPSeq/config.yaml` before running.

**Configuration (`ChIPSeq/config.yaml`):**
```yaml
threads: 8
directories:
  fastq: data/fastq
  aligned: results/aligned
  peaks: results/peaks
reference:
  bowtie2_index: /path/to/bowtie2/index/prefix  # MUST update
macs2:
  genome: hs       # hs=human, mm=mouse
  qvalue: 0.01
samples:           # Add SRA accession IDs
  SRR1234567:
    sra: SRR1234567
pairs:             # Define treatment/control relationships
  sample1:
    treatment: SRR1234567
    control: SRR1234568
```

**Running:**
```bash
cd ChIPSeq/
snakemake --cores 8
# Dry run first:
snakemake -n
```

**Pipeline steps:** `download_fastq → fastqc → align_bowtie2 → mark_duplicates → index_bam → macs2_callpeak`

**Output:** `.narrowPeak` files in `results/peaks/`

### WES (Exome) Visualization

Standalone R scripts for visualizing exome sequencing results. No pipeline orchestration — run individually after upstream variant analysis.

```bash
Rscript Exome/R/plotAnnotation.R annotation_counts.tsv 0
Rscript Exome/R/plot_chrom.r coverage_A.bed coverage_B.bed output.png
Rscript Exome/R/plot_chrom_diff.r coverage_A.bed coverage_B.bed output.png
```

---

## Development Conventions

### Languages and Style

**Bash scripts:**
- Shebang: `#!/bin/bash`
- Variables in `snake_case`
- Use `mkdir -p` for output directory creation
- Use `$(command)` syntax for command substitution
- Hard-coded reference paths — these must be updated per computing environment

**R scripts:**
- Use `ggplot2` for all visualizations
- Use `dplyr` for data manipulation
- Parse command-line args with `optparse` (complex) or `commandArgs(trailingOnly=TRUE)` (simple)
- Save plots with `ggsave()`, typically PNG at 1000 DPI, 12×8 inches
- Mix of CamelCase and snake_case naming — not strictly enforced

**Python scripts:**
- Python 3
- Use `pandas` for tabular data
- Validate command-line arguments manually via `sys.argv`
- Tab-separated files (`sep="\t"`) used throughout

**Snakemake:**
- YAML configuration loaded at top of Snakefile
- Use `expand()` for multi-sample output targets
- Use `lambda wildcards:` for dynamic parameter values
- Specify `threads:` in all compute-intensive rules

### Naming Conventions

| Context | Convention | Example |
|---------|-----------|---------|
| Bash variables | `snake_case` | `output_dir`, `n_threads` |
| R variables | Mixed (CamelCase or snake_case) | `MAplot`, `count_table` |
| Python variables | `snake_case` | `output_file`, `gene_col` |
| Script files | Numbered prefix for sequential steps | `03_star.sh` |
| Output directories | Numbered to match script | `outputs/03_star/` |

### File Formats

- Raw sequencing: `.fastq.gz` (paired-end: `_R1.fastq.gz`, `_R2.fastq.gz`)
- Alignments: `.bam` (coordinate-sorted)
- Count matrices: tab-separated `.txt` files
- DEG results: tab-separated with columns including `Gene`, `log2FoldChange`, `padj`
- Visualization output: `.png` (high-resolution, 1000 DPI)
- Configuration: `.yaml`

---

## Dependencies

There is no package manager configuration (no `requirements.txt`, `environment.yml`, or `Dockerfile`). Dependencies must be installed manually or via conda.

### Required External Tools

| Tool | Used in | Purpose |
|------|---------|---------|
| FastQ Screen | RNA-Seq | Contamination QC |
| FastQC | RNA-Seq, ChIP-Seq | Quality assessment |
| Trim Galore | RNA-Seq | Adapter trimming |
| STAR | RNA-Seq | Read alignment |
| featureCounts | RNA-Seq | Read counting |
| fasterq-dump | ChIP-Seq | SRA data download |
| Bowtie2 | ChIP-Seq | Read alignment |
| SAMtools | ChIP-Seq | BAM manipulation |
| MACS2 | ChIP-Seq | Peak calling |
| Snakemake | ChIP-Seq | Workflow orchestration |

### Required R Packages

| Package | Used for |
|---------|---------|
| ggplot2 | All visualization |
| dplyr | Data manipulation |
| ggrepel | Non-overlapping labels in plots |
| ggpmisc | Regression equation display |
| optparse | CLI argument parsing |
| DESeq2 | Differential expression |
| SARTools | DESeq2 wrapper with QC |
| biomaRt | Gene ID annotation |
| IdeoViz | Chromosome ideogram plots |
| RColorBrewer | Color palettes |

### Required Python Packages

| Package | Used for |
|---------|---------|
| pandas | Tabular data manipulation |

**Recommended setup:**
```bash
conda create -n bioinformatics -c bioconda -c conda-forge \
  snakemake star bowtie2 samtools macs2 fastqc trim-galore \
  fastq-screen subread sra-tools r-base r-ggplot2 r-dplyr \
  bioconductor-deseq2 bioconductor-biomart bioconductor-ideoviz \
  r-optparse r-ggrepel r-ggpmisc pandas
conda activate bioinformatics
```

---

## Environment-Specific Configuration

The following values are **hardcoded** in scripts and must be updated for each computing environment:

| Script | Hardcoded Value | What to Change |
|--------|----------------|----------------|
| `RNA/bash_pipeline/03_star.sh` | hg38 genome index path | Path to local STAR index |
| `RNA/bash_pipeline/03_star.sh` | GTF annotation path | Path to local GTF file |
| `RNA/bash_pipeline/04_featurecounts.sh` | GTF annotation path | Path to local GTF file |
| `ChIPSeq/config.yaml` | `bowtie2_index` | Path to local Bowtie2 index |
| `ChIPSeq/config.yaml` | Sample SRA IDs | Actual experiment accessions |

---

## Testing

There is no formal test suite. The `test/` directory contains demo scripts with synthetic data.

```bash
# Generate synthetic test data (100 genes, random SNPCount and FPKM)
Rscript test/codes/create_dummy.R
# Output: test/data/genes_data.csv

# Generate scatter plots with regression (quantile filtering)
Rscript test/codes/dotplot.R test/data/genes_data.csv test/outputs/output
# Output: test/outputs/output_2.png, output_3.png, output_4.png
```

Note: `test/data/` and `test/outputs/` are gitignored.

---

## Git Workflow

- Development follows a **pull request model** (PRs #1, #3 merged into `master`)
- Branch naming observed: `codex/<feature-description>` (e.g., `codex/create-snakemake-workflow-for-chip-seq`)
- The default branch is `master`
- Remote: `http://local_proxy@127.0.0.1:60622/git/CfierroR/BioinfoGenomicRepo`

**When making changes:**
1. Create a feature branch from `master`
2. Make focused, well-described commits
3. Open a pull request for review

---

## Key Things to Know for AI Assistance

1. **No automated test runner exists** — validate changes manually by checking script logic and expected outputs.

2. **Reference genome paths are hardcoded** — when suggesting edits to pipeline scripts, note that paths like `/path/to/hg38/...` are placeholders that must be replaced with real local paths.

3. **Scripts are meant to be run sequentially** — the numbered prefix on RNA-Seq bash scripts (00_, 01_, ...) indicates execution order. Do not reorder them.

4. **ChIP-Seq is the most automated** — it uses Snakemake and `config.yaml`. Prefer editing the config over modifying the Snakefile directly.

5. **No linting tools are configured** — there are no pre-commit hooks, no `.editorconfig`, and no style enforcement. Follow existing conventions in each language.

6. **Mixed comment languages** — Some scripts contain comments in Spanish, others in English. This is intentional (researcher preference). Do not change comment language unless asked.

7. **Visualization scripts are standalone** — R visualization scripts (`MAPlot.R`, `VolcanoPlot.R`, `plotAnnotation.R`, etc.) are independent tools, not part of any automated pipeline. They accept file paths as arguments.

8. **Dependency versions are unspecified** — if adding new functionality, prefer packages already in use (ggplot2, dplyr, pandas, etc.) to minimize new dependencies.
