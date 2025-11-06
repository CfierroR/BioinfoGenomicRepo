# ChIP-seq Snakemake Pipeline

This pipeline automates a standard ChIP-seq analysis workflow, starting with data retrieval from the SRA and ending with peak calling using MACS2. Update the configuration file before running the workflow.

## Requirements

Ensure the following tools are installed and available in your `$PATH`:

- [SRA Toolkit](https://github.com/ncbi/sra-tools) providing `fasterq-dump`
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [SAMtools](http://www.htslib.org/)
- [MACS2](https://github.com/macs3-project/MACS)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)

## Configuration

Edit `config.yaml` to match your experiment:

- `threads`: number of CPU threads to allocate for multithreaded steps.
- `directories`: destination folders for intermediate and final outputs. They will be created automatically.
- `reference.bowtie2_index`: prefix to the Bowtie2 genome index (built with `bowtie2-build`).
- `macs2`: parameters passed to MACS2 peak calling.
- `samples`: dictionary keyed by internal sample IDs with their corresponding SRA accessions.
- `pairs`: mapping between ChIP treatment and matched control samples. The keys are used to name the final peak outputs.

## Running the pipeline

1. Change into the `ChIPSeq` directory.
2. Review and edit `config.yaml` to include the desired SRA accessions, reference index, and output paths.
3. Execute Snakemake, specifying a desired level of parallelism:

   ```bash
   snakemake --use-conda --cores 8
   ```

   Remove `--use-conda` if you prefer to manage dependencies globally.

The workflow performs the following steps for each sample:

1. **Download SRA reads** with `fasterq-dump` (paired-end reads are split and compressed).
2. **Quality control** using FastQC.
3. **Alignment** of reads to the reference genome with Bowtie2, producing sorted BAM files.
4. **Duplicate removal** with `samtools markdup` and BAM indexing.
5. **Peak calling** with MACS2 using the ChIP/Input pairs defined in the configuration.

The final peaks are located in `results/peaks/<pair>/<pair>_peaks.narrowPeak`.
