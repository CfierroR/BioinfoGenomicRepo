#!/bin/bash
# pipeline_config.sh — Configuración central del pipeline RNA-Seq
# Source este archivo en cada script del pipeline:
#   source "$(dirname "$0")/pipeline_config.sh"

# ── Hilos ────────────────────────────────────────────────────────────────────
THREADS_HIGH=30   # FastQ Screen, FastQC, STAR
THREADS_MED=10    # featureCounts
THREADS_LOW=8     # Trim Galore

# ── Genoma de referencia (hg38) ───────────────────────────────────────────────
# Actualizar estas rutas según el entorno de cómputo local
REF_GENOME="/home/resources/genomes/Genomes/STAR_hg38"
GTF_FILE="/home/resources/genomes/Genomes/ensembl/hg38/Homo_sapiens.GRCh38.104.chr.gtf"
GENOME_VERSION="hg38"
