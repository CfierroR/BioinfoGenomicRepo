#!/usr/bin/env Rscript
# Thin wrapper — delega en plot_chrom.r con mode="diff".
# Usage: Rscript plot_chrom_diff.r ConditionA.bed ConditionB.bed output.png

args        <- commandArgs(trailingOnly = TRUE)
script_file <- commandArgs(trailingOnly = FALSE)
script_dir  <- dirname(normalizePath(sub("^--file=", "", script_file[grepl("^--file=", script_file)])))

system2("Rscript", c(file.path(script_dir, "plot_chrom.r"), args, "diff"))
