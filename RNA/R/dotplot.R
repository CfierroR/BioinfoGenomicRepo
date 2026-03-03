#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpmisc)
  library(dplyr)
})

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript dotplot.R <input_file> <output_file>")
}

input_file <- args[1]
output_file <- args[2]

for (n in 2:4) { 
# Load and filter data
x <- read.delim(input_file, header = TRUE, sep = "\t")
data <- x %>%
  filter(FPKM > quantile(x$FPKM)[n]) %>%
  filter(SNPCount > quantile(x$SNPCount)[n])

# Create the scatter plot with regression line and equation
p <- ggplot(data, aes(x = SNPCount, y = FPKM)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  geom_point()  +
theme_bw()

# Save the plot
ggsave(filename = paste0(output_file,"_",n,".png"), plot = p, height = 8, width = 10)
}
