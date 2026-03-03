library(IdeoViz)
library(RColorBrewer)
args = commandArgs(trailingOnly=TRUE)

# Usage:
#   Rscript plot_chrom.r ConditionA.bed ConditionB.bed output.png [mode]
# mode: "compare" (default) — dos condiciones superpuestas
#       "diff"              — diferencia B - A

fileA  = args[1]
fileB  = args[2]
output = args[3]
mode   = ifelse(length(args) >= 4, args[4], "compare")

A = read.table(fileA, as.is = T, sep="\t", header=F)
B = read.table(fileB, as.is = T, sep="\t", header=F)

# Set bin size to vertical barplot
windows_size = A$V3[1] - A$V2[1]

# Modify if conection to server fail. Need to work first or choose .rds in src folder
# Modify if genome version is wrong
ideo <- getIdeo("hg38")
saveRDS(ideo, "hg38_ideo.rds")
#ideo <- readRDS(args[4])

# Parse data
data = GRanges(A$V1, IRanges(start = A$V2, end = A$V3))

if (mode == "diff") {
  # Diferencia B - A en un solo canal
  mcols(data)$value = scale(B$V4 - A$V4)
  plot_col   = 'orange'
  val_range  = c(-5, 5)
  plot_type  = 'rect'
  plot_title = paste("Difference between",
                     gsub(".coverage", "", fileA),
                     gsub(".coverage", "", fileB),
                     "bin", windows_size)
} else {
  # Condición A izquierda (neg), condición B derecha (pos)
  mcols(data)$group1 = scale(B$V4)
  mcols(data)$group2 = scale(A$V4) * -1
  plot_col   = c('orange', 'blue')
  val_range  = c(-100, 100)
  plot_type  = 'lines'
  plot_title = NULL
}

png(output, width = 4000, height = 1200, units = "px")
par(mar = c(0, 0, 0, 0), mfrow = c(1,1))
plotOnIdeo(chrom         = seqlevels(data),
           ideoTable     = ideo,
           values_GR     = data,
           value_cols    = colnames(mcols(data)),
           col           = plot_col,
           addScale      = F,
           val_range     = val_range,
           plotType      = plot_type,
           plot_title    = plot_title,
           cex.axis      = 0.8,
           chromName_cex = 0.6,
           vertical      = T)

dev.off()
