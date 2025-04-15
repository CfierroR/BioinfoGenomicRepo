require(IdeoViz)
require(RColorBrewer)
args = commandArgs(trailingOnly=TRUE)

# Usage: Rscript plot_chrom.r ConditionA.bed ConditionB.bed output.png

fileA = args[1]
fileB = args[2]
output = args[3]

# fileA = Coverage in condition A
# fileB = Coverage in condition B

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
data = GRanges(A$V1, IRanges(start = A$V2, end = A$V3) )
#Condition A left with neg values, condition B right with pos values
mcols(data)$group1 = scale(B$V4) 
mcols(data)$group2 = scale(A$V4) * -1

png(output, width = 4000, height = 1200)
par(mar = c(0, 0, 0, 0), mfrow = c(1,1))
plotOnIdeo(chrom = seqlevels(data),
           ideoTable = ideo,
           values_GR = data,
           value_cols = colnames(mcols(data)),
	   #Define colors
           col = c('orange', 'blue'),
           addScale = F,
	   #Define ranges
           val_range=c(-100,100),
           plotType='lines',
           cex.axis = 0.8,
           chromName_cex = 0.6,
           vertical = T)

dev.off()
