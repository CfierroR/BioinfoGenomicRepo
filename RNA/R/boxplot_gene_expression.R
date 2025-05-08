# Argument receiver for the script

# Load necessary library
library(optparse)
library(ggplot2)
library(dplyr)
library(tidyr)

# Function to check if the gene column exists in the data and return the column name
checkGeneColumn <- function(data, gene) {
    # Check if the gene exists in any column
    for (i in 1:ncol(data)) {
        column_it <- colnames(data)[i]
        if (gene %in% data[[column_it]]) {
            # Extract the column with the gene name
            gene_data.col <- colnames(data)[i]
            return(gene_data.col)
        }
    }
    return(NULL)
}
#Function to check if normalized columns exist in the data
checkNormalizedColumns <- function(data) {
    norm_pos <- which(grepl("^norm\\.", colnames(data)))[1]
    return(norm_pos)
}

# Function to retrieve the group names from the data
getGroupNames <- function(data) {
    norm_pos <- checkNormalizedColumns(data)
    if(is.na(norm_pos)) {
        # Trim the column names to remove characters after the last underscore
        group_names <- sub("_(?!.*_).*", "", colnames(data), perl = TRUE)
        #Pop first column
        group_names <- group_names[-1]
        norm_pos = 0
    }
    else{
        
        norm_pos <- which(grepl("^norm\\.", colnames(data)))[1]
        group_names<-colnames(data)[2:(norm_pos - 1)]
        #Trim the column names to remove characters after the last underscore  
        group_names <- sub("_(?!.*_).*", "", group_names, perl = TRUE)
        
    }
    
    return(c(c(group_names),c(norm_pos)))
}

# Define the arguments
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                            help = "Path to the input file containing gene expression data", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                            help = "Path to save the output boxplot image [default: %default]", metavar = "character"),
    make_option(c("-g", "--gene"), type = "character", default = NULL,
                            help = "Gene name to plot the expression levels for", metavar = "character")
)

# Parse the arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if input file is provided
if (is.null(opt$input)) {
    stop("Input file is required. Use -i or --input to specify the file path.")
}

# Check if gene name is provided
if (is.null(opt$gene)) {
    stop("Gene name is required. Use -g or --gene to specify the gene.")
}
# Check if output file is provided
# If not provided, set a default output file name
if (is.null(opt$output)) {
    opt$output <- paste0(opt$gene,"_boxplot.png")  # Default output file name
}

# Read the input data
data <- read.csv(opt$input, header = TRUE, sep="\t", stringsAsFactors = FALSE)

gene_column <- checkGeneColumn(data, opt$gene)
# Check if the gene column was found
if(is.null(gene_column)){
    stop(paste("Gene", opt$gene, "not found in the input data."))
}

# Get group names in a vector
group_names <- getGroupNames(data)

# Get row with the gene name
gene_row <- data[data[[gene_column]] == opt$gene, ]

# Set data to long format for ggplot
if(as.numeric(tail(group_names, n=1)) == 0){
    print("No normalized columns found, using all columns")
    data_long <- gene_row %>%
    select(-gene_column) %>%
    pivot_longer(cols = everything(), names_to = "Sample", values_to = "Expression") %>%
    mutate(Group = group_names[1:length(group_names)-1])  # Set the order of samples
}else {
    print("Normalized columns found, using only unnormalized columns")
    data_long <- gene_row[2:(as.numeric(tail(group_names, n=1))-1)] %>%
    pivot_longer(cols = everything(), names_to = "Sample", values_to = "Expression") %>%
    mutate(Group = group_names[1:length(group_names)-1])  # Set the order of samples
}

print(data_long)
# Create the boxplot
boxplot <- ggplot(data_long, aes(x = Group, y = Expression)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(aes(color = Group), width = 0.2, alpha = 0.5) +
    labs(x = "Sample",
         y = "Expression Level") +
    theme_bw()
# Save the boxplot to a file
ggsave(filename = opt$output, plot = boxplot, width = 10, height = 6)
# Print a message indicating that the plot has been saved
cat("Boxplot saved to", opt$output, "\n")
