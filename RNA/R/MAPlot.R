# Cargar las librerías necesarias
library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

# Usage: Rscript MAPlot.R DEG_analysis.txt output.png

input_file = args[1]
output = args[2]

# Leer los datos desde archivo
data <- read.delim(input_file, header = TRUE, sep = "\t")

# Crear variable de significancia basada en el log fold change
data <- data %>%
  mutate(
    sig = case_when(
	  log2FoldChange <= -.5 ~ "Yes",
      log2FoldChange >=  .5 ~ "Yes",
      TRUE        ~ "No"
    )
  )

# Definir colores personalizados para el gráfico
custom_colors <- c("Yes" = "#FE4A49", "No" = "darkgrey")

FC_axis<-max(c(trunc(abs(max(data$log2FoldChange)))+1,trunc(abs(min(data$log2FoldChange))+1)))
baseMean_axis<-trunc(log10(max(data$baseMean,na.rm=TRUE)))+1

max(data$baseMean,na.rm=TRUE)
log10(max(data$baseMean,na.rm=TRUE))
trunc(log10(max(data$baseMean,na.rm=TRUE)))

# Crear el gráfico
plot <- ggplot(data) +
  geom_point(aes(x = log10(baseMean), y = log2FoldChange, color = sig), size = 1) +
  scale_colour_manual(values = custom_colors) +
  geom_hline(aes(yintercept = 0)) +
  theme_bw() +
  labs(
    x = "log baseMean",
    y = "log2(Fold Change)"
  ) +
  scale_x_continuous(limits = c(-1,baseMean_axis), breaks = seq(0, baseMean_axis, 1)) +
  scale_y_continuous(limits = c(-FC_axis, FC_axis), breaks = seq(-FC_axis, FC_axis, 1)) +
  guides(color = "none")

# Guardar el gráfico
ggsave(plot, filename = output, height = 8, width = 12, dpi = 1000)

