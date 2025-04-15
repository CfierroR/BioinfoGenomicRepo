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

# Asignar significancia basada en p-valor y log2FoldChange
data <- data %>%
  mutate(
    sig = case_when(
      pvalue <= 0.05 & log2FoldChange >= 1  ~ "Yes",
      pvalue <= 0.05 & log2FoldChange <= -1 ~ "Yes",
      TRUE                                  ~ "No"
    )
  )


FC_axis<-max(c(trunc(abs(max(data$log2FoldChange)))+1,trunc(abs(min(data$log2FoldChange))+1)))
padj_axis<--trunc(log10(min(data$padj,na.rm=TRUE)))+1


# Definir los colores para puntos significativos y no significativos
custom_colors <- c("Yes" = "#FE4A49", "No" = "darkgrey")

# Crear el gráfico tipo volcano
plot <- ggplot(data) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  scale_colour_manual(values = custom_colors) +
  theme_bw() +
  labs(
    x = "Fold Enrichment",
    y = "-log10(p-value)"
  ) +
  scale_x_continuous(limits=c(-FC_axis,FC_axis),breaks=seq(-FC_axis,FC_axis,1))+
  scale_y_continuous(limits=c(0,padj_axis),breaks=seq(0,padj_axis,1))+
  guides(color = "none")

# Guardar el gráfico
ggsave(plot, filename = output, width=12,height=8,dpi=1000)
