# Número de genes
n <- 100

# Crear el dataframe
df <- data.frame(
  gene = paste0("gene", 1:n),
  SNPCount = round(rnorm(n, mean = 0, sd = 100)),
  FPKM = round(rnorm(n, mean = 0, sd = 500), 2)
)

# Mostrar el dataframe
head(df)

# Exportar a archivo CSV

write.table(df, file="data/genes_data.csv", row.names=FALSE,sep="\t",quote=F)

