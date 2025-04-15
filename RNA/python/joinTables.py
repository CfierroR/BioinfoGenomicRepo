import pandas as pd
import sys
import os

#Usage: python joinTables.py archivo1.txt archivo2.txt salida_combinada.txt

def combine_gene_data(file1, file2, output_file):
    # Cargar los archivos en DataFrames con nombres genéricos
    df1 = pd.read_csv(file1, sep="\t", header=None, names=["Gene", "A_log2FoldChange", "A_pvalue"])
    df2 = pd.read_csv(file2, sep="\t", header=None, names=["Gene", "B_log2FoldChange", "B_pvalue"])

    # Fusionar ambos DataFrames en base a la columna 'Gene'
    merged_df = pd.merge(df1, df2, on="Gene", how="outer")

    # Guardar el resultado
    merged_df.to_csv(output_file, sep="\t", index=False)
    print(f"Archivo combinado guardado en: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python joinTables.py archivo_A archivo_B archivo_salida")
        sys.exit(1)

    archivo_A = sys.argv[1]
    archivo_B = sys.argv[2]
    archivo_salida = sys.argv[3]

    if not os.path.exists(archivo_A):
        print(f"Error: El archivo {archivo_A} no existe.")
        sys.exit(1)
    if not os.path.exists(archivo_B):
        print(f"Error: El archivo {archivo_B} no existe.")
        sys.exit(1)

    combine_gene_data(archivo_A, archivo_B, archivo_salida)

