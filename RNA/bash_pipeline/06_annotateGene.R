library('biomaRt')

#ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
#mouse_gene_ids <- 

#df <- read.table("../analyses/PCA/SARTools/run2/tables", sep="\t", header=T, stringsAsFactors=F, check.names=F)


mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#genes <- df$genes

#G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name","external_synonym"),values=genes, mart= mart)

# Merge data set
#merge(df,G_list,by.x="gene",by.y="ensembl_peptide_id")

# list tables
#files = list.files("../analyses/PCA/SARTools/run3/tables", full.names=T)
files = list.files("/home/cletelier/RNA_Seq_Barbara_31_01_24/outputs/05_SarTools/RNA_seq_Neuron/tables",full.names=T)
for ( file in files ) {
    df <- read.table(file, sep="\t", header=T, stringsAsFactors=F, check.names=F)
    if (nrow(df) > 0 ){
        genes <- df$Id
        G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=genes, mart= mart)
        # Merge data set
        new_df <- merge(df, G_list, by.x="Id", by.y="ensembl_gene_id", all.x=T)
        write.table(new_df, file.path(dirname(file), gsub(".txt",".annot.txt", basename(file))), sep="\t", col.names=T, row.names=F, quote=F)
    }
}


