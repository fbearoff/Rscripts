library(biomaRt)
library(data.table)

mart <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart) # select species here

tx2gene <- getBM(
  attributes = c(
    "ensembl_transcript_id_version",
    "ensembl_gene_id_version",
    "external_gene_name",
    "description",
    "chromosome_name",
    "start_position",
    "end_position"
  ),
  mart = ensembl
)

tx2gene <- as.data.table(tx2gene)[chromosome_name %in% c(1:22, "MT", "X", "Y")]
setnames(
  tx2gene,
  old = c(
    "ensembl_transcript_id_version",
    "ensembl_gene_id_version",
    "external_gene_name",
    "chromosome_name",
    "start_position",
    "end_position"
  ),
  new = c(
    "tx_id",
    "gene_id",
    "gene_symbol",
    "chr",
    "start",
    "end"
  )
)
