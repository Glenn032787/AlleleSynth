suppressMessages(library("tidyverse"))
suppressMessages(library("biomaRt"))
suppressMessages(library("optparse"))

option_list = list(
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="Sample name", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Gene expression data from RSEM", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default =NULL,
              help="Output directory name", metavar="character")
)

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$outdir
input <- opt$input
sample <- opt$sample


val <- read_delim(input, delim = "\t") %>%
  dplyr::select(gene_id, TPM)

mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="https://apr2020.archive.ensembl.org")
conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
             filters    = "ensembl_gene_id",
             values     = val$gene_id,
             mart       = mart)  %>%
	dplyr::filter(hgnc_symbol != "")

new_val <- left_join(conversion, val, by = c("ensembl_gene_id" = "gene_id")) %>%
  dplyr::select(hgnc_symbol, TPM) %>%
  `colnames<-`(c("gene", sample))

write.table(new_val, paste0(out, "/expression_matrix.tsv"), 
            col.names = T, row.names = F, quote = F, sep = "\t")
