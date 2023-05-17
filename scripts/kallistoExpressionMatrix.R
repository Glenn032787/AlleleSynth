#!/gsc/software/linux-x86_64-centos7/R-4.0.2/bin/Rscript --vanilla

suppressMessages(library("tidyverse"))
suppressMessages(library("biomaRt"))
suppressMessages(library("tximport"))

option_list = list(
  make_option(c("-p", "--sample"), type="character", default=NULL,
              help="Sample Name", metavar="character"),
  make_option(c("-a", "--abundance"), type="character", default=NULL,
              help="Kallisto abundance file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default = "mBASED",
              help="Output directory name", metavar="character"),
  make_option(c("-e", "--ensembl2hgnc"), type="character", default = NULL,
              help="Ensembl to hgnc symbol file", metavar="character")
)

# load in options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$outdir
sample <- opt$sample

ensembl2hgncPath <- opt$ensembl2hgnc
abundancefile <- opt$abundance

ensembl2hgnc <- read_delim(ensembl2hgncPath)
txi.kallisto <- tximport(abundancefile,
                         type = "kallisto",
                         txOut = TRUE,
                         ignoreTxVersion = TRUE)

transcript_names <-rownames(txi.kallisto$abundance)
transcript_names_convert <- gsub("\\..*","",transcript_names)
convert <-  data.frame(transcript=transcript_names,
                       ensembl_transcript_id = transcript_names_convert) %>%
   as_tibble()


txtgene <- left_join(convert, ensembl2hgnc, by = c("ensembl_transcript_id" = "transcript"), 
                     relationship = "many-to-many") %>%
   dplyr::filter(!is.na(hgnc_symbol)) %>%
   dplyr::select(transcript, hgnc_symbol)

txi.kallisto.gene <- summarizeToGene(txi.kallisto, as.data.frame(txtgene))
expression_matrix <- txi.kallisto.gene$abundance %>%
   as.data.frame() %>%
   rownames_to_column("gene") %>%
   dplyr::filter(gene != '')
colnames(expression_matrix) <- c('gene', sample)

write.table(expression_matrix, file=paste0(out, "/expression_matrix.tsv"), quote=FALSE, sep='\t', row.name = F)
