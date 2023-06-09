suppressMessages(library(polyester))
suppressMessages(library(Biostrings))
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))

option_list = list(
  make_option(c("-t", "--transcriptome"), type="character", default=NULL,
              help="Fasta of trnscriptome cDNA", metavar="character"),
  make_option(c("-n", "--numASE"), type="integer", default=100,
              help="Number of ASE Genes", metavar="character"),
  make_option(c("-f", "--fold"), type="integer", default=4,
              help="ASE fold change", metavar="character"),
  make_option(c("-e", "--expressed"), type="double", default=0.25,
              help="Proportion of transcript expressed", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default = "mBASED",
              help="Output directory name", metavar="character"),
  make_option(c("-d", "--depth"), type="integer", default = 10,
              help="Read Depth", metavar="character"),
  make_option(c("-c", "--convert"), type="character", default = NULL,
              help="Gene to Transcript TSV", metavar="character"),
  make_option(c("-s", "--snp"), type="character", default = NULL,
              help="Bed file of genes with SNP", metavar="character"),
  make_option(c("-u", "--tumorContent"), type="double", default = 1.0,
              help="Tumor Content", metavar="character"),
  make_option(c("-l", "--readLength"), type="integer", default = 150,
	      help="Read Length", metavar="character")
)

# load in options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$outdir
transcriptome <- opt$transcriptome
numASE <- opt$numASE
foldChange <- opt$fold
percentExpressed <- opt$expressed
depth <- opt$depth
convert <- opt$convert
snp_path <- opt$snp
tumorContent <- opt$tumorContent 
readLength <- opt$readLength

# Read in transcriptome
fasta <- readDNAStringSet(transcriptome)
snp <- read_delim(snp_path, col_names = F, delim="\t") %>%
  pull(X4)

# Generate baseline expression matrix (based on length)
fold_changes1 = matrix(1, nrow=length(fasta))
fold_changes2 = matrix(1, nrow=length(fasta))

readspertx <- round(depth * tumorContent * width(fasta) / readLength)

# Get transcript to gene table
transcriptGene <- tibble(original = names(fasta)) %>%
  separate(original, sep = " ", into=c("transcript", "cds"), fill = "right", remove = F) %>%
  dplyr::select(original, transcript)

convertList <- read_delim(convert, delim = "\t") %>%
  dplyr::filter(transcript %in% transcriptGene$transcript) %>%
  left_join(transcriptGene) %>%
  dplyr::filter((hgnc_symbol %in% snp))

# Randomly select genes to express and ASE
genelst <- pull(convertList, hgnc_symbol) %>% unique()
random_gene <- sample(genelst, length(genelst)*(percentExpressed))

ASE1_gene <- random_gene[1:(numASE %/% 2)]
ASE2_gene <- random_gene[(numASE %/% 2 + 1) : numASE]
Exp_gene <- random_gene[(numASE + 1):length(random_gene)]

# Convert randomly selected gene to transcript 
convert <- function(database, aseGene) {
  database %>%
    dplyr::filter(hgnc_symbol %in% aseGene) %>%
    pull(original)
}

ASE1_transcript <- convert(convertList, ASE1_gene)
ASE2_transcript <- convert(convertList, ASE2_gene)
Exp_transcript <- convert(convertList, Exp_gene)


ASE1_ID <- match(ASE1_transcript,fasta@ranges@NAMES)
ASE2_ID <- match(ASE2_transcript,fasta@ranges@NAMES)
Exp_ID <- match(Exp_transcript,fasta@ranges@NAMES)

# Change expression matrix for allele 1 and 2
fold_changes1[Exp_ID] = readspertx[Exp_ID] 
fold_changes1[ASE1_ID] = readspertx[ASE1_ID] * foldChange
fold_changes1[ASE2_ID] = readspertx[ASE2_ID]

fold_changes2[Exp_ID] = readspertx[Exp_ID] 
fold_changes2[ASE1_ID] = readspertx[ASE1_ID]
fold_changes2[ASE2_ID] = readspertx[ASE2_ID] * foldChange

# Save ASE list and expression matrix 
expressedGene <- convertList %>%
  dplyr::filter(hgnc_symbol %in% random_gene) %>%
  dplyr::mutate(ASE = hgnc_symbol %in% c(ASE1_gene, ASE2_gene)) %>%
  dplyr::select(-original)


saveRDS(c(fold_changes1), paste0(out, "/allele1/alt_expressionProfile.rds"))
saveRDS(c(fold_changes2), paste0(out, "/allele2/alt_expressionProfile.rds"))
write.table(expressedGene, file = paste0(out, "/expressedGene.tsv"), row.names = F, col.names = T, quote = F, sep="\t")

#### Ref genome
readspertx_ref <- round(depth * (1 - tumorContent) * width(fasta) / readLength)

allExp_transcript <- convert(convertList, Exp_gene)
allExp_ID <- match(allExp_transcript, fasta@ranges@NAMES) 

allExp <- matrix(1, nrow=length(fasta))
allExp[allExp_ID] <- readspertx_ref[allExp_ID]


saveRDS(c(allExp), paste0(out, "/allele1/ref_expressionProfile.rds"))
saveRDS(c(allExp), paste0(out, "/allele2/ref_expressionProfile.rds"))
