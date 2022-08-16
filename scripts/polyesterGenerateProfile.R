#!/gsc/software/linux-x86_64-centos7/R-4.1.3/bin/Rscript --vanilla
.libPaths("/home/glchang/R/x86_64-pc-linux-gnu-library/4.1")

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
              help="Read Depth", metavar="character")
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

fasta <- readDNAStringSet(transcriptome)

fold_changes1 = matrix(1, nrow=length(fasta))
fold_changes2 = matrix(1, nrow=length(fasta))

readspertx <- round(depth * width(fasta) / 100)

random_gene <- sample(fasta, length(fasta)*(percentExpressed))

random <- match(random_gene,fasta)
ASE1 <- random[1:(numASE %/% 2)]
ASE2 <- random[(numASE %/% 2 + 1) : numASE]
Exp <- random[(numASE + 1):length(random)]

fold_changes1[Exp] = readspertx[Exp] 
fold_changes1[ASE1] = readspertx[ASE1] * foldChange
fold_changes1[ASE2] = readspertx[ASE2]

fold_changes2[Exp] = readspertx[Exp] 
fold_changes2[ASE1] = readspertx[ASE1]
fold_changes2[ASE2] = readspertx[ASE2] * foldChange

ASEGenes <- tibble(expressedGene = random_gene[1:numASE]@ranges@NAMES) %>%
  separate(col = expressedGene, sep = ":", into = c("transcript", 'length'))

expressedGene <- random_gene[1:length(random)]@ranges@NAMES
expressGeneTibble <- tibble(expressedGene) %>%
  separate(col = expressedGene, sep = ":", into = c("transcript", 'length')) %>%
  dplyr::mutate(ASE = transcript %in% ASEGenes$transcript) 

saveRDS(c(fold_changes1), paste0(out, "/allele1/expressionProfile.rds"))
saveRDS(c(fold_changes2), paste0(out, "/allele2/expressionProfile.rds"))
write.table(expressGeneTibble, paste0(out, "/expressedGene.tsv"), row.names = F, col.names = T, quote = F)
