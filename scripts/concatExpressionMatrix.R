#!/gsc/software/linux-x86_64-centos7/R-4.1.3/bin/Rscript --vanilla
.libPaths("/home/glchang/R/x86_64-pc-linux-gnu-library/4.1")

suppressMessages(library(tidyverse))

basename.matches <- list.files(pattern='expression_matrix\\.tsv', recursive=TRUE,
                               full.names = TRUE)

input <- function(path) {
  file <- read_delim(path, show_col_types = FALSE) 
  file <- file[complete.cases(file), ] 
  file %>%
    distinct(gene, .keep_all = TRUE)
}

currPath <- basename.matches[1]
final <- input(currPath)


for (sample in basename.matches[-1]) {
  curr <- input(sample)
  final <- left_join(final, curr, by = "gene")
}

write.table(final, "final.expression.matrix.tsv", quote = F, row.names = F, sep = "\t")

