suppressMessages(library(polyester))
suppressMessages(library(Biostrings))
suppressMessages(library(optparse))

# Make help options
option_list = list(
  make_option(c("-t", "--transcriptome"), type="character", default=NULL,
              help="Fasta of transcriptome cDNA", metavar="character"),
  make_option(c("-e", "--expressionProfile"), type="character", default=NULL,
              help="Expression profile of fold change for each transcript", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default = "mBASED",
              help="Output directory name", metavar="character"),
  make_option(c("-r", "--readLength"), type="integer", default = NULL,
   	      help="Read Length", metavar="character")
)

# load in options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$outdir
readLength <- opt$readLength

expressionProfile <- readRDS(opt$expressionProfile)
transcriptome <- opt$transcriptome

fasta <- readDNAStringSet(transcriptome)

# Edit polyester script
sim <- function (fasta = NULL, gtf = NULL, seqpath = NULL, outdir = ".", 
                 num_reps = c(10, 10), reads_per_transcript = 300, size = NULL, 
                 fold_changes, paired = TRUE, reportCoverage = FALSE, ...) 
{
  extras = list(...)
  extras = polyester:::.check_extras(extras, paired, total.n = sum(num_reps))
  if ("seed" %in% names(extras)) {
    set.seed(extras$seed)
  }
  if (!is.null(fasta) & is.null(gtf) & is.null(seqpath)) {
    transcripts = readDNAStringSet(fasta)
  }
  else if (is.null(fasta) & !is.null(gtf) & !is.null(seqpath)) {
    message("parsing gtf and sequences...")
    if ("exononly" %in% names(extras)) {
      exononly = extras$exononly
    }
    else {
      exononly = TRUE
    }
    if ("idfield" %in% names(extras)) {
      idfield = extras$idfield
    }
    else {
      idfield = "transcript_id"
    }
    if ("attrsep" %in% names(extras)) {
      attrsep = extras$attrsep
    }
    else {
      attrsep = "; "
    }
    transcripts = seq_gtf(gtf, seqpath, feature = "transcript", 
                          exononly = exononly, idfield = idfield, attrsep = attrsep)
    message("done parsing")
  }
  else {
    stop("must provide either fasta or both gtf and seqpath")
  }
  polyester:::.check_fold_changes(fold_changes, num_reps, transcripts)
  if ("meanmodel" %in% names(extras)) {
    if (extras$meanmodel) {
      b0 = -3.0158
      b1 = 0.8688
      sigma = 4.152
      logmus = b0 + b1 * log2(width(transcripts)) + rnorm(length(transcripts), 
                                                          0, sigma)
      reads_per_transcript = 2^logmus - 1
      reads_per_transcript = pmax(reads_per_transcript, 
                                  0.000001)
    }
  }
  if (length(num_reps) == 1) {
    fold_changes = matrix(rep(1, length(transcripts)))
    basemeans = reads_per_transcript * fold_changes
  }
  else if (length(num_reps) == 2) {
    if (length(reads_per_transcript) == 1) {
      basemeans = matrix(reads_per_transcript, ncol = 2, 
                         nrow = length(transcripts))
      basemeans[, 2] = fold_changes * basemeans[, 1]
    }
    else if (length(reads_per_transcript) == length(transcripts)) {
      basemeans = matrix(c(reads_per_transcript, reads_per_transcript), 
                         nrow = length(reads_per_transcript))
      basemeans = basemeans * fold_changes
    }
    else {
      stop("reads_per_transcript is the wrong length.")
    }
  }
  else {
    basemeans = reads_per_transcript * fold_changes
  }
  if (is.null(size)) {
    size = basemeans/3
  }
  else if (class(size) == "numeric") {
    if (is.matrix(basemeans)) {
      num_rows = nrow(basemeans)
      num_cols = ncol(basemeans)
    }
    else {
      num_rows = length(basemeans)
      num_cols = 1
    }
    size = matrix(size, nrow = num_rows, ncol = num_cols)
  }
  else if (class(size) == "matrix") {
    if (!is.matrix(basemeans)) {
      stop("If you provide a matrix for size, you also need a matrix for reads_per_transcript.")
    }
    stopifnot(nrow(size) == nrow(basemeans))
    stopifnot(ncol(size) == ncol(basemeans))
  }
  else {
    stop("size must be a number, numeric vector, or matrix.")
  }
  if ("seed" %in% names(extras)) {
    set.seed(extras$seed)
  }
  group_ids = rep(1:length(num_reps), times = num_reps)
  numreadsList = vector("list", sum(num_reps))
  numreadsList = lapply(1:sum(num_reps), function(i) {
    group_id = group_ids[i]
    NB2(as.matrix(basemeans)[, group_id], as.matrix(size)[, 
                                                          group_id])
  })
  readmat = matrix(unlist(numreadsList), ncol = sum(num_reps))
  readmat = t(extras$lib_sizes * t(readmat))
  if ("gcbias" %in% names(extras)) {
    stopifnot(length(extras$gcbias) == sum(num_reps))
    gcclasses = sapply(extras$gcbias, class)
    if (!all(gcclasses %in% c("numeric", "loess"))) {
      stop(.makepretty("gc bias parameters must be integers 0 through 7\n                or loess objects"))
    }
    if (any(extras$gcbias != 0)) {
      readmat = add_gc_bias(readmat, extras$gcbias, transcripts)
    }
  }
  sysoutdir = gsub(" ", "\\\\ ", outdir)
  if (.Platform$OS.type == "windows") {
    shell(paste("mkdir", sysoutdir))
  }
  else {
    system(paste("mkdir -p", sysoutdir))
  }
  
  polyester:::sgseq(readmat, transcripts, paired, outdir, extras, reportCoverage)
  if (!("write_info" %in% names(extras))) {
    write_info = TRUE
  }
  if (write_info) {
    counts_matrix <- readmat
    polyester:::.write_info(extras, transcripts, num_reps, fold_changes, 
                            outdir, group_ids, counts_matrix)
  }
}


NB2 <- function (basemeans, size, seed = NULL) 
{
  if (!is.null(seed)) 
    set.seed(seed)
  
  numreads = basemeans
  numreads[numreads == 0] = 1
  return(numreads)
}

# Start simulation
print("Start simulation")
sim(transcriptome, reads_per_transcript=expressionProfile, 
                    num_reps=c(1), readlen=readLength, fold_changes=matrix(1, nrow=length(fasta)),outdir=out) 
