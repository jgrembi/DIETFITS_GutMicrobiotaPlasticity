## DADA2 pipline

## Use this on the server
# { time Rscript --no-save --no-restore --verbose dada2_pipeline.R > dietDADA2.Rout 2> errorDietDADA2.Rout ; } 2> timeDADA2.txt

## Load Package
.packages <- c("dada2", "doParallel", "foreach", "ggplot2", "dplyr", "phyloseq")
sapply(.packages, require, character.only = TRUE)

## SETTINGS ---------------------------
# DADA2 Paramers

## Filter and trim parameters
maxEE <- c(2, 2)
trimLeft <- c(10, 10) 
learnErrorsNreads <- 2e+06
ncharTabTrim <- seq(230, 260)

pool.list <- c("G1", "G2", "G3", "G4", "G5", "G9", "G10")
truncLen.lst <- list("G1" = c(240, 200), 
                     "G2" = c(240, 200), 
                     "G3" = c(200, 160),
                     "G4" = c(240, 200),
                     "G5" = c(240, 200),
                     "G9" = c(240, 200),
                     "G10" = c(230, 175))

## Directory to reference fasta files
path.to.reference <- "path/to/reference"

## Directory for split libraries 
main.dir <- path <- "path/to/processed_data"
## Directory to fastq raw files split into folders corresponding to batches 
# of sequencing run
path.to.split.lib <- file.path(main.dir, "split_libraries")

## Directory for trimmed files
dir.create(file.path(main.dir, "filtered_libraries"), showWarnings = FALSE)
filt.name <- paste0("maxEE_", paste0(maxEE, collapse = "_"))
path.to.filt <- file.path(main.dir, "filtered_libraries", filt.name)
dir.create(path.to.filt, showWarnings = FALSE)

## Directory to dada2 results
dir.create(file.path(main.dir, "dada2_results"), showWarnings = FALSE)
path.to.res <- file.path(main.dir, "dada2_results", filt.name)
dir.create(path.to.res, showWarnings = FALSE)

paramaters <- c("DADA2 parameters different from the default:", 
                paste0("[filterAndTrim()] maxEE: c(", maxEE[1], ",", maxEE[2], ")"),
                paste0("[filterAndTrim()] truncLengths: ", 
                       paste0(names(truncLen.lst)," = ", 
                              truncLen.lst, collapse = ", ")),
                paste0("[filterAndTrim()] trimLeft: ", trimLeft[1]),
                paste0("[learnErrors()] nreads = ", learnErrorsNreads), 
                "[learnErrors()] randomize = TRUE",
                paste0("[makeSequenceTab()] trimming seq by length: ", 
                       paste0("c(", ncharTabTrim[1], ", ", 
                              ncharTabTrim[length(ncharTabTrim)], ")"))
                
)
writeLines(paramaters, con = file.path(path.to.res, "dada2_params_used.txt"))

## RUN PIPELINE ---------------------------
registerDoParallel(length(pool.list))
foreach(i = seq_along(pool.list)) %dopar% {
  print("## Setting up file paths  ---------------------------")
  pool <- pool.list[i]
  path.to.filt.pool <- file.path(path.to.filt, pool)
  dir.create(path.to.filt.pool, showWarnings = FALSE)
  
  # Sort ensures forward/reverse reads samples are in same order
  path.to.pool.data <- file.path(path.to.split.lib, pool)
  fwdFQ <- sort(list.files(path.to.pool.data, pattern = "_R1_001"))
  revFQ <- sort(list.files(path.to.pool.data, pattern = "_R2_001"))
  # Extract sample names
  sample.names <- gsub(".fastq.gz", "", gsub(".*\\_", "", fwdFQ))
  sample.names[grepl(pool, sample.names)] <-  
    paste0("H2O_", sample.names[grepl(pool, sample.names)])
  # Specify the full path to the fwdFQ and revFQ
  fwdFQ <- file.path(path.to.pool.data, fwdFQ)
  revFQ <- file.path(path.to.pool.data, revFQ)
  
  # Names for the new filtered files
  filtFs <- file.path(path.to.filt.pool, paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(path.to.filt.pool, paste0(sample.names, "_R_filt.fastq.gz"))
  
  print("## Filtering and trimming ---------------------------")
  # CONSIDER increasing maxEE if too few reads are found !!!
  out <- filterAndTrim(fwdFQ, filtFs, revFQ, filtRs, 
                       truncLen = truncLen.lst[[pool]],
                       rm.phix = TRUE,
                       trimLeft = trimLeft,   # trim left end since q usually lower
                       maxEE = maxEE,   # THIS IS DIFFERENT THAN THE DEFAULT
                       matchIDs = TRUE, # Match ids between fwd and rev
                       multithread = TRUE) # Run in parallel
  
  filtered <- data.frame(row.names = sample.names,
                         sample.id = paste0(pool, "_", sample.names), 
                         file.name = rownames(out), out)
  # Save statistics
  write.csv(filtered, file.path(path.to.filt.pool, paste0(pool, "_track_filter.csv")))
}

registerDoParallel(length(pool.list))
foreach(i = seq_along(pool.list)) %dopar% {
  print("## Setting up file paths  ---------------------------")
  pool <- pool.list[i]
  path.to.filt.pool <- file.path(path.to.filt, pool)
  dir.create(path.to.filt.pool, showWarnings = FALSE)
    
  # Use only the filtered samples which were created for both
  # fwd and rev, as some of the samples might be discarded after 
  # filtering & trimming step due to low quality, and low read number
  filtFs <- list.files(path.to.filt.pool, pattern = "_F_filt.fastq.gz")
  filtRs <- list.files(path.to.filt.pool, pattern = "_R_filt.fastq.gz")
  filtFs.smp.name <- gsub("_F_filt.fastq.gz", "", filtFs)
  filtRs.smp.name <- gsub("_R_filt.fastq.gz", "", filtRs)
  non.zero.reads.smp <- as.character(intersect(filtFs.smp.name, filtRs.smp.name))
  filtFs <- file.path(path.to.filt.pool, paste0(non.zero.reads.smp, "_F_filt.fastq.gz"))
  filtRs <- file.path(path.to.filt.pool, paste0(non.zero.reads.smp, "_R_filt.fastq.gz"))
  
  print("## Learning error rates ---------------------------")
  # THIS PICKS only samples with numeric names and not the ones from baboon, chlamydia
  # etc from other studies. we also incorporate negative control H2O samples.
  keep.samples <- c(non.zero.reads.smp[!is.na(as.numeric(non.zero.reads.smp))],
                    non.zero.reads.smp[grepl("H2O", non.zero.reads.smp)])
  err.filtFs <- file.path(path.to.filt.pool, paste0(keep.samples, "_F_filt.fastq.gz"))
  err.filtRs <- file.path(path.to.filt.pool, paste0(keep.samples, "_R_filt.fastq.gz"))
  errF <- learnErrors(err.filtFs, nreads = 2e+06, randomize = TRUE, multithread = TRUE)
  errR <- learnErrors(err.filtRs, nreads = 2e+06, randomize = TRUE, multithread = TRUE)
  save(list = c("errF", "errR"), 
       file = file.path(path.to.filt.pool, paste0(pool, "_err.rda")))
  # We also save plots
  pdf(file = file.path(path.to.filt.pool, paste0(pool, "_err_rates.pdf")))
  print(plotErrors(errF, nominalQ = TRUE) + ggtitle(paste0(pool, ": forward Reads")))
  print(plotErrors(errR, nominalQ = TRUE) + ggtitle(paste0(pool, ": reverse Reads")))
  dev.off()
  
  print("## Dereplication ---------------------------")
  # Now dereplicate the samples that have filtered files for both fwd and rev
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  # Name the derep-class objects by the sample names
  names(derepFs) <- non.zero.reads.smp
  names(derepRs) <- non.zero.reads.smp
  
  print("## Sample Inference ---------------------------")
  # CONSIDER POOLING FOR BETTER RESULTS AT COST OF COMP TIME "pool = TRUE"
  dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
  dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
  save(list = c("dadaFs", "dadaRs"), 
       file = file.path(path.to.filt.pool, paste0(pool, "_dada.rda")))
  
  print("## Merge Reads ---------------------------")
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
  names(mergers) <- non.zero.reads.smp
  saveRDS(mergers, file.path(path.to.filt.pool, paste0(pool, "_mergers.rds")))
  
  print("## Track Reads ---------------------------")
  # Save the reads remaining after each step 
  getN <- function(x) sum(getUniques(x))
  filtered.file <- file.path(path.to.filt.pool, paste0(pool, "_track_filter.csv"))
  filtered <- read.csv(file = filtered.file, row.names = 1)
  filtered <- filtered[as.character(names(mergers)), ]
  track <- cbind(filtered, sapply(dadaFs, getN), sapply(mergers, getN))
  colnames(track) <- c("sample.id", "file.name", "input", "filtered", 
                       "denoised", "merged")
  rownames(track) <- paste0(pool, "_", rownames(track))
  write.csv(track, file.path(path.to.filt.pool, paste0(pool, "_track.csv")))
}

print("## Combine mergers ---------------------------")
track.list <- mergers.list <- list()
for (i in seq_along(pool.list)) {
  pool <- pool.list[i]
  path.to.filt.pool <- file.path(path.to.filt, pool)
  merged.reads <- file.path(path.to.filt.pool, paste0(pool, "_mergers.rds"))
  p.mergers <- readRDS(merged.reads)
  track.list[[pool]] <- 
    read.csv(file.path(path.to.filt.pool, paste0(pool, "_track.csv")),
             row.names = 1)
  names(p.mergers) <- paste0(pool, "_", names(p.mergers))
  mergers.list <- c(mergers.list, p.mergers)
}

print("## Save reads tracking data ---------------------------")
track.df <- plyr::ldply(track.list, .id = "pool")
rownames(track.df) <- track.df$sample.id
track.file <- file.path(path.to.res, "track_all_steps.csv")
write.csv(track.df, file = track.file)

print("## Save merger list ---------------------------")
mergers.file <- file.path(path.to.res, "all_mergers.rds")
saveRDS(mergers.list, file = mergers.file)
sample.id <- keep.samples <- names(mergers.list)
## Here you can remove some samples if you want
# keep.samples <- c(sample.id[!is.na(as.numeric(sample.id))],
#                   sample.id[grepl("H2O", sample.id)])
mergers.list <- mergers.list[keep.samples]

print("## Make table ---------------------------")
seqtab <- makeSequenceTable(mergers.list)
cat("Count table dimensions:", dim(seqtab), "\n")
seqtab.file <- file.path(path.to.res, "seqtab_init_transpose.csv")
write.csv(t(seqtab), seqtab.file)

# Inspect distribution of sequence lengths
cat("Distribution of sequence lengths: \n")
table(nchar(getSequences(seqtab)))
seqtab.rightlen <- seqtab[, nchar(colnames(seqtab)) %in% ncharTabTrim]
cat("Count table rightlen dim:", dim(seqtab.rightlen), "\n")
seqtab.file <- file.path(path.to.res, "seqtab_rightlen_transpose.csv")
write.csv(t(seqtab.rightlen), seqtab.file)

print("## Remove Chimeras ---------------------------")
seqtab.nochim <- removeBimeraDenovo(seqtab.rightlen, method="consensus", 
                                    multithread = TRUE, verbose = TRUE)
cat("Count table no chimera dim:", dim(seqtab.nochim), "\n")
# Fraction left with w/out chimera
cat("Fraction without chimera:", sum(seqtab.nochim)/sum(seqtab.rightlen), "\n")
seqtab.file <- file.path(path.to.res, "seqtab_nochim_transpose.csv")
write.csv(t(seqtab.nochim), file = seqtab.file)

print("## Track reads after seqtab ---------------------------")
track <- cbind(track.df[rownames(seqtab.rightlen), ], 
               rowSums(seqtab.rightlen), rowSums(seqtab.nochim))
colnames(track) <- c("pool", "sample.id", "file.name", "input", "filtered", 
                     "denoised", "merged", "tabled", "nonchim")
rownames(track) <- rownames(seqtab.rightlen)
track.file <- file.path(path.to.res, "track_nochim_reads.csv")
write.csv(track, file = track.file)

print("## Assign Taxonomy ---------------------------")
# Use RDP and Silva reference databases
refTaxa <- c(file.path(path.to.reference, "rdp_train_set_16.fa.gz"),
             file.path(path.to.reference, "silva_nr_v128_train_set.fa.gz"))

refSpecies <- c(file.path(path.to.reference, "rdp_species_assignment_16.fa.gz"),
                file.path(path.to.reference, "silva_species_assignment_v128.fa.gz"))
names(refTaxa) <- names(refSpecies) <- c("RDP", "Silva")

registerDoParallel(length(refTaxa))
foreach(ref = names(refTaxa)) %dopar% {
  print(paste0("## Running for ref:", ref))
  taxa <- assignTaxonomy(seqs = colnames(seqtab.nochim), verbose = TRUE,
                         refFasta = refTaxa[[ref]], multithread = TRUE)
  taxa <- addSpecies(taxa, refSpecies[[ref]], verbose=TRUE)
  saveRDS(taxa, file = file.path(path.to.res, paste0("taxonomy_", ref, ".rds")))
  ps <- phyloseq(otu_table(t(seqtab.nochim), taxa_are_rows = TRUE), 
                 tax_table(taxa))
  saveRDS(ps, file = file.path(path.to.res, paste0("phyloseq_", ref, ".rds")))
  
  multSpecies.file <- file.path(path.to.res, paste0("dietSpecies_", ref, ".rds"))
  species <- unname(
    assignSpecies(colnames(seqtab.nochim), refSpecies[[ref]], 
                  allowMultiple=TRUE)
  )
  saveRDS(species, file = multSpecies.file)
}