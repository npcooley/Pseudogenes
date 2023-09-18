###### -- read assembly pipeline ----------------------------------------------

###### -- libraries -----------------------------------------------------------

suppressMessages(library(SynExtend))

###### -- arguments -----------------------------------------------------------

ARGS <- commandArgs(trailingOnly = TRUE)

# SRA run to grab and manipulate
RUN <- ARGS[1L]
# logical, run trimDNA or not
TRIMDNA <- as.logical(ARGS[2L])
# int specifying quality bin to run
# 1 == < 20
# 2 == >= 20 & < 30
# 3 == >= 30
# 4 == all
QBIN <- as.integer(ARGS[3L])
# numeric specifying how aggressively to downsample remaining reads
# 1.0x, 0.5x, 0.25x, 0.05x
DSAMPLE <- as.numeric(ARGS[4L])
# unique ID to save off
OUTID <- ARGS[5L]

# set a seed for reproducibility --- the grid makes this a little funny though
set.seed(as.integer(OUTID))

###### -- code body part 1 ----------------------------------------------------
# given reads specified in RUN
# parse down to the reads to pass to each assembler

SCRIPTSTART <- Sys.time()

PREFETCH <- paste0("prefetch -O sra ",
                   RUN)
# grab READS from 
system(command = PREFETCH,
       intern = FALSE,
       timeout = 7200)

FASTQDUMP <- paste0("fastq-dump",
                    " --outdir fastq",
                    " --gzip",
                    " --skip-technical",
                    " --readids",
                    " --read-filter pass",
                    " --dumpbase",
                    " --split-3",
                    " --clip",
                    " sra/", # container sra directory
                    RUN,
                    "/",
                    RUN, # target
                    ".sra")

# should not take 10 hours, but who knows
system(command = FASTQDUMP,
       intern = FALSE,
       timeout = 36000)

# are reads paired
FILES01 <- list.files(path = "fastq",
                      full.names = TRUE)
print(FILES01)
if (length(FILES01) < 2L) {
  print(paste0("Skipping ",
               RUN,
               " because reads appear single ended"))
  system(command = "rm fastq/*")
  
  q(save = "no")
}

###### -- code body part 2 ----------------------------------------------------
# always trim with fastp to remove adaptor sequences when present

# trim with fastp
# mostly for adaptor trimming
FASTP <- paste0("fastp",
                " --in1",
                " ",
                FILES01[1L],
                " --in2",
                " ",
                FILES01[2L],
                " --out1",
                " fastq/res1a.fastq.gz",
                " --out2",
                " fastq/res2a.fastq.gz",
                " --detect_adapter_for_pe")
system(command = FASTP)

seqs1a <- readQualityScaledDNAStringSet(filepath = "fastq/res1a.fastq.gz")
seqs2a <- readQualityScaledDNAStringSet(filepath = "fastq/res2a.fastq.gz")

###### -- code body part 3 ----------------------------------------------------
# factorial requests from the job map

if (TRIMDNA) {
  # trim with trimDNA
  seqs1a <- TrimDNA(myDNAStringSet = seqs1a,
                    leftPatterns = "",
                    rightPatterns = "",
                    type = "sequences",
                    quality = seqs1a@quality)
  seqs2a <- TrimDNA(myDNAStringSet = seqs2a,
                    leftPatterns = "",
                    rightPatterns = "",
                    type = "sequences",
                    quality = seqs2a@quality)
  
  partner1 <- do.call(rbind,
                       strsplit(x = names(seqs1a),
                                split = " ",
                                fixed = TRUE))
  partner1 <- as.integer(partner1[, 2L])
  partner2 <- do.call(rbind,
                       strsplit(x = names(seqs2a),
                                split = " ",
                                fixed = TRUE))
  partner2 <- as.integer(partner2[, 2L])
  
  seqs1a <- seqs1a[partner1 %in% partner2]
  seqs2a <- seqs2a[partner2 %in% partner1]
  
  if (length(seqs1a) == 0 |
      length(seqs2a) == 0) {
    print("No reads remaining after quality filtering.")
    q(save = "no")
  }
  
  rm(list = c("partner1", "partner2"))
  
  quals1a <- mean(as(seqs1a@quality, "IntegerList"))
  quals2a <- mean(as(seqs2a@quality, "IntegerList"))
} else {
  # do not trim with trimDNA
  quals1a <- mean(as(seqs1a@quality, "IntegerList"))
  quals2a <- mean(as(seqs2a@quality, "IntegerList"))
  
}

# both reads are in the same bin
if (QBIN == 1L) {
  # only reads with quality scores below 20
  w1 <- quals1a < 20 & quals2a < 20
} else if (QBIN == 2L) {
  # reads with quality scores between 20 and 30
  w1 <- quals1a >= 20 & quals1a < 30 & quals2a >= 20 & quals2a < 30
} else if (QBIN == 3L) {
  # reads with quality scores above 30
  w1 <- quals1a > 30 & quals2a > 30
} else if (QBIN == 4L) {
  # all reads regardless of quality scores
  w1 <- rep(TRUE, length(quals1a))
}

seqs1a <- seqs1a[w1]
seqs2a <- seqs2a[w1]

if (length(seqs1a) == 0 |
    length(seqs2a) == 0) {
  print("No reads remaining after quality subsetting.")
  q(save = "no")
}

# downsample
SampleSize <- ceiling(length(seqs1a) * DSAMPLE)

if (SampleSize != length(seqs1a)) {
  SelectSample <- sample(x = seq(length(seqs1a)),
                         size = SampleSize,
                         replace = FALSE)
  seqs1a <- seqs1a[SelectSample]
  seqs2a <- seqs2a[SelectSample]
  
  if (length(seqs1a) == 0 |
      length(seqs2a) == 0) {
    print("No reads remaining after subsampling.")
    q(save = "no")
  }
}

writeQualityScaledXStringSet(x = seqs1a,
                             filepath = "fastq/res1b.fastq.gz",
                             compress = TRUE)
writeQualityScaledXStringSet(x = seqs2a,
                             filepath = "fastq/res2b.fastq.gz",
                             compress = TRUE)

# some testing stuff, no longer necessary
# writeQualityScaledXStringSet(x = seqs1a,
#                              filepath = paste0("readsa",
#                                                formatC(x = as.integer(OUTID),
#                                                        width = 6L,
#                                                        flag = 0,
#                                                        format = "d"),
#                                                ".fastq.gz"),
#                              compress = TRUE)
# writeQualityScaledXStringSet(x = seqs2a,
#                              filepath = paste0("readsb",
#                                                formatC(x = as.integer(OUTID),
#                                                        width = 6L,
#                                                        flag = 0,
#                                                        format = "d"),
#                                                ".fastq.gz"),
#                              compress = TRUE)
# 
# save(seqs1a,
#      seqs2a,
#      file = paste0("reads",
#                    formatC(x = as.integer(OUTID),
#                            width = 6L,
#                            flag = 0,
#                            format = "d"),
#                    ".RData"),
#      compress = "xz")

###### -- code body part 4 ----------------------------------------------------
# the 4 available assemblers

FILES02 <- c("fastq/res1b.fastq.gz",
             "fastq/res2b.fastq.gz")

###### -- SKESA ---------------------------------------------------------------
RUNSKESA <- paste0("skesa",
                   " ",
                   "--reads",
                   " ",
                   FILES02[1L],
                   ",",
                   FILES02[2L],
                   " ",
                   "--cores",
                   " ",
                   "4",
                   " ",
                   "--memory",
                   " ",
                   "32",
                   " > ",
                   "fastq/SKESA_ASSEMBLY.fa")

system(command = RUNSKESA,
       timeout = 18000)

SKESA <- try(readDNAStringSet("fastq/SKESA_ASSEMBLY.fa"))

if (!is(object = SKESA,
        class2 = "try-error")) {
  # remove Ns at start or end of contigs
  # for (m1 in seq_along(SKESA)) {
  #   N <- gregexpr(pattern = "N",
  #                 text = SKESA[[m1]])
  #   LAST <- width(SKESA[m1])
  #   if (1L %in% N[[1L]]) {
  #     SKESA[m1] <- subseq(x = SKESA[m1],
  #                         start = 2L,
  #                         end = LAST)
  #   }
  #   if (LAST %in% N[[1L]]) {
  #     SKESA[m1] <- subseq(x = SKESA[m1],
  #                         start = 1L,
  #                         end = LAST - 1L)
  #   }
  # }
  
  KEEPSEQ <- rep(TRUE, length(SKESA))
  NTFreq <- alphabetFrequency(x = SKESA)
  WIDTH <- width(SKESA)
  for (m1 in seq_along(SKESA)) {
    if (WIDTH[m1] < 200L) {
      KEEPSEQ[m1] <- FALSE
    }
    if (any(NTFreq[m1, ] == WIDTH[m1])) {
      KEEPSEQ[m1] <- FALSE
    }
  }
  
  SKESA <- SKESA[KEEPSEQ]
  
  writeXStringSet(x = SKESA,
                  filepath = paste0("Assembly_SKESA",
                                    formatC(x = as.integer(OUTID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
}

###### -- SPADES --------------------------------------------------------------

RUNSPADES <- paste0("spades.py",
                    " ",
                    "-1 ",
                    FILES02[1L],
                    " ",
                    "-2 ",
                    FILES02[2L],
                    " ",
                    "--isolate",
                    " ",
                    "--cov-cutoff auto",
                    " ",
                    "-t 4",
                    " ",
                    "-o SPADES_ASSEMBLY")
system(command = RUNSPADES,
       timeout = 7200,
       intern = TRUE)

system(command = "ls SPADES_ASSEMBLY")
system(command = "find . -name contigs.fasta", intern = TRUE) -> L
print(L)

SPADES <- try(readDNAStringSet(filepath = "SPADES_ASSEMBLY/contigs.fasta"))
if (!is(object = SPADES,
        class2 = "try-error")) {
  # remove Ns at start or end of contigs
  # for (m1 in seq_along(SPADES)) {
  #   N <- gregexpr(pattern = "N",
  #                 text = SPADES[[m1]])
  #   LAST <- width(SPADES[m1])
  #   if (1L %in% N[[1L]]) {
  #     SPADES[m1] <- subseq(x = SPADES[m1],
  #                          start = 2L,
  #                          end = LAST)
  #   }
  #   if (LAST %in% N[[1L]]) {
  #     SPADES[m1] <- subseq(x = SPADES[m1],
  #                          start = 1L,
  #                          end = LAST - 1L)
  #   }
  # }
  
  KEEPSEQ <- rep(TRUE, length(SPADES))
  NTFreq <- alphabetFrequency(x = SPADES)
  WIDTH <- width(SPADES)
  for (m1 in seq_along(SPADES)) {
    if (WIDTH[m1] < 200L) {
      KEEPSEQ[m1] <- FALSE
    }
    if (any(NTFreq[m1, ] == WIDTH[m1])) {
      KEEPSEQ[m1] <- FALSE
    }
  }
  
  SPADES <- SPADES[KEEPSEQ]
  
  writeXStringSet(x = SPADES,
                  filepath = paste0("Assembly_SPADES",
                                    formatC(x = as.integer(OUTID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
}

###### -- MEGAHIT -------------------------------------------------------------

RUNMEGAHIT <- paste0("megahit",
                     " -1 ",
                     FILES02[1L],
                     " -2 ",
                     FILES02[2L],
                     " -m 0.5",
                     " -t 4",
                     " -o MEGAHIT_ASSEMBLY")

MEGAHITOUT <- system(command = RUNMEGAHIT,
                     timeout = 7200,
                     intern = TRUE)

MEGAHIT <- try(readDNAStringSet("MEGAHIT_ASSEMBLY/final.contigs.fa"))
if (!is(object = MEGAHIT,
        class2 = "try-error")) {
  # remove Ns at start or end of contigs
  # for (m1 in seq_along(MEGAHIT)) {
  #   N <- gregexpr(pattern = "N",
  #                 text = MEGAHIT[[m1]])
  #   LAST <- width(MEGAHIT[m1])
  #   if (1L %in% N[[1L]]) {
  #     MEGAHIT[m1] <- subseq(x = MEGAHIT[m1],
  #                          start = 2L,
  #                          end = LAST)
  #   }
  #   if (LAST %in% N[[1L]]) {
  #     MEGAHIT[m1] <- subseq(x = MEGAHIT[m1],
  #                          start = 1L,
  #                          end = LAST - 1L)
  #   }
  # }
  
  KEEPSEQ <- rep(TRUE, length(MEGAHIT))
  NTFreq <- alphabetFrequency(x = MEGAHIT)
  WIDTH <- width(MEGAHIT)
  for (m1 in seq_along(MEGAHIT)) {
    if (WIDTH[m1] < 200L) {
      KEEPSEQ[m1] <- FALSE
    }
    if (any(NTFreq[m1, ] == WIDTH[m1])) {
      KEEPSEQ[m1] <- FALSE
    }
  }
  
  MEGAHIT <- MEGAHIT[KEEPSEQ]
  
  writeXStringSet(x = MEGAHIT,
                  filepath = paste0("Assembly_Megahit",
                                    formatC(x = as.integer(OUTID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
}

###### -- UNICYCLER -----------------------------------------------------------

RUNUNICYCLER <- paste0("unicycler",
                       " ",
                       "-1 ",
                       FILES02[1L],
                       " ",
                       "-2 ",
                       FILES02[2L],
                       " ",
                       "-o UNICYCLER_ASSEMBLY",
                       " ",
                       "--min_fasta_length 200",
                       " ",
                       "-t 4")

UNICYCLEROUT <- system(command = RUNUNICYCLER,
                       timeout = 18000,
                       intern = TRUE)
head(UNICYCLEROUT)
tail(UNICYCLEROUT)
system(command = "find . -name final.contigs.fa", intern = TRUE) -> x
print(x)

UNICYCLER <- try(readDNAStringSet("UNICYCLER_ASSEMBLY/assembly.fasta"))
if (!is(object = UNICYCLER,
        class2 = "try-error")) {
  # remove Ns at start or end of contigs
  # for (m1 in seq_along(UNICYCLER)) {
  #   N <- gregexpr(pattern = "N",
  #                 text = UNICYCLER[[m1]])
  #   LAST <- width(UNICYCLER[m1])
  #   if (1L %in% N[[1L]]) {
  #     UNICYCLER[m1] <- subseq(x = UNICYCLER[m1],
  #                          start = 2L,
  #                          end = LAST)
  #   }
  #   if (LAST %in% N[[1L]]) {
  #     UNICYCLER[m1] <- subseq(x = UNICYCLER[m1],
  #                          start = 1L,
  #                          end = LAST - 1L)
  #   }
  # }
  
  KEEPSEQ <- rep(TRUE, length(UNICYCLER))
  NTFreq <- alphabetFrequency(x = UNICYCLER)
  WIDTH <- width(UNICYCLER)
  for (m1 in seq_along(UNICYCLER)) {
    if (WIDTH[m1] < 200L) {
      KEEPSEQ[m1] <- FALSE
    }
    if (any(NTFreq[m1, ] == WIDTH[m1])) {
      KEEPSEQ[m1] <- FALSE
    }
  }
  
  UNICYCLER <- UNICYCLER[KEEPSEQ]
  
  writeXStringSet(x = UNICYCLER,
                  filepath = paste0("Assembly_Unicycler",
                                    formatC(x = as.integer(OUTID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
}

system(command = "ls -lh")

print("Script Completed!")

SCRIPTEND <- Sys.time()
print(SCRIPTEND - SCRIPTSTART)


