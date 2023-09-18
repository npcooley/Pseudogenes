###### -- assemble a single set of reads --------------------------------------

# simulate reads based on available read simulaters
# then assemble them in a uniform manner
suppressMessages(library(SynExtend))

ARGS <- commandArgs(trailingOnly = TRUE)

print(ARGS)
PersistentID <- as.integer(ARGS[1L])
QUAL <- ARGS[2L]
DEPTH <- ARGS[3L]
LENGTH <- ARGS[4L]
MODEL <- ARGS[5L]
FRAG <- if (LENGTH == 100) {
  400L
} else if (LENGTH == 150) {
  500L
} else if (LENGTH == 200) {
  600L
} else if (LENGTH == 250) {
  700L
}
DNA <- "ASSEMBLY01/assembly.fasta"

system(command = "gunzip Assembly.fna.gz")
REF3 <- "Assembly.fna"

# simulate reads
# fragment argument forces paired end simulations

# paired end simulations
system(paste("art_illumina",
             "-ss",
             MODEL,
             "-i",
             REF3,
             "-l",
             LENGTH,
             "-m",
             FRAG,
             "-s 10",
             "-o PAIRED_v01",
             "-f",
             DEPTH,
             "-rs",
             PersistentID,
             "-qs",
             QUAL,
             "-qs2",
             QUAL),
       intern = TRUE)

# single end simulations
system(paste("art_illumina",
             "-ss",
             MODEL,
             "-i",
             REF3,
             "-l",
             LENGTH,
             "-s 10",
             "-o SINGLE_v01",
             "-f",
             DEPTH,
             "-rs",
             PersistentID,
             "-qs",
             QUAL,
             "-qs2",
             QUAL),
       intern = TRUE)

# pb
system(command = paste("pbsim",
                       "--prefix pb",
                       "--strategy wgs",
                       "--method qshmm",
                       "--qshmm /pbsim3/data/QSHMM-RSII.model",
                       "--depth 20",
                       "--accuracy-mean .85",
                       "--seed",
                       PersistentID,
                       "--genome",
                       REF3),
       intern = TRUE)

# ont -- default error difference is for PB, set the recommended ONT
system(command = paste("pbsim",
                       "--prefix ont",
                       "--strategy wgs",
                       "--method qshmm",
                       "--qshmm /pbsim3/data/QSHMM-ONT.model",
                       "--depth 20",
                       "--accuracy-mean .85",
                       "--difference-ratio 39:24:36",
                       "--seed",
                       PersistentID,
                       "--genome",
                       REF3),
       intern = TRUE)

system("rm *.maf *.aln *.ref")

# try assemblers SKESA, Megahit, SPADES, Unicycler, Unicycler w/ONT, Unicycler w/PB
# skesa writes out a single file when successful
# spades writes outputs to a directory

FILES01 <- list.files()
PAIREDREADS <- FILES01[grepl(pattern = "PAIRED",
                             x = FILES01)]
SINGLEREADS <- FILES01[grepl(pattern = "SINGLE",
                             x = FILES01)]
PBREADS <- FILES01[grepl(pattern = "pb_[0-9]+",
                         x = FILES01)]
ONTREADS <- FILES01[grepl(pattern = "ont_[0-9]+",
                          x = FILES01)]

###### -- SKESA PAIRED --------------------------------------------------------
RUNSKESA <- paste0("skesa",
                   " ",
                   "--reads",
                   " ",
                   PAIREDREADS[1L],
                   ",",
                   PAIREDREADS[2L],
                   " ",
                   "--cores",
                   " ",
                   "4",
                   " ",
                   "--memory",
                   " ",
                   "32",
                   " > ",
                   "SKESA_ASSEMBLY_Paired.fa")

system(command = RUNSKESA,
       timeout = 18000)

SKESA <- try(readDNAStringSet("SKESA_ASSEMBLY_Paired.fa"))

if (!is(object = SKESA,
        class2 = "try-error")) {
  # if it's not a try error
  # remove any weird contigs and any contigs shorter than 200 nucs
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
                  filepath = paste0("Assembly_SKESA_Paired_",
                                    formatC(x = as.integer(PersistentID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
}

###### -- SKESA SINGLE --------------------------------------------------------
RUNSKESA <- paste0("skesa",
                   " ",
                   "--reads",
                   " ",
                   SINGLEREADS[1L],
                   " ",
                   "--cores",
                   " ",
                   "4",
                   " ",
                   "--memory",
                   " ",
                   "32",
                   " > ",
                   "SKESA_ASSEMBLY_Single.fa")

system(command = RUNSKESA,
       timeout = 18000)

SKESA <- try(readDNAStringSet("SKESA_ASSEMBLY_Single.fa"))

if (!is(object = SKESA,
        class2 = "try-error")) {
  # if it's not a try error
  # remove any weird contigs and any contigs shorter than 200 nucs
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
                  filepath = paste0("Assembly_SKESA_Single_",
                                    formatC(x = as.integer(PersistentID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
}

system(command = "rm -rf *.fa")

###### -- SPADES PAIRED -------------------------------------------------------

RUNSPADES <- paste0("spades.py",
                    " ",
                    "-1 ",
                    PAIREDREADS[1L],
                    " ",
                    "-2 ",
                    PAIREDREADS[2L],
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

SPADES <- try(readDNAStringSet("SPADES_ASSEMBLY/contigs.fasta"))

if (!is(object = SPADES,
        class2 = "try-error")) {
  # if it's not a try error
  # remove any weird contigs and any contigs shorter than 200 nucs
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
                  filepath = paste0("Assembly_SPADES_Paired_",
                                    formatC(x = as.integer(PersistentID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
}

system(command = "rm -rf SPADES_ASSEMBLY")

###### -- SPADES SINGLE -------------------------------------------------------

RUNSPADES <- paste0("spades.py",
                    " ",
                    "-s ",
                    SINGLEREADS[1L],
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

SPADES <- try(readDNAStringSet("SPADES_ASSEMBLY/contigs.fasta"))

if (!is(object = SPADES,
        class2 = "try-error")) {
  # if it's not a try error
  # remove any weird contigs and any contigs shorter than 200 nucs
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
                  filepath = paste0("Assembly_SPADES_Single_",
                                    formatC(x = as.integer(PersistentID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
}

system(command = "rm -rf SPADES_ASSEMBLY")

###### -- MEGAHIT PAIRED ------------------------------------------------------

RUNMEGAHIT <- paste0("megahit",
                     " -1 ",
                     PAIREDREADS[1L],
                     " -2 ",
                     PAIREDREADS[2L],
                     " -m 0.5",
                     " -t 4",
                     " -o MEGAHIT_ASSEMBLY")

system(command = RUNMEGAHIT,
       timeout = 7200,
       intern = TRUE)

MEGAHIT <- try(readDNAStringSet("MEGAHIT_ASSEMBLY/final.contigs.fa"))

if (!is(object = MEGAHIT,
        class2 = "try-error")) {
  # if it's not a try error
  # remove any weird contigs and any contigs shorter than 200 nucs
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
                  filepath = paste0("Assembly_MEGAHIT_Paired_",
                                    formatC(x = as.integer(PersistentID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
}

system(command = "rm -rf MEGAHIT_ASSEMBLY")

###### -- MEGAHIT SINGLE ------------------------------------------------------

RUNMEGAHIT <- paste0("megahit",
                     " -r ",
                     SINGLEREADS[1L],
                     " -m 0.5",
                     " -t 4",
                     " -o MEGAHIT_ASSEMBLY")

system(command = RUNMEGAHIT,
       timeout = 7200,
       intern = TRUE)

MEGAHIT <- try(readDNAStringSet("MEGAHIT_ASSEMBLY/final.contigs.fa"))

if (!is(object = MEGAHIT,
        class2 = "try-error")) {
  # if it's not a try error
  # remove any weird contigs and any contigs shorter than 200 nucs
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
                  filepath = paste0("Assembly_MEGAHIT_Single_",
                                    formatC(x = as.integer(PersistentID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
}

system(command = "rm -rf MEGAHIT_ASSEMBLY")

###### -- UNICYCLER PAIRED ----------------------------------------------------
# set spades kmers conservatively for better data capture

RUNUNICYCLER <- paste0("unicycler",
                       " ",
                       "-1 ",
                       PAIREDREADS[1L],
                       " ",
                       "-2 ",
                       PAIREDREADS[2L],
                       " ",
                       if (LENGTH == 100) {
                         "--kmers 21,39,53,65,75"
                       } else if (LENGTH == 150) {
                         "--kmers 21,39,53,65,75,83,89,95"
                       } else {
                         ""
                       },
                       " ",
                       "-o UNICYCLER_ASSEMBLY",
                       " ",
                       "--min_fasta_length 200",
                       " ",
                       "-t 4")

system(command = RUNUNICYCLER,
       timeout = 18000,
       intern = TRUE)


UNICYCLER <- try(readDNAStringSet("UNICYCLER_ASSEMBLY/assembly.fasta"))

if (!is(object = UNICYCLER,
        class2 = "try-error")) {
  # if it's not a try error
  # remove any weird contigs and any contigs shorter than 200 nucs
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
                  filepath = paste0("Assembly_UNICYCLER_Paired_",
                                    formatC(x = as.integer(PersistentID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
}

system(command = "rm -rf UNICYCLER_ASSEMBLY")

###### -- UNICYCLER SINGLE ----------------------------------------------------
# set spades kmers conservatively for better data capture

RUNUNICYCLER <- paste0("unicycler",
                       " ",
                       "-s ",
                       SINGLEREADS[1L],
                       " ",
                       if (LENGTH == 100) {
                         "--kmers 21,39,53,65,75"
                       } else if (LENGTH == 150) {
                         "--kmers 21,39,53,65,75,83,89,95"
                       } else {
                         ""
                       },
                       " ",
                       "-o UNICYCLER_ASSEMBLY",
                       " ",
                       "--min_fasta_length 200",
                       " ",
                       "-t 4")

system(command = RUNUNICYCLER,
       timeout = 18000,
       intern = TRUE)


UNICYCLER <- try(readDNAStringSet("UNICYCLER_ASSEMBLY/assembly.fasta"))

if (!is(object = UNICYCLER,
        class2 = "try-error")) {
  # if it's not a try error
  # remove any weird contigs and any contigs shorter than 200 nucs
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
                  filepath = paste0("Assembly_UNICYCLER_Single_",
                                    formatC(x = as.integer(PersistentID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
}

system(command = "rm -rf UNICYCLER_ASSEMBLY")

###### -- Paired End PB Hybrid ------------------------------------------------

RUNUNICYCLER <- paste0("unicycler",
                       " ",
                       "-1 ",
                       PAIREDREADS[1L],
                       " ",
                       "-2 ",
                       PAIREDREADS[2L],
                       " ",
                       "-l",
                       PBREADS[1L],
                       " ",
                       if (LENGTH == 100) {
                         "--kmers 21,39,53,65,75"
                       } else if (LENGTH == 150) {
                         "--kmers 21,39,53,65,75,83,89,95"
                       } else {
                         ""
                       },
                       " ",
                       "-o UNICYCLER_ASSEMBLY",
                       " ",
                       "--min_fasta_length 200",
                       " ",
                       "-t 4")

system(command = RUNUNICYCLER,
       timeout = 18000,
       intern = TRUE)


UNICYCLER <- try(readDNAStringSet("UNICYCLER_ASSEMBLY/assembly.fasta"))

if (!is(object = UNICYCLER,
        class2 = "try-error")) {
  # if it's not a try error
  # remove any weird contigs and any contigs shorter than 200 nucs
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
                  filepath = paste0("Assembly_UNICYCLER_Paired_PB_",
                                    formatC(x = as.integer(PersistentID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
}

system(command = "rm -rf UNICYCLER_ASSEMBLY")

###### -- Paired End ONT Hybrid ------------------------------------------------

RUNUNICYCLER <- paste0("unicycler",
                       " ",
                       "-1 ",
                       PAIREDREADS[1L],
                       " ",
                       "-2 ",
                       PAIREDREADS[2L],
                       " ",
                       "-l",
                       ONTREADS[1L],
                       " ",
                       if (LENGTH == 100) {
                         "--kmers 21,39,53,65,75"
                       } else if (LENGTH == 150) {
                         "--kmers 21,39,53,65,75,83,89,95"
                       } else {
                         ""
                       },
                       " ",
                       "-o UNICYCLER_ASSEMBLY",
                       " ",
                       "--min_fasta_length 200",
                       " ",
                       "-t 4")

system(command = RUNUNICYCLER,
       timeout = 18000,
       intern = TRUE)


UNICYCLER <- try(readDNAStringSet("UNICYCLER_ASSEMBLY/assembly.fasta"))

if (!is(object = UNICYCLER,
        class2 = "try-error")) {
  # if it's not a try error
  # remove any weird contigs and any contigs shorter than 200 nucs
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
                  filepath = paste0("Assembly_UNICYCLER_Paired_ONT_",
                                    formatC(x = as.integer(PersistentID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
}

system(command = "rm -rf UNICYCLER_ASSEMBLY")

###### -- Single Ended PB hybrid ----------------------------------------------

RUNUNICYCLER <- paste0("unicycler",
                       " ",
                       "-s ",
                       SINGLEREADS[1L],
                       " ",
                       "-l",
                       PBREADS[1L],
                       " ",
                       if (LENGTH == 100) {
                         "--kmers 21,39,53,65,75"
                       } else if (LENGTH == 150) {
                         "--kmers 21,39,53,65,75,83,89,95"
                       } else {
                         ""
                       },
                       " ",
                       "-o UNICYCLER_ASSEMBLY",
                       " ",
                       "--min_fasta_length 200",
                       " ",
                       "-t 4")

system(command = RUNUNICYCLER,
       timeout = 18000,
       intern = TRUE)


UNICYCLER <- try(readDNAStringSet("UNICYCLER_ASSEMBLY/assembly.fasta"))

if (!is(object = UNICYCLER,
        class2 = "try-error")) {
  # if it's not a try error
  # remove any weird contigs and any contigs shorter than 200 nucs
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
                  filepath = paste0("Assembly_UNICYCLER_Single_PB_",
                                    formatC(x = as.integer(PersistentID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
}

system(command = "rm -rf UNICYCLER_ASSEMBLY")

###### -- Single Ended ONT hybrid ---------------------------------------------

RUNUNICYCLER <- paste0("unicycler",
                       " ",
                       "-s ",
                       SINGLEREADS[1L],
                       " ",
                       "-l",
                       ONTREADS[1L],
                       " ",
                       if (LENGTH == 100) {
                         "--kmers 21,39,53,65,75"
                       } else if (LENGTH == 150) {
                         "--kmers 21,39,53,65,75,83,89,95"
                       } else {
                         ""
                       },
                       " ",
                       "-o UNICYCLER_ASSEMBLY",
                       " ",
                       "--min_fasta_length 200",
                       " ",
                       "-t 4")

system(command = RUNUNICYCLER,
       timeout = 18000,
       intern = TRUE)


UNICYCLER <- try(readDNAStringSet("UNICYCLER_ASSEMBLY/assembly.fasta"))

if (!is(object = UNICYCLER,
        class2 = "try-error")) {
  # if it's not a try error
  # remove any weird contigs and any contigs shorter than 200 nucs
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
                  filepath = paste0("Assembly_UNICYCLER_Single_ONT_",
                                    formatC(x = as.integer(PersistentID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
}

system(command = "rm -rf UNICYCLER_ASSEMBLY")
system(command = "rm -rf *.fastq *.fq Assembly.fna")


