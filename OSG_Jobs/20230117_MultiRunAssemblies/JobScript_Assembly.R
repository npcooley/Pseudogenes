###### -- assemble a single set of reads --------------------------------------

suppressMessages(library(SynExtend))

ARGS <- commandArgs(trailingOnly = TRUE)

ID <- ARGS[1L]
SRR <- ARGS[2L]
PLAT <- ARGS[4L]

PREFETCH <- paste("prefetch",
                  "-O sra",
                  SRR)


FASTQDUMP <- paste("fastq-dump",
                   "--outdir fastq",
                   "--gzip",
                   "--skip-technical",
                   "--readids",
                   "--read-filter pass",
                   "--dumpbase",
                   "--split-3",
                   "--clip",
                   paste0("sra/",
                          SRR,
                          "/",
                          SRR,
                          ".sra"))

print(PREFETCH)
system(command = PREFETCH,
       intern = FALSE,
       timeout = 7200)

# check to see that prefetch sent the .sra file to the correct place
list.files(path = "sra",
           full.names = TRUE)
list.files(path = paste0("sra/",
                         SRR),
           full.names = TRUE)

print(FASTQDUMP)
TIME01 <- Sys.time()
system(command = FASTQDUMP,
       intern = FALSE,
       timeout = 7200)
TIME02 <- Sys.time()
print(TIME02 - TIME01)

FILES01 <- list.files(path = "fastq",
                      full.names = TRUE)

if (PLAT %in% c("OXFORD_NANOPORE", "PACBIO_SMRT")) {
  # long reads
  RUNUNICYCLER <- paste("unicycler",
                        "-l",
                        FILES01[1L],
                        "-o ASSEMBLY",
                        "--min_fasta_length 200")
} else {
  if (length(FILES01) == 1L) {
    # SE short reads
    RUNUNICYCLER <- paste("unicycler",
                          "-s",
                          FILES01[1L],
                          "-o ASSEMBLY",
                          "--min_fasta_length 200")
  } else {
    # PE short reads
    RUNUNICYCLER <- paste("unicycler",
                          "-1",
                          FILES01[1L],
                          "-2",
                          FILES01[2L],
                          "-o ASSEMBLY",
                          "--min_fasta_length 200")
  }
}
print(RUNUNICYCLER)

system(command = RUNUNICYCLER,
       timeout = 18000,
       intern = TRUE)

dna <- try(readDNAStringSet(filepath = "UNICYCLER/assembly.fasta"))

if (!is(object = dna,
        class2 = "try-error")) {
  writeXStringSet(x = dna,
                  filepath = paste0("Assembly_",
                                    formatC(x = as.integer(ID),
                                            width = 6L,
                                            flag = 0,
                                            format = "d"),
                                    ".fna.gz"),
                  compress = TRUE)
} else {
  print("missing stringset")
}




