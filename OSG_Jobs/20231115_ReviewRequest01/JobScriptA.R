###### -- given a set of inputs, perform a reassembly given some things -------
# JobScriptA
# set out a starting point given what we assume the assembly size will be
# assemble with Unicycler at what should be approx 1000x coverage of short reads
# and when long reads are present, one of the three designated long read coverages
# check the actual coverage post assembly and report it

###### -- libraries -----------------------------------------------------------

suppressMessages(library(SynExtend))
suppressMessages(library(Rsamtools))

###### -- arguments -----------------------------------------------------------

ARGS <- commandArgs(trailingOnly = TRUE)
PersistentID <- as.integer(ARGS[1L])
ImpliedAssembly_Size <- ARGS[2L]
SR_Run <- ARGS[3L]
LR_Run <- ARGS[4L]
SR_Platform <- ARGS[5L]
LR_Platform <- ARGS[6L]
SR_Target <- ARGS[7L]
LR_Target <- ARGS[8L]
SR_Approx_Cov <- ARGS[9L]
LR_Approx_Cov <- ARGS[10L]

if (LR_Run == "x") {
  SR_ONLY <- TRUE
  SR_Target <- as.integer(SR_Target)
  SR_Approx_Cov <- as.numeric(SR_Approx_Cov)
} else {
  SR_ONLY <- FALSE
  SR_Target <- as.integer(SR_Target)
  LR_Target <- as.integer(LR_Target)
  SR_Approx_Cov <- as.numeric(SR_Approx_Cov)
  LR_Approx_Cov <- as.numeric(LR_Approx_Cov)
}

###### -- ad hoc functions ----------------------------------------------------

GetDepth_LR <- function(Ref,
                        Reads1,
                        MINIMAP_ARGS = list("x" = "map-ont")) {
  
  # tmp01 <- tempfile()
  SAMTOOLS_INDEX <- paste0("samtools faidx ",
                           Ref)
  
  system(command = SAMTOOLS_INDEX,
         intern = FALSE)
  
  tmp02 <- tempfile()
  MINIMAP2_COMMAND <- paste0("minimap2 -ax",
                             " ",
                             MINIMAP_ARGS[names(MINIMAP_ARGS) == "x"],
                             " ",
                             Ref,
                             " ",
                             Reads1,
                             " -o ",
                             tmp02)
  system(command = MINIMAP2_COMMAND,
         intern = FALSE)
  print(MINIMAP2_COMMAND)
  
  tmp03 <- tempfile()
  SAMTOOLS_CONVERT <- paste0("samtools view -S -b ",
                             tmp02,
                             " -o ",
                             tmp03)
  system(command = SAMTOOLS_CONVERT,
         intern = FALSE)
  
  tmp04 <- tempfile()
  SAMTOOLS_SORT <- paste0("samtools sort ",
                          tmp03,
                          " -o ",
                          tmp04)
  
  system(command = SAMTOOLS_SORT,
         intern = FALSE)
  
  tmp05 <- tempfile()
  SAMTOOLS_DEPTH <- paste0("samtools depth -aa ",
                           tmp04,
                           " -o ",
                           tmp05)
  
  system(command = SAMTOOLS_DEPTH,
         intern = FALSE)
  
  res <- read.table(file = tmp05,
                    sep = "\t",
                    header = FALSE)
  system(command = paste("rm",
                         paste(tmp03,
                               tmp04,
                               tmp05)))
  # system(command = paste0("rm ",
  #                         tmp02,
  #                         ".*"))
  return(res)
  
}

# get depth for short reads, paired end or otherwise
GetDepth_SR <- function(Ref,
                        Reads1,
                        Reads2,
                        BOWTIE_ARGS = list("-t" = NULL,
                                           "--no-unal" = NULL,
                                           "-p" = 6,
                                           "--ignore-quals" = NULL,
                                           "--end-to-end" = NULL,
                                           "-k" = 2000,
                                           "--xeq" = NULL,
                                           "--reorder" = NULL)) {
  
  # tmp01 <- tempfile()
  SAMTOOLS_INDEX <- paste0("samtools faidx ",
                           Ref)
  
  system(command = SAMTOOLS_INDEX,
         intern = FALSE)
  
  tmp02 <- tempfile()
  BOWTIE_INDEX <- paste0("bowtie2-build -q ",
                         '"',
                         Ref,
                         '"',
                         " ",
                         '"',
                         tmp02,
                         '"')
  
  print("Bowtie Index Command:")
  cat(paste0("\n",
             BOWTIE_INDEX,
             "\n"))
  system(command = BOWTIE_INDEX,
         intern = FALSE)
  
  tmp03 <- tempfile()
  w1 <- unname(sapply(X = BOWTIE_ARGS,
                      FUN = function(x) {
                        is.null(x)
                      },
                      USE.NAMES = FALSE))
  if (missing(Reads2)) {
    # no reads2, reads are single ended
    BOWTIE_MAP <- paste("bowtie2",
                        paste(names(BOWTIE_ARGS)[w1],
                              collapse = " "),
                        paste(names(BOWTIE_ARGS)[!w1],
                              unlist(BOWTIE_ARGS[!w1]),
                              collapse = " "),
                        "-x",
                        tmp02,
                        "-U",
                        Reads1,
                        "-S",
                        tmp03)
    
  } else {
    # reads are paired
    BOWTIE_MAP <- paste("bowtie2",
                        paste(names(BOWTIE_ARGS)[w1],
                              collapse = " "),
                        paste(names(BOWTIE_ARGS)[!w1],
                              unlist(BOWTIE_ARGS[!w1]),
                              collapse = " "),
                        "-x",
                        tmp02,
                        "-1",
                        Reads1,
                        "-2",
                        Reads2,
                        "-S",
                        tmp03)
  }
  
  print("Bowtie command:")
  print(BOWTIE_MAP)
  system(command = BOWTIE_MAP,
         intern = FALSE)
  
  tmp04 <- tempfile()
  SAMTOOLS_CONVERT <- paste0("samtools view -S -b ",
                             tmp03,
                             " -o ",
                             tmp04)
  
  system(command = SAMTOOLS_CONVERT,
         intern = FALSE)
  
  tmp05 <- tempfile()
  SAMTOOLS_SORT <- paste0("samtools sort ",
                          tmp04,
                          " -o ",
                          tmp05)
  
  system(command = SAMTOOLS_SORT,
         intern = FALSE)
  
  tmp06 <- tempfile()
  SAMTOOLS_DEPTH <- paste0("samtools depth -aa ",
                           tmp05,
                           " -o ",
                           tmp06)
  
  system(command = SAMTOOLS_DEPTH,
         intern = FALSE)
  
  res <- read.table(file = tmp06,
                    sep = "\t",
                    header = FALSE)
  # res <- split(x = res[, 3L],
  #              f = res[, 1L])
  
  system(command = paste("rm",
                         paste(tmp03,
                               tmp04,
                               tmp05,
                               tmp06)))
  system(command = paste0("rm ",
                          tmp02,
                          ".*"))
  return(res)
}

###### -- input data ----------------------------------------------------------

# not loading in any data directly from file transfers

###### -- code body -----------------------------------------------------------

# writeQualityScaledXStringSet reproducibly crashes when writing out large
# long read fastqs, so we're going to always convert them to quality score-less
# fastas -- not necessarily the most ideal solution, but we work with what we've got

# pull our fastqs into R
# when KEEPFILES is TRUE, the fastqs are saved into the working directory
# we don't want to do this on the grid because they will be sent back to our
# submit directory, also we're doing all our subsetting here in R

PullReadsStart <- Sys.time()
if (SR_ONLY) {
  t1 <- Sys.time()
  SR_Reads <- FastQFromSRR(SRR = SR_Run,
                           KEEPFILES = FALSE)
  t2 <- Sys.time()
  print("Short reads collected in:")
  print(t2 - t1)
  SR_Prop <- SR_Target / SR_Approx_Cov
  if (length(SR_Reads) == 1L) {
    # single ended
    Reads1 <- SR_Reads[[1]]
    Reads2 <- NULL
  } else {
    # paired end
    Reads1 <- SR_Reads[[1]]
    Reads2 <- SR_Reads[[2]]
  }
  s1 <- ceiling(length(Reads1) * SR_Prop)
  set.seed(PersistentID)
  s1 <- sample(x = length(Reads1),
               size = s1,
               replace = FALSE)
  
  Reads1 <- Reads1[s1]
  Reads2 <- Reads2[s1]
  Reads3 <- NULL
} else {
  t1 <- Sys.time()
  SR_Reads <- FastQFromSRR(SRR = SR_Run,
                           KEEPFILES = FALSE)
  t2 <- Sys.time()
  print("Short reads collected in:")
  print(t2 - t1)
  if (length(SR_Reads) == 1L) {
    # single ended
    Reads1 <- SR_Reads[[1]]
    Reads2 <- NULL
  } else {
    # paired end
    Reads1 <- SR_Reads[[1]]
    Reads2 <- SR_Reads[[2]]
  }
  t1 <- Sys.time()
  LR_Reads <- FastQFromSRR(SRR = LR_Run,
                           KEEPFILES = FALSE)
  t2 <- Sys.time()
  print("Long reads collected in:")
  print(t2 - t1)
  Reads3 <- LR_Reads[[1]]
  SR_Prop <- SR_Target / SR_Approx_Cov
  LR_Prop <- LR_Target / LR_Approx_Cov
  s1 <- ceiling(length(Reads1) * SR_Prop)
  s2 <- ceiling(length(Reads3) * LR_Prop)
  set.seed(PersistentID)
  s1 <- sample(x = length(Reads1),
               size = s1,
               replace = FALSE)
  s2 <- sample(x = length(Reads3),
               size = s2,
               replace = FALSE)
  
  Reads1 <- Reads1[s1]
  Reads2 <- Reads2[s1]
  Reads3 <- Reads3[s2]
}

PullReadsEnd <- Sys.time()
PullReadsTotalTime <- PullReadsEnd - PullReadsStart
print("Reads collected in:")
print(PullReadsTotalTime)

# clean up some things
rm(list = c("SR_Reads",
            "LR_Reads"))

tmp1 <- tempfile()
tmp1ext <- paste0(tmp1,
                  ".fastq.gz")
tmp2 <- tempfile()
tmp2ext <- paste0(tmp2,
                  ".fastq.gz")
tmp3 <- tempfile()
tmp3ext <- paste0(tmp3,
                  ".fasta.gz")
tmp4 <- tempfile()
writeQualityScaledXStringSet(x = Reads1,
                             filepath = tmp1ext,
                             compress = TRUE)
if (!is.null(Reads2)) {
  writeQualityScaledXStringSet(x = Reads2,
                               filepath = tmp2ext,
                               compress = TRUE)
}
if (!is.null(Reads3)) {
  writeXStringSet(x = Reads3,
                  filepath = tmp3ext,
                  compress = TRUE)
}

# assemble with unicycler
if (is.null(Reads2) & is.null(Reads3)) {
  # single ended short reads only
  RUN_UNICYCLER <- paste("unicycler",
                         "--unpaired",
                         tmp1ext,
                         "--out",
                         tmp4)
  
} else if (is.null(Reads2) & !is.null(Reads3)) {
  # single ended short reads + long reads
  RUN_UNICYCLER <- paste("unicycler",
                         "--unpaired",
                         tmp1ext,
                         "--long",
                         tmp3ext,
                         "--out",
                         tmp4)
  
} else if (!is.null(Reads2) & !is.null(Reads3)) {
  # paired end short reads + long reads
  RUN_UNICYCLER <- paste("unicycler",
                         "--short1",
                         tmp1ext,
                         "--short2",
                         tmp2ext,
                         "--long",
                         tmp3ext,
                         "--out",
                         tmp4)
  
} else if (!is.null(Reads2) & is.null(Reads3)) {
  # paired end short reads only
  RUN_UNICYCLER <- paste("unicycler",
                         "--short1",
                         tmp1ext,
                         "--short2",
                         tmp2ext,
                         "--out",
                         tmp4)
  
}

Unicycler_Start <- Sys.time()
x <- system(command = RUN_UNICYCLER,
            intern = TRUE,
            timeout = 14400)
Unicycler_End <- Sys.time()

UnicyclerTotalTime <- Unicycler_End - Unicycler_Start
print("Unicycler completed in:")
print(UnicyclerTotalTime)

res01 <- try(readDNAStringSet(filepath = paste0(tmp4,
                                                "/assembly.fasta")),
             silent = TRUE)

if (class(res01) != "try-error") {
  print("Successful assembly!")
  RefLocation <- paste0(tmp4,
                        "/assembly.fasta")
}

# map used reads back to the assembly
# get the 'TRUE' read coverage depth
# for any set of reads present
if (is.null(Reads2) & is.null(Reads3)) {
  # single ended short reads only
  dep1 <- GetDepth_SR(Ref = RefLocation,
                      Reads1 = tmp1ext)
  dep2 <- NULL
  
} else if (is.null(Reads2) & !is.null(Reads3)) {
  # single ended short reads + long reads
  dep1 <- GetDepth_SR(Ref = RefLocation,
                      Reads1 = tmp1ext)
  if (LR_Platform == "OXFORD_NANOPORE") {
    dep2 <- GetDepth_LR(Ref = RefLocation,
                        Reads1 = tmp3ext,
                        MINIMAP_ARGS = list("x" = "map-ont"))
  } else if (LR_Platform == "PACBIO_SMRT") {
    dep2 <- GetDepth_LR(Ref = RefLocation,
                        Reads1 = tmp3ext,
                        MINIMAP_ARGS = list("x" = "map-pb"))
  }
} else if (!is.null(Reads2) & !is.null(Reads3)) {
  # paired end short reads + long reads
  dep1 <- GetDepth_SR(Ref = RefLocation,
                      Reads1 = tmp1ext,
                      Reads2 = tmp2ext)
  if (LR_Platform == "OXFORD_NANOPORE") {
    dep2 <- GetDepth_LR(Ref = RefLocation,
                        Reads1 = tmp3ext,
                        MINIMAP_ARGS = list("x" = "map-ont"))
  } else if (LR_Platform == "PACBIO_SMRT") {
    dep2 <- GetDepth_LR(Ref = RefLocation,
                        Reads1 = tmp3ext,
                        MINIMAP_ARGS = list("x" = "map-pb"))
  }
} else if (!is.null(Reads2) & is.null(Reads3)) {
  # paired end short reads only
  dep1 <- GetDepth_SR(Ref = RefLocation,
                      Reads1 = tmp1ext,
                      Reads2 = tmp2ext)
  dep2 <- NULL
  
}


z1 <- rle(dep1$V1)
wtm1 <- weighted.mean(x = dep1$V3,
                      w = rep(x = (z1$lengths / nrow(dep1)),
                              times = z1$lengths))
sr_dep <- tapply(X = dep1$V3,
                 INDEX = dep1$V1,
                 FUN = function(x) {
                   mean(x)
                 })
circ_val <- sapply(X = strsplit(x = names(res01),
                                split = " ",
                                fixed = TRUE),
                   FUN = function(x) {
                     x[4]
                   })

if (!is.null(dep2)) {
  z1 <- rle(dep2$V1)
  wtm2 <- weighted.mean(x = dep2$V3,
                        w = rep(x = (z1$lengths / nrow(dep2)),
                                times = z1$lengths))
  lr_dep <- tapply(X = dep2$V3,
                   INDEX = dep2$V1,
                   FUN = function(x) {
                     mean(x)
                   })
  contig_names <- paste0("contig_",
                         seq(length(res01)),
                         " ",
                         sr_dep,
                         " ",
                         lr_dep,
                         " ",
                         circ_val)
} else {
  wtm2 <- NULL
  lr_dep <- NULL
  
  contig_names <- paste0("contig_",
                         seq(length(res01)),
                         " ",
                         sr_dep,
                         " ",
                         circ_val)
}

unicycler_names <- names(res01)
names(res01) <- contig_names

###### -- save results --------------------------------------------------------

writeXStringSet(x = res01,
                filepath = paste0("Result",
                                  formatC(x = PersistentID,
                                          width = 7,
                                          format = "d",
                                          flag = 0),
                                  ".fna.gz"),
                compress = TRUE)

save(dep1,
     dep2,
     unicycler_names,
     UnicyclerTotalTime,
     PullReadsTotalTime,
     file = paste0("Result",
                   formatC(x = PersistentID,
                           width = 7,
                           format = "d",
                           flag = 0),
                   ".RData"),
     compress = "xz")
