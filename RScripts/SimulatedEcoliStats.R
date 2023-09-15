###### -- plotting out the simulated data -------------------------------------

###### -- parsing simulated assemblies ----------------------------------------

suppressMessages(library(SynExtend))
library(stringr)

GRangeToDFrame <- function(GRangeObject,
                           FeaturesToCollect = c("gene",
                                                 "pseudogene"),
                           Verbose = FALSE) {
  # search until no more searches are necessary
  s1 <- as.character(GRangeObject$type)
  s2 <- GRangeObject$gene
  if (is.null(s2)) {
    s2 <- rep(NA, length(s1))
    s2[s1 %in% FeaturesToCollect] <- paste("unnamed_feature",
                                           seq(sum(s1 %in% FeaturesToCollect)))
  }
  s3 <- GRangeObject@strand
  s4 <- GRangeObject@ranges
  s5 <- as.character(GRangeObject@seqnames)
  s6 <- GRangeObject$Note
  s7 <- GRangeObject$Parent
  s8 <- GRangeObject$ID
  s9 <- !is.na(s2)
  
  # return(list(s1,
  #             s2,
  #             s3,
  #             s4,
  #             s5,
  #             s6,
  #             s7,
  #             s8,
  #             s9))
  
  TOTAL <- sum(table(s1[s1 %in% FeaturesToCollect]))
  
  if (TOTAL == 0) {
    if (Verbose) {
      print("No Features present to collect.")
    }
    return(NULL)
  }
  # print(TOTAL)
  CONTINUE <- TRUE
  KEEP <- vector(mode = "logical",
                 length = length(s1))
  START <- STOP <- vector(mode = "integer",
                          length = length(s1))
  NOTE <- CONTIG <- TYPE <- ID <- NAME <- vector(mode = "character",
                                                 length = length(s1))
  COUNT <- 1L
  FOUNDFEATURES <- 0L
  if (Verbose) {
    pBar <- txtProgressBar(style = 1L)
    TIMESTART <- Sys.time()
  }
  while (CONTINUE) {
    # is the line a line to evaluate
    # check its children
    if (s1[COUNT] %in% FeaturesToCollect) {
      w1 <- which(s7 == s8[COUNT])
      w1 <- which(lengths(w1) > 0L)
      # print(w1)
      # if the feature has any children
      if (length(w1) > 0L) {
        ph1 <- ""
        for (m2 in seq_along(w1)) {
          ph2 <- unlist(s6[w1[m2]])
          # print(ph2)
          if (length(ph2) > 0) {
            if (!is.na(ph2)) {
              # print(nchar(ph1))
              if (nchar(ph1) == 0) {
                ph1 <- ph2
              } else {
                ph1 <- paste(ph1, ph2, sep = "; ")
              }
            }
          } else {
            if (s1[COUNT] == "gene") {
              ph1 <- "normal feature"
            } else {
              ph1 <- "non-coding pseudofeature"
            }
          }
        }
        NOTE[COUNT] <- ph1
      } else {
        # feature has no children, what to do here?
        NOTE[COUNT] <- "child lines absent"
      }
      START[COUNT] <- s4@start[COUNT]
      STOP[COUNT] <- s4@start[COUNT] + s4@width[COUNT] - 1L
      CONTIG[COUNT] <- s5[COUNT]
      TYPE[COUNT] <- s1[COUNT]
      NAME[COUNT] <- s8[COUNT]
      
      if (s9[COUNT]) {
        ID[COUNT] <- s2[COUNT]
      } else {
        ID[COUNT] <- ""
      }
      KEEP[COUNT] <- TRUE
      FOUNDFEATURES <- FOUNDFEATURES + 1L
    }
    
    if (FOUNDFEATURES >= TOTAL) {
      CONTINUE <- FALSE
    } else {
      COUNT <- COUNT + 1L
    }
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = FOUNDFEATURES / TOTAL)
    }
    
  }
  if (Verbose) {
    close(pBar)
    cat("\n")
    TIMEEND <- Sys.time()
    print(TIMEEND - TIMESTART)
  }
  res <- DataFrame("Start" = START[KEEP],
                   "Stop" = STOP[KEEP],
                   "Type" = TYPE[KEEP],
                   "Contig" = CONTIG[KEEP],
                   "ID" = ID[KEEP],
                   "Note" = NOTE[KEEP],
                   "Name" = NAME[KEEP])
  return(res)
}

NoteCheck <- function(NoteVector,
                      CheckVector) {
  Res <- sapply(X = CheckVector,
                FUN = function(y) {
                  sum(grepl(pattern = y,
                            x = NoteVector))
                })
  return(Res)
}


FILES01 <- list.files(path = "~/Data/20230509_SimulatedEcoliReassemblies",
                      pattern = "Annotation")

FILES02 <- list.files(path = "~/Data/20230509_SimulatedEcoliReassemblies",
                      pattern = "Assembly")

FILES03 <- list.files(path = "~/Data/20230509_SimulatedEcoliReassemblies",
                      pattern = "Result")

JobMapB <- read.table("~/Data/20230509_SimulatedEcoliReassemblies/JobMapB.txt")
JobMapA <- read.table("~/Repos/20230404_SimulatedEColiReassemblies/JobMapA.txt")


PseudoTypeVector <- c("frameshifted",
                      "internal stop",
                      "partial abbutting assembly gap",
                      "partial in the middle",
                      "missing C-terminus",
                      "missing N-terminus")


RefAssembly <- readDNAStringSet(filepath = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz")
RefGC <- gffToDataFrame(GFF = "~/Data/20230509_SinglePGAP/annot.gff",
                        Verbose = TRUE)
RefGR <- rtracklayer::import("~/Data/20230509_SinglePGAP/annot.gff")
RefDF <- GRangeToDFrame(GRangeObject = RefGR,
                        Verbose = TRUE)
RefNotes <- NoteCheck(NoteVector = RefDF$Note,
                      CheckVector = PseudoTypeVector)
RefIS <- RefNotes[2L] / (sum(width(RefAssembly)) / 1000000)
RefFS <- RefNotes[1L] / (sum(width(RefAssembly)) / 1000000)

RES <- vector(mode = "list",
              length = length(FILES03))
MAP1 <- as.integer(unlist(regmatches(x = FILES03,
                                     m = gregexpr(pattern = "[0-9]+",
                                                  text = FILES03))))
MAP4 <- strsplit(x = JobMapB[, 2L],
                 split = "_",
                 fixed = TRUE)
MAP2 <- sapply(X = MAP4,
               FUN = function(x) {
                 x[2]
               })
MAP5 <- sapply(X = MAP4,
               FUN = function(x) {
                 x[3]
               })
MAP6 <- sapply(X = MAP4,
               FUN = function(x) {
                 if (length(x) == 5) {
                   x[4]
                 } else {
                   "None"
                 }
               })
# MAP2 <- unlist(regmatches(x = JobMapB[, 2L],
#                           m = gregexpr(pattern = "(?<=Assembly_)([A-Za-z]+)",
#                                        perl = TRUE,
#                                        text = JobMapB[, 2L])))
MAP3 <- as.integer(unlist(regmatches(x = JobMapB[, 2L],
                                     m = gregexpr(pattern = "[0-9]+",
                                                  text = JobMapB[, 2L]))))

Ref_FS_Names <- paste(RefDF$Name[grepl(pattern = "frameshifted", x = RefDF$Note)],
                      "B",
                      sep = "_")
Ref_IS_Names <- paste(RefDF$Name[grepl(pattern = "internal stop", x = RefDF$Note)],
                      "B",
                      sep = "_")

pBar <- txtProgressBar(style = 1L)
PBAR <- length(RES)
for (m1 in seq_along(FILES03)) {
  load(file = paste0("~/Data/20230509_SimulatedEcoliReassemblies",
                     "/",
                     FILES03[m1]),
       verbose = FALSE)
  # reference is B
  a <- c(res,
         list("assembler" = MAP2[MAP1[m1]]),
         list("arrangement" = MAP5[MAP1[m1]]),
         list("lr" = MAP6[MAP1[m1]]),
         list("total_genes" = sum(CurrentGC$Coding)),
         list("count1" = res$rel_fs * (res$total_length / 1000000)),
         list("count2" = res$rel_is * (res$total_length / 1000000)),
         list("found_fs" = sum(Ref_FS_Names %in% Cl1$p2)),
         list("found_is" = sum(Ref_IS_Names %in% Cl1$p2)),
         JobMapA[MAP3[MAP1[m1]], ])
  a <- data.frame(a)
  
  RES[[m1]] <- a
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)

dat1 <- do.call(rbind,
                RES)

dat2 <- data.frame("fs" = dat1$rel_fs - RefFS,
                   "is" = dat1$rel_is - RefIS,
                   "fragment" = dat1$n50 / dat1$total_length,
                   "total_length" = dat1$total_length,
                   "total_contigs" = dat1$tc,
                   # "cov" = log10(dat1$V3) - log10(100),
                   "cov" = log10(dat1$V3 / 50),
                   # "read_len" = log10(dat1$V4) - log10(100),
                   # "read_len" = log10(dat1$V4 / 100),
                   "read_len" = dat1$V4,
                   "assembler" = as.factor(dat1$assembler),
                   "qual" = dat1$V2,
                   "arrangement" = as.factor(paste(dat1$arrangement,
                                                   dat1$lr,
                                                   sep = "+")),
                   "sequencer" = as.factor(dat1$V5),
                   "total_genes" = dat1$total_genes,
                   "abs_fs" = dat1$count1,
                   "abs_is" = dat1$count2,
                   "found_ref_fs" = dat1$found_fs,
                   "found_ref_is" = dat1$found_is)
dat2 <- dat2[dat2$read_len == 150, ]
# subset out assemblies that are very wrong in length
dat3 <- dat2[dat2$total_length >= (sum(width(RefAssembly)) * 0.8) &
               dat2$total_length <= (sum(width(RefAssembly)) * 1.2) &
               dat2$total_genes >= (sum(RefGC$Coding) * 0.8) &
               dat2$total_genes <= (sum(RefGC$Coding) * 1.2), ]



is_dat01 <- by(data = dat3,
               INDICES = ~ assembler + sequencer + arrangement,
               FUN = function(x) {
                 # summary(lm(formula = is ~ cov * qual,
                 #            data = x))$coefficients
                 summary(glm(formula = (abs(abs_is - found_ref_is) / total_genes) ~ cov + qual,
                             data = x,
                             family = binomial(link = "logit"),
                             weights = total_genes))$coefficients
                 # length(x$fs)
               })
fs_dat01 <- by(data = dat3,
               INDICES = ~ assembler + sequencer + arrangement,
               FUN = function(x) {
                 # summary(lm(formula = fs ~ cov * qual,
                 #            data = x))$coefficients
                 summary(glm(formula = (abs(abs_fs - found_ref_fs) / total_genes) ~ cov + qual,
                             data = x,
                             family = binomial(link = "logit"),
                             weights = total_genes))$coefficients
                 # length(x$fs)
               })

is_vals01 <- apply(X = is_dat01,
                   MARGIN = 1:3,
                   FUN = function(x) {
                     x[[1]][, 1]
                   })
is_vals02 <- apply(X = is_dat01,
                   MARGIN = 1:3,
                   FUN = function(x) {
                     x[[1]][, 4]
                   })

is_vals03 <- do.call(rbind,
                     is_vals01[, , 1:6])
is_vals04 <- do.call(rbind,
                     is_vals02[, , 1:6])
is_vals04 <- matrix(data = p.adjust(p = is_vals04,
                                    method = "bonferroni"),
                    ncol = ncol(is_vals04))

fs_vals01 <- apply(X = fs_dat01,
                   MARGIN = 1:3,
                   FUN = function(x) {
                     x[[1]][, 1]
                   })
fs_vals02 <- apply(X = fs_dat01,
                   MARGIN = 1:3,
                   FUN = function(x) {
                     x[[1]][, 4]
                   })

fs_vals03 <- do.call(rbind,
                     fs_vals01[, , 1:6])
fs_vals04 <- do.call(rbind,
                     fs_vals02[, , 1:6])
fs_vals04 <- matrix(data = p.adjust(p = fs_vals04,
                                    method = "bonferroni"),
                    ncol = ncol(fs_vals04))

o1 <- c(1:4, 13, 16) # paired end
o2 <- c(5:8, 14, 17)
o3 <- c(9:12, 15, 18)
o4 <- c(19:22, 31, 34) # single end
o5 <- c(23:26, 32, 35)
o6 <- c(27:30, 33, 36)
fs_vals05 <- cbind(fs_vals03[o1, ],
                   fs_vals03[o2, ],
                   fs_vals03[o3, ])
fs_vals05 <- fs_vals05[, c(1,4,7,2,5,8,3,6,9)]
fs_vals06 <- cbind(fs_vals03[o4, ],
                   fs_vals03[o5, ],
                   fs_vals03[o6, ])
fs_vals06 <- fs_vals06[, c(1,4,7,2,5,8,3,6,9)]
fs_vals07 <- cbind(fs_vals04[o1, ],
                   fs_vals04[o2, ],
                   fs_vals04[o3, ])
fs_vals07 <- fs_vals07[, c(1,4,7,2,5,8,3,6,9)]
fs_vals08 <- cbind(fs_vals04[o4, ],
                   fs_vals04[o5, ],
                   fs_vals04[o6, ])
fs_vals08 <- fs_vals08[, c(1,4,7,2,5,8,3,6,9)]

is_vals05 <- cbind(is_vals03[o1, ],
                   is_vals03[o2, ],
                   is_vals03[o3, ])
is_vals05 <- is_vals05[, c(1,4,7,2,5,8,3,6,9)]
is_vals06 <- cbind(is_vals03[o4, ],
                   is_vals03[o5, ],
                   is_vals03[o6, ])
is_vals06 <- is_vals06[, c(1,4,7,2,5,8,3,6,9)]
is_vals07 <- cbind(is_vals04[o1, ],
                   is_vals04[o2, ],
                   is_vals04[o3, ])
is_vals07 <- is_vals07[, c(1,4,7,2,5,8,3,6,9)]
is_vals08 <- cbind(is_vals04[o4, ],
                   is_vals04[o5, ],
                   is_vals04[o6, ])
is_vals08 <- is_vals08[, c(1,4,7,2,5,8,3,6,9)]

# end of general data prep
save(fs_vals03,
     fs_vals04,
     is_vals03,
     is_vals04,
     dat3,
     file = "~/Repos/PseudogenePlots/InputData/SimulatedEcoliStats.RData",
     compress = "xz")




