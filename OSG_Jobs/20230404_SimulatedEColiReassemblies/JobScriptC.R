###### -- test of large scale reassembly comparisons --------------------------

# for each comparison save off:
# n50
# total contig count
# parsed clusterize results

###### -- libraries -----------------------------------------------------------

suppressMessages(library(SynExtend))

###### -- ad hoc functions ----------------------------------------------------

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

# names for seqs1 and seqs2 need to be unique and traceable
# back to whatever source df and feature variable you're interested in
# when BestMatch is true, returns the best matches for each unique name,
# meaning that a feature in either column can have multiple matches, but they must be
# at least someone's best match
PCOrtho <- function(Seqs1,
                    Seqs2,
                    ClusterCutoff = 0.1,
                    DistanceCutoff = 0.1,
                    BestMatch = FALSE,
                    Verbose = FALSE) {
  
  if (Verbose) {
    TimeStart <- Sys.time()
  }
  
  Seqs3 <- c(Seqs1,
             Seqs2)
  
  Clust01 <- Clusterize(myXStringSet = Seqs3,
                        cutoff = ClusterCutoff,
                        includeTerminalGaps = TRUE,
                        penalizeGapLetterMatches = TRUE,
                        verbose = Verbose)
  Clust02 <- tapply(X = seq_along(Clust01$cluster),
                    INDEX = Clust01$cluster,
                    FUN = c)
  Clust03 <- lapply(X = Clust02,
                    FUN = function(x) {
                      Seqs3[x]
                    })
  Clust04 <- lapply(X = Clust03,
                    FUN = function(x) {
                      if (length(x) >= 2L) {
                        AlignSeqs(myXStringSet = x,
                                  verbose = FALSE)
                      } else {
                        NULL
                      }
                    })
  Clust05 <- lapply(X = Clust04,
                    FUN = function(x) {
                      if (!is.null(x)) {
                        DistanceMatrix(myXStringSet = x,
                                       verbose = FALSE,
                                       includeTerminalGaps = TRUE)
                      }
                    })
  Clust07 <- lapply(X = Clust04,
                    FUN = function(x) {
                      if (!is.null(x)) {
                        DistanceMatrix(myXStringSet = x,
                                       verbose = FALSE,
                                       includeTerminalGaps = FALSE)
                      }
                    })
  
  Clust06 <- vector(mode = "list",
                    length = length(Clust05))
  for (m2 in seq_along(Clust05)) {
    if (is.null(Clust05[[m2]])) {
      next
    } else {
      l1 <- sum(upper.tri(Clust05[[m2]]))
      l2 <- nrow(Clust05[[m2]])
      v1 <- v2 <- vector(mode = "character",
                         length = l1)
      v3 <- v4 <- vector(mode = "numeric",
                         length = l1)
      count <- 1L
      for (m3 in seq_len(l2 - 1L)) {
        for (m4 in (m3 + 1L):l2) {
          v1[count] <- colnames(Clust05[[m2]])[m3]
          v2[count] <- colnames(Clust05[[m2]])[m4]
          v3[count] <- Clust05[[m2]][m3, m4]
          v4[count] <- Clust07[[m2]][m3, m4]
          count <- count + 1L
        }
      }
      Clust06[[m2]] <- data.frame("p1" = v1,
                                  "p2" = v2,
                                  "g_dist" = v3,
                                  "l_dist" = v4)
    }
    
  }
  DF <- do.call(rbind,
                Clust06)
  
  DF <- DF[DF$g_dist <= DistanceCutoff |
             DF$l_dist <= DistanceCutoff, ]
  # force arrangements so that p1 and p2 correspond to seqs1 and seqs2 names
  keep01 <- DF$p1 %in% names(Seqs1) & DF$p2 %in% names(Seqs2)
  keep02 <- DF$p1 %in% names(Seqs2) & DF$p2 %in% names(Seqs1)
  if (any(keep02)) {
    rep1 <- DF$p1[keep02]
    rep2 <- DF$p2[keep02]
    DF$p1[keep02] <- rep2
    DF$p2[keep02] <- rep1
  }
  DF <- DF[keep01 | keep02, ]
  
  if (BestMatch) {
    un01 <- unique(DF$p1)
    un02 <- unique(DF$p2)
    keep01 <- rep(FALSE,
                  nrow(DF))
    for (m3 in seq_along(un01)) {
      z1 <- which.min(DF$g_dist[DF$p1 == un01[m3]])
      keep01[which(DF$p1 == un01[m3])[z1]] <- TRUE
    }
    for (m3 in seq_along(un02)) {
      z1 <- which.min(DF$g_dist[DF$p2 == un02[m3]])
      keep01[which(DF$p2 == un02[m3])[z1]] <- TRUE
    }
    
    DF <- DF[keep01, ]
  }
  
  if (Verbose) {
    TimeEnd <- Sys.time()
    print("Sequences clustered in:")
    print(TimeEnd - TimeStart)
  }
  
  return(DF)
}

###### -- variables -----------------------------------------------------------

ARGS <- commandArgs(trailingOnly = TRUE)

PseudoTypeVector <- c("frameshifted",
                      "internal stop",
                      "partial abbutting assembly gap",
                      "partial in the middle",
                      "missing C-terminus",
                      "missing N-terminus")

###### -- code body -----------------------------------------------------------

print(ARGS)
list.files()

RefAssembly <- readDNAStringSet(filepath = "ref.fna.gz")
RefGC <- gffToDataFrame(GFF = "annot.gff",
                        Verbose = TRUE)
RefGR <- rtracklayer::import("annot.gff")
RefDF <- GRangeToDFrame(GRangeObject = RefGR,
                        Verbose = TRUE)
RefNotes <- NoteCheck(NoteVector = RefDF$Note,
                      CheckVector = PseudoTypeVector)

CurrentAssembly <- readDNAStringSet(filepath = ARGS[2L])
CurrentGC <- gffToDataFrame(GFF = ARGS[3L],
                            Verbose = TRUE)
CurrentGR <- rtracklayer::import(ARGS[3L])
CurrentDF <- GRangeToDFrame(GRangeObject = CurrentGR,
                            Verbose = TRUE)
CurrentNotes <- NoteCheck(NoteVector = CurrentDF$Note,
                          CheckVector = PseudoTypeVector)

head(CurrentGC)
head(RefGC)
ls()

Seqs1 <- ExtractBy(x = CurrentGC,
                   y = CurrentAssembly,
                   Verbose = TRUE)
names(Seqs1) <- paste(names(Seqs1),
                      "A",
                      sep = "_")
Seqs2 <- ExtractBy(x = RefGC,
                   y = RefAssembly,
                   Verbose = TRUE)
names(Seqs2) <- paste(names(Seqs2),
                      "B",
                      sep = "_")

Cl1 <- PCOrtho(Seqs1 = Seqs1, # the reassembly
               Seqs2 = Seqs2, # the reference
               BestMatch = TRUE,
               Verbose = TRUE)

ph1 <- paste(CurrentDF$Name[grepl(pattern = "internal stop",
                                  x = CurrentDF$Note)],
             "A",
             sep = "_")
ph2 <- paste(CurrentDF$Name[grepl(pattern = "frameshifted",
                                  x = CurrentDF$Note)],
             "A",
             sep = "_")
ph3 <- paste(RefDF$Name[grepl(pattern = "internal stop",
                              x = RefDF$Note)],
             "B",
             sep = "_")
ph4 <- paste(RefDF$Name[grepl(pattern = "frameshifted",
                              x = RefDF$Note)],
             "B",
             sep = "_")

is1 <- which(Cl1$p1 %in% ph1 &
               !(Cl1$p2 %in% ph3))
is2 <- which(!(Cl1$p1 %in% ph1) &
               Cl1$p2 %in% ph3)
fs1 <- which(Cl1$p1 %in% ph2 &
               !(Cl1$p2 %in% ph4))
fs2 <- which(!(Cl1$p1 %in% ph2) &
               Cl1$p2 %in% ph4)


w1 <- sort(width(CurrentAssembly), decreasing = TRUE)
TotalNucs <- sum(w1)
HalfNuc <- ceiling(TotalNucs / 2)
AllWidths <- cumsum(w1)
w2 <- which(AllWidths >= HalfNuc)[1L]
N50 <- w1[w2]

# frameshifts
res <- list("rel_fs" = unname(CurrentNotes[1]) / (sum(width(CurrentAssembly)) / 1000000),
            "rel_is" = unname(CurrentNotes[2]) / (sum(width(CurrentAssembly)) / 1000000),
            "fs1" = length(fs1) / nrow(Cl1),
            "is1" = length(is1) / nrow(Cl1),
            "fs2" = length(fs2) / nrow(Cl1),
            "is2" = length(is2) / nrow(Cl1),
            "n50" = N50,
            "total_length" = TotalNucs,
            "tc" = length(CurrentAssembly))

save(res,
     Cl1,
     is1,
     is2,
     fs1,
     fs2,
     CurrentGC,
     CurrentDF,
     file = paste0("Result",
                   formatC(x = as.integer(ARGS[1L]),
                           format = "d",
                           flag = 0,
                           width = 5L),
                   ".RData"),
     compress = "xz")

