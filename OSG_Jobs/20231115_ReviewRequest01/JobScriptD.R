###### -- run a battery of analyses on a pair of assemblies -------------------

###### -- libraries -----------------------------------------------------------

suppressMessages(library(SynExtend))

###### -- adhoc functions -----------------------------------------------------

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

GetN50 <- function(Assembly) {
  w1 <- sort(width(Assembly), decreasing = TRUE)
  TotalNucs <- sum(w1)
  HalfNuc <- ceiling(TotalNucs / 2)
  AllWidths <- cumsum(w1)
  w2 <- which(AllWidths >= HalfNuc)[1L]
  N50 <- w1[w2]
  return(N50)
}

###### -- arguments -----------------------------------------------------------

ARGS <- commandArgs(trailingOnly = TRUE)
print(ARGS)

PseudoTypeVector <- c("frameshifted",
                      "internal stop",
                      "partial abbutting assembly gap",
                      "partial in the middle",
                      "missing C-terminus",
                      "missing N-terminus")

PersistentID <- ARGS[1]
GFFs <- ARGS[c(3,5)]
Assemblies <- ARGS[c(2,4)]

###### -- code body part 1 ----------------------------------------------------
# load in GFFs and parse, get specific pseudogene notes

seqs01 <- seqs02 <- gc01 <- gc02 <- fs_ids <- is_ids <- N50List <- type_counts <- vector(mode = "list",
                                                                                         length = length(GFFs))
coding_lengths <- vector(mode = "list",
                         length = length(GFFs))

DBPATH <- tempfile()

for (m1 in seq_along(GFFs)) {
  seqs01[[m1]] <- readDNAStringSet(filepath = Assemblies[m1])
  Seqs2DB(seqs = seqs01[[m1]],
          dbFile = DBPATH,
          type = "XStringSet",
          identifier = as.character(m1),
          verbose = TRUE)
  y <- gffToDataFrame(GFF = GFFs[m1],
                      Verbose = TRUE)
  x <- rtracklayer::import(GFFs[m1])
  gc02[[m1]] <- GRangeToDFrame(GRangeObject = x,
                               Verbose = TRUE)
  
  gc01[[m1]] <- cbind(y,
                      "Note" = gc02[[m1]]$Note[match(x = gc02[[m1]]$Name,
                                                     table = y$ID)])
  seqs02[[m1]] <- ExtractBy(x = gc01[[m1]],
                            y = seqs01[[m1]],
                            Verbose = TRUE)
  names(seqs02[[m1]]) <- paste(names(seqs02[[m1]]),
                               LETTERS[m1],
                               sep = "_")
  fs_ids[[m1]] <- paste(gc01[[m1]]$ID[grepl(pattern = "frameshifted",
                                            x = gc01[[m1]]$Note)],
                        LETTERS[m1],
                        sep = "_")
  is_ids[[m1]] <- paste(gc01[[m1]]$ID[grepl(pattern = "internal stop",
                                            x = gc01[[m1]]$Note)],
                        LETTERS[m1],
                        sep = "_")
  
  N50List[[m1]] <- GetN50(Assembly = seqs01[[m1]])
  type_counts[[m1]] <- c(sum(gc01[[m1]]$Coding),
                         sum(gc01[[m1]]$Coding & gc01[[m1]]$Type != "pseudogene"),
                         sum(grepl(pattern = "frameshifted",
                                   x = gc01[[m1]]$Note)),
                         sum(grepl(pattern = "internal stop",
                                   x = gc01[[m1]]$Note)),
                         sum(gc01[[m1]]$Type == "pseudogene" &
                               !grepl(pattern = "frameshifted",
                                      x = gc01[[m1]]$Note) &
                               !grepl(pattern = "internal stop",
                                      x = gc01[[m1]]$Note)),
                         sum(width(seqs01[[m1]])))
  
  coding_lengths[[m1]] <- (gc01[[m1]]$Stop[gc01[[m1]]$Coding] - gc01[[m1]]$Start[gc01[[m1]]$Coding]) + 1L
  
}
names(gc01) <- seq(length(gc01))

###### -- code body part 2 ----------------------------------------------------
# run ID clusters

Cl1 <- PCOrtho(Seqs1 = seqs02[[1]], # the reassembly
               Seqs2 = seqs02[[2]], # the reference
               BestMatch = TRUE,
               Verbose = TRUE)
# internal stops in a but not b
is1 <- which(Cl1$p1 %in% is_ids[[1L]] &
               !(Cl1$p2 %in% is_ids[[2L]]))
# internal stops in b but not a
is2 <- which(!(Cl1$p1 %in% is_ids[[1L]]) &
               Cl1$p2 %in% is_ids[[2L]])
# shared internal stops == found reference
is3 <- which(Cl1$p1 %in% is_ids[[1L]] &
               Cl1$p2 %in% is_ids[[2L]])

# frameshifts in a but not in b
fs1 <- which(Cl1$p1 %in% fs_ids[[1L]] &
               !(Cl1$p2 %in% fs_ids[[2L]]))
# frameshifts in b but not in a
fs2 <- which(!(Cl1$p1 %in% fs_ids[[1L]]) &
               Cl1$p2 %in% fs_ids[[2L]])
# shared frameshifts == found reference
fs3 <- which(Cl1$p1 %in% fs_ids[[1L]] &
               Cl1$p2 %in% fs_ids[[2L]])

###### -- code body part 3 ----------------------------------------------------
# build the synteny map

syn <- FindSynteny(dbFile = DBPATH)
ali <- AlignSynteny(synteny = syn,
                    dbFile = DBPATH)

###### -- code body part 4 ----------------------------------------------------
# get ANI

# write out all the non-pseudogenes
writeXStringSet(x = seqs02[[1]][gc01[[1]]$Type != "pseudogene" & gc01[[1]]$Coding],
                filepath = "seqs1.fna")
writeXStringSet(x = seqs02[[2]][gc01[[2]]$Type != "pseudogene" & gc01[[2]]$Coding],
                filepath = "seqs2.fna")

ANICall <- "ANIcalculator -genome1fna seqs1.fna -genome2fna seqs2.fna -outfile RES.txt -outdir res"

system(command = ANICall,
       intern = FALSE,
       ignore.stderr = FALSE)

ani_res <- readLines("RES.txt")
file.remove(c("seqs1.fna",
              "seqs2.fna"))

###### -- gene length tests ---------------------------------------------------

z1 <- ks.test(x = coding_lengths[[1]],
              y = coding_lengths[[2]])

###### -- save data -----------------------------------------------------------

save(ani_res,
     syn,
     ali,
     Cl1,
     is1,
     is2,
     is3,
     fs1,
     fs2,
     fs3,
     gc01,
     N50List,
     type_counts,
     coding_lengths,
     z1,
     file = paste0("ResultD",
                   formatC(x = as.integer(PersistentID),
                           flag = 0,
                           width = 6,
                           format = "d"),
                   ".RData"),
     compress = "xz")

