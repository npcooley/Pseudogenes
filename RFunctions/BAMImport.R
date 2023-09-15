###### -- produce alignments from imported BAM --------------------------------

PairwiseBAMImport <- function(BAMImport,
                              Reference) {
  
  # tasks
  # match the reference names
  # assign regions to pull from references based on CIGAR stings
  # pull references in as few extracts as possible
  # split out to the sequences present in the BAM file which have already been
  # reverse complemented if necessary
  
  # do not assume suppressed records from unaligned sequences
  MatchedRecord <- !is.na(BAMImport$qwidth)
  RefNames <- unlist(regmatches(x = names(Reference),
                                m = gregexpr(pattern = "^[^ ]+",
                                             text = names(Reference))))
  BamRefNames <- as.character(BAMImport$rname)[MatchedRecord]
  BamRefCIGAR <- as.character(BAMImport$cigar)[MatchedRecord]
  BamRefPositions <- BAMImport$pos[MatchedRecord]
  BamRefWidth <- BAMImport$qwidth[MatchedRecord]
  BamReadNames <- as.character(BAMImport$qname)[MatchedRecord]
  SeqsPresent <- BAMImport$seq[MatchedRecord]
  # parse CIGAR strings for reference changes only
  # CIGAR is always zero base ? but that shouldn't matter for these parsings
  INS <- regmatches(x = BamRefCIGAR,
                    m = gregexpr(pattern = "([0-9]+)(?=I)",
                                 text = BamRefCIGAR,
                                 perl = TRUE))
  DEL <- regmatches(x = BamRefCIGAR,
                    m = gregexpr(pattern = "([0-9]+)(?=D)",
                                 text = BamRefCIGAR,
                                 perl = TRUE))
  SKIP <- regmatches(x = BamRefCIGAR,
                     m = gregexpr(pattern = "([0-9]+)(?=N)",
                                  text = BamRefCIGAR,
                                  perl = TRUE))
  for (m1 in seq_along(INS)) {
    if (length(INS[[m1]]) > 0) {
      INS[[m1]] <- sum(as.integer(INS[[m1]]))
    } else {
      INS[[m1]] <- 0L
    }
    if (length(DEL[[m1]]) > 0) {
      DEL[[m1]] <- sum(as.integer(DEL[[m1]]))
    } else {
      DEL[[m1]] <- 0L
    }
    if (length(SKIP[[m1]]) > 0) {
      SKIP[[m1]] <- sum(as.integer(SKIP[[m1]]))
    } else {
      SKIP[[m1]] <- 0L
    }
  }
  
  INS <- unlist(INS)
  DEL <- unlist(DEL)
  SKIP <- unlist(SKIP)
  
  BamRefWidth <- BamRefWidth + INS + SKIP - DEL
  
  u1 <- unique(BamRefNames)
  i1 <- match(x = BamRefNames,
              table = u1)
  InitialExtraction <- PairPartner <- vector(mode = "list",
                                             length = length(u1))
  # return(list(Reference,
  #             RefNames,
  #             u1,
  #             i1,
  #             BamRefPositions,
  #             BamRefWidth))
  for (m1 in seq_along(InitialExtraction)) {
    w1 <- i1 == m1
    InitialExtraction[[m1]] <- extractAt(x = Reference[RefNames == u1[m1]][[1L]],
                                         at = IRanges(start = BamRefPositions[w1],
                                                      width = BamRefWidth[w1]))
    names(InitialExtraction[[m1]]) <- rep(u1[m1],
                                          length(InitialExtraction[[m1]]))
    
    PairPartner[[m1]] <- SeqsPresent[w1]
    names(PairPartner[[m1]]) <- BamReadNames[w1]
  }
  
  InitialExtraction <- do.call(c,
                               InitialExtraction)
  PairPartner <- do.call(c,
                         PairPartner)
  
  PairedSeqs <- vector(mode = "list",
                       length = length(InitialExtraction))
  for (m1 in seq_along(PairedSeqs)) {
    PairedSeqs[[m1]] <- DNAStringSet(list(PairPartner[[m1]],
                                          InitialExtraction[[m1]]))
    names(PairedSeqs[[m1]]) <- c(names(PairPartner)[m1],
                                 names(InitialExtraction)[m1])
  }
  
  return(PairedSeqs)
  
}


