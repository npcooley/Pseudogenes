###### -- ask for syntenic blocks from the reads ------------------------------

SynMappingPipeline <- function(Reads,
                               Reference,
                               GeneCalls) {
  
  temp01 <- tempfile()
  # add reads
  Seqs2DB(seqs = Reads,
          type = "XStringSet",
          dbFile = temp01,
          identifier = "Reads",
          verbose = FALSE)
  Seqs2DB(seqs = Reference,
          type = "XStringSet",
          dbFile = temp01,
          identifier = "Reference",
          verbose = FALSE)
  
  res <- FindSynteny(dbFile = temp01,
                     useFrames = FALSE,
                     verbose = FALSE)
  # assume standard NCBI naming conventions
  # technically this will fail if a gene extends past the end of it's contig,
  # but that can be adjusted if necessary
  RefNames <- unlist(regmatches(x = names(Reference),
                                m = gregexpr(pattern = "^[^ ]+",
                                             text = names(Reference))))
  u1 <- unique(GeneCalls$Contig)
  seqs1 <- vector(mode = "list",
                  length = length(u1))
  for (m1 in seq_along(u1)) {
    # access reference string
    # print(m1)
    w1 <- which(RefNames == u1[m1])
    # grab IRangesList for the current contig
    z1 <- unname(GeneCalls$Range[GeneCalls$Contig == u1[m1]])
    # create a collapse map
    z2 <- lengths(z1)
    # unlist IRangesList to an IRanges object
    z1 <- unlist(z1,
                 recursive = FALSE)
    seqs1[[m1]] <- extractAt(x = Reference[w1][[1L]],
                             at = z1)
    # return(seqs1)
    collapsecount <- 0L
    w2 <- which(z2 > 1)
    # print(m1)
    # print(w2)
    if (length(w2) > 0L) {
      removevector <- vector(mode = "integer",
                             length = sum(z2[w2]) - length(w2))
      for (m2 in w2) {
        # print(m2)
        seqs1[[m1]][[m2 + collapsecount]] <- unlist(seqs1[[m1]][m2:(m2 + z2[m2] - 1L) + collapsecount])
        removevector[(collapsecount + 1L):(collapsecount + z2[m2] - 1L)] <- (m2 + 1L):(m2 + z2[m2] - 1L) + collapsecount
        collapsecount <- collapsecount + z2[m2] - 1L
      }
      seqs1[[m1]][removevector] <- NULL
    }
    # return(seqs1)
    names(seqs1[[m1]]) <- GeneCalls$ID[GeneCalls$Contig == u1[m1]]
  }
  
  seqs1 <- do.call(c,
                   seqs1)
  # return a list of stringsets of sanger reads and sequences
  w1 <- res[[2, 1]][, 1L]
  w2 <- res[[2, 1]][, 2L]
  s1 <- res[[2 ,1]][, 6L]
  s2 <- res[[2, 1]][, 8L]
  s3 <- res[[2, 1]][, 3L]
  
  # names of assembly by synteny object
  Ref2 <- unlist(regmatches(x = names(res[[2, 2]]),
                            m = gregexpr(pattern = "^[^ ]+",
                                         text = names(res[[2, 2]]))))
  # pull the sequences based on the query positions that
  w3 <- w4 <- w5 <- rep(NA_integer_,
                        length(w1))
  for (m1 in seq_along(w3)) {
    w3[m1] <- which(GeneCalls$Contig == Ref2[w2[m1]] &
                      GeneCalls$Start <= s1[m1] &
                      GeneCalls$Stop >= s2[m1])
  }
  # set w3 as the row of the GeneCalls object --- relatable to ID
  if (any(is.na(w3))) {
    w4 <- w1[!is.na(w3)]
    w5 <- s3[!is.na(w3)]
    w3 <- w3[!is.na(w3)]
  } else {
    w4 <- w1
    w5 <- s3
  }
  
  seqs2 <- seqs1[match(x = GeneCalls$ID[w3],
                       table = names(seqs1))]
  seqs3 <- Reads[w4]
  
  seqs4 <- vector(mode = "list",
                  length = length(s3))
  for (m1 in seq_along(w5)) {
    # if the synteny object
    # AND the genecalls tell you to flip
    # flip both?
    # 
    if (w5[m1] == 1L) {
      seqs4[[m1]] <- DNAStringSet(list(seqs2[[m1]],
                                       reverseComplement(seqs3[[m1]])))
    } else {
      seqs4[[m1]] <- DNAStringSet(list(seqs2[[m1]],
                                       seqs3[[m1]]))
    }
    names(seqs4[[m1]]) <- c(names(seqs2)[m1],
                            names(seqs3)[m1])
  }
  
  return(list(res,
              seqs4))
}

