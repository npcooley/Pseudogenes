###### -- describe pseudogenes ------------------------------------------------
# given 4 arguments
# perform a pairwise comparison
# and return a uniquely named data file

# goal is to return data objects that allow us to interrogate whether
# pseudogenizations are shared or not relative to the ANI

###### -- libraries -----------------------------------------------------------

suppressMessages(library(SynExtend))

###### -- data ----------------------------------------------------------------

# sent along to job by condor, will land in the working directory
load(file = "SearchResults.RData",
     verbose = TRUE)

# debugging the ANI calculator -- it has odd permissions which make it slightly
# less accessible on the node than when just spinning docker up locally
# system(command = "echo $PATH")
# system(command = "pwd")
# system(command = "ls -lh")
# system(command = "ls -lh /")
# system(command = "ANIcalculator -help")
# system(command = "ls -lh ANIcalculator_v1")

###### -- arguments -----------------------------------------------------------

ARGS <- commandArgs(trailingOnly = TRUE)
SELECT01 <- ARGS[1L] # genus / taxID selector
SELECT02 <- ARGS[2L] # pair partner 1
SELECT03 <- ARGS[3L] # pair partner 2
SELECT04 <- ARGS[4L] # unique queue ID

###### -- functions -----------------------------------------------------------

# adhoc function for pseudogene interrogation 
GRangeToDFrame <- function(GRangeObject,
                           FeaturesToCollect = c("gene",
                                                 "pseudogene")) {
  # search until no more searches are necessary
  s1 <- as.character(GRangeObject$type)
  s2 <- GRangeObject$gene
  s3 <- GRangeObject@strand
  s4 <- GRangeObject@ranges
  s5 <- as.character(GRangeObject@seqnames)
  s6 <- GRangeObject$Note
  s7 <- GRangeObject$Parent
  s8 <- GRangeObject$ID
  s9 <- !is.na(s2)
  
  TOTAL <- sum(table(s1[s1 %in% FeaturesToCollect]))
  # print(TOTAL)
  CONTINUE <- TRUE
  KEEP <- vector(mode = "logical",
                 length = length(s1))
  START <- STOP <- vector(mode = "integer",
                          length = length(s1))
  NOTE <- CONTIG <- TYPE <- ID <- GENE <- vector(mode = "character",
                                                 length = length(s1))
  COUNT <- 1L
  FOUNDFEATURES <- 0L
  pBar <- txtProgressBar(style = 1L)
  while (CONTINUE) {
    # is the line a line to evaluate
    # check its children
    if (s1[COUNT] %in% FeaturesToCollect) {
      if (s1[COUNT] == "pseudogene") {
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
              ph1 <- "pseudofeature with absent note"
            }
          }
          NOTE[COUNT] <- ph1
        } else {
          # feature has no children, what to do here?
          NOTE[COUNT] <- "child lines absent"
        }
      } else {
        # ph1 <- "normal feature"
        NOTE[COUNT] <- "normal feature"
      }
      START[COUNT] <- s4@start[COUNT]
      STOP[COUNT] <- s4@start[COUNT] + s4@width[COUNT] - 1L
      CONTIG[COUNT] <- s5[COUNT]
      TYPE[COUNT] <- s1[COUNT]
      ID[COUNT] <- s8[COUNT]
      
      if (s9[COUNT]) {
        GENE[COUNT] <- s2[COUNT]
      } else {
        GENE[COUNT] <- ""
      }
      KEEP[COUNT] <- TRUE
      FOUNDFEATURES <- FOUNDFEATURES + 1L
      
    } # end if s2 is a feature to collect or not
    
    if (FOUNDFEATURES >= TOTAL) {
      CONTINUE <- FALSE
    } else {
      COUNT <- COUNT + 1L
    }
    
    setTxtProgressBar(pb = pBar,
                      value = FOUNDFEATURES / TOTAL)
  }
  close(pBar)
  cat("\n")
  # return(list(START,
  #             STOP,
  #             TYPE,
  #             CONTIG,
  #             ID,
  #             NOTE,
  #             KEEP))
  res <- DataFrame("Start" = START[KEEP],
                   "Stop" = STOP[KEEP],
                   "Type" = TYPE[KEEP],
                   "Contig" = CONTIG[KEEP],
                   "ID" = ID[KEEP],
                   "Gene" = GENE[KEEP],
                   "Note" = NOTE[KEEP])
  return(res)
}

###### -- code body part 1: data collection -----------------------------------

# load in everything that matters:
DBPATH01 <- tempfile() # CDSs
DBPATH02 <- tempfile() # genomes

w1 <- c(SELECT02,
        SELECT03)
w2 <- as.integer(w1)
w3 <- EntrezResults$FTP[w2]

system(command = "mkdir seqdir")

GC01 <- GC02 <- GC03 <- vector(mode = "list",
                               length = length(w1))

for (m1 in seq_along(w1)) {
  
  add1 <- paste0(w3[m1],
                 "/",
                 strsplit(w3[m1],
                          split = "/",
                          fixed = TRUE)[[1]][10],
                 "_cds_from_genomic.fna.gz")
  add2 <- paste0(w3[m1],
                 "/",
                 strsplit(w3[m1],
                          split = "/",
                          fixed = TRUE)[[1]][10],
                 "_genomic.gff.gz")
  add3 <- paste0(w3[m1],
                 "/",
                 strsplit(w3[m1],
                          split = "/",
                          fixed = TRUE)[[1]][10],
                 "_genomic.fna.gz")
  
  # the CDSs
  RETRY01 <- try(Seqs2DB(seqs = add1,
                         type = "FASTA",
                         dbFile = DBPATH01,
                         identifier = w1[m1],
                         verbose = TRUE))
  RETRY02 <- 1L
  while (is(object = RETRY01,
            class2 = "try-error")) {
    Sys.sleep(5L)
    RETRY01 <- try(Seqs2DB(seqs = add1,
                           type = "FASTA",
                           dbFile = DBPATH01,
                           identifier = w1[m1],
                           verbose = TRUE))
    RETRY02 <- RETRY02 + 1L
    
    if (RETRY02 >= 5L) {
      print("Huh")
      quit(save = "no")
    }
  }
  print(paste0("CDSs for ",
               w1[m1],
               " collected."))
  
  DB2Seqs(file = paste0("seqdir/seqs",
                        w1[m1],
                        ".fna"),
          dbFile = DBPATH01,
          identifier = w1[m1],
          verbose = TRUE)
  
  # the full genomic seqs
  RETRY01 <- try(Seqs2DB(seqs = add3,
                         type = "FASTA",
                         dbFile = DBPATH02,
                         identifier = w1[m1],
                         verbose = TRUE))
  RETRY02 <- 1L
  while (is(object = RETRY01,
            class2 = "try-error")) {
    Sys.sleep(5L)
    RETRY01 <- try(Seqs2DB(seqs = add3,
                           type = "FASTA",
                           dbFile = DBPATH02,
                           identifier = w1[m1],
                           verbose = TRUE))
    RETRY02 <- RETRY02 + 1L
    
    if (RETRY02 >= 5L) {
      print("Huh")
      quit(save = "no")
    }
  }
  print(paste0("genomic seqs for ",
               w1[m1],
               " collected."))
  
  # the gffs
  RETRY01 <- try(gffToDataFrame(GFF = add2,
                                Verbose = TRUE))
  RETRY02 <- 1L
  while (is(object = RETRY01,
            class2 = "try-error")) {
    Sys.sleep(5L)
    RETRY01 <- try(gffToDataFrame(GFF = add2,
                                  Verbose = TRUE))
    RETRY02 <- RETRY02 + 1L
    
    if (RETRY02 >= 5L) {
      print("Huh")
      quit(save = "no")
    }
  }
  GC01[[m1]] <- RETRY01
  print(paste0("GFFs part 1 for ",
               w1[m1],
               " collected."))
  
  # the gffs again
  RETRY01 <- try(rtracklayer::import(add2))
  RETRY02 <- 1L
  while (is(object = RETRY01,
            class2 = "try-error")) {
    Sys.sleep(5L)
    RETRY01 <- try(rtracklayer::import(add2))
    RETRY02 <- RETRY02 + 1L
    
    if (RETRY02 >= 5L) {
      print("Huh")
      quit(save = "no")
    }
  }
  GC02[[m1]] <- RETRY01
  RETRY03 <- GRangeToDFrame(GRangeObject = RETRY01)
  GC03[[m1]] <- RETRY03
  print(paste0("GFFs part 2 for ",
               w1[m1],
               " collected."))
  
}
names(GC01) <- names(GC02) <- names(GC03) <- w1

GC04 <- vector(mode = "list",
               length = length(GC01))
for (m1 in seq_along(GC01)) {
  GC04[[m1]] <- cbind(GC01[[m1]],
                      "FeatureNote" = GC03[[m1]][match(x = GC01[[m1]]$ID,
                                                       table = GC03[[m1]]$ID), "Note"])
}
names(GC04) <- w1

# z1 <- gffToDataFrame(GFF = x,
#                      Verbose = FALSE)
# z2 <- rtracklayer::import(x)
# z3 <- GRangeToDFrame(GRangeObject = z2,
#                      FeaturesToCollect = c("gene",
#                                            "pseudogene"))
# genecalls <- cbind(z1,
#                    "Note" = z3$Note[match(x = z1$ID,
#                                           table = z3$ID)])

###### -- code body part 2: JGI calculator ------------------------------------

CBT01 <- Sys.time()

ANICALL <- paste0("ANIcalculator",
                  " -genome1fna ",
                  paste0("seqdir/seqs",
                         w1[1],
                         ".fna"),
                  " -genome2fna ",
                  paste0("seqdir/seqs",
                         w1[2],
                         ".fna"),
                  " -outfile CURRENTRESULT.txt")

system(command = ANICALL,
       intern = FALSE)

ANIRes <- readLines("CURRENTRESULT.txt")
ANIRes <- as.numeric(strsplit(x = ANIRes,
                             split = "\t",
                             fixed = TRUE)[[2]][3:6])

system(command = "rm CURRENTRESULT.txt ANIcalculator.log")

CBT02 <- Sys.time()
ANITime <- CBT02 - CBT01

###### -- code body part 3: IdClusters / Clusterize ---------------------------

CBT01 <- Sys.time()

genome1 <- SearchDB(dbFile = DBPATH02,
                    identifier = w1[1L],
                    nameBy = "description")

genes1 <- ExtractBy(x = GC01[[w1[1L]]],
                    y = genome1)

genome2 <- SearchDB(dbFile = DBPATH02,
                    identifier = w1[2L],
                    nameBy = "description")

genes2 <- ExtractBy(x = GC01[[w1[2L]]],
                    y = genome2)
# reshape to only coding genes
geneskey1 <- GC01[[w1[1L]]]$Type
j1 <- grepl(pattern = "internal stop",
            x = GC04[[w1[1L]]]$FeatureNote)
j2 <- grepl(pattern = "frameshifted",
            x = GC04[[w1[1L]]]$FeatureNote)
geneskey1[geneskey1 == "pseudogene"] <- "other"
geneskey1[j1 & j2] <- "both"
geneskey1[j1 & !j2] <- "internal_stop"
geneskey1[!j1 & j2] <- "frameshift"

geneskey2 <- GC01[[w1[2L]]]$Type
j1 <- grepl(pattern = "internal stop",
            x = GC04[[w1[2L]]]$FeatureNote)
j2 <- grepl(pattern = "frameshifted",
            x = GC04[[w1[2L]]]$FeatureNote)
geneskey2[geneskey2 == "pseudogene"] <- "other"
geneskey2[j1 & j2] <- "both"
geneskey2[j1 & !j2] <- "internal_stop"
geneskey2[!j1 & j2] <- "frameshift"

genes1 <- genes1[GC01[[w1[1L]]]$Coding]
genes2 <- genes2[GC01[[w1[2L]]]$Coding]
geneskey1 <- geneskey1[GC01[[w1[1L]]]$Coding]
geneskey2 <- geneskey2[GC01[[w1[2L]]]$Coding]
# rename
gk01 <- paste("genome_a", geneskey1, which(GC01[[w1[1L]]]$Coding))
gk02 <- paste("genome_b", geneskey2, which(GC01[[w1[2L]]]$Coding))
names(genes1) <- gk01
names(genes2) <- gk02

genes <- c(genes1, genes2)

clust <- Clusterize(genes,
                    cutoff = 0.1)
clust2 <- tapply(X = rownames(clust),
                 INDEX = clust$cluster,
                 FUN = c)

clust3 <- clust2[lengths(clust2) > 1]
prefixes <- c("genome_a",
              "genome_b")
clust4 <- vector(mode = "list",
                 length = length(w2))
for (m1 in seq_along(clust4)) {
  clust4[[m1]] <- lapply(X = clust3,
                         FUN = function(x) {
                           x[grepl(pattern = prefixes[m1],
                                   x = x)]
                         })
}
clust5 <- clust4[[1]][lengths(clust4[[1]]) > 0 & lengths(clust4[[2]]) > 0]
clust6 <- clust4[[2]][lengths(clust4[[1]]) > 0 & lengths(clust4[[2]]) > 0]
clust7 <- clust3[lengths(clust4[[1]]) > 0 & lengths(clust4[[2]]) > 0]

g1 <- unlist(clust5)
g2 <- unlist(clust6)

featuretypes <- c("gene",
                  "both",
                  "internal_stop",
                  "frameshift",
                  "other")

l1 <- l2 <- vector(mode = "list",
                   length = length(featuretypes))
for (m1 in seq_along(featuretypes)) {
  l1[[m1]] <- grepl(pattern = featuretypes[m1],
                    x = g1)
  l1[[m1]] <- split(x = l1[[m1]],
                    f = rep(seq(length(clust5)),
                            lengths(clust5)))
  l1[[m1]] <- sapply(X = l1[[m1]],
                     FUN = function(x) {
                       any(x)
                     })
  l2[[m1]] <- grepl(pattern = featuretypes[m1],
                    x = g2)
  l2[[m1]] <- split(x = l2[[m1]],
                    f = rep(seq(length(clust6)),
                            lengths(clust6)))
  l2[[m1]] <- sapply(X = l2[[m1]],
                     FUN = function(x) {
                       any(x)
                     })
}

l3 <- do.call(cbind,
              c(l1, l2))
colnames(l3) <- paste(rep(c("a", "b"),
                          c(5, 5)),
                      featuretypes,
                      sep = "_")

# genome is represented solely by features with this specific pseudogenization
# in this homolog group -- if it contains features with the other pseudogenization
# that doesn't matter from this focal point
# side a internal stops
side_a_is <- (!l3[, 1L] & !l3[, 5L]) & (l3[, 2] | l3[, 3])
# side b internal stops
side_b_is <- (!l3[, 6L] & !l3[, 10L]) & (l3[, 7] | l3[, 8])
# side a frameshifts
side_a_fs <- (!l3[, 1L] & !l3[, 5L]) & (l3[, 2] | l3[, 4])
# side b frameshifts
side_b_fs <- (!l3[, 6L] & !l3[, 10L]) & (l3[, 7] | l3[, 9])
# genome is represented solely by features without this specific pseudogenization
# in this homolog group == both must be false
# it doesn't matter if the set contains pseudos of another type
# side a non-internal stops
side_a_non_is <- (!l3[, 2] & !l3[, 3])
# side b non-internal stops
side_b_non_is <- (!l3[, 7] & !l3[, 8])
# side a non-frameshifts
side_a_non_fs <- (!l3[, 2] & !l3[, 4])
# side b non-frameshifts
side_b_non_fs <- (!l3[, 7] & !l3[, 9])

incongruent_is <- (side_a_is & side_b_non_is) | (side_a_non_is & side_b_is)
incongruent_fs <- (side_a_fs & side_b_non_fs) | (side_a_non_fs & side_b_fs)
congruent_is <- side_a_is & side_b_is
congruent_fs <- side_a_fs & side_b_fs

CBT02 <- Sys.time()
IDCTime <- CBT02 - CBT01

###### -- code body part 4: synextend -----------------------------------------

# CBT01 <- Sys.time()
# 
# syn <- FindSynteny(dbFile = DBPATH02,
#                    verbose = TRUE)
# 
# links <- NucleotideOverlap(SyntenyObject = syn,
#                            GeneCalls = GC04,
#                            Verbose = TRUE)
# p01 <- PairSummaries(SyntenyLinks = links,
#                      DBPATH = DBPATH02,
#                      PIDs = TRUE,
#                      Score = TRUE,
#                      Verbose = TRUE)
# 
# CBT02 <- Sys.time()
# SYNTime <- CBT02 - CBT01

###### -- save data -----------------------------------------------------------

save(clust7,
     l3,
     congruent_is,
     congruent_fs,
     incongruent_is,
     incongruent_fs,
     geneskey1,
     geneskey2,
     ANIRes,
     w2,
     file = paste0("Result_",
                   SELECT01,
                   "_",
                   SELECT04,
                   ".RData"),
     compress = "xz")

