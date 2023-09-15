###### -- redo pairwise comparison with gene type info ------------------------

suppressMessages(library(SynExtend))

path01 <- "~/Data/20220920_SetComparison"

files01 <- list.files(path = path01)

load(file = "~/Repos/20220920_SetComparison/SearchResults.RData",
     verbose = TRUE)
jm01 <- read.table(file = "~/Repos/20220920_SetComparison/JobMap01.txt")

###### -- adhoc function for parsing GFFs -------------------------------------

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

# build out the list of all possible gffs that need to get parsed

u1 <- unique(unlist(jm01[,2:3]))

gff_locs <- paste0(EntrezResults$FTP[u1],
                   "/",
                   sapply(X = EntrezResults$FTP[u1],
                          FUN = function(x) {
                            strsplit(x = x,
                                     split = "/",
                                     fixed = TRUE)[[1]][10]
                          }),
                   "_genomic.gff.gz")

gff_loc2 <- split(x = gff_locs,
                  f = 1:100)
u2 <- split(u1,
            f = 1:100)

gffs1 <- vector(mode = "list",
                length = length(gff_loc2))

for (m1 in seq_along(gff_loc2)) {
  # it's better to not try and store all these in memory at the same time ...
  # MCTIMESTART <- Sys.time()
  # gffs1[[m1]] <- mclapply(X = gff_loc2[[m1]],
  #                         FUN = function(x) {
  #                           z1 <- gffToDataFrame(GFF = x,
  #                                                Verbose = FALSE)
  #                           # return(z1)
  #                           z2 <- rtracklayer::import(x)
  #                           z3 <- GRangeToDFrame(GRangeObject = z2,
  #                                                FeaturesToCollect = c("gene",
  #                                                                      "pseudogene"))
  #                           z4 <- cbind(z1,
  #                                       "Note" = z3$Note[match(x = z1$ID,
  #                                                              table = z3$ID)])
  #                           return(z4)
  #                         },
  #                         mc.cores = 10)
  # MCTIMEEND <- Sys.time()
  # print(MCTIMEEND - MCTIMESTART)
  # print(m1)
  
  MCTIMESTART <- Sys.time()
  ph1 <- mcmapply(FUN = function(x, y) {
    z1 <- gffToDataFrame(GFF = x,
                         Verbose = FALSE)
    z2 <- rtracklayer::import(x)
    z3 <- GRangeToDFrame(GRangeObject = z2,
                         FeaturesToCollect = c("gene",
                                               "pseudogene"))
    genecalls <- cbind(z1,
                       "Note" = z3$Note[match(x = z1$ID,
                                              table = z3$ID)])
    save(genecalls,
         file = paste0("~/Data/20230811_GeneCalls/GeneCalls",
                       formatC(x = y,
                               width = 7,
                               flag = 0,
                               format = "d"),
                       ".RData"),
         compress = "xz")
    return(NULL)
  },
  mc.cores = 10,
  x = gff_loc2[[m1]],
  y = u2[[m1]],
  USE.NAMES = FALSE)
  
  MCTIMEEND <- Sys.time()
  print(MCTIMEEND - MCTIMESTART)
  print(m1)
}


pBar <- txtProgressBar(style = 1)
PBAR <- length(files01)
res01 <- vector(mode = "list",
                length = PBAR)

# in blocks of 1000 ask MCL to load and process the set comparison results

