###### -- parse a lot of GFFs from factorial read subsetting ------------------

suppressMessages(library(SynExtend))

###### -- ad hoc functions ----------------------------------------------------

# erik's Split Violin plotting function
split_violin <- function(X, # a factor
                         Y1, # series of values
                         Y2, # a second series of values
                         colors,
                         YMAX = max(Y1,
                                    Y2,
                                    na.rm = TRUE),
                         adjust = 1,
                         las = 2,
                         CustomX = FALSE,
                         ...) {
  u_x <- sort(unique(X))
  m <- match(X, u_x)
  
  plot(NA,
       xlim=c(0.5, length(u_x) + 0.5),
       ylim=c(0, YMAX),
       xaxt="n",
       ...)
  if (CustomX) {
    # do not plot x labels, user must now supply an `axis` call post function
  } else {
    axis(1, seq_along(u_x), u_x, las=las, ...)
  }
  #abline(v=seq_along(u_x))
  
  for (i in seq_along(u_x)) {
    w <- which(m==i & !is.na(Y1))
    if (length(w) > 1) {
      d1 <- density(Y1[w], adjust=adjust, from=min(Y1[w]), to=max(Y1[w]))
      x1 <- d1$x
      y1 <- d1$y
      y1 <- y1/sum(y1)*length(y1)
      y1 <- y1/max(y1)/2 + i
      polygon(c(y1, i, i),
              c(x1, max(x1), min(x1)),
              col=colors[1],
              border=substring(colors[1], 1, 7))
      segments(i,
               median(Y1[w]),
               y1[which.min(abs(median(Y1[w]) - x1))],
               median(Y1[w]),
               col=substring(colors[1], 1, 7))
    }
    
    w <- which(m==i & !is.na(Y2))
    if (length(w) > 1) {
      d2 <- density(Y2[w], adjust=adjust, from=min(Y2[w]), to=max(Y2[w]))
      x2 <- d2$x
      y2 <- d2$y
      y2 <- y2/sum(y2)*length(y2)
      y2 <- -y2/max(y2)/2 + i
      polygon(c(y2, i, i),
              c(x2, max(x2), min(x2)),
              col=colors[2],
              border=substring(colors[2], 1, 7))
      segments(i,
               median(Y2[w]),
               y2[which.min(abs(median(Y2[w]) - x2))],
               median(Y2[w]),
               col=substring(colors[2], 1, 7))
    }
  }
}

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
  NOTE <- CONTIG <- TYPE <- ID <- vector(mode = "character",
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
                   "Note" = NOTE[KEEP])
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

ExtractParents <- function(GRangeObject,
                           AcceptableParents = NULL,
                           Verbose = FALSE) {
  
  
  p1 <- GRangeObject$ID
  p2 <- GRangeObject$Parent
  p4 <- as.character(GRangeObject@seqnames)
  p5 <- as.character(GRangeObject$type)
  
  if (is.null(p1) | is.null(p2)) {
    print("GRange object does not appear to have the expected ID and Parent relationships.")
    return(NULL)
  }
  if (length(p1) != length(p2)) {
    print("Feature IDs and Parents have an unexpected relationship.")
    return(NULL)
  }
  
  
  # loop through IDs,
  # if the ID has a parent(s)
  # create the directional edges
  # if it does not, skip
  p3 <- vector(mode = "list",
               length = length(p1))
  l1 <- lengths(p2)
  if (Verbose) {
    TIMESTART <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
    PBAR <- length(p3)
  }
  for (m1 in seq_along(p3)) {
    if (l1[m1] == 1) {
      # the ID has an associated parent, insert the edges
      p3[[m1]] <- c(p1[m1], p2[[m1]], p4[m1])
    } else if (l1[m1] > 1) {
      # the ID has multiple associated parents
      p3[[m1]] <- matrix(data = c(rep(p1[m1], l1[m1]),
                                  p2[[m1]],
                                  rep(p4[m1], l1[m1])),
                         ncol = 3L)
    } # else do nothing, leave the list position as NULL
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = m1 / PBAR)
    }
  }
  if (Verbose) {
    close(pBar)
    TIMEEND <- Sys.time()
    print("Parent IDs identified in:")
    print(TIMEEND - TIMESTART)
  }
  
  p3 <- do.call(rbind,
                p3)
  ParentFeatures <- p3[!(p3[, 1L] %in% p3[, 2L]), 2L]
  
  # AcceptableParents is null by default, if it's left at the default, match
  # a series of `normal-ish` acceptable parents,
  # else use user supplied vector to select out parent types
  if (is.null(AcceptableParents)) {
    AcceptableParents <- c("gene",
                           "pseudogene",
                           "tRNA",
                           "rRNA")
  }
  
  w1 <- p3[, 2L] %in% ParentFeatures
  Contigs <- unique(p3[w1, 3L])
  Parents2 <- which(p1 %in% ParentFeatures & p5 %in% AcceptableParents)
  
  ChildMatch <- findMatches(x = p3[w1, 1L],
                            table = p1)
  ParentMatch <- findMatches(x = p3[w1, 2L],
                             table = p1)
  # return(list(p3,
  #             w1,
  #             ParentMatch,
  #             ChildMatch))
  Relationships <- Hits(from = ParentMatch@to,
                        to = ChildMatch@to,
                        nLnode = max(ParentMatch@to),
                        nRnode = max(ChildMatch@to))
  
  # from DisjointSets
  p1 <- Relationships@from
  p2 <- Relationships@to
  P1Ids <- unique(p1)
  P2Ids <- unique(p2)
  AllIds <- unique(c(P1Ids,
                     P2Ids))
  F1 <- factor(p1,
               levels = AllIds)
  F2 <- factor(p2,
               levels = AllIds)
  
  FI1 <- as.integer(F1)
  FI2 <- as.integer(F2)
  FC1 <- as.character(F1)
  FC2 <- as.character(F2)
  
  FInts <- c(FI1, FI2)
  FChars <- c(FC1, FC2)
  Key1 <- !duplicated(FInts)
  IntMap <- data.frame("FactorRep" = FInts[Key1],
                       "UniqueIDs" = FChars[Key1])
  # run FindSets
  IntRes <- FindSets(p1 = FI1,
                     p2 = FI2,
                     Verbose = FALSE)
  rm(list = c("FI1",
              "FI2",
              "FC1",
              "FC2",
              "FInts",
              "FChars",
              "Key1"))
  # all nodes
  Nodes <- IntRes[, 1L, drop = TRUE]
  # all representatives
  Reps <- IntRes[, 2L, drop = TRUE]
  # set the order to the unique reps
  AllIds <- as.integer(IntMap[match(x = IntMap[, 1L],
                                    table = Nodes), 2L])
  
  UClusts <- split(x = AllIds,
                   f = Reps)
  # return(list(ParentFeatures,
  #             Relationships,
  #             UClusts,
  #             AllIds,
  #             Reps,
  #             IntMap))
  
  PHASE <- vector(mode = "integer",
                  length = length(UClusts))
  CONTIG <- ID <- NOTE <- vector(mode = "character",
                                 length = length(UClusts))
  PARENTSBYPRIORITY <- sapply(UClusts, min)
  START <- GRangeObject@ranges@start[PARENTSBYPRIORITY]
  STOP <- START + GRangeObject@ranges@width[PARENTSBYPRIORITY] - 1L
  STRAND <- ifelse(test = GRangeObject@strand[PARENTSBYPRIORITY] == "+",
                   yes = 0L,
                   no = 1L)
  ID <- GRangeObject$ID[PARENTSBYPRIORITY]
  TYPE <- as.character(GRangeObject$type[PARENTSBYPRIORITY])
  CONTIG <- p4[PARENTSBYPRIORITY]
  INFERENCE <- PRODUCT <- NOTE <- vector(mode = "character",
                                         length = length(UClusts))
  if (Verbose) {
    TIMESTART <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
    PBAR <- length(UClusts)
  }
  for (m1 in seq_along(UClusts)) {
    ph1 <- !is.na(GRangeObject$inference[UClusts[[m1]]])
    # if(sum(ph1) > 1) {
    #   return(list(m1,
    #               UClusts,
    #               ph1))
    # }
    INFERENCE[m1] <- paste(unique(GRangeObject$inference[UClusts[[m1]]][ph1]),
                           collapse = "::")
    ph1 <- !is.na(GRangeObject$product[UClusts[[m1]]])
    PRODUCT[m1] <- paste(unique(GRangeObject$product[UClusts[[m1]]][ph1]),
                         collapse = "::")
    ph1 <- lengths(GRangeObject$Note[UClusts[[m1]]]) > 0
    if (sum(ph1) > 0) {
      NOTE[m1] <- paste(unique(unlist(GRangeObject$Note[UClusts[[m1]]][ph1])),
                        collapse = "::")
    }
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = m1 / PBAR)
    }
  }
  if (Verbose) {
    close(pBar)
    TIMEEND <- Sys.time()
    print("Parent IDs identified in:")
    print(TIMEEND - TIMESTART)
  }
  RES <- DataFrame("Start" = START,
                   "Stop" = STOP,
                   "Type" = TYPE,
                   "Contig" = CONTIG,
                   "ID" = ID,
                   "Note" = NOTE,
                   "Inference" = INFERENCE,
                   "Product" = PRODUCT)
  
  return(RES)
  
  # for every ph1, record:
  # Index
  # Strand
  # Start
  # Stop
  # Type
  # ID
  # Range
  # 
  # Product
  # Coding
  # Translation Table
  # Contig name
  
}

###### -- some hard coded variables -------------------------------------------

PseudoTypeVector <- c("frameshifted",
                      "internal stop",
                      "partial abbutting assembly gap",
                      "partial in the middle",
                      "missing C-terminus",
                      "missing N-terminus")

# erik's original colors
ColVector <- c(total = "#00000033",
               lines = "green4",
               frame = "#B0130033",
               stop = "#00169933",
               points = "gray",
               transposase = "#FF0000",
               transporter = "#40E0D0",
               `hypothetical protein` = "#9932CC",
               other = "#FFDEAD") -> ColVec01

PLOTINLINE <- FALSE

###### -- code body -----------------------------------------------------------

# each assembly and annotation is an observation
# total genes, total pseudogenes, pseudogenes by split i.e.
# internal stop vs frameshift
# 
# subset total data to assemblies with more than half a million bp (neisseria should be ~2.8)

# hard coded paths:
# completed GFF
PATH01 <- "~/Data/20230111_Factorial_Annotations"
# initial data pulls, rebuild the job map inside this R session
# to match completed jobs against recorded NCBI calls if desired
PATH02 <- "~/Repos/Pseudogenes/SearchResults.RData"
PATH03 <- "~/Repos/Pseudogenes/SearchResults2.RData"
# completed fnas
PATH04 <- "~/Data/20230111_Factorial_Assemblies"
PATH05 <- "~/Repos/20221207_FactorialAssembly"

FILES01 <- list.files(path = PATH01,
                      full.names = TRUE)
FILES02 <- list.files(path = PATH04,
                      full.names = TRUE)
load(file = PATH02,
     verbose = TRUE)
load(file = PATH03,
     verbose = TRUE)

# wrangle pseudogenes into pg per kb
frameshiftperkb <- PseudosCondensed[, 1L] / (EntrezResults$Total_Length / 1000)
internalstopperkb <- PseudosCondensed[, 4] / (EntrezResults$Total_Length / 1000)
partialperkb <- PseudosCondensed[, 3L] / (EntrezResults$Total_Length / 1000)


b1 <- EntrezResults$Biosample[EntrezResults$Coverage >= 100]

SRAResults <- SRAResults[SRAResults$BioSample %in% b1, ]

# wrangle genus names
g1 <- unlist(regmatches(x = SRAResults$ScientificName,
                        m = gregexpr(pattern = "^[^ ]+",
                                     text = SRAResults$ScientificName)))
g2 <- grep(pattern = "candidatus",
           ignore.case = TRUE,
           x = g1)
g1[g2] <- unlist(regmatches(x = SRAResults$ScientificName[g2],
                            m = gregexpr(pattern = "(?<=Candidatus )([^ ]+)",
                                         text = SRAResults$ScientificName[g2],
                                         perl = TRUE)))
g2 <- grep(pattern = "[^A-Za-z]",
           x = g1)
g1[g2] <- unlist(regmatches(x = g1[g2],
                            m = gregexpr(pattern = "[A-Za-z]+",
                                         text = g1[g2])))

w1 <- SRAResults$Platform == "ILLUMINA"

g1 <- g1[w1]
SRAResults <- SRAResults[w1, ]

# head(sort(table(g1), decreasing = TRUE), n = 100)

pin <- "Neisseria"
w2 <- which(g1 == pin)
w2 <- unique(SRAResults$BioSample[w2])
w3 <- internalstopperkb[EntrezResults$Biosample %in% w2]
w4 <- frameshiftperkb[EntrezResults$Biosample %in% w2]
w5 <- names(table(EntrezResults$SpeciesName[EntrezResults$Biosample %in% w2]))
w5 <- match(x = EntrezResults$SpeciesName[EntrezResults$Biosample %in% w2],
            table = w5)

if (PLOTINLINE) {
  
  plot(y = w3,
       x = w4,
       ylab = "internal stops / kb",
       xlab = "frameshifts / kb",
       # xlim = c(0, 0.1),
       # ylim = c(0, 0.1),
       # xaxs="i",
       # yaxs="i",
       main = paste(pin, "pseudogenes"),
       col = w5)
}

SelectSRA <- SRAResults[SRAResults$BioSample %in% w2, ]

JobMap <- expand.grid(list("Run" = unique(SelectSRA$Run),
                           "TRIM" = c(TRUE, FALSE),
                           "Quality" = 1:4,
                           "DownSample" = c(1, 0.5, 0.25, 0.05)))
scn1 <- unlist(regmatches(x = SelectSRA$ScientificName,
                          m = gregexpr(pattern = "^[A-Za-z]+ [A-Za-z]+",
                                       text = SelectSRA$ScientificName)))


# split FILES01 apart to compartmentalize data parsing
s1 <- split(FILES01,
            1:700)

for (m1 in seq_along(s1)) {
  TIMESTART <- Sys.time()
  Res4 <- mclapply(X = s1[[m1]],
                   FUN = function(x) {
                     y <- rtracklayer::import(x)
                     z1 <- GRangeToDFrame(GRangeObject = y,
                                          Verbose = FALSE)
                     # z1 <- ExtractParents(GRangeObject = y,
                     #                      Verbose = FALSE)
                     return(z1)
                     # if (is.null(z1)) {
                     #   return(list(NULL, NULL))
                     # } else {
                     #   z2 <- gffToDataFrame(GFF = x,
                     #                        Verbose = FALSE)
                     #   return(list(z1, z2))
                     # }
                   },
                   mc.cores = 6L)
  TIMEEND <- Sys.time()
  print(m1)
  print(TIMEEND - TIMESTART)
  CurrentFiles <- s1[[m1]]
  save(Res4,
       CurrentFiles,
       file = paste0("~/tempcounts",
                     formatC(x = m1,
                             width = 3,
                             flag = 0,
                             format = "d"),
                     ".RData"),
       compress = "xz")
}

j1 <- list.files(path = "~",
                 pattern = "tempcounts",
                 full.names = TRUE)

j2 <- vector(mode = "list",
             length = length(j1))
pBar <- txtProgressBar(style = 1L)
PBAR <- length(j2)
LOOPTIMESTART <- Sys.time()
for (m1 in seq_along(j1)) {
  load(file = j1[m1],
       verbose = FALSE)
  j2[[m1]] <- vector(mode = "list",
                     length = length(Res4))
  
  for (m2 in seq_along(Res4)) {
    z1 <- Res4[[m2]]$Note[Res4[[m2]]$Type == "pseudogene"]
    if (length(z1) == 0) {
      # PGAP reported no pseudogenes
      j2[[m1]][[m2]] <- rep(0, length(PseudoTypeVector))
    } else {
      # PGAP reported any pseudogenes
      j2[[m1]][[m2]] <- NoteCheck(NoteVector = z1,
                                  CheckVector = PseudoTypeVector)
    }
  }
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)
LOOPTIMEEND <- Sys.time()
print(LOOPTIMEEND - LOOPTIMESTART)

j3 <- unlist(j2, recursive = FALSE)
j3 <- do.call(rbind,
              j3)

# save(Res4,
#      file = "~/tempdata.RData",
#      compress = "xz")

ID1 <- unlist(regmatches(x = unname(unlist(s1)),
                         m = gregexpr(pattern = "(?<=Annotation_)([A-Za-z0-9]+)(?=\\.gff)",
                                      text = unname(unlist(s1)),
                                      perl = TRUE)))
ID2 <- unlist(regmatches(x = FILES02,
                         m = gregexpr(pattern = "(?<=Assembly_)([A-Za-z0-9]+)(?=\\.fna\\.gz)",
                                      text = FILES02,
                                      perl = TRUE)))

j4 <- match(x = ID1,
            table = ID2)
j5 <- vector(mode = "list",
             length = length(j4))

pBar <- txtProgressBar(style = 1)
PBAR <- length(j4)
LOOPTIMESTART <- Sys.time()
for (m1 in seq_along(j4)) {
  dna <- readDNAStringSet(filepath = FILES02[j4[m1]])
  w6 <- sort(width(dna), decreasing = TRUE)
  TotalNucs <- sum(w6)
  HalfNuc <- ceiling(TotalNucs / 2)
  AllWidths <- cumsum(w6)
  w7 <- which(AllWidths >= HalfNuc)[1L]
  N50 <- w6[w7]
  L50 <- w7
  TC <- length(AllWidths)
  
  j5[[m1]] <- c(TotalNucs,
                N50,
                L50,
                TC)
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)
LOOPTIMEEND <- Sys.time()
print(LOOPTIMEEND - LOOPTIMESTART)
j6 <- do.call(rbind,
              j5)

ID3 <- as.integer(unlist(regmatches(x = ID1,
                                    m = gregexpr(pattern = "[0-9]+",
                                                 text = ID1))))
ID4 <- unlist(regmatches(x = ID1,
                         m = gregexpr(pattern = "[^0-9]+",
                                      text = ID1)))
FactorSet <- JobMap[ID3, 2:4]
SRADat <- SRAResults[match(x = JobMap$Run[ID3],
                           table = SRAResults$Run), 2:3]

ResTable <- data.frame("IS" = j3[, 2L],
                       "FR" = j3[, 1L],
                       "AssemblySize" = j6[, 1L],
                       "N50" = j6[, 2L],
                       "L50" = j6[, 3L],
                       "TC" = j6[, 4L],
                       "Assembler" = ID4,
                       "Spots" = SRADat[, 1L],
                       "TotalBases" = SRADat[, 2L],
                       "Trim" = FactorSet[, 1L],
                       "QualBin" = FactorSet[, 2L],
                       "DownSample" = FactorSet[, 3L],
                       "Run" = JobMap$Run[ID3])

save(ResTable,
     file = paste0(PATH05,
                   "/ParsedPseudogenes.RData"),
     compress = "xz")
