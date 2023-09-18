###### -- parse GFFs for the multi-run re-assemblies --------------------------

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

PseudoTypeVector <- c("frameshifted",
                      "internal stop",
                      "partial abbutting assembly gap",
                      "partial in the middle",
                      "missing C-terminus",
                      "missing N-terminus")

PATH01 <- "~/Data/20230117_MultiRunAssemblies"
PATH02 <- "~/Data/20230117_MultiRunAnnotations"
PATH03 <- "~/Data/20230117_MultiRunParsedGFFs"

RefTable1 <- "~/Repos/Pseudogenes/SearchResults.RData"
RefTable2 <- "~/Repos/Pseudogenes/SearchResults2.RData"

load(file = RefTable1,
     verbose = TRUE)
load(file = RefTable2,
     verbose = TRUE)

# keep only illumina and ONT, pacbio doesn't produce meaningful quality scores
Keep1 <- c("ILLUMINA",
           "OXFORD_NANOPORE",
           "PACBIO_SMRT",
           "LS454")
SRA <- SRAResults[SRAResults$Platform %in% Keep1 &
                    SRAResults$LibraryStrategy == "WGS" &
                    SRAResults$LibrarySource == "GENOMIC", ]
Keep2 <- table(SRA$BioSample)
Keep2 <- names(Keep2[Keep2 >= 2L])
SRA <- SRA[SRA$BioSample %in% Keep2, ]

w1 <- match(x = SRA$BioSample,
            table = EntrezResults$Biosample)
JobMap <- data.frame("PersistentID" = seq(nrow(SRA)),
                     "SRR" = SRA$Run,
                     "BioSample" = SRA$BioSample,
                     "Platform" = SRA$Platform,
                     "Model" = SRA$Model,
                     "SpeciesName" = EntrezResults$SpeciesName[w1],
                     "NucsTotal" = EntrezResults$Total_Length[w1],
                     "N50" = EntrezResults$ScaffoldN50[w1],
                     "Frameshifts" = PseudosCondensed[w1, 1L],
                     "InternalStops" = PseudosCondensed[w1, 4L])

FILES01 <- list.files(path = PATH01,
                      full.names = TRUE)
FILES02 <- list.files(path = PATH02,
                      full.names = TRUE)

f1 <- as.integer(unlist(regmatches(x = FILES01,
                                   m = gregexpr(pattern = "(?<=Assembly_)([0-9]+)(?=\\.fna)",
                                                text = FILES01,
                                                perl = TRUE))))
f2 <- as.integer(unlist(regmatches(x = FILES02,
                                   m = gregexpr(pattern = "(?<=Annotation_)([0-9]+)(?=\\.gff)",
                                                text = FILES02,
                                                perl = TRUE))))

s1 <- split(FILES02,
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
                   mc.cores = 10L)
  TIMEEND <- Sys.time()
  print(m1)
  print(TIMEEND - TIMESTART)
  CurrentFiles <- s1[[m1]]
  print(paste0(PATH03,
               "/tempcounts",
               formatC(x = m1,
                       width = 3,
                       flag = 0,
                       format = "d"),
               ".RData"))
  save(Res4,
       CurrentFiles,
       file = paste0(PATH03,
                     "/tempcounts",
                     formatC(x = m1,
                             width = 3,
                             flag = 0,
                             format = "d"),
                     ".RData"),
       compress = "xz")
}

FILES03 <- list.files(path = PATH03,
                      full.names = TRUE)

j2 <- vector(mode = "list",
             length = length(FILES03))

FILES04 <- vector(mode = "list",
                  length = length(j2))
pBar <- txtProgressBar(style = 1L)
PBAR <- length(j2)
LOOPTIMESTART <- Sys.time()
for (m1 in seq_along(FILES03)) {
  load(file = FILES03[m1],
       verbose = FALSE)
  j2[[m1]] <- vector(mode = "list",
                     length = length(Res4))
  FILES04[[m1]] <- CurrentFiles
  
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
FILES04 <- unlist(FILES04)

ID1 <- unlist(regmatches(x = unname(unlist(s1)),
                         m = gregexpr(pattern = "(?<=Annotation_)([A-Za-z0-9]+)(?=\\.gff)",
                                      text = unname(unlist(s1)),
                                      perl = TRUE)))
ID2 <- unlist(regmatches(x = FILES01,
                         m = gregexpr(pattern = "(?<=Assembly_)([A-Za-z0-9]+)(?=\\.fna\\.gz)",
                                      text = FILES01,
                                      perl = TRUE)))

j4 <- match(x = ID1,
            table = ID2)
j5 <- vector(mode = "list",
             length = length(j4))

pBar <- txtProgressBar(style = 1)
PBAR <- length(j4)
LOOPTIMESTART <- Sys.time()
for (m1 in seq_along(j4)) {
  dna <- readDNAStringSet(filepath = FILES01[j4[m1]])
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

ID3 <- as.integer(ID1)
ID4 <- as.integer(ID2)

ID5 <- as.integer(unlist(regmatches(x = FILES04,
                                    m = gregexpr(pattern = "(?<=Annotation_)([A-Za-z0-9]+)(?=\\.gff)",
                                                 text = FILES04,
                                                 perl = TRUE))))

ph1 <- match(x = ID5,
             table = ID3)
ph2 <- match(x = ID5,
             table = JobMap$PersistentID)

ResTable <- data.frame("IS" = j3[, 2L],
                       "FS" = j3[, 1L],
                       "AssemblySize" = j6[ph1, 1L],
                       "N50" = j6[ph1, 2L],
                       "L50" = j6[ph1, 3L],
                       "TC" = j6[ph1, 4L],
                       "Source_Size" = JobMap$NucsTotal[ph2],
                       "Source_FS" = JobMap$Frameshifts[ph2],
                       "Source_IS" = JobMap$InternalStops[ph2],
                       "Platform" = JobMap$Platform[ph2],
                       "Model" = JobMap$Model[ph2],
                       "BioSample" = JobMap$BioSample[ph2],
                       "Species" = JobMap$SpeciesName[ph2],
                       "Run" = JobMap$SRR[ph2])
# drop entries in the table that don't seem to have a matched entrez result
ResTable <- ResTable[!is.na(ResTable$Source_Size), ]

# content control
ph1 <- ResTable$Source_Size * 0.75
ph2 <- ResTable$Source_Size * 1.25
ph3 <- ResTable$AssemblySize >= ph1 & ResTable$AssemblySize <= ph2
# ph4 <- ResTable$Platform == "ILLUMINA"

ResTable <- ResTable[ph3, ]
w2 <- table(ResTable$BioSample)
w2 <- names(w2[w2 >= 2L])
ResTable <- ResTable[ResTable$BioSample %in% w2, ]

dir.create(path = "~/Data/20230221_MultiRunIntermediateData")
TestTable <- ResTable[ResTable$BioSample == "SAMN11786966", ]
save(TestTable,
     file = "~/Data/20230221_MultiRunIntermediateData/TestTable.RData",
     compress = "xz")
save(ResTable,
     file = "~/Data/20230221_MultiRunIntermediateData/ResTable.RData",
     compress = "xz")

# PseudosPerMB <- data.frame("Run_IS" = (ResTable$IS / ResTable$AssemblySize) * 1000000,
#                            "Run_FS" = (ResTable$FS / ResTable$AssemblySize) * 1000000,
#                            "Source_IS" = (ResTable$Source_IS / ResTable$Source_Size) * 1000000,
#                            "Source_FS" = (ResTable$Source_IS / ResTable$Source_Size) * 1000000)
# 
# SampleSet <- unique(ResTable$BioSample)
# 
# for (m1 in seq_along(SampleSet)) {
#   
#   ResSubSet <- ResTable$BioSample == SampleSet[m1]
#   
#   if (sum(ResSubSet) >= 5L) {
#     pdf(file = paste0("~/Plots/20230220_FactorialInitialPlots/BioSample_",
#                       SampleSet[m1],
#                       ".pdf"))
#     
#     plot(x = PseudosPerMB$Run_IS[ResSubSet],
#          y = PseudosPerMB$Run_FS[ResSubSet],
#          pch = match(x = ResTable$Platform[ResSubSet],
#                      table = unique(ResTable$Platform[ResSubSet])),
#          xlab = "IS per MB",
#          ylab = "FS per MB",
#          main = unique(ResTable$Species[ResSubSet])[1L])
#     
#     dev.off()
#   }
#   
#   
# }



