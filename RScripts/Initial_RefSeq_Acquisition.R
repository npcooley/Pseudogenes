###### -- Query pseudogenes based on available meta data ----------------------
# usage
# Rscript <this_script.R> <Arg1> <Arg2 - optional>
# Arg1 == directory location to place results and plots
# Arg2 == optional results already generated in the format that the plotting functions
# here expect

# script will attempt to create a directory if the target directory does not exist
# but it will not be able to create multiple levels at once

# requirements
# Bioconductor package SynExtend
# R (base?) package parallel
# NCBI's eutils

###### -- libraries -----------------------------------------------------------

suppressMessages(library(parallel))
suppressMessages(library(SynExtend))

###### -- arguments -----------------------------------------------------------

ARGS <- commandArgs(trailingOnly = TRUE)

DIR <- ARGS[1L]
if (!dir.exists(paths = DIR)) {
  dir.create(DIR)
}
RES <- ARGS[2L]
setwd(DIR)

###### -- adhoc functions -----------------------------------------------------

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

# erik's Split Violin plotting function
split_violin <- function(X,
                         Y1,
                         Y2,
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

###### -- code body 1 ---------------------------------------------------------
# use esearch to go out and grab some data
# scrape the counts of pseudogenes and available assembly and assembler information
# from the NCBI ftp site

# user supplied only a directory
# generate some files
if (length(ARGS) == 1L) {
  
  # part 1 a
  # use esearch to construct an initial query
  # this is a rather large query, so even though esearch is very fast,
  # this will still take a few seconds / minutes
  EntrezQuery <- paste("esearch -db assembly ",
                       "-query '",
                       '("Bacteria"[Organism] OR "Archaea"[Organism]) ',
                       'AND "latest refseq"[properties] ',
                       'AND "refseq has annotation"[properties] ',
                       'AND "taxonomy check ok"[filter] ',
                       "NOT partial[filter] ",
                       "NOT anomalous[filter]' ",
                       '| ',
                       'esummary ',
                       '| ',
                       'xtract -pattern DocumentSummary -element ',
                       'FtpPath_RefSeq ',
                       'BioSampleAccn ',
                       'AssemblyStatus ',
                       'SubmitterOrganization ',
                       'SubmissionDate ',
                       'Organism ',
                       'Taxid ',
                       'SpeciesName ',
                       'ContigN50 ',
                       'ScaffoldN50 ',
                       'Coverage ',
                       '-block Stat -if "@category" -equals total_length -element Stat',
                       sep = "")
  EntrezReplies1 <- system(command = EntrezQuery,
                           intern = TRUE,
                           timeout = 2000L)
  EntrezReplies2 <- strsplit(x = EntrezReplies1,
                             split = "\t",
                             fixed = TRUE)
  # EntrezReplies2 <- do.call(rbind,
  #                           strsplit(x = EntrezReplies1,
  #                                    split = "\t",
  #                                    fixed = TRUE))
  # expect a length of 12
  EntrezReplies2 <- do.call(rbind,
                            EntrezReplies2[lengths(EntrezReplies2) == 12])
  EntrezReplies3 <- data.frame("FTP" = EntrezReplies2[, 1L],
                               "Biosample" = EntrezReplies2[, 2L],
                               "Assembly_Status" = EntrezReplies2[, 3L],
                               "Submitter_Org" = EntrezReplies2[, 4L],
                               "Submission_Date" = EntrezReplies2[, 5L],
                               "Organism" = EntrezReplies2[, 6L],
                               "TaxID" = EntrezReplies2[, 7L],
                               "SpeciesName" = EntrezReplies2[, 8L],
                               "ContigN50" = as.integer(EntrezReplies2[, 9L]),
                               "ScaffoldN50" = as.integer(EntrezReplies2[, 10L]),
                               "Coverage" = as.numeric(gsub(pattern = "([^0-9]*)([0-9.]+)([^0-9]*)",
                                                            replacement = "\\2",
                                                            x = EntrezReplies2[, 11L])),
                               "Total_Length" = as.integer(EntrezReplies2[, 12L]))
  
  # part 1 b
  # given our e-search results scrape data from NCBI assembly reports and GFFs
  # if a second argument is given to this script this section will be skipped

  # scrape data associated with table
  # as of the writing of this script, the initial query returns approximately
  # 225,000 ftp folders to grab data from
  # we want to grab data from these as efficiently as possible
  
  ScrapeResults <- vector(mode = "list",
                          length = nrow(EntrezReplies3))
  TestList <- c("frameshifted",
                "partial abutting assembly gap",
                "partial in the middle",
                "internal stop",
                "missing C-terminus",
                "missing N-terminus")
  
  BatchSize <- 1000L
  MaxSearches <- 300 # if the first query returns 
  CurrentSearches <- 1L
  Continue <- TRUE
  Missing <- lengths(ScrapeResults) == 0L
  
  
  while (CurrentSearches <= MaxSearches &
         Continue) {
    w <- which(Missing)[seq(BatchSize)]
    if (any(is.na(w))) {
      w <- w[!is.na(w)]
    }
    T01 <- Sys.time()
    ScrapeResults[w] <- mclapply(X = EntrezReplies3$FTP[w],
                                 FUN = function(x) {
                                   
                                   # grab the GFF and parse out pseudogene annotations
                                   GFF <- paste0(x,
                                                 "/",
                                                 strsplit(x = x,
                                                          split = "/",
                                                          fixed = TRUE)[[1]][10],
                                                 "_genomic.gff.gz")
                                   gc01 <- try(rtracklayer::import(GFF),
                                               silent = TRUE)
                                   if (!is(object = gc01,
                                           class2 = "try-error")) {
                                     gc02 <- GRangeToDFrame(gc01)
                                     
                                     t1 <- grepl(pattern = TestList[1L],
                                                 x = gc02$Note[gc02$Type == "pseudogene"])
                                     t2 <- grepl(pattern = TestList[2L],
                                                 x = gc02$Note[gc02$Type == "pseudogene"])
                                     t3 <- grepl(pattern = TestList[3L],
                                                 x = gc02$Note[gc02$Type == "pseudogene"])
                                     t4 <- grepl(pattern = TestList[4L],
                                                 x = gc02$Note[gc02$Type == "pseudogene"])
                                     t5 <- grepl(pattern = TestList[5L],
                                                 x = gc02$Note[gc02$Type == "pseudogene"])
                                     t6 <- grepl(pattern = TestList[6L],
                                                 x = gc02$Note[gc02$Type == "pseudogene"])
                                     t7 <- cbind("fr" = t1,
                                                 "paag" = t2,
                                                 "pim" = t3,
                                                 "is" = t4,
                                                 "mct" = t5,
                                                 "mnt" = t6)
                                     t8 <- table(apply(X = t7,
                                                       MARGIN = 1,
                                                       FUN = function(x) {
                                                         paste(TestList[x],
                                                               collapse = " + ")
                                                       }))
                                   }
                                   
                                   # grab the assembly report parse out the sequencing technology
                                   # and the assembler
                                   REPORT <- paste0(x,
                                                    "/",
                                                    strsplit(x = x,
                                                             split = "/",
                                                             fixed = TRUE)[[1]][10],
                                                    "_assembly_report.txt")
                                   
                                   assembly_report_01 <- try(readLines(REPORT),
                                                             silent = TRUE)
                                   if (!is(object = assembly_report_01,
                                           class2 = "try-error")) {
                                     Assembly_Method <- grep(x = assembly_report_01,
                                                             pattern = "# Assembly method:")
                                     Assembly_Method <- assembly_report_01[Assembly_Method]
                                     Assembly_Technology <- grep(x = assembly_report_01,
                                                                 pattern = "# Sequencing technology:")
                                     Assembly_Technology <- assembly_report_01[Assembly_Technology]
                                     Assembly_Name <- grep(x = assembly_report_01,
                                                           pattern = "# Assembly name:")
                                     Assembly_Name <- assembly_report_01[Assembly_Name]
                                   }
                                   if (!is(object = gc01,
                                           class2 = "try-error") &
                                       !is(object = assembly_report_01,
                                           class2 = "try-error")) {
                                     list(t8,
                                          Assembly_Method,
                                          Assembly_Technology,
                                          Assembly_Name)
                                   } else {
                                     NULL
                                   }
                                 },
                                 mc.cores = detectCores())
    T02 <- Sys.time()
    print(T02 - T01)
    print(paste0(sum(!Missing),
                 " completed on iteration ",
                 CurrentSearches,
                 "!"))
    CurrentSearches <- CurrentSearches + 1L
    Missing <- lengths(ScrapeResults) == 0L
    if (sum(!Missing) == length(ScrapeResults)) {
      Continue <- FALSE
    }
    # setTxtProgressBar(pb = pBar,
    #                   value = CurrentSearches / MaxSearches)
  }
  
  meta1 <- lapply(X = ScrapeResults,
                  FUN = function(x) {
                    if (length(x) > 1L) {
                      x[[1]]
                    } else {
                      NULL
                    }
                  })
  possiblecats <- unique(names(unlist(meta1)))
  meta2 <- lapply(X = meta1,
                  FUN = function(x) {
                    if (length(x) > 0L) {
                      y <- rep(0, length(possiblecats))
                      y[match(x = names(x),
                              table = possiblecats)] <- unname(x)
                      y
                    } else {
                      rep(0, length(possiblecats))
                    }
                  })
  meta2 <- do.call(rbind, meta2)
  colnames(meta2) <- possiblecats
  
  meta3 <- lapply(X = ScrapeResults,
                  FUN = function(x) {
                    if (length(x) > 1) {
                      unlist(x[2:length(x)])
                    } else {
                      NULL
                    }
                  })
  meta4 <- sapply(X = meta3,
                  FUN = function(x) {
                    unlist(regmatches(x = x,
                                      m = gregexpr(pattern = "(?<=# Sequencing technology: )(.+)",
                                                   text = x,
                                                   perl = TRUE)))
                  })
  meta4 <- sapply(X = meta4,
                  FUN = function(x) {
                    if (length(x) > 0) {
                      paste(x,
                            sep = "; ")
                    } else {
                      ""
                    }
                  })
  # sequencing keywords sets:
  # incudes others added by checking what remained uncaptured originally
  
  # Illumina
  # HiSeq
  # MiSeq
  # NextSeq
  # NovaSeq
  # Nextera
  
  ## fisher technology
  # IonTorrent
  # Ion Torrent
  # PGM
  
  # Sanger
  
  # Roche
  # 454
  
  # Pacific
  # PacBio
  # SMRT
  
  # ABI
  # SOLID
  
  # Oxford
  # Nanopore
  # minION
  # gridION
  # PromethION
  
  meta5 <- cbind("Illumina" = grepl(pattern = "Illumina|HiSeq|MiSeq|NextSeq|NovaSeq|iSeq|Nextera|Solexa",
                                    x = meta4,
                                    ignore.case = TRUE),
                 "Fisher" = grepl(pattern = "IonTorrent|Ion Torrent|PGM|Genome Analyzer|3730xl|Ion Personal Genome|Ion GeneStudio",
                                  x = meta4,
                                  ignore.case = TRUE),
                 "Sanger" = grepl(pattern = "Sanger",
                                  x = meta4,
                                  ignore.case = TRUE),
                 "Roche" = grepl(pattern = "Roche|454|Titanium|GS-FLX",
                                 x = meta4,
                                 ignore.case = TRUE),
                 "ABI" = grepl(pattern = "ABI|SOLID",
                               x = meta4,
                               ignore.case = TRUE),
                 "Oxford_Nanopore" = grepl(pattern = "Oxford|Nanopore|minION|gridION|PromethION",
                                           x = meta4,
                                           ignore.case = TRUE),
                 "PacBio" = grepl(pattern = "Pacific|PacBio|SMRT",
                                  x = meta4,
                                  ignore.case = TRUE),
                 "MGI" = grepl(pattern = "Complete Genomics|MGI|DBN",
                               x = meta4,
                               ignore.case = TRUE))
  meta5 <- cbind(meta5,
                 "Unknown_Other_Uncaptured" = apply(X = meta5,
                                                    MARGIN = 1,
                                                    FUN = function(x) {
                                                      !any(x)
                                                    }))
  
  
  meta6 <- sapply(X = meta3,
                  FUN = function(x) {
                    if (length(x) > 0 &
                        any(grepl(pattern = "# Assembly method: ",
                                  x = x))) {
                      paste(unlist(regmatches(x = x,
                                        m = gregexpr(pattern = "(?<=# Assembly method: )(.+)",
                                                     text = x,
                                                     perl = TRUE))),
                            sep = "; ")
                    } else {
                      ""
                    }
                    
                  })
  
  # Assembler keyword list:
  # SPADES
  # Unicycler
  # HGAP
  # miniasm
  # Racon
  # abyss
  # flye
  # SOAPdenovo
  # SRMT
  # Megahit
  # PATRIC
  # CANU
  # RAST
  # Trycycler
  # Newbler
  # Barapost
  # Pilon
  # skesa
  # CLC
  # Geneious
  # at some point it is not my problem if these aren't spelled correctly
  
  meta7 <- cbind("SPADES" = grepl(pattern = "SPADES",
                                  x = meta6,
                                  ignore.case = TRUE),
                 "Unicycler" = grepl(pattern = "Unicycler|Unicycle|Unicyler|Unicylcer|Unlcycler|Uniciclyer|Uniclycler|Uncycler|Unycicler|Unicyckler",
                                     x = meta6,
                                     ignore.case = TRUE),
                 "HGAP" = grepl(pattern = "HGAP",
                                x = meta6,
                                ignore.case = TRUE),
                 "miniasm" = grepl(pattern = "miniasm",
                                   x = meta6,
                                   ignore.case = TRUE),
                 "Racon" = grepl(pattern = "Racon",
                                 x = meta6,
                                 ignore.case = TRUE),
                 "Pilon" = grepl(pattern = "Pilon",
                                 x = meta6,
                                 ignore.case = TRUE),
                 "Abyss" = grepl(pattern = "Abyss",
                                 x = meta6,
                                 ignore.case = TRUE),
                 "Flye" = grepl(pattern = "Flye|Fyle|Fly",
                                x = meta6,
                                ignore.case = TRUE),
                 "SOAPDenovo" = grepl(pattern = "SOAPDenovo|Soap|S0ap|NextDenovo",
                                      x = meta6,
                                      ignore.case = TRUE),
                 "SMRT" = grepl(pattern = "SMRT|SMART",
                                x = meta6,
                                ignore.case = TRUE),
                 "Megahit" = grepl(pattern = "Megahit",
                                   x = meta6,
                                   ignore.case = TRUE),
                 "PATRIC" = grepl(pattern = "PATRIC",
                                  x = meta6,
                                  ignore.case = TRUE),
                 "CANU" = grepl(pattern = "CANU",
                                x = meta6,
                                ignore.case = TRUE),
                 "RAST" = grepl(pattern = "RAST",
                                x = meta6,
                                ignore.case = TRUE),
                 "Trycycler" = grepl(pattern = "Trycycler|Tycycler",
                                     x = meta6,
                                     ignore.case = TRUE),
                 "Newbler" = grepl(pattern = "Newbler",
                                   x = meta6,
                                   ignore.case = TRUE),
                 "Barapost" = grepl(pattern = "Barapost",
                                    x = meta6,
                                    ignore.case = TRUE),
                 "Pilon" = grepl(pattern = "Pilon",
                                 x = meta6,
                                 ignore.case = TRUE),
                 "Skesa" = grepl(pattern = "Skesa",
                                 x = meta6,
                                 ignore.case = TRUE),
                 "CLC" = grepl(pattern = "CLC|CLGenomics| CL genomics|CLS|workbench",
                               x = meta6,
                               ignore.case = TRUE),
                 "Geneious" = grepl(pattern = "Geneious|Genious",
                                    x = meta6,
                                    ignore.case = TRUE),
                 "FALCON" = grepl(pattern = "FALCON",
                                  x = meta6,
                                  ignore.case = TRUE),
                 "DNASTAR" = grepl(pattern = "DNASTAR",
                                   x = meta6,
                                   ignore.case = TRUE),
                 "Velvet" = grepl(pattern = "Velvet",
                                  x = meta6,
                                  ignore.case = TRUE),
                 "Shovill" = grepl(pattern = "Shovill|shovel",
                                   x = meta6,
                                   ignore.case = TRUE),
                 "AllPaths" = grepl(pattern = "AllPaths",
                                    x = meta6,
                                    ignore.case = TRUE),
                 "GS_De_Novo" = grepl(pattern = "GS De Novo",
                                      x = meta6,
                                      ignore.case = TRUE),
                 "Celera" = grepl(pattern = "Celera",
                                  x = meta6,
                                  ignore.case = TRUE),
                 "MaSuRCa" = grepl(pattern = "Masurca",
                                   x = meta6,
                                   ignore.case = TRUE),
                 "IBDA" = grepl(pattern = "IBDA|IDBA",
                                x = meta6,
                                ignore.case = TRUE),
                 "A5" = grepl(pattern = "A5|miseq",
                              x = meta6,
                              ignore.case = TRUE),
                 "Platanus" = grepl(pattern = "Platanus|platunus",
                                    x = meta6,
                                    ignore.case = TRUE),
                 "BWA" = grepl(pattern = "BWA",
                               x = meta6,
                               ignore.case = TRUE),
                 "MIRA" = grepl(pattern = "MIRA",
                                x = meta6,
                                ignore.case = TRUE),
                 "SNIPPY" = grepl(pattern = "SNIPPY",
                                  x = meta6,
                                  ignore.case = TRUE),
                 "bacass" = grepl(pattern = "bacass",
                                  x = meta6,
                                  ignore.case = TRUE),
                 "Ridom" = grepl(pattern = "Ridom",
                                  x = meta6,
                                  ignore.case = TRUE),
                 "Edena" = grepl(pattern = "Edena",
                                 x = meta6,
                                 ignore.case = TRUE),
                 "Raven" = grepl(pattern = "Raven",
                                 x = meta6,
                                 ignore.case = TRUE),
                 "Galaxy" = grepl(pattern = "Galaxy",
                                 x = meta6,
                                 ignore.case = TRUE),
                 "gcType" = grepl(pattern = "gcType",
                                  x = meta6,
                                  ignore.case = TRUE),
                 "BGI" = grepl(pattern = "BGI",
                                  x = meta6,
                                  ignore.case = TRUE))
  
  meta7 <- cbind(meta7,
                 "Unknown_Other_Uncaptured" = apply(X = meta7,
                                                    MARGIN = 1,
                                                    FUN = function(x) {
                                                      !any(x)
                                                    }))
  
  EntrezResults <- EntrezReplies3
  AssemblersLogical <- meta7
  TechnologyLogical <- meta5
  # ProksData <- cbind(EntrezReplies3,
  #                    meta5,
  #                    meta7)
  PseudosAll <- meta2
  # capture total counts of all annotation occurrences 
  CaptureVector <- c("frameshift",
                     "partial abutting assembly gap",
                     "partial in the middle",
                     "internal stop",
                     "missing")
  CurrentPseudoOpts <- colnames(PseudosAll)
  PseudosCondensed <- t(apply(X = PseudosAll,
                              MARGIN = 1,
                              FUN = function(x) {
                                c(sum(x[grepl(pattern = CaptureVector[1L],
                                              x = CurrentPseudoOpts)]),
                                  sum(x[grepl(pattern = CaptureVector[2L],
                                              x = CurrentPseudoOpts)]),
                                  sum(x[grepl(pattern = CaptureVector[3L],
                                              x = CurrentPseudoOpts)]),
                                  sum(x[grepl(pattern = CaptureVector[4L],
                                              x = CurrentPseudoOpts)]),
                                  sum(x[grepl(pattern = CaptureVector[5L],
                                              x = CurrentPseudoOpts) &
                                          !grepl(pattern = CaptureVector[2L],
                                                 x = CurrentPseudoOpts) &
                                          !grepl(pattern = CaptureVector[3L],
                                                 x = CurrentPseudoOpts)]))
                              }))
  colnames(PseudosCondensed) <- CaptureVector
  
  save(ScrapeResults,
       file = paste0(getwd(),
                     "/ScrapedData.RData"),
       compress = "xz")
  save(EntrezResults,
       AssemblersLogical,
       TechnologyLogical,
       PseudosAll,
       PseudosCondensed,
       file = paste0(getwd(),
                     "/SearchResults.RData"),
       compress = "xz")
  
} else if (length(ARGS) == 2L) {
  print("Loading pre-specified data.")
  load(file = RES,
       verbose = TRUE)
  LOADEDFILES <- ls()
  # expect a data frame of esearch results
  # a logical matrix of specified assemblers
  # a logical matrix of specified assembly technology
  # an integer matrix of pseudogene counts all unique coding assignment possibilities
  # an integer matrix of pseudogene counts condensed counts
  EXPECTEDFILES <- c("EntrezResults",
                     "AssemblersLogical",
                     "TechnologyLogical",
                     "PseudosAll",
                     "PseudosCondensed")
  if (!all(EXPECTEDFILES %in% LOADEDFILES)) {
    stop ("An expected file is missing or not named correctly.")
  }
}


###### -- code body part 2 ----------------------------------------------------

if (dir.exists(paths = past0(getwd(),
                             "/Figures"))) {
  # dir already exists, do nothing
} else {
  dir.create(path = past0(getwd(),
                          "/Figures"))
}

# plots
###### -- convenience variables -----------------------------------------------
# erik's original colors
ColVector <- c(total = "#00000033",
               lines = "green4",
               frame = "#B0130033",
               stop = "#00169933",
               points = "gray",
               transposase = "#FF0000",
               transporter = "#40E0D0",
               `hypothetical protein` = "#9932CC",
               other = "#FFDEAD")
# pseudogenization events (limit of 1 per gene) per kb
z1 <- PseudosCondensed[, 1L] / (EntrezResults$Total_Length / 1000) # frame shifts
z2 <- PseudosCondensed[, 2L] / (EntrezResults$Total_Length / 1000) # assembly gap
z3 <- PseudosCondensed[, 3L] / (EntrezResults$Total_Length / 1000) # partial
z4 <- PseudosCondensed[, 4L] / (EntrezResults$Total_Length / 1000) # internal stop
z5 <- PseudosCondensed[, 5L] / (EntrezResults$Total_Length / 1000) # truncated but not partial

# submission date
# year
z6 <- strsplit(x = EntrezResults$Submission_Date,
               split = "/",
               fixed = TRUE)
z6 <- as.integer(sapply(X = z6,
                        FUN = function(x) {
                          x[1]
                        },
                        simplify = TRUE))
# month as a percentage
z7 <- strsplit(x = EntrezResults$Submission_Date,
               split = "/",
               fixed = TRUE)
z7 <- as.numeric(sapply(X = z7,
                        FUN = function(x) {
                          x[2]
                        },
                        simplify = TRUE)) / 12
# day as a percentage
z8 <- strsplit(x = EntrezResults$Submission_Date,
               split = "/| ")
z8 <- as.numeric(sapply(X = z8,
                        FUN = function(x) {
                          x[3]
                        },
                        simplify = TRUE)) / 365

# named genus
g1 <- unlist(regmatches(x = EntrezResults$SpeciesName,
                        m = gregexpr(pattern = "^[^ ]+",
                                     text = EntrezResults$SpeciesName)))
g2 <- grep(pattern = "candidatus",
           ignore.case = TRUE,
           x = g1)
g1[g2] <- unlist(regmatches(x = EntrezResults$SpeciesName[g2],
                            m = gregexpr(pattern = "(?<=Candidatus )([^ ]+)",
                                         text = EntrezResults$SpeciesName[g2],
                                         perl = TRUE)))
g2 <- grep(pattern = "[^A-Za-z]",
           x = g1)
g1[g2] <- unlist(regmatches(x = g1[g2],
                            m = gregexpr(pattern = "[A-Za-z]+",
                                         text = g1[g2])))

# assembler breakdown
a1 <- apply(X = AssemblersLogical,
            MARGIN = 1,
            FUN = function(x) {
              paste(colnames(AssemblersLogical)[x],
                    collapse = "\n+ ")
            })

# technology breakdown
a2 <- apply(X = TechnologyLogical,
            MARGIN = 1,
            FUN = function(x) {
              paste(colnames(TechnologyLogical)[x],
                    collapse = "\n+ ")
            })
# assembly level
a3 <- factor(paste0(ifelse(test = nchar(EntrezResults$Assembly_Status) > 6,
                           yes = substring(text = EntrezResults$Assembly_Status,
                                           first = 1,
                                           last = 5),
                           no = EntrezResults$Assembly_Status),
                    ifelse(test = nchar(EntrezResults$Assembly_Status) > 6,
                           yes = ".",
                           no = "")))
###### -- plot 1 --------------------------------------------------------------
# plots 1! general information scraped from the 
# meta data
pdf(file = paste0(getwd(),
                  "/Figures/GeneralPlots01.pdf"),
    height = 7,
    width = 7)
layout(mat = matrix(data = 1:4,
                    nrow = 2,
                    byrow = TRUE))
par(mar = c(4, 4, 1, 1),
    mgp = c(2.3,1, 0))
# panel a
# pseudogenes by type
w <- z6 >= 2000
split_violin(X = z6[w],
             Y1 = z1[w],
             Y2 = z4[w],
             colors = ColVector[c("frame", "stop")],
             YMAX = 0.1,
             xlab = "",
             ylab = "Pseudogenes per kbp")

# panel b
# T1 <- tapply((z1 + z4), g1, c)
T1 <- tapply(X = (z1 + z2)[w],
             INDEX = g1[w],
             FUN = function(x) {
               c(x)
             })
GX <- sapply(X = T1,
             FUN = function(x) {
               sum(!is.na(x))
             })
GY <- sapply(X = T1,
             FUN = function(x) {
               median(x)
             })
UQ <- sapply(X = T1,
             FUN = function(x) {
               quantile(x, 0.75, na.rm = TRUE)
             })
LQ <- sapply(X = T1,
             FUN = function(x) {
               quantile(x, 0.25, na.rm = TRUE)
             })
# arrow gets mad about very short segments...
MQ <- abs(UQ - LQ) > 0.0005
plot(x = GX,
     y = GY,
     xlab = "Genomes per genus",
     ylab = "Average pseudogenes per kbp",
     log = "x",
     ylim = c(0, 0.2),
     xlim = c(min(GX),
              max(GX*10)),
     col = ColVector[5])
# arrows(GX, UQ, GY, LQ, code = 3, length = 0.02, angle = 90, col = ColVector[5])
arrows(x0 = GX[MQ],
       y0 = UQ[MQ],
       x1 = GX[MQ],
       y1 = LQ[MQ],
       code = 3,
       length = 0.02,
       angle = 90,
       col = ColVector[5])
# w1 <- names(head(sort(table(g1), decreasing = TRUE), n = 10))
# w1 <- names(GY[head(order(GY, decreasing = TRUE), n = 10)])
# text(x = GX[names(GX) %in% w1],
#      y = GY[names(GY) %in% w1],
#      labels = w1)
# points(x = GX[names(GX) %in% w1],
#        y = GY[names(GY) %in% w1],
#        col = ColVector[2L],
#        pch = 16)

# panel c & d
genera <- c("Shigella", "Francisella")
for (m1 in genera) {
  W <- which(g1 == m1)
  X <- (z6 + z7 + z8)[W]
  t_frame <- tapply(z1[W], X, range)
  t_stop <- tapply(z4[W], X, range)
  # t_frame <- tapply((pseudo_frame[w]/as.numeric(r$Size)[w]/1000)[W], X, range)
  # t_stop <- tapply((pseudo_stop[w]/as.numeric(r$Size)[w]/1000)[W], X, range)
  plot(x = 0,
       y = 0,
       type = "n",
       xlab = "Release date",
       ylab = "Pseudogenes per kbp",
       xlim = c(2000, 2021),
       ylim = range(unlist(t_frame),
                    unlist(t_stop),
                    na.rm=TRUE))
  segments(as.numeric(names(t_frame)),
           sapply(t_frame, `[`, 1L),
           as.numeric(names(t_frame)),
           sapply(t_frame, `[`, 2L),
           col = paste0(substring(ColVector["frame"], 1, 7), "66"),
           lwd = 2)
  segments(as.numeric(names(t_stop)),
           sapply(t_stop, `[`, 1L),
           as.numeric(names(t_stop)),
           sapply(t_stop, `[`, 2L),
           col = paste0(substring(ColVector["stop"], 1, 7), "66"),
           lwd = 2)
  legend("topleft", m1, bty = "n", text.font = 3)
}
dev.off()

###### -- plot 2 --------------------------------------------------------------

pdf(file = paste0(getwd(),
                  "/Figures/GeneralPlots02.pdf"),
    height = 7,
    width = 7)
layout(mat = matrix(data = 1:4,
                    nrow = 2,
                    byrow = TRUE))
par(mar = c(5, 4, 4, 1),
    mgp = c(2.3, 1, 0))
# plot 2! more meta data, submission completeness
# panel a
# violin of pseudogenes by assembly level
split_violin(X = a3,
             Y1 = z1,
             Y2 = z4,
             colors = ColVector[c("frame", "stop")],
             YMAX = 0.1,
             xlab = "Assembly Status",
             ylab = "Pseudogenes per kbp",
             las = 1)

# panel b
# reported coverage
CovPresent <- EntrezResults$Coverage > 0.1
plot(x = 1,
     y = 1,
     type = "n",
     xlab = "Reported Coverage",
     ylab = "Pseudogenes per kb",
     xlim = c(1, 10000),
     ylim = c(0, 0.4),
     log = "x")
points(x = EntrezResults$Coverage[CovPresent],
       y = z1[CovPresent],
       pch = 46,
       col = paste0(substring(ColVector["frame"], 1, 7), "66"))
points(x = EntrezResults$Coverage[CovPresent],
       y = z4[CovPresent],
       pch = 46,
       col = paste0(substring(ColVector["stop"], 1, 7), "66"))

# panel c
# technology
TechSubSet01 <- table(a2)
TechSubSet01 <- names(TechSubSet01[TechSubSet01 > 1000])
TechSubSet02 <- a2 %in% TechSubSet01
# a4 <- gsub(pattern = "Illumina",
#            replacement = "Ill.",
#            x = a2)
a4 <- gsub(pattern = "Oxford_Nanopore",
           replacement = "Ox. Np.",
           x = a2)
a4 <- gsub(pattern = "Unknown_Other_Uncaptured",
           replacement = "Other",
           x = a4)


split_violin(X = as.factor(a4[TechSubSet02]),
             Y1 = z1[TechSubSet02],
             Y2 = z4[TechSubSet02],
             colors = ColVector[c("frame", "stop")],
             YMAX = 0.2,
             xlab = "",
             ylab = "Pseudogenes per kbp",
             CustomX = TRUE)
AXISLABELS01 <- sort(unique(as.factor(a4[TechSubSet02])))
axis(1, seq_along(AXISLABELS01), AXISLABELS01, cex.axis = 0.75, las = 2)

# panel d
# assembler
AssemblerSubSet01 <- table(a1)
AssemblerSubSet01 <- names(AssemblerSubSet01[AssemblerSubSet01 > 1000])
AssemblerSubSet02 <- a1 %in% AssemblerSubSet01
a5 <- gsub(pattern = "Unknown_Other_Uncaptured",
           replacement = "Other",
           x = a1)
split_violin(X = as.factor(a5[AssemblerSubSet02]),
             Y1 = z1[AssemblerSubSet02],
             Y2 = z4[AssemblerSubSet02],
             colors = ColVector[c("frame", "stop")],
             YMAX = 0.2,
             xlab = "",
             ylab = "Pseudogenes per kbp",
             CustomX = TRUE)
AXISLABELS01 <- sort(unique(as.factor(a5[AssemblerSubSet02])))
axis(1, seq_along(AXISLABELS01), AXISLABELS01, cex.axis = 0.5, las = 2)
dev.off()

###### -- plot 3, 4, and 5 ----------------------------------------------------
# plot pseudogenes per kb for each assembler / technology
# pair that produced more than 100 mb
# at the genus level, and again at the taxID level
# if all points are tightly clustered, the combinations are indistinct
# if they are not, some combos are producing different outcomes
# the context necessary to determine whether this is directly a result of that
# combination and not experimental effects from the sources themselves does not
# appear to be present in this data as it is collected so far

pdf(file = paste0(getwd(),
                  "/Figures/GeneralPlots03.pdf"),
    height = 7,
    width = 7)
layout(mat = matrix(data = 1:6,
                    nrow = 2,
                    ncol = 3))
par(mar = c(3.5, 2.5, 2.5, 1),
    mgp = c(1.35, 0.45, 0))
zz1 <- names(table(g1)[table(g1) > 10000])
zzz1 <- unname(table(g1)[table(g1) > 10000])
for (m1 in seq_along(zz1)) {
  w1 <- g1 == zz1[m1]
  w2 <- names(table(a1[w1])[table(a1[w1]) > 100])
  w3 <- names(table(a2[w1])[table(a2[w1]) > 100])
  w4 <- w5 <- w6 <- w9 <-  matrix(data = 0,
                                  nrow = length(w2),
                                  ncol = length(w3))
  for (m2 in seq_along(w2)) {
    for (m3 in seq_along(w3)) {
      w2a <- a1 == w2[m2]
      w3a <- a2 == w3[m3]
      w4[m2, m3] <- sum(PseudosCondensed[w1 & w2a & w3a, 1L])
      w5[m2, m3] <- sum(PseudosCondensed[w1 & w2a & w3a, 4L])
      w6[m2, m3] <- sum(EntrezResults$Total_Length[w1 & w2a & w3a])
      w9a <- gsub(x = w2[m2],
                  pattern = "\n+ ",
                  replacement = "; ")
      w9a <- gsub(x = w9a,
                  pattern = "Unknown_Other_Uncaptured",
                  replacement = "Other")
      w9b <- gsub(x = w3[m3],
                  pattern = "\n+ ",
                  replacement = "; ")
      w9b <- gsub(x = w9b,
                  pattern = "Unknown_Other_Uncaptured",
                  replacement = "Other")
      w9[m2, m3] <- paste(w9a,
                          w9b,
                          sep = "\n+ ")
    }
  }
  # limit to 100 million bp combos
  w7 <- w4[w6 > 100000000] / w6[w6 > 100000000] # frameshifts
  # shift to pseudogenes per kbp -- was pseudogenes per bp
  w7 <- w7 * 1000
  w8 <- w5[w6 > 100000000] / w6[w6 > 100000000] # internal stops
  # shift to pseudogenes per kbp -- was pseudogenes per bp
  w8 <- w8 * 1000
  w9 <- w9[w6 > 100000000]
  plot(x = w7,
       y = w8,
       xlim = c(0, 0.04),
       ylim = c(0, 0.04),
       xlab = if (m1 %in% c(2,4,6)) {
         "frameshifts per kb"
       } else {
         NA
       },
       ylab = if (m1 %in% 1:2) {
         "internal stops per kb"
       } else {
         NA
       },
       main = paste0(zz1[m1],
                     " (",
                     zzz1[m1],
                     ")"))
  # text(x = w7,
  #      y = w8,
  #      labels = w9)
}
dev.off()

# repeat the above plot for the top 6 taxids - our usual suspects
pdf(file = paste0(getwd(),
                  "/Figures/GeneralPlots04.pdf"),
    height = 7,
    width = 7)
layout(mat = matrix(data = 1:6,
                    nrow = 2,
                    ncol = 3))
par(mar = c(3.5, 2.5, 2.5, 1),
    mgp = c(1.35, 0.45, 0))
# grab the top 6 appearing taxIDs
zz1 <- names(head(sort(table(EntrezResults$TaxID), decreasing = TRUE)))
zzz1 <- unname(head(sort(table(EntrezResults$TaxID), decreasing = TRUE)))
TitleVector <- c("E. coli",
                 "K pneumonia",
                 "S. aureus",
                 "S. pneumonia",
                 "P. aeruginosa",
                 "A. baumanii")
for (m1 in seq_along(zz1)) {
  w1 <- EntrezResults$TaxID == zz1[m1]
  w2 <- names(table(a1[w1])[table(a1[w1]) > 100])
  w3 <- names(table(a2[w1])[table(a2[w1]) > 100])
  w4 <- w5 <- w6 <- w9 <- matrix(data = 0,
                                 nrow = length(w2),
                                 ncol = length(w3))
  for (m2 in seq_along(w2)) {
    for (m3 in seq_along(w3)) {
      w2a <- a1 == w2[m2]
      w3a <- a2 == w3[m3]
      w4[m2, m3] <- sum(PseudosCondensed[w1 & w2a & w3a, 1L])
      w5[m2, m3] <- sum(PseudosCondensed[w1 & w2a & w3a, 4L])
      w6[m2, m3] <- sum(EntrezResults$Total_Length[w1 & w2a & w3a])
      w9a <- gsub(x = w2[m2],
                  pattern = "\n+ ",
                  replacement = "; ")
      w9a <- gsub(x = w9a,
                  pattern = "Unknown_Other_Uncaptured",
                  replacement = "Other")
      w9b <- gsub(x = w3[m3],
                  pattern = "\n+ ",
                  replacement = "; ")
      w9b <- gsub(x = w9b,
                  pattern = "Unknown_Other_Uncaptured",
                  replacement = "Other")
      w9[m2, m3] <- paste(w9a,
                          w9b,
                          sep = "\n+ ")
    }
  }
  # limit to 100 million bp combos
  w7 <- w4[w6 > 100000000] / w6[w6 > 100000000] # frameshifts
  # shift to pseudogenes per kbp -- was pseudogenes per bp
  w7 <- w7 * 1000
  w8 <- w5[w6 > 100000000] / w6[w6 > 100000000] # internal stops
  # shift to pseudogenes per kbp -- was pseudogenes per bp
  w8 <- w8 * 1000
  w9 <- w9[w6 > 100000000]
  plot(x = w7,
       y = w8,
       xlim = c(0, 0.04),
       ylim = c(0, 0.04),
       xlab = if (m1 %in% c(2,4,6)) {
         "frameshifts per kb"
       } else {
         NA
       },
       ylab = if (m1 %in% 1:2) {
         "internal stops per kb"
       } else {
         NA
       },
       main = paste0(TitleVector[m1],
                     " (",
                     zzz1[m1],
                     ")"))
  # text(x = w7,
  #      y = w8,
  #      labels = w9)
}
dev.off()

# second 6 taxids, is there enough data to evaluate here
# if so, is there anything interesting
pdf(file = paste0(getwd(),
                  "/Figures/GeneralPlots05.pdf"),
    height = 7,
    width = 7)
layout(mat = matrix(data = 1:6,
                    nrow = 2,
                    ncol = 3))
par(mar = c(3.5, 2.5, 2.5, 1),
    mgp = c(1.35, 0.45, 0))
zz1 <- names(head(sort(table(EntrezResults$TaxID), decreasing = TRUE), n = 12))[7:12]
zzz1 <- unname(head(sort(table(EntrezResults$TaxID), decreasing = TRUE), n = 12))[7:12]
TitleVector <- c("M. tuberculosis",
                 "L. monocytogenes",
                 "E. faecium",
                 "S. enterica",
                 "C. difficile",
                 "N. meningitidis")
for (m1 in seq_along(zz1)) {
  w1 <- EntrezResults$TaxID == zz1[m1]
  w2 <- names(table(a1[w1])[table(a1[w1]) > 100])
  w3 <- names(table(a2[w1])[table(a2[w1]) > 100])
  w4 <- w5 <- w6 <- w9 <- matrix(data = 0,
                                 nrow = length(w2),
                                 ncol = length(w3))
  for (m2 in seq_along(w2)) {
    for (m3 in seq_along(w3)) {
      w2a <- a1 == w2[m2]
      w3a <- a2 == w3[m3]
      w4[m2, m3] <- sum(PseudosCondensed[w1 & w2a & w3a, 1L])
      w5[m2, m3] <- sum(PseudosCondensed[w1 & w2a & w3a, 4L])
      w6[m2, m3] <- sum(EntrezResults$Total_Length[w1 & w2a & w3a])
      w9a <- gsub(x = w2[m2],
                  pattern = "\n+ ",
                  replacement = "; ")
      w9a <- gsub(x = w9a,
                  pattern = "Unknown_Other_Uncaptured",
                  replacement = "Other")
      w9b <- gsub(x = w3[m3],
                  pattern = "\n+ ",
                  replacement = "; ")
      w9b <- gsub(x = w9b,
                  pattern = "Unknown_Other_Uncaptured",
                  replacement = "Other")
      w9[m2, m3] <- paste(w9a,
                          w9b,
                          sep = "\n+ ")
    }
  }
  # limit to 100 million bp combos
  w7 <- w4[w6 > 100000000] / w6[w6 > 100000000] # frameshifts
  # shift to pseudogenes per kbp -- was pseudogenes per bp
  w7 <- w7 * 1000
  w8 <- w5[w6 > 100000000] / w6[w6 > 100000000] # internal stops
  # shift to pseudogenes per kbp -- was pseudogenes per bp
  w8 <- w8 * 1000
  w9 <- w9[w6 > 100000000]
  plot(x = w7,
       y = w8,
       xlim = c(0, 0.04),
       ylim = c(0, 0.04),
       xlab = if (m1 %in% c(2,4,6)) {
         "frameshifts per kb"
       } else {
         NA
       },
       ylab = if (m1 %in% 1:2) {
         "internal stops per kb"
       } else {
         NA
       },
       main = paste0(TitleVector[m1],
                     " (",
                     zzz1[m1],
                     ")"))
  # text(x = w7,
  #      y = w8,
  #      labels = w9)
}
dev.off()


# plot(x = EntrezResults$Coverage,
#      y = z1,
#      xlim = c(0, 5000),
#      pch = 46,
#      col = "#00000001",
#      main = "frameshifts")
# abline(lm(z1 ~ EntrezResults$Coverage))
# plot(x = EntrezResults$Coverage,
#      y = z3,
#      xlim = c(0, 5000),
#      pch = 46,
#      main = "partial")
# plot(x = EntrezResults$Coverage,
#      y = z4,
#      xlim = c(0, 5000),
#      pch = 46,
#      main = "internal stop")
# plot(x = EntrezResults$Coverage[EntrezResults$Coverage > 0],
#      y = (z1 + z4 + z3)[EntrezResults$Coverage > 0],
#      xlim = c(10, 5000),
#      pch = 46,
#      main = "all",
#      log = "x")
# 
# plot(x = EntrezResults$Coverage,
#      y = z2,
#      xlim = c(0, 5000),
#      pch = 46,
#      main = "assembly gaps")
# 
# plot(x = EntrezResults$Total_Length[EntrezResults$Coverage > 0],
#      y = (z1 + z4 + z3)[EntrezResults$Coverage > 0],
#      pch = 46)
