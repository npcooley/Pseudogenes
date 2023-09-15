


load(file = "~/Repos/PseudogenePlots/InputData/MonoResults.RData",
     verbose = TRUE)

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

###### -- code body -----------------------------------------------------------

PseudoTypeVector <- c("frameshifted",
                      "internal stop",
                      "partial abbutting assembly gap",
                      "partial in the middle",
                      "missing C-terminus",
                      "missing N-terminus")


REFDNA <- readDNAStringSet("~/Repos/20230302_SimulatedMonoReassemblies/GCF_000005845.2_ASM584v2_genomic.fna")
REFGFF <- gffToDataFrame(GFF = "~/Repos/20230302_SimulatedMonoReassemblies/reference.gff",
                         Verbose = TRUE)
REFGR <- rtracklayer::import("~/Repos/20230302_SimulatedMonoReassemblies/reference.gff")

REFNotes <- GRangeToDFrame(GRangeObject = REFGR,
                           Verbose = TRUE)
REFCounts <- NoteCheck(NoteVector = REFNotes$Note,
                       CheckVector = PseudoTypeVector)


datIS <- data.frame("ISperMB" = (res[, 1L] / res[, 4L]) * 1000000,
                    # "FSperMB" = (res[, 2L] / res[, 4L]) * 1000000,
                    "Platform" = chemvec,
                    "Cov" = coveragevec,
                    "Quality" = qualityvec)
datIS <- datIS[res[, 4] > 4e6, ]
datIS <- datIS[datIS$Platform %in% c("HS25",
                                     "HSXt",
                                     "MSv3"), ]
datIS$Cov <- log10(datIS$Cov)
datIS$Platform <- factor(datIS$Platform)
datIS$ISperMB <- datIS$ISperMB - ((REFCounts[2] / sum(width(REFDNA))) * 1000000)

datFS <- data.frame("FSperMB" = (res[, 2L] / res[, 4L]) * 1000000,
                    # "FSperMB" = (res[, 2L] / res[, 4L]) * 1000000,
                    "Platform" = chemvec,
                    "Cov" = coveragevec,
                    "Quality" = qualityvec)
datFS <- datFS[res[, 4] > 4e6, ]
datFS <- datFS[datFS$Platform %in% c("HS25",
                                     "HSXt",
                                     "MSv3"), ]
datFS$Cov <- log10(datFS$Cov)
datFS$Platform <- factor(datFS$Platform)
datFS$FSperMB <- datFS$FSperMB - ((REFCounts[1] / sum(width(REFDNA))) * 1000000)

# datIS <- data.frame("ISperMB" = datIS$ISperMB - ((REFCounts[2] / sum(width(REFDNA))) * 1000000),
#                     "Plat_MSv3" = as.integer(grepl(pattern = "MSv3",
#                                                    x = datIS$Platform)),
#                     "Plat_HS25" = as.integer(grepl(pattern = "HS25",
#                                                    x = datIS$Platform)),
#                     "Plat_HSXt" = as.integer(grepl(pattern = "HSXt",
#                                                    x = datIS$Platform)),
#                     "Cov" = datIS$Cov,
#                     "Quality" = datIS$Quality)


save(datIS,
     datFS,
     file = "~/Repos/PseudogenePlots/InputData/ModelTables.RData",
     compress = "xz")

mod1 <- lm(formula = ISperMB ~ . + 0,
            data = datIS)
z1 <- coefficients(mod1)
z2 <- confint(mod1)
# confint(mod1)
# coef(summary(mod1))
# symnum(x = coef(summary(mod1))[, "Pr(>|t|)"], cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))

mod2 <- lm(formula = FSperMB ~ . + 0,
           data = datFS)
z3 <- coefficients(mod2)
z4 <- confint(mod2)

# plot coefficients for variables with a bar for the confidence interval

par(mar = c(3.5,2.5,2.5,1),
    mgp = c(1.5,0.5,0))
plot(x = 0,
     y = 0,
     type = "n",
     xlim = c(1, 10),
     ylim = c(-.175, .175),
     ylab = "Coefficient",
     xlab = "",
     xaxt = "n",
     yaxt = "n")
mtext(side = 3,
      at = 3,
      line = 0.25,
      cex = 1,
      text = "Internal Stops")
mtext(side = 3,
      at = 8,
      line = 0.25,
      cex = 1,
      text = "Frameshifts")
axis(side = 1,
     at = 1:10,
     labels = rep(c("HS25",
                    "HSXt",
                    "MSv3",
                    "Cov.",
                    "Qual."),
                  2),
     las = 2)
axis(side = 2,
     at = c(-.1, 0, .1),
     labels = c(-.1, 0, .1))
points(x = seq(length(z1)),
       y = z1,
       pch = 20,
       col = ColVec1[1])
segments(x0 = seq(nrow(z2)),
         x1 = seq(nrow(z2)),
         y0 = z2[, 1L],
         y1 = z2[, 2L],
         lty = 2)
points(x = seq(length(z3)) + 5,
       y = z3,
       pch = 20,
       col = ColVec1[2])
segments(x0 = seq(nrow(z4)) + 5,
         x1 = seq(nrow(z4)) + 5,
         y0 = z4[, 1L],
         y1 = z4[, 2L],
         lty = 2)
# 
# plot(z3[!is.na(z3)],
#      ylim = c(-.2, .1))
# segments(x0 = seq(nrow(z2)),
#          x1 = seq(nrow(z2)),
#          y0 = z2[rownames(z2) %in% z4, 1L],
#          y1 = z2[rownames(z2) %in% z4, 2L],
#          lty = 2)
# points(z2[rownames(z2) %in% z4, 1L], pch = 20)
# points(z2[rownames(z2) %in% z4, 2L], pch = 20)
