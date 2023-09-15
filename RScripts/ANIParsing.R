###### -- parse ANI results ---------------------------------------------------

# requires parsing the old results together as well because IDs weren't carried through
# in a meaningful way

suppressMessages(library(SynExtend))
library(stringr)

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

FILES01 <- list.files(path = "~/Data/20230509_SimulatedEcoliReassemblies",
                      pattern = "Annotation")
FILES02 <- list.files(path = "~/Data/20230509_SimulatedEcoliReassemblies",
                      pattern = "Assembly")
FILES03 <- list.files(path = "~/Data/20230509_SimulatedEcoliReassemblies",
                      pattern = "Result")
FILES04 <- list.files(path = "~/Data/20230509_SimulatedEcoliReassemblies",
                      pattern = "ANI")

JobMapB <- read.table("~/Data/20230509_SimulatedEcoliReassemblies/JobMapB.txt")
JobMapA <- read.table("~/Repos/20230404_SimulatedEColiReassemblies/JobMapA.txt")

PseudoTypeVector <- c("frameshifted",
                      "internal stop",
                      "partial abbutting assembly gap",
                      "partial in the middle",
                      "missing C-terminus",
                      "missing N-terminus")

RefAssembly <- readDNAStringSet(filepath = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz")
RefGC <- gffToDataFrame(GFF = "~/Data/20230509_SinglePGAP/annot.gff",
                        Verbose = TRUE)
RefGR <- rtracklayer::import("~/Data/20230509_SinglePGAP/annot.gff")
RefDF <- GRangeToDFrame(GRangeObject = RefGR,
                        Verbose = TRUE)
RefNotes <- NoteCheck(NoteVector = RefDF$Note,
                      CheckVector = PseudoTypeVector)
RefIS <- RefNotes[2L] / (sum(width(RefAssembly)) / 1000000)
RefFS <- RefNotes[1L] / (sum(width(RefAssembly)) / 1000000)

w1 <- as.integer(unlist(regmatches(x = FILES03,
                                   m = gregexpr(pattern = "[0-9]+",
                                                text = FILES03))))
w2 <- as.integer(unlist(regmatches(x = FILES04,
                                   m = gregexpr(pattern = "[0-9]+",
                                                text = FILES04))))
w3 <- w1 %in% w2
w4 <- w2 %in% w1

FILES05 <- FILES03[w3]
FILES06 <- FILES04[w4]

RES <- vector(mode = "list",
              length = length(FILES05))
MAP1 <- as.integer(unlist(regmatches(x = FILES05,
                                     m = gregexpr(pattern = "[0-9]+",
                                                  text = FILES05))))
MAP4 <- strsplit(x = JobMapB[, 2L],
                 split = "_",
                 fixed = TRUE)
MAP2 <- sapply(X = MAP4,
               FUN = function(x) {
                 x[2]
               })
MAP5 <- sapply(X = MAP4,
               FUN = function(x) {
                 x[3]
               })
MAP6 <- sapply(X = MAP4,
               FUN = function(x) {
                 if (length(x) == 5) {
                   x[4]
                 } else {
                   "None"
                 }
               })
MAP3 <- as.integer(unlist(regmatches(x = JobMapB[, 2L],
                                     m = gregexpr(pattern = "[0-9]+",
                                                  text = JobMapB[, 2L]))))

Ref_FS_Names <- paste(RefDF$Name[grepl(pattern = "frameshifted", x = RefDF$Note)],
                      "B",
                      sep = "_")
Ref_IS_Names <- paste(RefDF$Name[grepl(pattern = "internal stop", x = RefDF$Note)],
                      "B",
                      sep = "_")

pBar <- txtProgressBar(style = 1L)
PBAR <- length(RES)
for (m1 in seq_along(FILES05)) {
  load(file = paste0("~/Data/20230509_SimulatedEcoliReassemblies",
                     "/",
                     FILES05[m1]),
       verbose = FALSE)
  load(file = paste0("~/Data/20230509_SimulatedEcoliReassemblies",
                     "/",
                     FILES06[m1]),
       verbose = FALSE)
  z1 <- as.numeric(strsplit(x = x,
                            split = "\t",
                            fixed = TRUE)[[2]][3:6])
  z2 <- mean(z1[1:2])
  z3 <- mean(z1[3:4])
  
  # reference is B
  a <- c(res,
         list("assembler" = MAP2[MAP1[m1]]),
         list("arrangement" = MAP5[MAP1[m1]]),
         list("lr" = MAP6[MAP1[m1]]),
         list("total_genes" = sum(CurrentGC$Coding)),
         list("count1" = res$rel_fs * (res$total_length / 1000000)),
         list("count2" = res$rel_is * (res$total_length / 1000000)),
         list("found_fs" = sum(Ref_FS_Names %in% Cl1$p2)),
         list("found_is" = sum(Ref_IS_Names %in% Cl1$p2)),
         JobMapA[MAP3[MAP1[m1]], ],
         list("ANI" = z2),
         list("AF" = z3))
  a <- data.frame(a)
  
  RES[[m1]] <- a
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)

dat1 <- do.call(rbind,
                RES)

dat2 <- data.frame("fs" = dat1$rel_fs - RefFS,
                   "is" = dat1$rel_is - RefIS,
                   "fragment" = dat1$n50 / dat1$total_length,
                   "total_length" = dat1$total_length,
                   "total_contigs" = dat1$tc,
                   # "cov" = log10(dat1$V3) - log10(100),
                   "cov" = log10(dat1$V3 / 50),
                   # "read_len" = log10(dat1$V4) - log10(100),
                   # "read_len" = log10(dat1$V4 / 100),
                   "read_len" = dat1$V4,
                   "assembler" = as.factor(dat1$assembler),
                   "qual" = dat1$V2,
                   "arrangement" = as.factor(paste(dat1$arrangement,
                                                   dat1$lr,
                                                   sep = "+")),
                   "sequencer" = as.factor(dat1$V5),
                   "total_genes" = dat1$total_genes,
                   "abs_fs" = dat1$count1,
                   "abs_is" = dat1$count2,
                   "found_ref_fs" = dat1$found_fs,
                   "found_ref_is" = dat1$found_is,
                   "ANI" = dat1$ANI,
                   "AF" = dat1$AF)
dat2 <- dat2[dat2$read_len == 150, ]
# subset out assemblies that are very wrong in length
dat3 <- dat2[dat2$total_length >= (sum(width(RefAssembly)) * 0.8) &
               dat2$total_length <= (sum(width(RefAssembly)) * 1.2) &
               dat2$total_genes >= (sum(RefGC$Coding) * 0.8) &
               dat2$total_genes <= (sum(RefGC$Coding) * 1.2), ]

# some vals
dat4 <- data.frame("ANI" = c(cor(x = dat3$ANI,
                                 y = abs(dat3$abs_fs - dat3$found_ref_fs) / dat3$total_genes)^2,
                             cor(x = dat3$ANI,
                                 y = abs(dat3$abs_is - dat3$found_ref_is) / dat3$total_genes)^2,
                             cor(x = dat3$ANI,
                                 y = dat3$total_contigs)^2,
                             cor(x = dat3$ANI,
                                 y = dat3$fragment)^2),
                   "AF" = c(cor(x = dat3$AF,
                                y = abs(dat3$abs_fs - dat3$found_ref_fs) / dat3$total_genes)^2,
                            cor(x = dat3$AF,
                                y = abs(dat3$abs_is - dat3$found_ref_is) / dat3$total_genes)^2,
                            cor(x = dat3$AF,
                                y = dat3$total_contigs)^2,
                            cor(x = dat3$AF,
                                y = dat3$fragment)^2),
                   row.names = c("Delta_Rel_FS",
                                 "Delta_Rel_IS",
                                 "Total_Contigs",
                                 "Norm_N50"))

save(dat3,
     dat4,
     file = "~/Repos/PseudogenePlots/InputData/SimulatedEcoliANIResults.RData",
     compress = "xz")
# 
# pdf(file = "~/tempplot01.pdf",
#     height = 7,
#     width = 7)
# layout(mat = matrix(data = 1:4,
#                     ncol = 2))
# plot(x = dat3$ANI,
#      y = dat3$fs,
#      ylim = c(-1, 10),
#      pch = match(table = unique(dat3$assembler),
#                  x = dat3$assembler))
# plot(x = dat3$ANI,
#      y = dat3$is,
#      ylim = c(-1, 10),
#      pch = match(table = unique(dat3$assembler),
#                  x = dat3$assembler))
# plot(x = dat3$AF,
#      y = dat3$fs,
#      ylim = c(-1, 10),
#      pch = match(table = unique(dat3$assembler),
#                  x = dat3$assembler))
# plot(x = dat3$AF,
#      y = dat3$is,
#      ylim = c(-1, 10),
#      pch = match(table = unique(dat3$assembler),
#                  x = dat3$assembler))
# dev.off()
# 
# pdf(file = "~/tempplot02.pdf",
#     height = 7,
#     width = 7)
# layout(mat = matrix(data = 1:4,
#                     ncol = 2))
# plot(x = dat3$ANI,
#      y = dat3$fs,
#      pch = match(table = unique(dat3$assembler),
#                  x = dat3$assembler))
# plot(x = dat3$ANI,
#      y = dat3$is,
#      pch = match(table = unique(dat3$assembler),
#                  x = dat3$assembler))
# plot(x = dat3$AF,
#      y = dat3$fs,
#      pch = match(table = unique(dat3$assembler),
#                  x = dat3$assembler))
# legend("topright",
#        legend = unique(dat3$assembler),
#        pch = seq(length(dat3$assembler)),
#        cex = 0.75)
# plot(x = dat3$AF,
#      y = dat3$is,
#      pch = match(table = unique(dat3$assembler),
#                  x = dat3$assembler))
# dev.off()
# 
# pdf(file = "~/tempplot03.pdf",
#     height = 3.5,
#     width = 7)
# layout(mat = matrix(data = 1:2,
#                     nrow = 1))
# plot(density(dat3$ANI))
# plot(density(dat3$AF, adjust = 4))
# dev.off()
# 
# plot(x = dat3$ANI,
#      y = dat3$AF)
# 
# plot(x = 0,
#      y = 0,
#      xlim = range(dat3$ANI),
#      ylim = range(dat3$fs),
#      type = "n")
# points(x = dat3$ANI[as.character(dat3$assembler) == "UNICYCLER"],
#        y = dat3$fs[as.character(dat3$assembler) == "UNICYCLER"],
#        pch = 20,
#        col = ifelse(test = grepl(pattern = "PB|ONT",
#                                  x = dat3$arrangement[as.character(dat3$assembler) == "UNICYCLER"]),
#                     yes = "#22222250",
#                     no = "#99550050"))
# points(x = dat3$ANI[as.character(dat3$assembler) == "MEGAHIT"],
#        y = dat3$fs[as.character(dat3$assembler) == "MEGAHIT"],
#        pch = 20,
#        col = "#77755533")
# points(x = dat3$ANI[as.character(dat3$assembler) == "SPADES"],
#        y = dat3$fs[as.character(dat3$assembler) == "SPADES"])
# points(x = dat3$ANI[as.character(dat3$assembler) == "SKESA"],
#        y = dat3$fs[as.character(dat3$assembler) == "SKESA"])
# 
# plot(x = 0,
#      y = 0,
#      xlim = range(dat3$ANI),
#      ylim = range(dat3$fs),
#      # ylim = c(0, 50),
#      xlab = "ANI",
#      ylab = "PG per Mbp",
#      type = "n")
# points(x = dat3$ANI[as.character(dat3$assembler) == "UNICYCLER"],
#        y = dat3$fs[as.character(dat3$assembler) == "UNICYCLER"],
#        pch = ifelse(test = grepl(pattern = "PB|ONT",
#                                  x = dat3$arrangement[as.character(dat3$assembler) == "UNICYCLER"]),
#                     yes = 1,
#                     no = 4),
#        col = ifelse(test = grepl(pattern = "PB|ONT",
#                                  x = dat3$arrangement[as.character(dat3$assembler) == "UNICYCLER"]),
#                     yes = "#112299",
#                     no = "#995500"))
# points(x = dat3$ANI[as.character(dat3$assembler) == "UNICYCLER"],
#        y = dat3$is[as.character(dat3$assembler) == "UNICYCLER"],
#        pch = ifelse(test = grepl(pattern = "PB|ONT",
#                                  x = dat3$arrangement[as.character(dat3$assembler) == "UNICYCLER"]),
#                     yes = 3,
#                     no = 2),
#        cex = 1.5,
#        col = ifelse(test = grepl(pattern = "PB|ONT",
#                                  x = dat3$arrangement[as.character(dat3$assembler) == "UNICYCLER"]),
#                     yes = "#112299",
#                     no = "#995500"))
# 
# 
# plot(x = 0,
#      y = 0,
#      xlim = range(dat3$AF),
#      # ylim = range(dat3$fs),
#      # ylim = c(0, 50),
#      ylim = c(0, 20),
#      xlab = "AF",
#      ylab = "PG per Mbp",
#      type = "n")
# points(x = dat3$AF[as.character(dat3$assembler) != "UNICYCLER"],
#        y = abs(dat3$fs[as.character(dat3$assembler) != "UNICYCLER"]),
#        pch = c(5,6,7)[match(x = dat3$assembler[as.character(dat3$assembler) != "UNICYCLER"],
#                             table = unique(dat3$assembler[as.character(dat3$assembler) != "UNICYCLER"]))],
#        col = c("green", "violet", "tomato")[match(x = dat3$assembler[as.character(dat3$assembler) != "UNICYCLER"],
#                                                   table = unique(dat3$assembler[as.character(dat3$assembler) != "UNICYCLER"]))])
# points(x = dat3$AF[as.character(dat3$assembler) == "UNICYCLER"],
#        y = dat3$fs[as.character(dat3$assembler) == "UNICYCLER"],
#        pch = ifelse(test = grepl(pattern = "PB|ONT",
#                                  x = dat3$arrangement[as.character(dat3$assembler) == "UNICYCLER"]),
#                     yes = 1,
#                     no = 4),
#        col = ifelse(test = grepl(pattern = "PB|ONT",
#                                  x = dat3$arrangement[as.character(dat3$assembler) == "UNICYCLER"]),
#                     yes = "#112299",
#                     no = "#995500"))
# points(x = dat3$AF[as.character(dat3$assembler) == "UNICYCLER"],
#        y = dat3$is[as.character(dat3$assembler) == "UNICYCLER"],
#        pch = ifelse(test = grepl(pattern = "PB|ONT",
#                                  x = dat3$arrangement[as.character(dat3$assembler) == "UNICYCLER"]),
#                     yes = 3,
#                     no = 2),
#        cex = 1.5,
#        col = ifelse(test = grepl(pattern = "PB|ONT",
#                                  x = dat3$arrangement[as.character(dat3$assembler) == "UNICYCLER"]),
#                     yes = "#112299",
#                     no = "#995500"))
# 
# PBAR <- length(files01)
# pBar <- txtProgressBar(style = 1)
# ANIRes <- vector(mode = "list",
#                  length = PBAR)
# for (m1 in seq_along(files01)) {
#   load(file = files01[m1],
#        verbose = FALSE)
#   
#   z1 <- as.numeric(strsplit(x = x,
#                             split = "\t",
#                             fixed = TRUE)[[2]][3:6])
#   z2 <- mean(z1[1:2])
#   z3 <- mean(z1[3:4])
#   ANIRes[[m1]] <- c(z2, z3)
#   
#   setTxtProgressBar(pb = pBar,
#                     value = m1 / PBAR)
# }
# 
# ANIRes <- do.call(rbind,
#                   ANIRes)
# 
# plot(x = ANIRes[, 1L],
#      y = ANIRes[, 2L],
#      pch = 46,
#      xlim = c(95, 100),
#      ylim = c(0.3, 1.0),
#      xlab = "ANI",
#      ylab = "AF")
# 
# # something else
# ColVec1 <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000')
# ColVec2 <- paste0(ColVec1,
#                   "33")
# ColVec3 <- paste0(ColVec1,
#                   "15")
# 
# 
# # sorthing orders
# o1 <- c(1:4, 13, 16) # paired end
# o2 <- c(5:8, 14, 17)
# o3 <- c(9:12, 15, 18)
# o4 <- c(19:22, 31, 34) # single end
# o5 <- c(23:26, 32, 35)
# o6 <- c(27:30, 33, 36)
# # reorder columns
# o7 <- c(1,4,7,2,5,8,3,6,9)
# 
# # ani data parsing
# val1 <- by(data = dat3,
#            INDICES = ~ assembler + sequencer + arrangement,
#            FUN = function(x) {
#              summary(glm(formula = (100 - ANI) ~ cov + qual,
#                          data = x,
#                          family = binomial(link = "logit"),
#                          weights = total_genes))$coefficients
#            })
# #coefficients
# ani_vals01a <- apply(X = val1,
#                      MARGIN = 1:3,
#                      FUN = function(x) {
#                        x[[1]][, 1]
#                      })
# ani_vals01b <- do.call(rbind,
#                        ani_vals01a[, , 1:6])
# # paired end ani coefficients
# ani_vals01c <- cbind(ani_vals01b[o1, ],
#                      ani_vals01b[o2, ],
#                      ani_vals01b[o3, ])
# ani_vals01c <- ani_vals01c[, o7]
# # single end ani coefficients
# ani_vals01d <- cbind(ani_vals01b[o4, ],
#                      ani_vals01b[o5, ],
#                      ani_vals01b[o6, ])
# ani_vals01d <- ani_vals01d[, o7]
# # p-vals
# ani_vals02a <- apply(X = val1,
#                      MARGIN = 1:3,
#                      FUN = function(x) {
#                        x[[1]][, 4]
#                      })
# ani_vals02b <- do.call(rbind,
#                        ani_vals02a[, , 1:6])
# ani_vals02b <- matrix(data = p.adjust(p = ani_vals02b,
#                                       method = "bonferroni"),
#                       ncol = ncol(ani_vals02b))
# # paired end p-vals
# ani_vals02c <- cbind(ani_vals02b[o1, ],
#                      ani_vals02b[o2, ],
#                      ani_vals02b[o3, ])
# ani_vals02c <- ani_vals02c[, o7]
# # single end p-vals
# ani_vals02d <- cbind(ani_vals02b[o4, ],
#                      ani_vals02b[o5, ],
#                      ani_vals02b[o6, ])
# ani_vals02d <- ani_vals02d[, o7]
# 
# val2 <- by(data = dat3,
#            INDICES = ~ assembler + sequencer + arrangement,
#            FUN = function(x) {
#              summary(glm(formula = (1 - AF) ~ cov + qual,
#                          data = x,
#                          family = binomial(link = "logit"),
#                          weights = total_genes))$coefficients
#            })
# #coefficients
# af_vals01a <- apply(X = val2,
#                     MARGIN = 1:3,
#                     FUN = function(x) {
#                       x[[1]][, 1]
#                     })
# af_vals01b <- do.call(rbind,
#                       af_vals01a[, , 1:6])
# # paired end ani coefficients
# af_vals01c <- cbind(af_vals01b[o1, ],
#                     af_vals01b[o2, ],
#                     af_vals01b[o3, ])
# af_vals01c <- af_vals01c[, o7]
# # single end ani coefficients
# af_vals01d <- cbind(af_vals01b[o4, ],
#                     af_vals01b[o5, ],
#                     af_vals01b[o6, ])
# af_vals01d <- af_vals01d[, o7]
# # p-vals
# af_vals02a <- apply(X = val2,
#                     MARGIN = 1:3,
#                     FUN = function(x) {
#                       x[[1]][, 4]
#                     })
# af_vals02b <- do.call(rbind,
#                       af_vals02a[, , 1:6])
# af_vals02b <- matrix(data = p.adjust(p = af_vals02b,
#                                      method = "bonferroni"),
#                      ncol = ncol(af_vals02b))
# # paired end p-vals
# af_vals02c <- cbind(af_vals02b[o1, ],
#                     af_vals02b[o2, ],
#                     af_vals02b[o3, ])
# af_vals02c <- af_vals02c[, o7]
# # single end p-vals
# af_vals02d <- cbind(af_vals02b[o4, ],
#                     af_vals02b[o5, ],
#                     af_vals02b[o6, ])
# af_vals02d <- af_vals02d[, o7]
# 
# # simulated pins
# tot_genes <- 4379
# tot_nts <- 4641652
# cov_set <- seq(-1, 1, 0.01)
# qual_set <- seq(-4, 6, 0.01)
# 
# # given coefficients, predict |delta pseudogenes| per Mbp
# predictResponse <- function(coef,
#                             COV,
#                             QUAL,
#                             TOT_GENE,
#                             TOT_NT) {
#   Y <- coef[1] + coef[2]*COV + coef[3]*QUAL
#   Y <- exp(Y)
#   Y <- Y/(1 + Y) # |delta pseudogenes| / (total genes)
#   Y*TOT_GENE/TOT_NT*1e6 # |delta pseudogenes| / Mbp
# }
# convertCoverage <- function(cov) {
#   50*(10^cov) ### NEED TO VERIFY ### verified
# }
# convertQuality <- function(qual) {
#   qual + 34 ### NEED TO VERIFY ### still need to verify -- build a quick plot for this ...
# }
# 
# # ANI vs Coverage
# plot(x = 0,
#      y = 0,
#      type = "n",
#      xlim = c(5, 500),
#      ylim = c(0, 5),
#      xaxs = "i",
#      # xaxt = "n",
#      # yaxs = "i",
#      # xlab = "fold coverage",
#      main = "Paired End Reads",
#      xlab = "Coverage",
#      ylab = "100 - ANI",
#      log = "x")
# for (m1 in seq_len(nrow(ani_vals01c))) {
#   for (m2 in seq_len(ncol(ani_vals01c) / 3L)) {
#     c1 <- predictResponse(coef = ani_vals01c[m1, c(m2, m2 + 3L, m2 + 6L)],
#                           COV = cov_set,
#                           QUAL = 0,
#                           TOT_GENE = tot_genes,
#                           TOT_NT = tot_nts)
#     c2 <- convertCoverage(cov = cov_set)
#     lines(x = c2,
#           y = c1,
#           col = m1,
#           lty = m2)
#   }
# }
# 
# # ANI vs quality
# plot(x = 0,
#      y = 0,
#      type = "n",
#      xlim = c(30, 40),
#      ylim = c(0, 5),
#      xaxs = "i",
#      # xaxt = "n",
#      # yaxs = "i",
#      # xlab = "fold coverage",
#      main = "Paired End Reads",
#      xlab = "Quality",
#      ylab = "100 - ANI",
#      log = "x")
# for (m1 in seq_len(nrow(ani_vals01c))) {
#   for (m2 in seq_len(ncol(ani_vals01c) / 3L)) {
#     c1 <- predictResponse(coef = ani_vals01c[m1, c(m2, m2 + 3L, m2 + 6L)],
#                           COV = 0,
#                           QUAL = qual_set,
#                           TOT_GENE = tot_genes,
#                           TOT_NT = tot_nts)
#     c2 <- convertQuality(qual = qual_set)
#     lines(x = c2,
#           y = c1,
#           col = m1,
#           lty = m2)
#   }
# }
# L <- legend(x = 34,
#             y = 4.5,
#             legend = rep(NA, 18),
#             col = rep(ColVec1[seq_len(6)],
#                       3L),
#             lty = c(rep(1, 6),
#                     rep(2, 6),
#                     rep(3, 6)),
#             pch = rep(NA, 18),
#             ncol = 3L,
#             bty = 'n',
#             x.intersp = 0.5,
#             inset = 0.02,
#             cex = 0.75)
# legend(x = 34,
#        y = 4.5,
#        legend = c("MEGAHIT",
#                   "SKESA",
#                   "SPAdes",
#                   "Unicycler",
#                   "Unicycler + ONT",
#                   "Unicycler + PB"),
#        # col = rep(NA, 6),
#        # lty = rep(NA, 6),
#        ncol = 1,
#        x.intersp = 7.5,
#        bg = NA,
#        bty = "n",
#        cex = 0.75)
# text(x = c(34, 35, 36),
#      y = 4.525,
#      pos = 4,
#      xpd = TRUE,
#      cex = 0.75,
#      labels = rep(c("HS25",
#                     "HSXt",
#                     "MSv3"),
#                   3L),
#      srt = 45)
# 
# plot(x = 0,
#      y = 0,
#      type = "n",
#      xlim = c(5, 500),
#      ylim = c(0, 100),
#      xaxs = "i",
#      # xaxt = "n",
#      # yaxs = "i",
#      # xlab = "fold coverage",
#      xlab = "Quality",
#      ylab = "1 - AF",
#      log = "x")
# for (m1 in seq_len(nrow(af_vals01c))) {
#   for (m2 in seq_len(ncol(af_vals01c) / 3L)) {
#     c1 <- predictResponse(coef = af_vals01c[m1, c(m2, m2 + 3L, m2 + 6L)],
#                           COV = cov_set,
#                           QUAL = 0,
#                           TOT_GENE = tot_genes,
#                           TOT_NT = tot_nts)
#     c2 <- convertCoverage(cov = cov_set)
#     lines(x = c2,
#           y = c1,
#           col = m1,
#           lty = m2)
#   }
# }

###### -- just plot out ANI / AF vs N50,IS,FS,TC ------------------------------

pdf(file = "tempplot01.pdf",
    height = 7,
    width = 7)
layout(mat = matrix(data = 1:4,
                    ncol = 2L))
par(mar = c(2,4,4,1.5),
    mgp = c(2.25, 1, 0))
plot(x = dat3$ANI,
     y = abs(dat3$abs_fs - dat3$found_ref_fs) / dat3$total_genes,
     pch = 16,
     ylab = expression("|"*Delta*" Rel. Frameshifts from Reference|"),
     xlab = "",
     xaxt = "n",
     # yaxs = "i",
     col = "#11111133")
par(mar = c(5,4,1,1.5),
    mgp = c(2.25, 1, 0))
plot(x = dat3$ANI,
     y = dat3$total_contigs,
     pch = 16,
     ylab = "Total Contigs",
     xlab = "ANI",
     # yaxs = "i",
     col = "#11111133")
par(mar = c(2,3.5,4,2),
    mgp = c(2.25, 1, 0))
plot(x = dat3$ANI,
     y = abs(dat3$abs_is - dat3$found_ref_is) / dat3$total_genes,
     pch = 16,
     ylab = expression("|"*Delta*" Rel. Internal Stops from Reference|"),
     xlab = "",
     xaxt = "n",
     # yaxs = "i",
     col = "#11111133")
par(mar = c(5,3.5,1,2),
    mgp = c(2.25, 1, 0))
plot(x = dat3$ANI,
     y = dat3$fragment,
     pch = 16,
     ylab = "N50 Normalized to Total NTs",
     xlab = "ANI",
     # yaxs = "i",
     col = "#11111133")
dev.off()

