###### -- parse set comparison results ----------------------------------------

suppressMessages(library(SynExtend))


ColVec1 <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000')
ColVec2 <- paste0(ColVec1,
                  "33")
ColVec3 <- paste0(ColVec1,
                  "15")

# expects arguments in the format of 
# folder location of set comparison results
# file location of job map
# file location of original entrez search results
# file location for the output
ARGS <- commandArgs(trailingOnly = TRUE)

if (length(ARGS) == 0L) {
  stop("when run from the terminal this script expects trailing arguments")
}

# ARGS <- c("~/Data/20230902_SetComparison",
#           "~/Repos/20230814_SetComparison/JobMap01.txt",
#           "~/Repos/20230814_SetComparison/SearchResults.RData",
#           "~/Repos/PseudogenePlots/InputData/Counts_Orthos_v02.RData")

load(file = ARGS[3L],
     verbose = TRUE)
jobmap <- read.table(file = ARGS[2L])

path01 <- ARGS[1L]
files01 <- list.files(path = path01)
completedjobs <- as.integer(unlist(regmatches(x = files01,
                                              m = gregexpr(pattern = "[0-9]+",
                                                           text = files01))))
# jobs don't always complete neatly on the grid, accountting for that,
# build out a table of the newer refseq gene counts for the genomes that are present in the completed set
completedjobmap <- jobmap[completedjobs, ]
u1 <- unique(c(completedjobmap$V2, completedjobmap$V3))
newcounts <- vector(mode = "list",
                    length = length(u1))
names(newcounts) <- u1

current_lengths <- lengths(newcounts)
currently_missing <- current_lengths == 0
count <- 1L
while(any(currently_missing)) {
  u2 <- which(currently_missing)
  if (length(u2) > 0) {
    pin1 <- names(u2[1])
    pin2 <- as.integer(pin1)
    # print(pin1)
    # print(pin2)
    u3 <- which(completedjobmap$V2 == pin2)
    if (length(u3) > 0) {
      u3 <- u3[1L]
    } else {
      u3 <- which(completedjobmap$V3 == pin2)
      if (length(u3) > 0) {
        u3 <- u3[1L]
      } else {
        stop("what?")
      }
    }
    # print(u3)
    load(file = paste0(path01,
                       "/",
                       files01[match(x = completedjobmap$V4[u3],
                                     table = completedjobs)]),
         verbose = FALSE)
    
    print(w2)
    input_val1 <- which(names(newcounts) == as.character(w2[1]))
    input_val2 <- which(names(newcounts) == as.character(w2[2]))
    nt1 <- EntrezResults[w2[1], 12]
    nt2 <- EntrezResults[w2[2], 12]
    if (length(newcounts[[input_val1]]) == 0) {
      t1 <- table(geneskey1)
      if ("gene" %in% names(t1)) {
        c1 <- unname(t1[names(t1) == "gene"])
      } else {
        c1 <- 0
      }
      if ("frameshift" %in% names(t1) & "both" %in% names(t1)) {
        c2 <- unname(t1[names(t1) == "frameshift"] + t1[names(t1) == "both"])
      } else if ("frameshift" %in% names(t1) & !("both" %in% names(t1))) {
        c2 <- unname(t1[names(t1) == "frameshift"])
      } else if (!("frameshift" %in% names(t1)) & "both" %in% names(t1)) {
        c2 <- unname(t1[names(t1) == "both"])
      } else {
        c2 <- 0
      }
      if ("internal_stop" %in% names(t1) & "both" %in% names(t1)) {
        c3 <- unname(t1[names(t1) == "internal_stop"] + t1[names(t1) == "both"])
      } else if ("internal_stop" %in% names(t1) & !("both" %in% names(t1))) {
        c3 <- unname(t1[names(t1) == "internal_stop"])
      } else if (!("internal_stop" %in% names(t1)) & "both" %in% names(t1)) {
        c3 <- unname(t1[names(t1) == "both"])
      } else {
        c3 <- 0
      }
      if ("other" %in% names(t1)) {
        c4 <- unname(t1[names(t1) == "other"])
      } else {
        c4 <- 0
      }
      c5 <- sum(t1)
      newcounts[[input_val1]] <- c(c1, c2, c3, c4, c5, nt1)
    }
    
    if (length(newcounts[[input_val2]]) == 0) {
      t1 <- table(geneskey2)
      if ("gene" %in% names(t1)) {
        c1 <- unname(t1[names(t1) == "gene"])
      } else {
        c1 <- 0
      }
      if ("frameshift" %in% names(t1) & "both" %in% names(t1)) {
        c2 <- unname(t1[names(t1) == "frameshift"] + t1[names(t1) == "both"])
      } else if ("frameshift" %in% names(t1) & !("both" %in% names(t1))) {
        c2 <- unname(t1[names(t1) == "frameshift"])
      } else if (!("frameshift" %in% names(t1)) & "both" %in% names(t1)) {
        c2 <- unname(t1[names(t1) == "both"])
      } else {
        c2 <- 0
      }
      if ("internal_stop" %in% names(t1) & "both" %in% names(t1)) {
        c3 <- unname(t1[names(t1) == "internal_stop"] + t1[names(t1) == "both"])
      } else if ("internal_stop" %in% names(t1) & !("both" %in% names(t1))) {
        c3 <- unname(t1[names(t1) == "internal_stop"])
      } else if (!("internal_stop" %in% names(t1)) & "both" %in% names(t1)) {
        c3 <- unname(t1[names(t1) == "both"])
      } else {
        c3 <- 0
      }
      if ("other" %in% names(t1)) {
        c4 <- unname(t1[names(t1) == "other"])
      } else {
        c4 <- 0
      }
      c5 <- sum(t1)
      newcounts[[input_val2]] <- c(c1, c2, c3, c4, c5, nt2)
    }
    
  } else {
    break
  }
  current_lengths <- lengths(newcounts)
  currently_missing <- current_lengths == 0
  print(count)
  count <- count + 1L
}

adjusted_counts <- do.call(rbind,
                           newcounts)
adjusted_counts <- data.frame("regular_coding" = adjusted_counts[, 1L],
                              "frameshifts" = adjusted_counts[, 2L],
                              "internal_stops" = adjusted_counts[, 3L],
                              "other" = adjusted_counts[, 4L],
                              "all_coding" = adjusted_counts[, 5L],
                              "total_length" = adjusted_counts[, 6L],
                              row.names = rownames(adjusted_counts))

PBAR <- length(files01)
pBar <- txtProgressBar(style = 1)
res01 <- vector(mode = "list",
                length = PBAR)

TSTART <- Sys.time()
for (m1 in seq_along(files01)) {
  load(file = paste0(path01,
                     "/",
                     files01[m1]),
       verbose = FALSE)
  
  # build out a table of:
  # ani results
  # pair identifiers
  # congruent groups/pairs
  # incongruent groups/pairs
  # total groups
  # total pairs
  pairedonly <- lengths(clust7) == 2L
  res01[[m1]] <- c(ANIRes,
                   w2,
                   sum(congruent_fs),
                   sum(congruent_is),
                   sum(congruent_fs & pairedonly),
                   sum(congruent_is & pairedonly),
                   sum(incongruent_fs),
                   sum(incongruent_is),
                   sum(incongruent_fs & pairedonly),
                   sum(incongruent_is & pairedonly),
                   length(clust7),
                   sum(pairedonly))
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)
TEND <- Sys.time()
print(TEND - TSTART)

res02 <- do.call(rbind,
                 res01)

dat1 <- data.frame("ANI" = apply(X = res02[, 1:2],
                                 MARGIN = 1,
                                 FUN = function(x) {
                                   mean(x)
                                 }),
                   "AF" = apply(X = res02[, 3:4],
                                MARGIN = 1,
                                FUN = function(x) {
                                  mean(x)
                                }),
                   "Congruent_FS" = as.integer(res02[, 7]),
                   "Congruent_IS" = as.integer(res02[, 8]),
                   "Congruent_FS_Pairs" = as.integer(res02[, 9]),
                   "Congruent_IS_Pairs" = as.integer(res02[, 10]),
                   "Incongruent_FS" = as.integer(res02[, 11]),
                   "Incongruent_IS" = as.integer(res02[, 12]),
                   "Incongruent_FS_Pairs" = as.integer(res02[, 13]),
                   "Incongruent_IS_Pairs" = as.integer(res02[, 14]),
                   "AllGroups" = as.integer(res02[, 15]),
                   "AllPairs" = as.integer(res02[, 16]),
                   "id1" = as.integer(res02[, 5]),
                   "id2" = as.integer(res02[, 6]))

# subset to best matches
u4 <- sort(unique(c(dat1$id1, dat1$id2)))
pBar <- txtProgressBar(style = 1)
PBAR <- length(u4)
keepset <- vector(mode = "integer",
                  length = PBAR)
TSTART <- Sys.time()
for (m1 in seq_along(u4)) {
  u5 <- which(dat1$id1 == u4[m1] | dat1$id2 == u4[m1])
  u6 <- which.max(dat1$ANI[u5])
  keepset[m1] <- u5[u6]
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)
TEND <- Sys.time()
print(TEND - TSTART)

dat2 <- dat1[keepset, ]

z1 <- dat2$Incongruent_FS / dat2$Congruent_FS
z1 <- z1[!is.infinite(z1) & !is.na(z1)]
z2 <- dat2$Incongruent_IS / dat2$Congruent_IS
z2 <- z2[!is.infinite(z2) & !is.na(z2)]
hist(x = z1[z1 <= 11],
     breaks = seq(from = 0,
                  to = 11,
                  by = 1),
     right = FALSE,
     xlim = c(0, 10))
hist(x = z2[z2 <= 11],
     breaks = seq(from = 0,
                  to = 11,
                  by = 1),
     right = FALSE,
     xlim = c(0, 10))

plot(x = dat2$ANI,
     xlab = "Average Nucleotide Identity (%)",
     xlim = c(90, 100),
     y = dat2$Incongruent_FS_Pairs / dat2$AllPairs,
     ylab = "Incongruent Frameshifts / all pairs",
     ylim = c(0, 0.01),
     pch = 16,
     xaxs = "i",
     cex = 0.5,
     col = "#11227733")
spfit1 <- smooth.spline(x = dat2$ANI[dat2$ANI >= 90],
                        y = (dat2$Incongruent_FS_Pairs / dat2$AllPairs)[dat2$ANI >= 90],
                        df = 5)
spfit3 <- predict(spfit1)
lines(x = spfit3$x,
      y = spfit3$y,
      col = "orange")

plot(x = dat2$ANI,
     xlab = "Average Nucleotide Identity (%)",
     xlim = c(90, 100),
     y = dat2$Incongruent_IS_Pairs / dat2$AllPairs,
     ylab = "Incongruent Internal Stops / all pairs",
     ylim = c(0, 0.01),
     pch = 16,
     xaxs = "i",
     cex = 0.5,
     col = "#77221133")
spfit1 <- smooth.spline(x = dat2$ANI[dat2$ANI >= 90],
                        y = (dat2$Incongruent_IS_Pairs / dat2$AllPairs)[dat2$ANI >= 90],
                        df = 5)
spfit3 <- predict(spfit1)
lines(x = spfit3$x,
      y = spfit3$y,
      col = "orange")

hist((dat2$Incongruent_IS_Pairs / dat2$AllPairs)[dat2$ANI >= 90],
     breaks = 1000,
     xlim = c(0, 0.01))

plot(x = (dat2$Congruent_FS_Pairs / dat2$AllPairs)[dat2$ANI >= 90],
     xlab = "Relative congruent frameshift pairs",
     y = (dat2$Incongruent_FS_Pairs / dat2$AllPairs)[dat2$ANI >= 90],
     ylab = "Relative incongruent frameshift pairs",
     pch = 46,
     xaxs = "i",
     yaxs = "i")
abline(a = 0, b = 1, col = "orange")

plot(x = (dat2$Congruent_IS_Pairs / dat2$AllPairs)[dat2$ANI >= 90],
     xlab = "Relative congruent internal stop pairs",
     y = (dat2$Incongruent_IS_Pairs / dat2$AllPairs)[dat2$ANI >= 90],
     ylab = "Relative incongruent internal stop pairs",
     pch = 46,
     xaxs = "i",
     yaxs = "i")
abline(a = 0, b = 1, col = "orange")

# plot(x = dat2$ANI,
#      xlab = "Average Nucleotide Identity (%)",
#      y = dat2$Incongruent_FS / dat2$Congruent_FS,
#      ylab = "Incongruent / Congruent",
#      pch = 16,
#      # ylim = c(0, 10),
#      # xlim = c(90, 100),
#      xaxs = "i",
#      col = "#11227733",
#      cex = 0.5)
# abline(h = 1, lty = 2, col = "orange")
# 
# 
# plot(x = dat2$ANI,
#      xlab = "Average Nucleotide Identity (%)",
#      y = dat2$Incongruent_IS / dat2$Congruent_IS,
#      ylab = "Incongruent / Congruent",
#      pch = 16,
#      ylim = c(0, 10),
#      xlim = c(90, 100),
#      xaxs = "i",
#      col = "#77221133",
#      cex = 0.5)
# abline(h = 1, lty = 2, col = "orange")
# hist(dat2$Incongruent_IS / dat2$Congruent_IS, xlim = c(0, 10), breaks = 100)

# map new IS and FS counts to original lengths, and get the differences in PG per Mbp
mat1 <- match(x = dat2$id1,
              table = as.integer(rownames(adjusted_counts)))
mat2 <- match(x = dat2$id2,
              table = as.integer(rownames(adjusted_counts)))
dat3 <- cbind(dat2,
              "l1" = adjusted_counts$total_length[mat1],
              "fs1" = adjusted_counts$frameshifts[mat1],
              "is1" = adjusted_counts$internal_stops[mat1],
              "l2" = adjusted_counts$total_length[mat2],
              "fs2" = adjusted_counts$frameshifts[mat2],
              "is2" = adjusted_counts$internal_stops[mat2])

pdf(file = "testplot.pdf",
    width = 7,
    height = 7)

layout(mat = matrix(data = 1:6,
                    nrow = 3,
                    byrow = TRUE))

plot(x = (dat3$Congruent_FS_Pairs / dat3$AllPairs)[dat2$ANI >= 90],
     xlab = "Relative congruent frameshift pairs",
     y = (dat3$Incongruent_FS_Pairs / dat3$AllPairs)[dat2$ANI >= 90],
     ylab = "Relative incongruent frameshift pairs",
     pch = 16,
     cex = 0.5,
     xaxs = "i",
     yaxs = "i",
     col = ColVec2[1L])
abline(a = 0, b = 1, col = "black")

plot(x = (dat2$Congruent_IS_Pairs / dat2$AllPairs)[dat2$ANI >= 90],
     xlab = "Relative congruent internal stop pairs",
     y = (dat2$Incongruent_IS_Pairs / dat2$AllPairs)[dat2$ANI >= 90],
     ylab = "Relative incongruent internal stop pairs",
     pch = 16,
     cex = 0.5,
     xaxs = "i",
     yaxs = "i",
     col = ColVec2[2L])
abline(a = 0, b = 1, col = "black")

plot(x = dat3$ANI,
     xlab = "Average Nucleotide Identity (%)",
     xlim = c(90, 100),
     y = dat3$Incongruent_FS_Pairs / dat3$AllPairs,
     ylab = "Incongruent Frameshifts / all pairs",
     ylim = c(0, 0.02),
     pch = 16,
     xaxs = "i",
     yaxs = "i",
     cex = 0.5,
     col = ColVec2[1L])
spfit1 <- smooth.spline(x = dat3$ANI[dat2$ANI >= 90],
                        y = (dat3$Incongruent_FS_Pairs / dat3$AllPairs)[dat3$ANI >= 90],
                        df = 5)
spfit3 <- predict(spfit1)
lines(x = spfit3$x,
      y = spfit3$y,
      col = "black")

plot(x = dat3$ANI,
     xlab = "Average Nucleotide Identity (%)",
     xlim = c(90, 100),
     y = dat3$Incongruent_IS_Pairs / dat3$AllPairs,
     ylab = "Incongruent Internal Stops / all pairs",
     ylim = c(0, 0.005),
     pch = 16,
     xaxs = "i",
     yaxs = "i",
     cex = 0.5,
     col = ColVec2[2L])
spfit1 <- smooth.spline(x = dat3$ANI[dat3$ANI >= 90],
                        y = (dat3$Incongruent_IS_Pairs / dat3$AllPairs)[dat2$ANI >= 90],
                        df = 5)
spfit3 <- predict(spfit1)
lines(x = spfit3$x,
      y = spfit3$y,
      col = "black")

plot(x = dat3$ANI[dat3$ANI >= 90],
     xlab = "Average Nucleotide Identity (%)",
     y = (abs((dat3$fs1 / dat3$l1) - (dat3$fs2 / dat3$l2)) * 1000000)[dat3$ANI >= 90],
     ylab = expression("|"*Delta*" Frameshifts per Mbp|"),
     ylim = c(0, 10),
     xaxs = "i",
     yaxs = "i",
     pch = 16,
     cex = 0.5,
     col = ColVec2[1L])
spfit1 <- smooth.spline(x = dat3$ANI[dat3$ANI >= 90],
                        y = (abs((dat3$fs1 / dat3$l1) - (dat3$fs2 / dat3$l2)) * 1000000)[dat3$ANI >= 90],
                        df = 5)
spfit3 <- predict(spfit1)
lines(x = spfit3$x,
      y = spfit3$y,
      col = "black")

plot(x = dat3$ANI[dat3$ANI >= 90],
     xlab = "Average Nucleotide Identity (%)",
     y = (abs((dat3$is1 / dat3$l1) - (dat3$is2 / dat3$l2)) * 1000000)[dat3$ANI >= 90],
     ylab = expression("|"*Delta*" Internal stops per Mbp|"),
     ylim = c(0, 4),
     xaxs = "i",
     yaxs = "i",
     pch = 16,
     cex = 0.5,
     col = ColVec2[2L])
spfit1 <- smooth.spline(x = dat3$ANI[dat3$ANI >= 90],
                        y = (abs((dat3$is1 / dat3$l1) - (dat3$is2 / dat3$l2)) * 1000000)[dat3$ANI >= 90],
                        df = 5)
spfit3 <- predict(spfit1)
lines(x = spfit3$x,
      y = spfit3$y,
      col = "black")
dev.off()

save(adjusted_counts,
     dat3,
     file = ARGS[4],
     compress = "xz")


