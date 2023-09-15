###### -- is there anything interesting in the metadata about technology ------

load(file = "~/Repos/Pseudogenes/ScrapedData.RData",
     verbose = TRUE)
load(file = "~/Repos/Pseudogenes/SearchResults.RData",
     verbose = TRUE)
load(file = "~/Repos/Pseudogenes/SearchResults2.RData",
     verbose = TRUE)

ColVec1 <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000')
ColVec2 <- paste0(ColVec1,
                  "33")

# create a CDF showing the relative shift in pseudogenizations per MB for technologies
# vs 1x ILLUMINA generations within taxIDs

# create a series of vectors for platform mixes associated with RefSeq sequences
# ILLUMINA 1x
# >= 2x ILLUMINA
# ILLUMINA + ONT
# ILLUMINA + PacBio
# ILLUMINA + Other
# etc

w1 <- EntrezResults$Biosample[duplicated(EntrezResults$Biosample)]
w2 <- EntrezResults$Biosample %in% SRAResults$BioSample
w3 <- EntrezResults$Coverage > 1

SET1 <- EntrezResults[!(EntrezResults$Biosample %in% w1) &
                        w2 &
                        w3, ]
SET2 <- PseudosCondensed[!(EntrezResults$Biosample %in% w1) &
                           w2 &
                           w3, ]
SET6 <- AssemblersLogical[!(EntrezResults$Biosample %in% w1) &
                            w2 &
                            w3, ]
SET7 <- apply(X = SET6,
              MARGIN = 1L,
              FUN = function(x) {
                paste(colnames(AssemblersLogical)[x],
                      collapse = " + ")
              },
              simplify = TRUE)

SET3 <- unname(sapply(X = SET1$Biosample,
                      FUN = function(x) {
                        paste(sort(unique(SRAResults$Platform[SRAResults$BioSample == x])),
                              collapse = " ")
                      }))
SET4 <- unname(sapply(X = SET1$Biosample,
                      FUN = function(x) {
                        sum(SRAResults$BioSample == x)
                      }))
SET5 <- unname(sapply(X = SET1$Biosample,
                      FUN = function(x) {
                        SRAResults$TaxID[SRAResults$BioSample == x][1]
                      }))

SET8 <- strsplit(SET1$Submission_Date,
                 split = "/",
                 fixed = TRUE)
SET8 <- sapply(X = SET8,
               FUN = function(x) {
                 as.integer(x[1])
               },
               simplify = TRUE)


dat1 <- data.frame("BioSample" = SET1$Biosample,
                   "IS_per_MB" = (SET2[, 4L] / SET1$Total_Length) * 1000000,
                   "FS_per_MB" = (SET2[, 1L] / SET1$Total_Length) * 1000000,
                   "Technology" = SET3,
                   "Total_Runs" = SET4,
                   "TaxID" = SET5,
                   "Assembler" = SET7,
                   "Submission_Year" = SET8,
                   "Cov" = SET1$Coverage)
for (m1 in seq_along(dat1$Technology)) {
  if (dat1$Technology[m1] == "ILLUMINA" & 
      dat1$Total_Runs[m1] == 1L) {
    dat1$Technology[m1] <- "ILLUMINA 1x"
  } else if (dat1$Technology[m1] == "ILLUMINA" & 
             dat1$Total_Runs[m1] > 1L) {
    dat1$Technology[m1] <- "ILLUMINA 2x <="
  }
}

utax <- unique(dat1$TaxID)
res1 <- res2 <- vector(mode = "list",
                       length = length(utax))
for (m1 in seq_along(utax)) {
  res1[[m1]] <- tapply(X = dat1$IS_per_MB[dat1$TaxID == utax[m1]],
                       INDEX = dat1$Technology[dat1$TaxID == utax[m1]],
                       FUN = mean)
  res2[[m1]] <- tapply(X = dat1$FS_per_MB[dat1$TaxID == utax[m1]],
                       INDEX = dat1$Technology[dat1$TaxID == utax[m1]],
                       FUN = mean)
}

w1 <- which(lengths(res1) > 1)
res3 <- res1[w1]
res4 <- res2[w1]
w2 <- sapply(X = res3,
             FUN = function(x) {
               "ILLUMINA 1x" %in% names(x)
             })
res3 <- res3[w2]
res4 <- res4[w2]

res3 <- lapply(X = res3,
               FUN = function(x) {
                 x - x[names(x) == "ILLUMINA 1x"]
               })
res3 <- unlist(res3)
res3 <- tapply(X = unname(res3),
               INDEX = names(res3),
               FUN = c)
res3 <- res3[order(lengths(res3), decreasing = TRUE)]
res3 <- res3[!(names(res3) == "ILLUMINA 1x")]
res4 <- lapply(X = res4,
               FUN = function(x) {
                 x - x[names(x) == "ILLUMINA 1x"]
               })
res4 <- unlist(res4)
res4 <- tapply(X = unname(res4),
               INDEX = names(res4),
               FUN = c)
res4 <- res4[order(lengths(res4), decreasing = TRUE)]
res4 <- res4[!(names(res4) == "ILLUMINA 1x")]

plot(x = 0,
     y = 0,
     type = "n",
     xlab = "dev. from Ill. mean",
     ylab = "frequency",
     xlim = c(-5, 5),
     ylim = c(0, 1))
for (m1 in 1:8) {
  points(x = sort(res4[[m1]]),
         y = seq_along(res4[[m1]]) / length(res4[[m1]]),
         col = m1,
         pch = 46)
}
legend("topleft",
       legend = names(res4)[1:8],
       pch = 20,
       col = 1:8,
       cex = 0.5)

PIN <- 562
t1 <- tapply(X = dat1$IS_per_MB[dat1$TaxID == PIN],
             INDEX = dat1$Assembler[dat1$TaxID == PIN],
             FUN = c)
t1 <- t1[order(lengths(t1), decreasing = TRUE)]
t2 <- tapply(X = dat1$FS_per_MB[dat1$TaxID == PIN],
             INDEX = dat1$Assembler[dat1$TaxID == PIN],
             FUN = c)
t2 <- t2[order(lengths(t2), decreasing = TRUE)]

PIN <- 562 
t1 <- tapply(X = dat1$IS_per_MB[dat1$TaxID == PIN],
             INDEX = dat1$Submission_Year[dat1$TaxID == PIN],
             FUN = c)
# t1 <- t1[order(lengths(t1), decreasing = TRUE)]
t2 <- tapply(X = dat1$FS_per_MB[dat1$TaxID == PIN],
             INDEX = dat1$Submission_Year[dat1$TaxID == PIN],
             FUN = c)
# t2 <- t2[order(lengths(t2), decreasing = TRUE)]

plot(x = 0,
     y = 0,
     type = "n",
     xlab = "IS per MB",
     ylab = "frequency",
     xlim = c(0, 30),
     ylim = c(0, 1),
     main = "",
     xaxs = "i",
     yaxs = "i")
for (m1 in seq_len(length(t1))) {
  points(x = sort(t1[[m1]]),
         y = seq_along(t1[[m1]]) / length(t1[[m1]]),
         col = ColVec1[m1],
         pch = 46)
}
legend("bottomright",
       legend = names(t1)[seq_len(length(t1))],
       pch = 20,
       col = ColVec1[seq_len(length(t2))],
       cex = 0.5)

plot(x = 0,
     y = 0,
     type = "n",
     xlab = "FS per MB",
     ylab = "frequency",
     xlim = c(0, 30),
     ylim = c(0, 1),
     main = "",
     xaxs = "i",
     yaxs = "i")
for (m1 in seq_len(length(t2))) {
  points(x = sort(t2[[m1]]),
         y = seq_along(t2[[m1]]) / length(t2[[m1]]),
         col = ColVec1[m1],
         pch = 46)
}
legend("topleft",
       legend = names(t2)[seq_len(length(t2))],
       pch = 20,
       col = ColVec1[seq_len(length(t2))],
       cex = 0.5)


PIN <- 562
# BINS <- .bincode(x = dat1$Cov, breaks = sort(c(range(dat1$Cov), 1e1, 1e2, 1e3)))
# right = TRUE ~~ <= right bound
BINS <- .bincode(x = dat1$Cov, breaks = c(0, 10, 100, 1000, max(dat1$Cov)))
t1 <- tapply(X = dat1$IS_per_MB[dat1$TaxID == PIN],
             INDEX = BINS[dat1$TaxID == PIN],
             FUN = c)
# t1 <- t1[order(lengths(t1), decreasing = TRUE)]
t2 <- tapply(X = dat1$FS_per_MB[dat1$TaxID == PIN],
             INDEX = BINS[dat1$TaxID == PIN],
             FUN = c)
# t2 <- t2[order(lengths(t2), decreasing = TRUE)]


plot(x = 0,
     y = 0,
     type = "n",
     xlab = "IS per MB",
     ylab = "frequency",
     xlim = c(0, 30),
     ylim = c(0, 1),
     main = "",
     xaxs = "i",
     yaxs = "i")
for (m1 in seq_len(length(t1))) {
  points(x = sort(t1[[m1]]),
         y = seq_along(t1[[m1]]) / length(t1[[m1]]),
         col = ColVec1[m1],
         pch = 46)
}
legend("bottomright",
       legend = names(t1)[seq_len(length(t1))],
       pch = 20,
       col = ColVec1[seq_len(length(t2))],
       cex = 0.5)

plot(x = 0,
     y = 0,
     type = "n",
     xlab = "FS per MB",
     ylab = "frequency",
     xlim = c(0, 30),
     ylim = c(0, 1),
     main = "",
     xaxs = "i",
     yaxs = "i")
for (m1 in seq_len(length(t2))) {
  points(x = sort(t2[[m1]]),
         y = seq_along(t2[[m1]]) / length(t2[[m1]]),
         col = ColVec1[m1],
         pch = 46)
}
legend("topleft",
       legend = names(t2)[seq_len(length(t2))],
       pch = 20,
       col = ColVec1[seq_len(length(t2))],
       cex = 0.5)


save(dat1,
     res3,
     res4,
     file = "~/Repos/PseudogenePlots/InputData/Technology_Observations_v01.RData",
     compress = "xz")
