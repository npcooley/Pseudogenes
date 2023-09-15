###### -- plot distributions by assembler / technology choices ----------------

suppressMessages(library(SynExtend))

load(file = "~/Repos/PseudogenePlots/InputData/Technology_Observations_v01.RData",
     verbose = TRUE)

ColVec1 <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000')
ColVec2 <- paste0(ColVec1,
                  "55")
TOPTAX <- head(names(sort(table(dat1$TaxID), decreasing = TRUE)))
dat2 <- dat1[dat1$TaxID %in% as.integer(TOPTAX), ]
pvals1 <- pvals2 <- submitterchoices <- vector(mode = "list",
                                               length = 6L)

for (m2 in seq_along(TOPTAX)) {
  
  pvals1[[m2]] <- pvals2[[m2]] <- vector(mode = "list",
                                         length = 5L)
  
  dat3 <- by(data = dat2[dat2$TaxID %in% as.integer(TOPTAX)[m2], ],
             INDICES = ~ Technology + Assembler,
             FUN = function(x) {
               x$FS_per_MB
             })
  dat4 <- by(data = dat2[dat2$TaxID %in% as.integer(TOPTAX)[m2], ],
             INDICES = ~ Technology + Assembler,
             FUN = function(x) {
               x$IS_per_MB
             })
  w1 <- lengths(dat3)
  w2 <- which(w1 >= sort(w1, decreasing = TRUE)[6], arr.ind = TRUE)
  o1 <- order(w1, decreasing = TRUE)
  head(sort(w1, decreasing = TRUE), n = 15)
  
  # find the most populated distribution
  # ask how different all of the smaller top six are from the largest
  w3 <- w4 <- quants1 <- quants2 <- vector(mode = "list",
                                           length = 6L)
  combonames <- vector(mode = "character",
                       length = 6L)
  combosizes <- vector(mode = "integer",
                       length = 6L)
  
  for (m1 in seq_len(nrow(w2))) {
    w3[[m1]] <- dat3[[w2[m1, 1], w2[m1, 2]]]
    w4[[m1]] <- dat4[[w2[m1, 1], w2[m1, 2]]]
    combonames[m1] <- paste(dimnames(dat3)[[1]][w2[m1, 1]],
                            dimnames(dat3)[[2]][w2[m1, 2]],
                            sep = " + ")
    combosizes[m1] <- length(w3[[m1]])
    quants1[[m1]] <- quantile(x = w3[[m1]],
                              probs = c(0.01, 0.99))
    quants2[[m1]] <- quantile(x = w4[[m1]],
                              probs = c(0.01, 0.99))
  }
  
  o2 <- order(combosizes, decreasing = TRUE)
  w3 <- w3[o2]
  w4 <- w4[o2]
  combonames <- combonames[o2]
  combosizes <- combosizes[o2]
  submitterchoices[[m2]] <- combonames
  
  plot(x = 0,
       y = 0,
       type = "n",
       # xlim = c(min(mean(sapply(X = quants1,
       #                          FUN = function(x) {
       #                            x[1]
       #                          })),
       #              mean(sapply(X = quants2,
       #                          FUN = function(x) {
       #                            x[1]
       #                          }))),
       #          max(mean(sapply(X = quants1,
       #                          FUN = function(x) {
       #                            x[2]
       #                          })),
       #              mean(sapply(X = quants2,
       #                          FUN = function(x) {
       #                            x[2]
       #                          })))),
       xlim = c(min(c(sapply(X = quants1,
                             FUN = function(x) {
                               x[1]
                             }),
                      sapply(X = quants2,
                             FUN = function(x) {
                               x[1]
                             }))),
                max(mean(sapply(X = quants1,
                                FUN = function(x) {
                                  x[2]
                                })),
                    mean(sapply(X = quants2,
                                FUN = function(x) {
                                  x[2]
                                })))),
       ylim = c(0, 1),
       xaxs = "i",
       yaxs = "i",
       xlab = "PG per Mbp",
       ylab = "Cumulative Frequency",
       main = TOPTAX[m2])
  
  for (m1 in seq_len(6L)) {
    
    lines(x = sort(w3[[m1]]),
          y = seq(length(w3[[m1]])) / length(w3[[m1]]),
          lty = 1,
          lwd = if (m1 == 1L) {
            3
          } else {
            1
          },
          col = if (m1 == 1) {
            ColVec1[m1]
          } else {
            ColVec2[m1]
          })
    
    if (m1 > 1) {
      pvals1[[m2]][[m1 - 1L]] <- suppressWarnings(ks.test(x = w3[[1L]],
                                                          y = w3[[m1]]))
    }
    
    lines(x = sort(w4[[m1]]),
          y = seq(length(w4[[m1]])) / length(w4[[m1]]),
          lty = 3,
          lwd = if (m1 == 1L) {
            3
          } else {
            1
          },
          col = if (m1 == 1) {
            ColVec1[m1]
          } else {
            ColVec2[m1]
          })
    if (m1 > 1) {
      pvals2[[m2]][[m1 - 1L]] <- suppressWarnings(ks.test(x = w4[[1L]],
                                                          y = w4[[m1]]))
    }
  }
}



