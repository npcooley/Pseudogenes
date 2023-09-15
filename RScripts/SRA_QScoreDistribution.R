###### -- parse the very condensed SRA stats ----------------------------------

suppressMessages(library(SynExtend))

files01 <- list.files(path = "~/Data/20230605_ReadStatCollection/",
                      full.names = TRUE)
jobmap <- read.table(file = "~/Repos/20230605_ReadStatCollection/JobMap.txt")
currentcompletedjobs <- as.integer(unlist(regmatches(x = files01,
                                                     m = gregexpr(pattern = "(?<=ReadStats)([0-9]+)(?=\\.RData)",
                                                                  perl = TRUE,
                                                                  text = files01))))


resa <- resb <- vector(mode = "numeric",
                       length = length(files01))
resc <- vector(mode = "character",
               length = length(files01))

resd <- rese <- vector(mode = "list",
                       length = length(files01))

pBar <- txtProgressBar(style = 1)
PBAR <- length(files01)
for (m1 in seq_along(resa)) {
  
  load(file = files01[m1],
       verbose = FALSE)
  
  if (length(RES) == 0) {
    next
  }
  
  if (length(RES[[1]]) == 3L) {
    resa[m1] <- RES[[1]][[3]]
    
    if (length(RES) > 1L) {
      resb[m1] <- RES[[2]][[3]]
    }
  } else if (length(RES[[1]]) == 5L) {
    resa[m1] <- RES[[1]][[3]]
    
    if (length(RES) > 1L) {
      resb[m1] <- RES[[2]][[3]]
    }
    resd[[m1]] <- RES[[1]][[4]]
    rese[[m1]] <- RES[[1]][[5]]
  }
  
  resc[m1] <- jobmap[currentcompletedjobs[m1], 3L]
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}

save(resa,
     resb,
     resc,
     file = "~/Repos/PseudogenePlots/InputData/SRA_Stats01.RData",
     compress = "xz")

# mean mean qscore to center simulated data plots.
# > mean(resa[resc == "ILLUMINA"])
# [1] 35.33467
# > mean(resb[resc == "ILLUMINA"])
# [1] 33.92757

plot(x = 0,
     y = 0,
     ylim = c(0, 1),
     xlim = c(0, 50),
     xlab = "mean Q score",
     ylab = "Cumulative Frequency",
     type = "n")
z <- resa[resc == "ILLUMINA"]
points(x = sort(z),
       y = seq_along(z) / length(z),
       pch = 46,
       col = 1)
z <- resb[resc == "ILLUMINA"]
points(x = sort(z),
       y = seq_along(z) / length(z),
       pch = 46,
       col = 2)
z <- resa[resc == "LS454"]
points(x = sort(z),
       y = seq_along(z) / length(z),
       pch = 46,
       col = 3)
z <- resa[resc == "PACBIO_SMRT"]
points(x = sort(z),
       y = seq_along(z) / length(z),
       pch = 46,
       col = 4)
z <- resa[resc == "OXFORD_NANOPORE"]
points(x = sort(z),
       y = seq_along(z) / length(z),
       pch = 46,
       col = 5)
legend("topleft",
       legend = c("Illumina p1",
                  "Illumina p2",
                  "454",
                  "PB",
                  "ONT"),
       lty = 1,
       col = 1:5,
       bty = "n")




