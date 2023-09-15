###### -- set comparison post script ------------------------------------------

suppressMessages(library(SynExtend))

path01 <- "~/Data/20220920_SetComparison"

files01 <- list.files(path = path01)

pBar <- txtProgressBar(style = 1)
PBAR <- length(files01)
res01 <- vector(mode = "list",
                length = PBAR)

for (m1 in seq_along(files01)) {
  load(file = paste0(path01,
                     "/",
                     files01[m1]),
       verbose = FALSE)
  
  # z <- mean(ANIRes[1:2])
  # if (z < 90) {
  #   next
  # }
  
  z1 <- strsplit(x = rownames(clust),
                 split = " ",
                 fixed = TRUE)
  z1 <- sapply(X = z1,
               FUN = function(x) {
                 x[2]
               })
  z1 <- grepl(pattern = "pseudo",
              x = z1)
  
  z2 <- tapply(X = rownames(clust),
               INDEX = clust$cluster,
               FUN = c)
  z3 <- z2[lengths(z2) > 1]
  
  # any case where a pseudogene's sole partner in the other genome is a regular gene
  res2 <- vector(mode = "logical",
                 length = length(z3))
  for (m2 in seq_along(z3)) {
    t1 <- any(grepl(pattern = paste(w2[1],
                                    "pseudogene"),
                    x = z3[m2]))
    t2 <- any(grepl(pattern = paste(w2[1],
                                    "gene"),
                    x = z3[m2]))
    t3 <- any(grepl(pattern = paste(w2[2],
                                    "pseudogene"),
                    x = z3[m2]))
    t4 <- any(grepl(pattern = paste(w2[2],
                                    "gene"),
                    x = z3[m2]))
    
    # IF:
    # T F F T
    # OR:
    # F T T F
    if ((t1 & !t2 & !t3 & t4) | (!t1 & t2 & t3 & !t4)) {
      res2[m2] <- TRUE
      # no else, logical vectors are FALSE by default
    }
  } # end m2 loop
  res01[[m1]] <- c(length(res2),
                   sum(res2),
                   ANIRes,
                   PseudosFrameShift,
                   PseudosIntStop,
                   w2)
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}

res02 <- do.call(rbind,
                 res01)
res03 <- data.frame("pairs" = as.integer(res02[, 1]),
                    "incongruent" = as.integer(res02[, 2]),
                    "ani" = apply(X = res02,
                                  MARGIN = 1L,
                                  FUN = function(x) {
                                    mean(x[3:4])
                                  }),
                    "af" = apply(X = res02,
                                 MARGIN = 1L,
                                 FUN = function(x) {
                                   mean(x[5:6])
                                 }),
                    "fs_kbp_1" = res02[, 7],
                    "fs_kbp_2" = res02[, 8],
                    "is_kbp_1" = res02[, 9],
                    "is_kbp_2" = res02[, 10],
                    "partner1" = as.integer(res02[, 11]),
                    "partner2" = as.integer(res02[, 12]))

# by best partner for ANI and AF
u1 <- unique(unlist(res03[, 9:10]))

pBar <- txtProgressBar(style = 1)
PBAR <- length(u1)
res04 <- res05 <- vector(mode = "integer",
                         length = PBAR)
for (m1 in seq_along(u1)) {
  
  w1 <- unique(which(res03$partner1 == u1[m1] | res03$partner2 == u1[m1]))
  w2 <- which.max(res03$ani[w1])
  res04[m1] <- w1[w2]
  
  w2 <- which.max(res03$af[w1])
  res05[m1] <- w1[w2]
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}

ani_dat <- res03[res04, ]
af_dat <- res03[res05, ]

ani_dat <- ani_dat[!duplicated(ani_dat[, 9:10]), ]
af_dat <- af_dat[!duplicated(af_dat[, 9:10]), ]

save(ani_dat,
     af_dat,
     res03,
     file = "~/Repos/PseudogenePlots/InputData/ClusteredPseudos.RData",
     compress = "xz")

# pdf(file = "~/testplotpairs01.pdf",
#     height = 3.5,
#     width = 3.5)
# par(mar = c(3,3,2,2),
#     mgp = c(2, 0.85, 0),
#     cex.lab = 0.95)
# plot(x = ani_dat$ani,
#      y = ani_dat$incongruent / ani_dat$pairs,
#      pch = 20,
#      cex = 0.5,
#      col = "#11558833",
#      xlim = c(90, 100),
#      ylab = "incongruent pairs / all pairs",
#      xlab = "ANI",
#      ylim = c(0, 0.1))
# spfit01 <- smooth.spline(x = ani_dat$ani,
#                          y = ani_dat$incongruent / ani_dat$pairs,
#                          df = 5)
# spfit02 <- predict(spfit01)
# lines(x = spfit02$x,
#       y = spfit02$y)
# dev.off()
# 
# pdf(file = "~/testplotpairs02.pdf",
#     height = 3.5,
#     width = 3.5)
# par(mar = c(3,3,2,2),
#     mgp = c(2, 0.85, 0),
#     cex.lab = 0.95)
# plot(x = af_dat$af,
#      y = af_dat$incongruent / af_dat$pairs,
#      pch = 20,
#      cex = 0.5,
#      col = "#11558833",
#      xlim = c(.6, 1),
#      ylab = "incongruent pairs / all pairs",
#      xlab = "af",
#      ylim = c(0, 0.1))
# spfit01 <- smooth.spline(x = af_dat$af,
#                          y = af_dat$incongruent / af_dat$pairs,
#                          df = 5)
# spfit02 <- predict(spfit01)
# lines(x = spfit02$x,
#       y = spfit02$y)
# dev.off()

# collect 

