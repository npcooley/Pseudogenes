###### -- MultiRunAssemblyDataParsing -----------------------------------------

load(file = "~/Data/20230221_MultiRunIntermediateData/ResTable.RData",
     verbose = TRUE)

ColVec1 <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000')
ColVec2 <- paste0(ColVec1,
                  "33")

plot(ResTable$AssemblySize / ResTable$Source_Size,
     pch = 46)
abline(h = 0.95, lty = 2, lwd = 2, col = "orange")
abline(h = 1.05, lty = 2, lwd = 2, col = "orange")

SubSetTable <- ResTable[(ResTable$AssemblySize >= (ResTable$Source_Size * 0.95)) &
                          (ResTable$AssemblySize <= (ResTable$Source_Size * 1.05)), ]

ph1 <- table(SubSetTable$BioSample)
ph1 <- names(ph1[ph1 == 2])
SubSetTable <- SubSetTable[SubSetTable$BioSample %in% ph1, ]

Cat1 <- vector(mode = "character",
               length = length(ph1))
Diff1 <- Diff2 <- Diff3 <- vector(mode = "numeric",
                                  length = length(ph1))

pBar <- txtProgressBar(style = 1)
PBAR <- length(ph1)
for (m1 in seq_along(ph1)) {
  ph2 <- SubSetTable$BioSample == ph1[m1]
  CurrentSet <- SubSetTable[ph2, ]
  o1 <- order(CurrentSet$Platform)
  CurrentSet <- CurrentSet[o1, ]
  Cat1[m1] <- paste(CurrentSet$Platform,
                    collapse = "_")
  
  ISpMB <- (CurrentSet$IS / CurrentSet$AssemblySize) * 1000000
  FSpMB <- (CurrentSet$FS / CurrentSet$AssemblySize) * 1000000
  
  Diff1[m1] <- ISpMB[2] - ISpMB[1]
  Diff2[m1] <- FSpMB[2] - FSpMB[1]
  Diff3[m1] <- min(CurrentSet$AssemblySize) / max(CurrentSet$AssemblySize)
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)

# size restriction has already been imposed
# ph3 <- Diff3 >= 0.9
# Diff1 <- Diff1[ph3]
# Diff2 <- Diff2[ph3]
# Cat1 <- Cat1[ph3]


Cat2 <- names(table(Cat1)[table(Cat1) > 50])
ph3 <- Cat1 %in% Cat2
ISDiffs <- Diff1[ph3]
FSDiffs <- Diff2[ph3]
CatSubSet <- Cat1[ph3]
table(CatSubSet)

ISDiffs1 <- tapply(X = ISDiffs,
                   INDEX = CatSubSet,
                   FUN = function(x) {
                     sort(x)
                   })
FSDiffs1 <- tapply(X = FSDiffs,
                   INDEX = CatSubSet,
                   FUN = function(x) {
                     sort(x)
                   })

save(FSDiffs1,
     ISDiffs1,
     ISDiffs,
     FSDiffs,
     Cat2,
     CatSubSet,
     file = "~/Repos/PseudogenePlots/InputData/MultRunReassemblies_v02.RData",
     compress = "xz")


layout(mat = matrix(data = c(1,2),
                    nrow = 1))
par(mar = c(4,3,2,0.75),
    mgp = c(2, 0.75, 0))
plot(x = 0,
     y = 0,
     type = "n",
     xlab = "Diff. in IS per MB",
     ylab = "Frequency",
     ylim = c(0, 1),
     xlim = c(-20, 80))
for (m1 in seq_along(Cat2)) {
  # points(x = ISDiffs1[[m1]],
  #        y = seq_along(ISDiffs1[[m1]]) / length(ISDiffs1[[m1]]),
  #        col = ColVec1[m1],
  #        pch = 20)
  yvals <- seq_along(ISDiffs1[[m1]]) / length(ISDiffs1[[m1]])
  p1 <- mean(ISDiffs1[[m1]])
  p2 <- yvals[which.min(abs(ISDiffs1[[m1]] - p1))]
  lines(x = ISDiffs1[[m1]],
        y = yvals,
        col = ColVec1[m1],
        lty = 1)
  points(x = p1,
         y = p2,
         pch = 2,
         col = ColVec1[m1])
}
abline(v = 0, lty = 2)
# legend("topleft",
#        legend = names(ISDiffs1),
#        lty = 1,
#        col = ColVec1,
#        cex = 0.5,
#        bg = NA,
#        bty = "n")
par(mar = c(4,2.75,2,1))
plot(x = 0,
     y = 0,
     type = "n",
     xlab = "Diff. in FS per MB",
     ylab = "Frequency",
     ylim = c(0, 1),
     xlim = c(-100, 500))
for (m1 in seq_along(Cat2)) {
  # points(x = ISDiffs1[[m1]],
  #        y = seq_along(ISDiffs1[[m1]]) / length(ISDiffs1[[m1]]),
  #        col = ColVec1[m1],
  #        pch = 20)
  yvals <- seq_along(FSDiffs1[[m1]]) / length(FSDiffs1[[m1]])
  p1 <- mean(FSDiffs1[[m1]])
  p2 <- yvals[which.min(abs(FSDiffs1[[m1]] - p1))]
  lines(x = FSDiffs1[[m1]],
        y = yvals,
        col = ColVec1[m1],
        lty = 1)
  points(x = p1,
         y = p2,
         pch = 2,
         col = ColVec1[m1])
}
abline(v = 0, lty = 2)
# legend("topleft",
#        legend = Cat2,
#        lty = 1,
#        col = ColVec1,
#        cex = 0.5,
#        bg = NA,
#        bty = "n")

dev.off()
