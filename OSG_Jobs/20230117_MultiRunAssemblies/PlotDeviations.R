###### -- plot out results from the paired assemblies -------------------------

suppressMessages(library(SynExtend))

load(file = "~/Data/20230221_MultiRunIntermediateData/ResTable.RData",
     verbose = TRUE)

ph1 <- table(ResTable$BioSample)
ph1 <- names(ph1[ph1 == 2])
SubSetTable <- ResTable[ResTable$BioSample %in% ph1, ]

Cat1 <- vector(mode = "character",
               length = length(ph1))
Diff1 <- Diff2 <- vector(mode = "numeric",
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
  
  Diff1[m1] <- ISpMB[1] - ISpMB[2]
  Diff2[m1] <- FSpMB[1] - FSpMB[2]
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}

Cat2 <- names(table(Cat1)[table(Cat1) > 100])
ph3 <- Cat1 %in% Cat2

ISDiffs <- Diff1[ph3]
FSDiffs <- Diff2[ph3]
CatSubSet <- Cat1[ph3]

# ColVec <- c("olivedrab", "springgreen1", "orchid", "tomato", "blue", "navy")
#ColVec <- c('#ffe119', '#4363d8', '#f58231', '#dcbeff', '#800000', '#000075')
ColVec <- c('#ffe119', '#4363d8', '#f58231', '#dcbeff', '#32CD32', '#000075')
#ColVec2 <- gplots::col2hex(ColVec)
ColVec2 <- paste0(ColVec,
                  "33")
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

pdf(file = "~/Repos/20230117_MultiRunAssemblies/PairedRunDifferences.pdf",
    height = 7,
    width = 7)

layout(mat = matrix(data = c(1,1,2,2,
                             1,1,2,2,
                             1,1,3,3,
                             1,1,3,3),
                    nrow = 4,
                    byrow = TRUE))
par(mar = c(4,3,2,1),
    mgp = c(2, 0.75, 0))
plot(x = Diff1[Cat1 %in% Cat2],
     y = Diff2[Cat1 %in% Cat2],
     col = ColVec2[match(x = Cat1[Cat1 %in% Cat2],
                         table = Cat2)],
     xlab = "Diff. IS per MB",
     ylab = "Diff. FS per MB",
     pch = ifelse(test = Cat1[Cat1 %in% Cat2] %in% Cat2[2:4],
                  yes = 17,
                  no = 20),
     # pch.cex = 0.75,
     xlim = c(-200, 100),
     xaxs = "i")
legend("topleft",
       legend = c("Illumina - Illumina, n = 3203",
                  "Illumina - 454, n = 276",
                  "Illumina - ONT, n = 1836",
                  "Illumina - PacBio, n = 1017",
                  "454 - 454, n = 193",
                  "PacBio - PacBio, n = 293"),
       pch = ifelse(test = Cat2 %in% Cat2[2:4],
                    yes = 17,
                    no = 20),
       col = ColVec,
       lty = 1,
       cex = 1,
       bg = NA,
       bty = "n")

plot(x = 0,
     y = 0,
     type = "n",
     xlab = "Diff. in IS per MB",
     ylab = "Frequency",
     ylim = c(0, 1),
     xlim = c(-80, 20))
for (m1 in seq_along(ColVec)) {
  # points(x = ISDiffs1[[m1]],
  #        y = seq_along(ISDiffs1[[m1]]) / length(ISDiffs1[[m1]]),
  #        col = ColVec[m1],
  #        pch = 20)
  yvals <- seq_along(ISDiffs1[[m1]]) / length(ISDiffs1[[m1]])
  p1 <- mean(ISDiffs1[[m1]])
  p2 <- yvals[which.min(abs(ISDiffs1[[m1]] - p1))]
  lines(x = ISDiffs1[[m1]],
        y = yvals,
        col = ColVec[m1],
        lty = 1)
  points(x = p1,
         y = p2,
         pch = 2,
         col = ColVec[m1])
}
abline(v = 0, lty = 2)
# legend("topleft",
#        legend = names(ISDiffs1),
#        lty = 1,
#        col = ColVec,
#        cex = 0.5,
#        bg = NA,
#        bty = "n")

plot(x = 0,
     y = 0,
     type = "n",
     xlab = "Diff. in FS per MB",
     ylab = "Frequency",
     ylim = c(0, 1),
     xlim = c(-500, 100))
for (m1 in seq_along(ColVec)) {
  # points(x = ISDiffs1[[m1]],
  #        y = seq_along(ISDiffs1[[m1]]) / length(ISDiffs1[[m1]]),
  #        col = ColVec[m1],
  #        pch = 20)
  yvals <- seq_along(FSDiffs1[[m1]]) / length(FSDiffs1[[m1]])
  p1 <- mean(FSDiffs1[[m1]])
  p2 <- yvals[which.min(abs(FSDiffs1[[m1]] - p1))]
  lines(x = FSDiffs1[[m1]],
        y = yvals,
        col = ColVec[m1],
        lty = 1)
  points(x = p1,
         y = p2,
         pch = 2,
         col = ColVec[m1])
}
abline(v = 0, lty = 2)
# legend("topleft",
#        legend = Cat2,
#        lty = 1,
#        col = ColVec,
#        cex = 0.5,
#        bg = NA,
#        bty = "n")

dev.off()

# ISDiffs2 <- tapply(X = ISDiffs,
#                    INDEX = CatSubSet,
#                    FUN = function(x) {
#                      density(x)
#                    })
# FSDiffs2 <- tapply(X = FSDiffs,
#                    INDEX = CatSubSet,
#                    FUN = function(x) {
#                      density(x)
#                    })
# 
# plot(x = 0,
#      y = 0,
#      type = "n",
#      ylim = c(0, 0.4),
#      xlim = c(-20, 20),
#      xlab = "Diff. in I.S.",
#      ylab = "Density")
# for (m1 in seq_along(ISDiffs)) {
#   lines(ISDiffs[[m1]],
#         col = m1)
# }
# legend("topleft",
#        legend = Cat2,
#        lty = 1,
#        col = seq(length(Cat2)),
#        cex = 0.5,
#        bg = NA,
#        bty = "n")
# abline(v = 0, lty = 2)
# 
# plot(x = 0,
#      y = 0,
#      type = "n",
#      ylim = c(0, 0.0175),
#      xlim = c(-400, 400),
#      xlab = "Diff. in F.S.",
#      ylab = "Density")
# for (m1 in seq_along(ISDiffs)) {
#   lines(FSDiffs[[m1]],
#         col = m1)
# }
# legend("topleft",
#        legend = Cat2,
#        lty = 1,
#        col = seq(length(Cat2)),
#        cex = 0.5,
#        bg = NA,
#        bty = "n")
# abline(v = 0, lty = 2)

ISDevFromSource <- ((SubSetTable$IS / SubSetTable$AssemblySize) - (SubSetTable$Source_IS / SubSetTable$Source_Size)) * 1000000
FSDevFromSource <- ((SubSetTable$FS / SubSetTable$AssemblySize) - (SubSetTable$Source_FS / SubSetTable$Source_Size)) * 1000000

plot(y = FSDevFromSource,
     x = ISDevFromSource,
     pch = 46,
     ylab = "Diff. in FS",
     xlab = "Diff. in IS",
     col = match(x = SubSetTable$Platform,
                 table = unique(SubSetTable$Platform)),
     # ylim = c(-50, 350),
     # xlim = c(-100, 100),
     xaxs = "i",
     yaxs = "i",
     main = "Deviation in Reassembly")
abline(h = -10, lty = 2)
abline(v = -10, lty = 2)
legend("right",
       legend = unique(SubSetTable$Platform),
       pch = 20,
       col = seq_along(unique(SubSetTable$Platform)),
       cex = 0.5)


plot(y = ((SubSetTable$FS[SubSetTable$Platform == "ILLUMINA"] / SubSetTable$AssemblySize[SubSetTable$Platform == "ILLUMINA"]) -
       (SubSetTable$Source_FS[SubSetTable$Platform == "ILLUMINA"] / SubSetTable$Source_Size[SubSetTable$Platform == "ILLUMINA"])) * 1000000,
     x = ((SubSetTable$IS[SubSetTable$Platform == "ILLUMINA"] / SubSetTable$AssemblySize[SubSetTable$Platform == "ILLUMINA"]) -
       (SubSetTable$Source_IS[SubSetTable$Platform == "ILLUMINA"] / SubSetTable$Source_Size[SubSetTable$Platform == "ILLUMINA"])) * 1000000,
     pch = 46,
     ylab = "FS",
     xlab = "IS",
     xlim = c(-30, 10),
     ylim = c(-50, 10))

plot(y = ((SubSetTable$FS / SubSetTable$AssemblySize) -
            (SubSetTable$Source_FS / SubSetTable$Source_Size)) * 1000000,
     x = ((SubSetTable$IS / SubSetTable$AssemblySize) -
            (SubSetTable$Source_IS / SubSetTable$Source_Size)) * 1000000,
     pch = 46,
     ylab = "FS",
     xlab = "IS",
     col = match(x = SubSetTable$Platform,
                 table = unique(SubSetTable$Platform)),
     ylim = c(-50, 350),
     xlim = c(-100, 100),
     xaxs = "i",
     yaxs = "i")
abline(h = -10, lty = 2)
abline(v = -10, lty = 2)
legend("topleft",
       legend = unique(SubSetTable$Platform),
       pch = 20,
       col = seq_along(unique(SubSetTable$Platform)),
       cex = 0.65)




