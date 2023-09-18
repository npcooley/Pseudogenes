###### -- plot out GFFs ------------------------------------------------------

suppressMessages(library(SynExtend))

PATH05 <- "~/Repos/20221207_FactorialAssembly"

load(file = paste0(PATH05,
                   "ParsedPseudogenes.RData"),
     verbose = TRUE)

# select out anomalously sized assemblies and only trimmed assemblies
pv2 <- median(ResTable$AssemblySize)
pv1 <- ResTable$AssemblySize > (pv2 * 0.75) & ResTable$AssemblySize < (pv2 * 1.25)  # assembly is not anomalously sized
pv2 <- ResTable$Trim
TableForPlots <- ResTable[pv1 & pv2, ]

FRperKB <- (TableForPlots$FR / TableForPlots$AssemblySize) * 1000
ISperKB <- (TableForPlots$IS / TableForPlots$AssemblySize) * 1000

ColVec02 <- c(a = "firebrick",
              b = "plum1",
              c = "orchid",
              d = "tomato")
ColVec02 <- c(a = "purple",
              b = "firebrick",
              c = "turquoise",
              d = "goldenrod")

matmatch <- function(z1) {
  
  rac <- apply(X = z1,
               MARGIN = 1L,
               FUN = function(x) {
                 paste(x, collapse = "_")
               })
  
  rtab <- table(rac)
  ans <- as.integer(unname(rtab[match(x = rac,
                                      table = names(rtab))]))
  return(ans)
  
}


KEY01 <- c(NA,
           NA,
           "Pseudogenes per Mb")
KEY03 <- list(c("Megahit", "Unicycler", "SKESA", "SPADES"), # all four assemblers
              c(2, 3, 4), # quality bins
              c(0.05, 0.25, 0.5, 1)) # down sampling
KEY06 <- c("Assembler",
           "QualBin",
           "DownSample")
KEY07 <- list(c("Megahit", "Unicycler", "SKESA", "SPADES"), # all four assemblers
              c("20-30", "30+", "All"), # quality bins
              c("0.05x", "0.25x", "0.50x", "1.00x")) # sparsification

pdf(file = paste0(PATH05,
                  "/Neisseria_Factorial_Distributions.pdf"),
    height = 10.5,
    width = 3.5)
StatsTables <- vector(mode = "list",
                      length = 3L)

# plot margins and spacing
par(mar = c(4, 4, 1, 1),
    mgp = c(2.3,1, 0))
layout(mat = matrix(data = 1:3,
                    nrow = 3))

for (m1 in seq(3L)) {
  
  CurrentTarget <- KEY06[m1]
  SelectCols <- unique(KEY06)
  SelectCols <- SelectCols[!(SelectCols == CurrentTarget)]
  SelectCols <- c(SelectCols, "Run")
  
  CurrentCounts <- matmatch(z1 = TableForPlots[, SelectCols])
  CurrentSelection <- CurrentCounts == max(CurrentCounts)
  print(sum(CurrentSelection) / max(CurrentCounts))
  
  StatsTables[[m1]] <- matrix(data = 0,
                              nrow = length(KEY03[[m1]]),
                              ncol = 4L)
  plot(x = 0,
       y = 0,
       type = "n",
       xlab = KEY01[m1],
       ylab = "Frequency",
       xlim = c(0,0.07) * 1000,
       ylim = c(0, 1),
       yaxs = "i")
  # in each plot window, plot IS as solid
  # plot FS as dashed
  for (m2 in seq_along(KEY03[[m1]])) {
    pv3 <- TableForPlots[, KEY06[m1]] == KEY03[[m1]][m2]
    
    StatsTables[[m1]][m2, ] <- c(mean(TableForPlots[pv3 & CurrentSelection, "AssemblySize"]),
                                 mean(TableForPlots[pv3 & CurrentSelection, "N50"]),
                                 mean(TableForPlots[pv3 & CurrentSelection, "L50"]),
                                 mean(TableForPlots[pv3 & CurrentSelection, "TC"]))
    lines(x = sort(ISperKB[pv3 & CurrentSelection]) * 1000,
          y = seq(sum(pv3 & CurrentSelection)) / sum(pv3 & CurrentSelection),
          col = ColVec02[m2],
          lty = 1)
    
    lines(x = sort(FRperKB[pv3 & CurrentSelection]) * 1000,
          y = seq(sum(pv3 & CurrentSelection)) / sum(pv3 & CurrentSelection),
          col = ColVec02[m2],
          lty = 3)
  }
  
  L <- legend(x = 45,
              y = 0.3,
              legend = rep(NA, length(KEY07[[m1]]) * 2),
              col = rep(ColVec02[seq(length(KEY07[[m1]]))],
                        2),
              lty = c(rep(1, length(KEY07[[m1]])),
                      rep(3, length(KEY07[[m1]]))),
              # pch = rep(NA_integer_, length(length(KEY07[[m1]])) * 2L),
              ncol = 2,
              bty = 'n',
              x.intersp = 0.5,
              inset = 0.02,
              cex = 0.75)
  legend(x = L$rect$left,
         y = L$rect$top,
         legend = KEY07[[m1]],
         col = rep(NA,2),
         lty = c(1,3),
         ncol = 1,
         x.intersp = 5,
         bg = NA,
         bty = "n",
         cex = 0.75)
  mtext("IS",
        cex = 0.55,
        at = 48.5,
        line = -15)
  mtext("FS",
        cex = 0.55,
        at = 54.75,
        line = -15)
  
  # 1st n = 3877
  # 2nd n = 857
  # 3rd n = 1242
  # simple legend
  # legend("bottomright",
  #        legend = KEY07[[m1]],
  #        lty = 1,
  #        pch = 20,
  #        col = ColVec02[seq_along(KEY03[[m1]])])
}

dev.off()




