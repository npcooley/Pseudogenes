###### -- subset neisseria data slightly more intelligently -------------------

ColVec1 <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000')
ColVec2 <- paste0(ColVec1,
                  "33")

load(file = "~/Repos/Pseudogenes/SearchResults.RData",
     verbose = TRUE)
load(file = "~/Repos/Pseudogenes/SearchResults2.RData",
     verbose = TRUE)

load(file = "~/Repos/20221207_FactorialAssembly/ParsedPseudogenes.RData",
     verbose = TRUE)
# Subset to N. gonorrhoeae only
snvec <- SRAResults$ScientificName[match(x = as.character(ResTable$Run),
                                         table = SRAResults$Run)]
w1 <- grepl(pattern = "gonorrhoeae",
            x = snvec)

ResTable <- ResTable[w1, ]
ResTable$Run <- as.character(ResTable$Run)

# select out anomalously sized assemblies and retain:
# trimmed assemblies that didn't quality bin and didn't downsample
pv2 <- median(ResTable$AssemblySize)
pv1 <- ResTable$AssemblySize > (pv2 * 0.85) & ResTable$AssemblySize < (pv2 * 1.15)  # assembly is not anomalously sized
pv2 <- ResTable$Trim
pv3 <- ResTable$DownSample == 1
pv4 <- ResTable$QualBin == 4
dat1 <- ResTable[pv1 &
                   pv2 &
                   pv3 &
                   pv4, ]


# relate sets back to the reported vals in RefSeq via a source biosample
SBS <- SRAResults$BioSample[match(x = dat1$Run,
                                  table = SRAResults$Run)]
UBS <- unique(SBS)
w1 <- match(x = UBS,
            table = EntrezResults$Biosample)
AssemblerCode <- apply(X = AssemblersLogical[w1, ],
                       MARGIN = 1L,
                       FUN = function(x) {
                         paste0(colnames(AssemblersLogical)[x],
                                collapse = " + ")
                       })
dat2 <- data.frame("BioSample" = UBS,
                   "ISperMB" = (PseudosCondensed[w1, 4L] / EntrezResults$Total_Length[w1]) * 1000000,
                   "FSperMB" = (PseudosCondensed[w1, 1L] / EntrezResults$Total_Length[w1]) * 1000000,
                   "Contig_N50" = EntrezResults$ContigN50[w1],
                   "Total_Nucs" = EntrezResults$Total_Length[w1],
                   "SPADES" = AssemblersLogical[w1, "SPADES"],
                   "Megahit" = AssemblersLogical[w1, "Megahit"],
                   "Unicycler" = AssemblersLogical[w1, "Unicycler"],
                   "SKESA" = AssemblersLogical[w1, "Skesa"])

dat1 <- cbind(dat1,
              "BioSample" = SBS)

U_Assembler <- unique(dat1$Assembler)

ISperMB <- (dat1$IS / dat1$AssemblySize) * 1000000
FSperMB <- (dat1$FR / dat1$AssemblySize) * 1000000
# negative vals mean reported counts are higher than reassembled counts
ISdev <- ISperMB - dat2$ISperMB[match(x = dat1$BioSample,
                                      table = dat2$BioSample)]
FSdev <- FSperMB - dat2$FSperMB[match(x = dat1$BioSample,
                                      table = dat2$BioSample)]

plot(x = 0,
     y = 0,
     type = "n",
     xlab = "PG per MB",
     ylab = "Frequency",
     ylim = c(0, 1),
     xlim = c(10, 60))
for (m1 in seq_along(U_Assembler)) {
  pv3 <- dat1$Assembler == U_Assembler[m1]
  
  lines(x = sort(ISperMB[pv3]),
        y = seq(sum(pv3)) / sum(pv3),
        col = ColVec1[m1],
        lty = 1)
  
  lines(x = sort(FSperMB[pv3]),
        y = seq(sum(pv3)) / sum(pv3),
        col = ColVec1[m1],
        lty = 3)
}

plot(x = 0,
     y = 0,
     type = "n",
     xlab = "PG per MB",
     ylab = "Frequency",
     ylim = c(0, 1),
     xlim = c(-10, 5))
for (m1 in seq_along(U_Assembler)) {
  # drop points that were originally assembled by the current assembler
  keep <- !dat2[, U_Assembler[m1]]
  keep <- dat2$BioSample[keep]
  keep <- (dat1$BioSample %in% keep) & dat1$Assembler == U_Assembler[m1]
  lines(x = sort(ISdev[keep]),
        y = seq_len(sum(keep)) / sum(keep),
        lty = 1,
        col = ColVec1[m1])
  
  lines(x = sort(FSdev[keep]),
        y = seq_len(sum(keep)) / sum(keep),
        lty = 2,
        col = ColVec1[m1])
}
legend("topleft",
       legend = U_Assembler,
       lty = 1,
       col = ColVec1[seq_along(U_Assembler)])

save(dat1,
     dat2,
     AssemblerCode,
     file = "~/Repos/PseudogenePlots/InputData/Neisseria_v02.RData",
     compress = "xz")

plot(x = dat2$FSperMB,
     y = dat2$Contig_N50 / dat2$Total_Nucs,
     pch = 16,
     col = "red", ylim = c(0, 0.1))
points(x = (dat1$FR / dat1$AssemblySize) * 1000000,
       y = dat1$N50 / dat1$AssemblySize,
       pch = 16,
       col = match(x = dat1$Assembler,
                   table = unique(dat1$Assembler)))
legend("topright",
       legend = unique(dat1$Assembler),
       pch = 20,
       col = seq_along(unique(dat1$Assembler)))

plot(x = dat2$ISperMB,
     y = dat2$Contig_N50 / dat2$Total_Nucs,
     pch = 16,
     col = "red", ylim = c(0, 0.1))
points(x = (dat1$IS / dat1$AssemblySize) * 1000000,
       y = dat1$N50 / dat1$AssemblySize,
       pch = 16,
       col = match(x = dat1$Assembler,
                   table = unique(dat1$Assembler)))
legend("topright",
       legend = unique(dat1$Assembler),
       pch = 20,
       col = seq_along(unique(dat1$Assembler)))







