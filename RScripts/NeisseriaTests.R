###### -- Data from neisseria reassemblies ------------------------------------


ColVec1 <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000')
ColVec2 <- paste0(ColVec1,
                  "33")

load(file = "~/Repos/20221207_FactorialAssembly/ParsedPseudogenes.RData",
     verbose = TRUE)

# select out anomalously sized assemblies and retain only trimmed assemblies
pv2 <- median(ResTable$AssemblySize)
pv1 <- ResTable$AssemblySize > (pv2 * 0.75) & ResTable$AssemblySize < (pv2 * 1.25)  # assembly is not anomalously sized
pv2 <- ResTable$Trim
TableForPlots <- ResTable[pv1 & pv2, ]

FSperMB <- (TableForPlots$FR / TableForPlots$AssemblySize) * 1000000
ISperMB <- (TableForPlots$IS / TableForPlots$AssemblySize) * 1000000


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

# 
# KEY01 <- c(NA,
#            NA,
#            "Pseudogenes per Mb")
# KEY03 <- list(c("Megahit", "Unicycler", "SKESA", "SPADES"), # all four assemblers
#               c(2, 3, 4), # quality bins
#               c(0.05, 0.25, 0.5, 1)) # down sampling
# KEY06 <- c("Assembler",
#            "QualBin",
#            "DownSample")
# KEY07 <- list(c("Megahit", "Unicycler", "SKESA", "SPADES"), # all four assemblers
#               c("20-30", "30+", "All"), # quality bins
#               c("0.05x", "0.25x", "0.50x", "1.00x")) # sparsification

# vals for assembler
# pin the columns that aren't the target
A_counts <- matmatch(z1 = TableForPlots[, c("QualBin", "DownSample", "Run")])
A_select <- A_counts == max(A_counts)
U_Assembler <- unique(TableForPlots$Assembler)
# vals for downsampling
# pin the columns that aren't the target
D_counts <- matmatch(z1 = TableForPlots[, c("Assembler", "QualBin", "Run")])
D_select <- D_counts == max(D_counts)
U_Down <- unique(TableForPlots$DownSample)

load(file = "~/Repos/Pseudogenes/SearchResults.RData",
     verbose = TRUE)
load(file = "~/Repos/Pseudogenes/SearchResults2.RData",
     verbose = TRUE)

g1 <- unlist(regmatches(x = EntrezResults$SpeciesName,
                        m = gregexpr(pattern = "^[^ ]+",
                                     text = EntrezResults$SpeciesName)))
g2 <- grep(pattern = "candidatus",
           ignore.case = TRUE,
           x = g1)
g1[g2] <- unlist(regmatches(x = EntrezResults$SpeciesName[g2],
                            m = gregexpr(pattern = "(?<=Candidatus )([^ ]+)",
                                         text = EntrezResults$SpeciesName[g2],
                                         perl = TRUE)))
g2 <- grep(pattern = "[^A-Za-z]",
           x = g1)
g1[g2] <- unlist(regmatches(x = g1[g2],
                            m = gregexpr(pattern = "[A-Za-z]+",
                                         text = g1[g2])))

w2 <- g1 == "Neisseria"
w3 <- grepl(pattern = "gonorrhoeae",
            x = EntrezResults$SpeciesName)
w4 <- grepl(pattern = "meningitidis",
            x = EntrezResults$SpeciesName)
w5 <- match(x = TableForPlots$Run,
            table = SRAResults$Run)
w6 <- grepl(pattern = "gonorrhoeae",
            x = SRAResults$ScientificName[w5])
w7 <- grepl(pattern = "meningitidis",
            x = SRAResults$ScientificName[w5])

TableForPlots <- cbind(TableForPlots,
                       "species" = ifelse(test = w6 & !w7,
                                          yes = "gonorrhoeae",
                                          no = ifelse(test = !w6 & w7,
                                                      yes = "meningitidis",
                                                      no = "other")))


dat2 <- data.frame("BioSample" = EntrezResults$Biosample[(w2 & w3) | (w2 & w4)],
                   "SpeciesName" = EntrezResults$SpeciesName[(w2 & w3) | (w2 & w4)],
                   "TaxID" = EntrezResults$TaxID[(w2 & w3) | (w2 & w4)],
                   "Submitter" = EntrezResults$Submitter_Org[(w2 & w3) | (w2 & w4)],
                   "ISperMB" = (PseudosCondensed[(w2 & w3) | (w2 & w4), 4L] / EntrezResults$Total_Length[(w2 & w3) | (w2 & w4)]) * 1000000,
                   "FSperMB" = (PseudosCondensed[(w2 & w3) | (w2 & w4), 1L] / EntrezResults$Total_Length[(w2 & w3) | (w2 & w4)]) * 1000000)


save(TableForPlots,
     A_select,
     U_Assembler,
     D_select,
     U_Down,
     ISperMB,
     FSperMB,
     dat2,
     file = "~/Repos/PseudogenePlots/InputData/Neisseria_v01.RData",
     compress = "xz")

layout(mat = matrix(data = c(1,2),
                    nrow = 1))
par(mar = c(4,3,2,0.75),
    mgp = c(2, 0.75, 0))
plot(x = 0,
     y = 0,
     type = "n",
     xlab = "PG per MB",
     ylab = "Frequency",
     ylim = c(0, 1),
     xlim = c(0, 60))
for (m1 in seq_along(U_Assembler)) {
  pv3 <- TableForPlots$Assembler == U_Assembler[m1]
  
  lines(x = sort(ISperMB[pv3 & A_select]),
        y = seq(sum(pv3 & A_select)) / sum(pv3 & A_select),
        col = ColVec1[m1],
        lty = 1)
  
  lines(x = sort(FSperMB[pv3 & A_select]),
        y = seq(sum(pv3 & A_select)) / sum(pv3 & A_select),
        col = ColVec1[m1],
        lty = 3)
}

plot(x = 0,
     y = 0,
     type = "n",
     xlab = "PG per MB",
     ylab = "Frequency",
     ylim = c(0, 1),
     xlim = c(0, 60))
for (m1 in seq_along(U_Down)) {
  pv3 <- TableForPlots$DownSample == U_Down[m1]
  
  lines(x = sort(ISperMB[pv3 & D_select]),
        y = seq(sum(pv3 & D_select)) / sum(pv3 & D_select),
        col = ColVec1[m1],
        lty = 1)
  
  lines(x = sort(FSperMB[pv3 & D_select]),
        y = seq(sum(pv3 & D_select)) / sum(pv3 & D_select),
        col = ColVec1[m1],
        lty = 3)
}




