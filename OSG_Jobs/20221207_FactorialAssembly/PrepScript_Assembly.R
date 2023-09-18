###### -- prepscript ----------------------------------------------------------

library(SynExtend)
load(file = "~/Repos/Pseudogenes/SearchResults2.RData",
     verbose = TRUE)
load(file = "~/Repos/Pseudogenes/SearchResults.RData",
     verbose = TRUE)



# wrangle pseudogenes into pg per kb
frameshiftperkb <- PseudosCondensed[, 1L] / (EntrezResults$Total_Length / 1000)
internalstopperkb <- PseudosCondensed[, 4] / (EntrezResults$Total_Length / 1000)
partialperkb <- PseudosCondensed[, 3L] / (EntrezResults$Total_Length / 1000)

# wrangle genus names in the assembly results, select for biosamples with
# high coverage
# g1 <- unlist(regmatches(x = EntrezResults$SpeciesName,
#                         m = gregexpr(pattern = "^[^ ]+",
#                                      text = EntrezResults$SpeciesName)))
# g2 <- grep(pattern = "candidatus",
#            ignore.case = TRUE,
#            x = g1)
# g1[g2] <- unlist(regmatches(x = EntrezResults$SpeciesName[g2],
#                             m = gregexpr(pattern = "(?<=Candidatus )([^ ]+)",
#                                          text = EntrezResults$SpeciesName[g2],
#                                          perl = TRUE,
#                                          ignore.case = TRUE)))
# g2 <- grep(pattern = "[^A-Za-z]",
#            x = g1)
# g1[g2] <- unlist(regmatches(x = g1[g2],
#                             m = gregexpr(pattern = "[A-Za-z]+",
#                                          text = g1[g2])))

b1 <- EntrezResults$Biosample[EntrezResults$Coverage >= 100]

SRAResults <- SRAResults[SRAResults$BioSample %in% b1, ]

# wrangle genus names
g1 <- unlist(regmatches(x = SRAResults$ScientificName,
                        m = gregexpr(pattern = "^[^ ]+",
                                     text = SRAResults$ScientificName)))
g2 <- grep(pattern = "candidatus",
           ignore.case = TRUE,
           x = g1)
g1[g2] <- unlist(regmatches(x = SRAResults$ScientificName[g2],
                            m = gregexpr(pattern = "(?<=Candidatus )([^ ]+)",
                                         text = SRAResults$ScientificName[g2],
                                         perl = TRUE)))
g2 <- grep(pattern = "[^A-Za-z]",
           x = g1)
g1[g2] <- unlist(regmatches(x = g1[g2],
                            m = gregexpr(pattern = "[A-Za-z]+",
                                         text = g1[g2])))

w1 <- SRAResults$Platform == "ILLUMINA"

g1 <- g1[w1]
SRAResults <- SRAResults[w1, ]

head(sort(table(g1), decreasing = TRUE), n = 100)

pin <- "Neisseria"
w2 <- which(g1 == pin)
w2 <- unique(SRAResults$BioSample[w2])
w3 <- internalstopperkb[EntrezResults$Biosample %in% w2]
w4 <- frameshiftperkb[EntrezResults$Biosample %in% w2]
w5 <- names(table(EntrezResults$SpeciesName[EntrezResults$Biosample %in% w2]))
w5 <- match(x = EntrezResults$SpeciesName[EntrezResults$Biosample %in% w2],
            table = w5)
plot(y = w3,
     x = w4,
     ylab = "internal stops / kb",
     xlab = "frameshifts / kb",
     # xlim = c(0, 0.1),
     # ylim = c(0, 0.1),
     # xaxs="i",
     # yaxs="i",
     main = paste(pin, "pseudogenes"),
     col = w5)

SelectSRA <- SRAResults[SRAResults$BioSample %in% w2, ]

JobMap <- expand.grid(list("Run" = unique(SelectSRA$Run),
                           "TRIM" = c(TRUE, FALSE),
                           "Quality" = 1:4,
                           "DownSample" = c(1, 0.5, 0.25, 0.05)))
scn1 <- unlist(regmatches(x = SelectSRA$ScientificName,
                          m = gregexpr(pattern = "^[A-Za-z]+ [A-Za-z]+",
                                       text = SelectSRA$ScientificName)))

JobMap <- cbind(JobMap,
                "ID" = seq(nrow(JobMap)),
                "NAME" = gsub(x = scn1[match(x = JobMap$Run,
                                             table = SelectSRA$Run)],
                              pattern = " ",
                              replacement = "_"))

# split jobmap into 3 groups to maybe alleviate storage issues on the OSG
chunks <- split(x = sample(x = nrow(JobMap),
                           size = nrow(JobMap),
                           replace = FALSE),
                f = factor(sort(rank(seq(nrow(JobMap)) %% 3L))))

for (m1 in seq_along(chunks)) {
  write.table(x = JobMap[chunks[[m1]], ],
              file = paste0("~/Repos/20221207_FactorialAssembly/JobMap_Assembly_",
                            m1,
                            ".txt"),
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              append = FALSE)
}

# write.table(x = JobMap,
#             file = "~/Repos/20221207_FactorialAssembly/JobMap.txt",
#             quote = FALSE,
#             row.names = FALSE,
#             col.names = FALSE,
#             append = FALSE)

TestMap <- JobMap[sample(x = seq(nrow(JobMap)),
                         size = 10,
                         replace = FALSE), ]

write.table(x = TestMap,
            file = "~/Repos/20221207_FactorialAssembly/TestMap.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = FALSE)


