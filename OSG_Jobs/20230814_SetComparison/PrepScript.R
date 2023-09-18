###### -- build jobmap --------------------------------------------------------

suppressMessages(library(SynExtend))

load(file = "~/Repos/20230814_SetComparison/SearchResults.RData",
     verbose = TRUE)


# make genus table as JobMap01
g1 <- unlist(regmatches(x = EntrezResults$SpeciesName,
                        m = gregexpr(pattern = "^[^ ]+",
                                     text = EntrezResults$SpeciesName)))
g2 <- grep(pattern = "candidatus",
           ignore.case = TRUE,
           x = g1)
g1[g2] <- unlist(regmatches(x = EntrezResults$SpeciesName[g2],
                            m = gregexpr(pattern = "(?<=Candidatus )([^ ]+)",
                                         text = EntrezResults$SpeciesName[g2],
                                         ignore.case = TRUE,
                                         perl = TRUE)))
g2 <- grep(pattern = "[^A-Za-z]",
           x = g1)
g1[g2] <- unlist(regmatches(x = g1[g2],
                            m = gregexpr(pattern = "[A-Za-z]+",
                                         text = g1[g2])))

t1 <- table(g1)
t2 <- names(t1[t1 >= 100])
t3 <- data.frame("a" = rep("Genus",
                           length(t2)),
                 "b" = t2,
                 "c" = seq(length(t2)))
t4 <- vector(mode = "list",
             length = nrow(t3))
pBar <- txtProgressBar(style = 1L)
PBAR <- nrow(t3)
# build a large job map,
# for every genus that meets the criteria, build out each pairwise comparison
for (m1 in seq_len(PBAR)) {
  w1 <- which(g1 == t3[m1, 2L])
  w2 <- sort(sample(w1,
                    size = 100L,
                    replace = FALSE))
  l01 <- length(w2)
  l02 <- l01 * (l01 - 1L) / 2L
  count <- 1L
  w3 <- vector(mode = "list",
               length = l02)
  
  for (m2 in seq_len(length(w2) - 1L)) {
    for (m3 in (m2 + 1L):length(w2)) {
      w3[[count]] <- w2[c(m2, m3)]
      count <- count + 1L
    }
  }
  
  w3 <- do.call(rbind,
                w3)
  
  t4[[m1]] <- data.frame("genus" = t3[m1, 2L],
                         "p1" = w3[, 1L],
                         "p2" = w3[, 2L])
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}

t4 <- do.call(rbind,
              t4)
t4 <- cbind(t4,
            "queue" = seq(nrow(t4)))
write.table(x = t4,
            file = "~/Repos/20230814_SetComparison/JobMap01.txt",
            quote = FALSE,
            append = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(x = t4[sample(x = nrow(t4),
                          size = 10L,
                          replace = FALSE),],
            file = "~/Repos/20230814_SetComparison/TestMap.txt",
            quote = FALSE,
            append = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# build from the taxIDs
t1 <- table(EntrezResults$TaxID)
t2 <- names(t1[t1 >= 100])
t3 <- data.frame("a" = rep("TaxID",
                           length(t2)),
                 "b" = t2,
                 "c" = seq(length(t2)))
t4 <- vector(mode = "list",
             length = nrow(t3))
pBar <- txtProgressBar(style = 1L)
PBAR <- nrow(t3)
# build a large job map,
# for every genus that meets the criteria, build out each pairwise comparison
for (m1 in seq_len(PBAR)) {
  w1 <- which(EntrezResults$TaxID == t3[m1, 2L])
  w2 <- sort(sample(w1,
                    size = 100L,
                    replace = FALSE))
  l01 <- length(w2)
  l02 <- l01 * (l01 - 1L) / 2L
  count <- 1L
  w3 <- vector(mode = "list",
               length = l02)
  
  for (m2 in seq_len(length(w2) - 1L)) {
    for (m3 in (m2 + 1L):length(w2)) {
      w3[[count]] <- w2[c(m2, m3)]
      count <- count + 1L
    }
  }
  
  w3 <- do.call(rbind,
                w3)
  
  t4[[m1]] <- data.frame("TaxID" = t3[m1, 2L],
                         "p1" = w3[, 1L],
                         "p2" = w3[, 2L])
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}

t4 <- do.call(rbind,
              t4)
t4 <- cbind(t4,
            "queue" = seq(nrow(t4)))

write.table(x = t4,
            file = "~/Repos/20220920_SetComparison/JobMap02.txt",
            quote = FALSE,
            append = FALSE,
            row.names = FALSE,
            col.names = FALSE)
