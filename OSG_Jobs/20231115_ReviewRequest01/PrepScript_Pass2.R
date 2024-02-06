###### -- re-cast the jobs, try and fill in gaps ------------------------------

suppressMessages(library(SynExtend))

file01 <- "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/JobMap_A_B.RData"
file02 <- "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/JobMapA.txt"
file03 <- "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/JobMapB.txt"
files01 <- list.files(path = "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01",
                      pattern = "Result[0-9]+\\.RData")
files02 <- list.files(path = "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01",
                      pattern = "ResultB[0-9]+\\.RData")

load(file = file01,
     verbose = TRUE)
tab1 <- read.table(file02)
tab2 <- read.table(file03)

completedpins <- as.integer(unlist(regmatches(x = files01,
                                              m = gregexpr(pattern = "[0-9]+",
                                                           text = files01))))
completedtests <- as.integer(unlist(regmatches(x = files02,
                                               m = gregexpr(pattern = "[0-9]+",
                                                            text = files02))))

missingpins <- tab1[!(tab1$V1 %in% completedpins), ]
missingtests <- res03[!(as.integer(rownames(res03)) %in% completedtests), 1:9]

write.table(x = missingpins,
            file = "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/JobMapA_v2.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = FALSE)
write.table(x = missingtests,
            file = "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/JobMapA_v2.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = FALSE,
            append = TRUE)
tab3 <- read.table("~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/JobMapA_v2.txt")
write.table(x = tab3[sample(x = nrow(tab3),
                            size = 10,
                            replace = FALSE), ],
            file = "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/TestMapA_v2.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = FALSE)







