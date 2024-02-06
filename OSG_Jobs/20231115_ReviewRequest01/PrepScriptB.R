###### -- prep job b ----------------------------------------------------------
# given successfully completed job a jobs subset the planned job b jobs
# to only those whose 'pin' SRA run was originally successful AND hit close enough
# to the target original read coverage

suppressMessages(library(SynExtend))

TargetDir <- "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01"

load(file = "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/JobMap_A_B.RData",
     verbose = TRUE)
load(file = paste0(TargetDir,
                   "/JobMap_A_B.RData"),
     verbose = TRUE)

files01 <- list.files(path = TargetDir,
                      pattern = "Result[0-9]+\\.fna\\.gz")
files02 <- list.files(path = TargetDir,
                      pattern = "Result[0-9]+\\.RData")

w1 <- as.integer(unlist(regmatches(x = files01,
                                   m = gregexpr(pattern = "(?<=Result)([0-9]+)(?=\\.fna)",
                                                text = files01,
                                                perl = TRUE))))
w2 <- as.integer(rownames(res04)) %in% w1
w3 <- res04$sr_r[w2]

pBar <- txtProgressBar(style = 1)
PBAR <- length(files01)
sr_vec <- lr_vec <- vector(mode = "numeric",
                           length = PBAR)

for (m1 in seq_len(PBAR)) {
  load(file = paste0(TargetDir,
                     "/",
                     files02[m1]),
       verbose = FALSE)
  
  z1 <- rle(dep1$V1)
  wtm1 <- weighted.mean(x = dep1$V3,
                        w = rep(x = (z1$lengths / nrow(dep1)),
                                times = z1$lengths))
  sr_dep <- tapply(X = dep1$V3,
                   INDEX = dep1$V1,
                   FUN = function(x) {
                     mean(x)
                   })
  
  if (!is.null(dep2)) {
    z1 <- rle(dep2$V1)
    wtm2 <- weighted.mean(x = dep2$V3,
                          w = rep(x = (z1$lengths / nrow(dep2)),
                                  times = z1$lengths))
    lr_dep <- tapply(X = dep2$V3,
                     INDEX = dep2$V1,
                     FUN = function(x) {
                       mean(x)
                     })
  } else {
    wtm2 <- NULL
    lr_dep <- NULL
    
  }
  rm(list = c("dep1",
              "dep2"))
  
  sr_vec[m1] <- wtm1
  if (!is.null(wtm2)) {
    lr_vec[m1] <- wtm2
  }
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}

res05 <- res04[res04$sr_r %in% w3, ]
res05 <- res05[sr_vec >= 950 &
                 lr_vec <= 1200, ]

res06 <- res03[res03$sr_r %in% res05$sr_r, ]


write.table(x = res06[, 1:9],
            row.names = TRUE,
            col.names = FALSE,
            append = FALSE,
            quote = FALSE,
            file = "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/JobMapB.txt")

write.table(x = res06[sample(x = nrow(res06),
                             size = 10,
                             replace = FALSE), 1:9],
            row.names = TRUE,
            col.names = FALSE,
            append = FALSE,
            quote = FALSE,
            file = "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/TestMapB.txt")

