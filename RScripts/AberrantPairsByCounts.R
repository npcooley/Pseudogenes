###### -- how many close assemblies are weird? --------------------------------

suppressMessages(library(SynExtend))

load("~/Repos/PseudogenePlots/InputData/Counts_Orthos_v02.RData",
     verbose = TRUE)

mat1 <- match(x = dat3$id1,
              table = as.integer(rownames(adjusted_counts)))
mat2 <- match(x = dat3$id2,
              table = as.integer(rownames(adjusted_counts)))

dat3$features1 <- adjusted_counts$all_coding[mat1]
dat3$features2 <- adjusted_counts$all_coding[mat2]

size <- ifelse(test = dat3$is1 > dat3$is2,
               yes = dat3$features1,
               no = dat3$features2) # features in assembly with more IS
rate <- ifelse(test = dat3$is1 < dat3$is2,
               yes = dat3$is1/dat3$features1,
               no = dat3$is2/dat3$features2) # rate of IS in assembly with fewer IS
num <- ifelse(test = dat3$is1 > dat3$is2,
              yes = dat3$is1,
              no = dat3$is2) # number of IS in assembly with more IS

reject <- num > qbinom(0.9975, size, rate)
sum(dat3$ANI >= 99.9)
mean(reject[dat3$ANI >= 99.9])


size <- ifelse(test = dat3$fs1 > dat3$fs2,
               yes = dat3$features1,
               no = dat3$features2) # features in assembly with more FS
rate <- ifelse(test = dat3$fs1 < dat3$fs2,
               yes = dat3$fs1/dat3$features1,
               no = dat3$is2/dat3$features2) # rate of FS in assembly with fewer FS
num <- ifelse(test = dat3$fs1 > dat3$fs2,
              yes = dat3$fs1,
              no = dat3$fs2) # number of FS in assembly with more FS

reject <- num > qbinom(0.9975, size, rate)
sum(dat3$ANI >= 99.9)
mean(reject[dat3$ANI >= 99.9])


