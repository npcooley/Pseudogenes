###### -- re-assembling components of assemblies ------------------------------

suppressMessages(library(SynExtend))


RefTable1 <- "~/Repos/Pseudogenes/SearchResults.RData"
RefTable2 <- "~/Repos/Pseudogenes/SearchResults2.RData"

load(file = RefTable1,
     verbose = TRUE)
load(file = RefTable2,
     verbose = TRUE)

# keep only illumina and ONT, pacbio doesn't produce meaningful quality scores
Keep1 <- c("ILLUMINA",
           "OXFORD_NANOPORE",
           "PACBIO_SMRT",
           "LS454")
SRA <- SRAResults[SRAResults$Platform %in% Keep1 &
                    SRAResults$LibraryStrategy == "WGS" &
                    SRAResults$LibrarySource == "GENOMIC", ]
Keep2 <- table(SRA$BioSample)
Keep2 <- names(Keep2[Keep2 >= 2L])
SRA <- SRA[SRA$BioSample %in% Keep2, ]

# build a job map
JobMap <- data.frame("PersistentID" = seq(nrow(SRA)),
                     "SRR" = SRA$Run,
                     "BioSample" = SRA$BioSample,
                     "Platform" = SRA$Platform)


write.table(x = JobMap,
            file = "~/Repos/20230117_MultiRunAssemblies/JobMap_Assembly.txt",
            append = FALSE,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


write.table(x = JobMap[sample(x = seq(nrow(JobMap)),
                              size = 10L,
                              replace = FALSE), ],
            file = "~/Repos/20230117_MultiRunAssemblies/TestMap_Assembly.txt",
            append = FALSE,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


# how many assemblies from paired biosamples are short/long vs short/short vs long/long

# z1 <- table(SRA$BioSample)
# z1 <- names(z1[z1 == 2L])
# z1 <- unique(SRA$BioSample)
# 
# z2 <- sapply(X = z1,
#              FUN = function(x) {
#                z2 <- unique(SRA$Platform[SRA$BioSample == x])
#                if (length(z2) == 1L) {
#                  if (z2 == "ILLUMINA") {
#                    z2 <- 1L
#                  } else if (z2 == "OXFORD_NANOPORE") {
#                    z2 <- 2L
#                  }
#                } else {
#                  z2 <- 3L
#                }
#                return(z2)
#              },
#              simplify = TRUE,
#              USE.NAMES = FALSE)

# z2 <- sapply(X = z1,
#              FUN = function(x) {
#                z2 <- unique(SRA$Platform[SRA$BioSample == x])
#                if (length(z2) == 1L) {
#                  if (z2 == "ILLUMINA") {
#                    z2 <- 1L
#                  } else if (z2 == "PACBIO_SMRT") {
#                    z2 <- 2L
#                  } else {
#                    z2 <- 3L
#                  }
#                } else {
#                  z2 <- 4L
#                }
#                return(z2)
#              },
#              simplify = TRUE,
#              USE.NAMES = FALSE)

# z2
# 1    2    3    4
# 3708  325   12 4199









