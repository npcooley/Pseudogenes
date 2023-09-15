###### -- reported coverages in genbank ---------------------------------------

# no required packages here, just base, though edirect needs to be accessible

query1 <- paste("esearch",
                "-db assembly",
                "-query '",
                '"Archaea"[Organism]',
                'AND "latest genbank"[properties]',
                "'",
                "|",
                "esummary",
                "|",
                "xtract -pattern DocumentSummary",
                "-element Coverage SubmissionDate Organism Taxid")

SEARCHTIMESTART <- Sys.time()
replies01 <- system(command = query1,
                    intern = TRUE,
                    timeout = 6000L)
SEARCHTIMEEND <- Sys.time()
GENBANKSEARCHTIME <- SEARCHTIMEEND - SEARCHTIMESTART
print(GENBANKSEARCHTIME)

query2 <- paste("esearch",
                "-db assembly",
                "-query '",
                '"Bacteria"[Organism]',
                'AND "latest genbank"[properties]',
                "'",
                "|",
                "esummary",
                "|",
                "xtract -pattern DocumentSummary",
                "-element Coverage SubmissionDate Organism Taxid")

SEARCHTIMESTART <- Sys.time()
replies02 <- system(command = query2,
                    intern = TRUE,
                    timeout = 12000L)
SEARCHTIMEEND <- Sys.time()
GENBANKSEARCHTIME <- SEARCHTIMEEND - SEARCHTIMESTART
print(GENBANKSEARCHTIME)

# combine results
res01 <- strsplit(x = replies01,
                  split = "\t",
                  fixed = TRUE)
res01 <- do.call(rbind,
                 res01[lengths(res01) == 4])
res02 <- data.frame("cov" = as.numeric(gsub(x = res01[, 1L],
                                            pattern = "[^0-9.]",
                                            replacement = "")),
                    "submissiondate" = res01[, 2L],
                    "sp" = res01[, 3L],
                    "taxid" = as.integer(res01[, 4L]))

res03 <- strsplit(x = replies02,
                  split = "\t",
                  fixed = TRUE)
res03 <- do.call(rbind,
                 res03[lengths(res03) == 4])
res04 <- data.frame("cov" = as.numeric(gsub(x = res03[, 1L],
                                            pattern = "[^0-9.]",
                                            replacement = "")),
                    "submissiondate" = res03[, 2L],
                    "sp" = res03[, 3L],
                    "taxid" = as.integer(res03[, 4L]))

res05 <- rbind(res02,
               res04)

# sum(res05$cov < 10) / nrow(res05)
# [1] 0.05439376
# sum(res05$cov <= 10) / nrow(res05)
# [1] 0.0760108
# sum(res05$cov < 20) / nrow(res05)
# [1] 0.1206815
# sum(res05$cov <= 20) / nrow(res05)
# [1] 0.1263315

save(replies01,
     replies02,
     res05,
     file = "~/Repos/PseudogenePlots/InputData/Genbank_Summary.RData",
     compress = "xz")


