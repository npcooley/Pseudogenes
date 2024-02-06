###### -- submit the second pass jobs for annotation --------------------------
# this'll take a hot second

suppressMessages(library(SynExtend))
suppressMessages(library(yaml))

TargetDir <- "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01"
TargetDir2 <- "~/Repos/PseudogenePlots/InputData"
YAMLDir <- paste0(TargetDir,
                  "/YAMLFiles_v2")
load(file = paste0(TargetDir,
                   "/JobMap_A_B.RData"),
     verbose = TRUE)

# we need the names for these assemblies, and forgot to save them initially
load(file = paste0(TargetDir2,
                   "/GenBank_Assembly_Summary.RData"),
     verbose = TRUE)
load(file = paste0(TargetDir2,
                   "/GenBank_SRA_Summary.RData"),
     verbose = TRUE)


files01 <- list.files(path = TargetDir,
                      pattern = "ResultC[0-9]+\\.RData")
files02 <- list.files(path = TargetDir,
                      pattern = "ResultC[0-9]+\\.fna\\.gz")
files03 <- list.files(path = TargetDir,
                      pattern = "anno[0-9]+\\.gff\\.gz")

# create a job map that builds out jobs for all the persistent IDs that don't have annotations

completedids01 <- as.integer(unlist(regmatches(x = files01,
                                               m = gregexpr(pattern = "(?<=ResultC)([0-9]+)(?=\\.RData)",
                                                            text = files01,
                                                            perl = TRUE))))
completedids02 <- as.integer(unlist(regmatches(x = files02,
                                               m = gregexpr(pattern = "(?<=ResultC)([0-9]+)(?=\\.fna\\.gz)",
                                                            text = files02,
                                                            perl = TRUE))))
completedids03 <- as.integer(unlist(regmatches(x = files03,
                                               m = gregexpr(pattern = "(?<=anno)([0-9]+)(?=\\.gff\\.gz)",
                                                            text = files03,
                                                            perl = TRUE))))

totalexpectedids <- max(unique(c(completedids01,
                                 completedids02)))

missingannotationjobs <- completedids01[!(completedids01 %in% completedids03)]
files04 <- files01[!(completedids01 %in% completedids03)]
files05 <- files02[!(completedids01 %in% completedids03)]

# still need to read in and get the weighted means


pBar <- txtProgressBar(style = 1)
PBAR <- length(missingannotationjobs)
val3 <- val4 <- vector(mode = "numeric",
                       length = PBAR)
TimeStart <- Sys.time()
for (m1 in seq_along(missingannotationjobs)) {
  x <- try(load(file = paste0(TargetDir,
                              "/",
                              files04[m1]),
                verbose = FALSE),
           silent = TRUE)
  if (is(object = x,
         class2 = "try-error")) {
    next
  }
  
  z1 <- rle(dep3$V1)
  wtm1 <- weighted.mean(x = dep3$V3,
                        w = rep(x = (z1$lengths / nrow(dep3)),
                                times = z1$lengths))
  val3[m1] <- wtm1
  # sr_dep <- tapply(X = dep3$V3,
  #                  INDEX = dep3$V1,
  #                  FUN = function(x) {
  #                    mean(x)
  #                  })
  
  if (!is.null(dep2)) {
    z1 <- rle(dep4$V1)
    wtm2 <- weighted.mean(x = dep4$V3,
                          w = rep(x = (z1$lengths / nrow(dep4)),
                                  times = z1$lengths))
    # lr_dep <- tapply(X = dep4$V3,
    #                  INDEX = dep4$V1,
    #                  FUN = function(x) {
    #                    mean(x)
    #                  })
    val4[m1] <- wtm2
  } else {
    wtm2 <- NULL
    # lr_dep <- NULL
    val4[m1] <- -1
  }
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)
TimeStop <- Sys.time()
print(TimeStop - TimeStart)

# build the jobmap
pgapjobmap <- data.frame("PersistentID" = missingannotationjobs,
                         "fn1" = files05,
                         "fn2" = substr(x = files05,
                                        start = 1,
                                        stop = 18),
                         "wtm1" = val3,
                         "wtm2" = val4)

save(pgapjobmap,
     file = paste0(TargetDir,
                   "/PGAPJobTables_v2.RData"),
     compress = "xz")


# write a single controller file and multiple submol files
# that specify the genus and species
# after the input assembly is unzipped,
# it needs to be renamed to 'assembly.fna'
# the submol file needs to be similarly renamed on the node

if (!dir.exists(YAMLDir)) {
  dir.create(YAMLDir)
}

Y1 <- as.yaml(list("fasta" = list("class" = "File",
                                  "location" = "assembly.fna"),
                   "submol" = list("class" = "File",
                                   "location" = "/srv/submol.yaml"),
                   "supplemental_data" = "{ class: Directory, location: /srv/input-2022-04-14.build6021 }",
                   "report_usage" = "true",
                   "ignore_all_errors" = "true"))
# Y1 <- str_replace_all(string = Y1,
#                       pattern = "'",
#                       replacement = "")
Y1 <- gsub(pattern = "'",
           replacement = "",
           x = Y1)

cat(Y1,
    file = paste0(YAMLDir,
                  "/controller.yaml"))

# RData files that successfully loaded
w2 <- which(pgapjobmap$wtm1 != 0)
pgapjobmap <- pgapjobmap[w2, ]

runvec <- vector(mode = "character",
                 length = nrow(pgapjobmap))
rowvals1 <- as.integer(rownames(res03))
rowvals2 <- as.integer(rownames(res04))
for (m1 in seq_along(runvec)) {
  if (pgapjobmap$PersistentID[m1] %in% rowvals1) {
    w1 <- match(x = pgapjobmap$PersistentID[m1],
                table = rowvals1)
    if (length(w1) != 1 |
        is.na(w1)) {
      stop("check here")
    }
    runvec[m1] <- res03$sr_r[w1]
  } else {
    w1 <- match(x = pgapjobmap$PersistentID[m1],
                table = rowvals2)
    if (length(w1) != 1) {
      stop("check here")
    }
    runvec[m1] <- res04$sr_r[w1]
  }
}

u1 <- match(x = runvec,
            table = SRA_Meta$Run)
u2 <- SRA_Meta$ScientificName[u1]
g1 <- grepl(pattern = "^candidatus",
            ignore.case = TRUE,
            x = u2)
if (any(g1)) {
  u2[g1] <- gsub(pattern = "^candidatus ",
                 x = u2[g1],
                 replacement = "",
                 ignore.case = TRUE)
}
u2 <- gsub(pattern = "'",
           replacement = "",
           x = u2)
u2 <- gsub(pattern = '"',
           replacement = "",
           x = u2)
u2 <- gsub(pattern = "[",
           replacement = "",
           x = u2,
           fixed = TRUE)
u2 <- gsub(pattern = "]",
           replacement = "",
           x = u2,
           fixed = TRUE)
u2 <- unlist(regmatches(m = gregexpr(pattern = "^[^ ]+ [^ ]+",
                                     text = u2),
                        x = u2))

u3 <- unique(u2)
submolfiles <- vector(mode = "character",
                      length = length(u3))

for (m1 in seq_along(u3)) {
  Y2 <- as.yaml(list("topology" = "linear",
                     "organism" = list("genus_species" = u3[m1]),
                     "authors" = list(list("author" = list("last_name" = "Cooley",
                                                           "first_name" = "Nicholas",
                                                           "middle_initial" = "P"))),
                     "contact_info" = list("first_name" = "Nicholas",
                                           "last_name" = "Cooley",
                                           "email" = "npc19@pitt.edu",
                                           "organization" = "University of Pittsburgh",
                                           "department" = "Department of Biomedical Informatics",
                                           "phone" = "412-624-5100",
                                           "street" = "100 Technology Dr.",
                                           "city" = "Pittsburgh",
                                           "postal_code" = "15219",
                                           "country" = "USA")))
  Y2 <- gsub(pattern = "(?<=: )([^'\n]+)",
             replacement = "'\\1'",
             x = Y2,
             perl = TRUE)
  
  Y3 <- paste0("submol",
               formatC(x = m1,
                       width = 5,
                       format = "d",
                       flag = 0),
               ".yaml")
  submolfiles[m1] <- Y3
  
  cat(Y2,
      file = paste0(YAMLDir,
                    "/",
                    Y3))
}

# match the unique names back

targetfile <- submolfiles[match(x = u2,
                                table = u3)]
outfile <- paste0("anno",
                  formatC(x = pgapjobmap$PersistentID,
                          format = "d",
                          flag = 0,
                          width = 6),
                  ".gff")

pgapjobmap2 <- cbind(pgapjobmap,
                     "submolfile" = targetfile,
                     "outfile" = outfile,
                     "sci_name" = u2)

save(pgapjobmap2,
     file = paste0(TargetDir,
                   "/PGAPJobTables_v2.RData"),
     compress = "xz")


write.table(x = pgapjobmap2[, -ncol(pgapjobmap2)],
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = FALSE,
            file = paste0(TargetDir,
                          "/JobMapC_v2.txt"))
write.table(x = pgapjobmap2[sample(x = nrow(pgapjobmap2),
                                   size = 10,
                                   replace = FALSE),
                            -ncol(pgapjobmap2)],
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = FALSE,
            file = paste0(TargetDir,
                          "/TestMapC_v2.txt"))
