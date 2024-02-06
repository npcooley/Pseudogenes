###### -- prepare jobs for PGAP -----------------------------------------------
# take in the results of both the 'pin assemblies'
# and the subset assemblies
# build a job map to annotate all of them
# this includes building out the yaml files for PGAP

suppressMessages(library(SynExtend))
suppressMessages(library(yaml))

TargetDir <- "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01"
TargetDir2 <- "~/Repos/PseudogenePlots/InputData"
YAMLDir <- paste0(TargetDir,
                  "/YAMLFiles")
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
                      pattern = "ResultB[0-9]+\\.RData")
files02 <- list.files(path = TargetDir,
                      pattern = "ResultB[0-9]+\\.fna\\.gz")

completedids01 <- as.integer(unlist(regmatches(x = files01,
                                               m = gregexpr(pattern = "(?<=ResultB)([0-9]+)(?=\\.RData)",
                                                            text = files01,
                                                            perl = TRUE))))
completedids02 <- as.integer(unlist(regmatches(x = files02,
                                               m = gregexpr(pattern = "(?<=ResultB)([0-9]+)(?=\\.fna\\.gz)",
                                                            text = files02,
                                                            perl = TRUE))))

totalexpectedids <- max(unique(c(completedids01,
                                 completedids02)))
completedids03 <- seq(totalexpectedids)[seq(totalexpectedids) %in% completedids01 & seq(totalexpectedids) %in% completedids02]

files03 <- files01[completedids01 %in% completedids03]
files04 <- files02[completedids02 %in% completedids03]

pBar <- txtProgressBar(style = 1)
PBAR <- length(completedids03)
val1 <- val2 <- val6 <- vector(mode = "character",
                               length = PBAR)
val3 <- val4 <- vector(mode = "numeric",
                       length = PBAR)
val5 <- vector(mode = "integer",
               length = PBAR)

# for (m1 in 5216:length(files03)) {
for (m1 in seq_along(files03)) {
  # i haven't had this problem in a while, but sometimes files get corrupted on the grid
  # this is relatively rare, but annoying,
  # wrap our load call inside try and skip files that error out when loaded
  x <- try(load(file = paste0(TargetDir,
                              "/",
                              files03[m1]),
                verbose = FALSE),
           silent = TRUE)
  if (is(object = x,
         class2 = "try-error")) {
    next
  }
  fna1 <- readDNAStringSet(filepath = paste0(TargetDir,
                                             "/",
                                             files04[m1]))
  
  w2 <- width(fna1) >= 200
  fna2 <- fna1[w2]
  
  dep3 <- dep1[dep1$V1 %in% which(w2), ]
  if (!is.null(dep2)) {
    dep4 <- dep2[dep2$V1 %in% which(w2), ]
  } else {
    dep4 <- NULL
  }
  
  z1 <- rle(dep3$V1)
  wtm1 <- weighted.mean(x = dep3$V3,
                        w = rep(x = (z1$lengths / nrow(dep3)),
                                times = z1$lengths))
  sr_dep <- tapply(X = dep3$V3,
                   INDEX = dep3$V1,
                   FUN = function(x) {
                     mean(x)
                   })
  
  if (!is.null(dep2)) {
    z1 <- rle(dep4$V1)
    wtm2 <- weighted.mean(x = dep4$V3,
                          w = rep(x = (z1$lengths / nrow(dep4)),
                                  times = z1$lengths))
    lr_dep <- tapply(X = dep4$V3,
                     INDEX = dep4$V1,
                     FUN = function(x) {
                       mean(x)
                     })
    val4[m1] <- wtm2
  } else {
    wtm2 <- NULL
    lr_dep <- NULL
    val4[m1] <- -1
  }
  val1[m1] <- gsub(pattern = "(?<=Result)(B)(?=0)",
                   replacement = "C",
                   perl = TRUE,
                   x = files02[m1])
  val2[m1] <- gsub(pattern = "\\.gz",
                   replacement = "",
                   x = val1[m1])
  val3[m1] <- wtm1
  val5[m1] <- completedids03[m1]
  val6[m1] <- gsub(pattern = "(?<=Result)(B)(?=0)",
                   replacement = "C",
                   perl = TRUE,
                   x = files03[m1])
  
  save(dep1,
       dep2,
       dep3,
       dep4,
       fna1,
       w2,
       unicycler_names,
       UnicyclerTotalTime,
       PullReadsTotalTime,
       file = paste0(TargetDir,
                     "/",
                     val6[m1]),
       compress = "xz")
  
  writeXStringSet(x = fna2,
                  filepath = paste0(TargetDir,
                                    "/",
                                    val1[m1]),
                  compress = TRUE)
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)

pgapjobmap <- data.frame("PersistentID" = val5,
                         "fn1" = val1,
                         "fn2" = val2,
                         "wtm1" = val3,
                         "wtm2" = val4)

save(pgapjobmap,
     file = paste0(TargetDir,
                   "/PGAPJobTables.RData"),
     compress = "xz")

# do the same thing with the "pin assemblies"
files01 <- list.files(path = TargetDir,
                      pattern = "Result[0-9]+\\.RData")
files02 <- list.files(path = TargetDir,
                      pattern = "Result[0-9]+\\.fna\\.gz")

completedids01 <- as.integer(unlist(regmatches(x = files01,
                                               m = gregexpr(pattern = "(?<=Result)([0-9]+)(?=\\.RData)",
                                                            text = files01,
                                                            perl = TRUE))))
completedids02 <- as.integer(unlist(regmatches(x = files02,
                                               m = gregexpr(pattern = "(?<=Result)([0-9]+)(?=\\.fna\\.gz)",
                                                            text = files02,
                                                            perl = TRUE))))

totalexpectedids <- max(unique(c(completedids01,
                                 completedids02)))
completedids03 <- seq(totalexpectedids)[seq(totalexpectedids) %in% completedids01 & seq(totalexpectedids) %in% completedids02]

files03 <- files01[completedids01 %in% completedids03]
files04 <- files02[completedids02 %in% completedids03]

pBar <- txtProgressBar(style = 1)
PBAR <- length(completedids03)
val1 <- val2 <- val6 <- vector(mode = "character",
                               length = length(completedids03))
val3 <- val4 <- vector(mode = "numeric",
                       length = length(completedids03))
val5 <- vector(mode = "integer",
               length = length(completedids03))

for (m1 in seq_along(files03)) {
  x <- try(load(file = paste0(TargetDir,
                              "/",
                              files03[m1]),
                verbose = FALSE),
           silent = TRUE)
  
  if (is(object = x,
         class2 = "try-error")) {
    next
  }
  
  fna1 <- readDNAStringSet(filepath = paste0(TargetDir,
                                             "/",
                                             files04[m1]))
  
  w2 <- width(fna1) >= 200
  fna2 <- fna1[w2]
  
  dep3 <- dep1[dep1$V1 %in% which(w2), ]
  if (!is.null(dep2)) {
    dep4 <- dep2[dep2$V1 %in% which(w2), ]
  } else {
    dep4 <- NULL
  }
  
  z1 <- rle(dep3$V1)
  wtm1 <- weighted.mean(x = dep3$V3,
                        w = rep(x = (z1$lengths / nrow(dep3)),
                                times = z1$lengths))
  sr_dep <- tapply(X = dep3$V3,
                   INDEX = dep3$V1,
                   FUN = function(x) {
                     mean(x)
                   })
  
  if (!is.null(dep2)) {
    z1 <- rle(dep4$V1)
    wtm2 <- weighted.mean(x = dep4$V3,
                          w = rep(x = (z1$lengths / nrow(dep4)),
                                  times = z1$lengths))
    lr_dep <- tapply(X = dep4$V3,
                     INDEX = dep4$V1,
                     FUN = function(x) {
                       mean(x)
                     })
    val4[m1] <- wtm2
  } else {
    wtm2 <- NULL
    lr_dep <- NULL
    val4[m1] <- -1
  }
  val1[m1] <- gsub(pattern = "Result",
                   replacement = "ResultC",
                   perl = TRUE,
                   x = files02[m1])
  val2[m1] <- gsub(pattern = "\\.gz",
                   replacement = "",
                   x = val1[m1])
  val3[m1] <- wtm1
  val5[m1] <- completedids03[m1]
  val6[m1] <- gsub(pattern = "Result",
                   replacement = "ResultC",
                   perl = TRUE,
                   x = files03[m1])
  
  save(dep1,
       dep2,
       dep3,
       dep4,
       fna1,
       w2,
       unicycler_names,
       UnicyclerTotalTime,
       PullReadsTotalTime,
       file = paste0(TargetDir,
                     "/",
                     val6[m1]),
       compress = "xz")
  
  writeXStringSet(x = fna2,
                  filepath = paste0(TargetDir,
                                    "/",
                                    val1[m1]),
                  compress = TRUE)
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
  
}
close(pBar)

pinassemblies <- data.frame("PersistentID" = val5,
                            "fn1" = val1,
                            "fn2" = val2,
                            "wtm1" = val3,
                            "wtm2" = val4)

save(pinassemblies,
     pgapjobmap,
     file = paste0(TargetDir,
                   "/PGAPJobTables.RData"),
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

# drop anyone that was weird

assemblies1 <- pgapjobmap[nchar(pgapjobmap$fn1) > 0, ]
assemblies2 <- pinassemblies[nchar(pinassemblies$fn1) > 0, ]

assemblies3 <- rbind(assemblies1,
                     assemblies2)
runvec <- vector(mode = "character",
                 length = nrow(assemblies3))
rowvals1 <- as.integer(rownames(res03))
rowvals2 <- as.integer(rownames(res04))
for (m1 in seq_along(runvec)) {
  if (assemblies3$PersistentID[m1] %in% rowvals1) {
    w1 <- match(x = assemblies3$PersistentID[m1],
                table = rowvals1)
    if (length(w1) != 1 |
        is.na(w1)) {
      stop("check here")
    }
    runvec[m1] <- res03$sr_r[w1]
  } else {
    w1 <- match(x = assemblies3$PersistentID[m1],
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
                  formatC(x = assemblies4$PersistentID,
                          format = "d",
                          flag = 0,
                          width = 6),
                  ".gff")

assemblies4 <- cbind(assemblies3,
                     "submolfile" = targetfile,
                     "outfile" = outfile,
                     "sci_name" = u2)

save(assemblies4,
     file = paste0(TargetDir,
                   "/PGAPJobTables.RData"),
     compress = "xz")

write.table(x = assemblies4[, -ncol(assemblies4)],
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = FALSE,
            file = paste0(TargetDir,
                          "/JobMapC.txt"))
write.table(x = assemblies4[sample(x = nrow(assemblies4),
                                   size = 10,
                                   replace = FALSE),
                            -ncol(assemblies4)],
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = FALSE,
            file = paste0(TargetDir,
                          "/TestMapC.txt"))


