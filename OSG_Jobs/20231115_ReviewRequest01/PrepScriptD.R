###### -- build job d map and things ------------------------------------------

suppressMessages(library(SynExtend))

targetdir <- "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01"


load(file = paste0(targetdir,
                   "/PGAPJobTables.RData"),
     verbose = TRUE)
load(file = paste0(targetdir,
                   "/JobMap_A_B.RData"),
     verbose = TRUE)

files01 <- list.files(path = targetdir,
                      pattern = "ResultC[0-9]+\\.fna\\.gz")
files02 <- list.files(path = targetdir,
                      pattern = "anno[0-9]+\\.gff\\.gz")

# res03 and res04 rownames are related to the persistentid column in assemblies4
# return a structure that allows us to run all the correct comparisons of the
# "test" assemblies vs their "pin" assemblies

completed_assemblies <- as.integer(unlist(regmatches(x = files01,
                                                     m = gregexpr(pattern = "[0-9]+",
                                                                  text = files01))))
completed_annotations <- as.integer(unlist(regmatches(x = files02,
                                                      m = gregexpr(pattern = "[0-9]+",
                                                                   text = files02))))
pin_assemblies01 <- as.integer(rownames(res04))
test_assemblies01 <- as.integer(rownames(res03))

pin_assemblies02 <- pin_assemblies01[pin_assemblies01 %in% completed_assemblies &
                                       pin_assemblies01 %in% completed_annotations]
test_assemblies02 <- test_assemblies01[test_assemblies01 %in% completed_assemblies &
                                         test_assemblies01 %in% completed_annotations]
res05 <- res03[as.integer(rownames(res03)) %in% test_assemblies02, ]
test_assemblies03 <- as.integer(rownames(res05))

res06 <- vector(mode = "list",
                length = length(pin_assemblies02))
pBar <- txtProgressBar(style = 1)
PBAR <- length(pin_assemblies02)

for (m1 in seq_along(pin_assemblies02)) {
  # reference back to the original table, not all of which produced the expected files
  # on the grid
  w1 <- which(pin_assemblies01 == pin_assemblies02[m1])
  o1 <- pin_assemblies02[m1]
  pin_run <- res04$sr_r[w1]
  pin_cov_1 <- res04$sr_cov_target[w1]
  pin_cov_2 <- res04$lr_cov_target[w1]
  
  w2 <- res05$sr_r == pin_run & res05$lr_cov_target == pin_cov_2
  w3 <- res05[w2, ]
  rel_assemblies <- as.integer(rownames(w3))
  
  if (length(rel_assemblies) > 0) {
    # build out the a DF in the form of:
    # filename for reference GFF, filename for reassembly GFF, filename for ref assembly fna,
    # filename for ressembly fna
    refassemblyfile <- files01[completed_assemblies == o1]
    refannotationfile <- files02[completed_annotations == o1]
    
    files03 <- files01[match(x = rel_assemblies,
                             table = completed_assemblies)]
    files04 <- files02[match(x = rel_assemblies,
                             table = completed_annotations)]
    
    res06[[m1]] <- data.frame("ref_assembly" = rep(refassemblyfile,
                                                   length(files03)),
                              "ref_annotation" = rep(refannotationfile,
                                                     length(files03)),
                              "test_assemblies" = files03,
                              "test_annotations" = files04)
  }
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}

res06 <- do.call(rbind,
                 res06)

res06 <- cbind("PersistentID" = seq(nrow(res06)),
               res06)

save(res06,
     file = paste0(targetdir,
                   "/ComparisonTables.RData"),
     compress = "xz")

write.table(x = res06,
            row.names = FALSE,
            col.names = FALSE,
            append = FALSE,
            quote = FALSE,
            file = paste0(targetdir,
                          "/JobMapD.txt"))
write.table(x = res06[sample(x = nrow(res06),
                             size = 10,
                             replace = FALSE), ],
            row.names = FALSE,
            col.names = FALSE,
            append = FALSE,
            quote = FALSE,
            file = paste0(targetdir,
                          "/TestMapD.txt"))

