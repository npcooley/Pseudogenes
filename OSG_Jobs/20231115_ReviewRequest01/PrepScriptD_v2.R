###### -- build job d map and things ------------------------------------------

suppressMessages(library(SynExtend))

targetdir01 <- "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01"
targetdir02 <- "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/Annotations"

load(file = paste0(targetdir01,
                   "/PGAPJobTables.RData"),
     verbose = TRUE)
load(file = paste0(targetdir01,
                   "/PGAPJobTables_v2.RData"),
     verbose = TRUE)
load(file = paste0(targetdir01,
                   "/JobMap_A_B.RData"),
     verbose = TRUE)
load(file = paste0(targetdir01,
                   "/ComparisonTables.RData"),
     verbose = TRUE)

files01 <- list.files(path = targetdir01,
                      pattern = "ResultC[0-9]+\\.fna\\.gz")
files02 <- list.files(path = targetdir02,
                      pattern = "anno[0-9]+\\.gff\\.gz")
files03 <- list.files(path = targetdir01,
                      pattern = "ResultD[0-9]+\\.RData")

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

all_requested_annotations <- rbind(assemblies4,
                                   pgapjobmap2)
rownames(all_requested_annotations) <- NULL


# these comparisons are relatively cheap, so i can just ask for all the missing ones

res07 <- vector(mode = "list",
                length = length(pin_assemblies02))
pBar <- txtProgressBar(style = 1)
PBAR <- length(pin_assemblies02)

for (m1 in seq_along(pin_assemblies01)) {
  
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
    
    res07[[m1]] <- data.frame("ref_assembly" = rep(refassemblyfile,
                                                   length(files03)),
                              "ref_annotation" = rep(refannotationfile,
                                                     length(files03)),
                              "test_assemblies" = files03,
                              "test_annotations" = files04)
  }
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}

res07 <- do.call(rbind,
                 res07)
res07 <- cbind("PersistentID" = seq(nrow(res07)),
               res07)

save(res07,
     all_requested_annotations,
     file = paste0(targetdir01,
                   "/ComparisonTables_v2.RData"),
     compress = "xz")

write.table(x = res07,
            row.names = FALSE,
            col.names = FALSE,
            append = FALSE,
            quote = FALSE,
            file = paste0(targetdir01,
                          "/JobMapD_v2.txt"))
write.table(x = res06[sample(x = nrow(res07),
                             size = 10,
                             replace = FALSE), ],
            row.names = FALSE,
            col.names = FALSE,
            append = FALSE,
            quote = FALSE,
            file = paste0(targetdir01,
                          "/TestMapD_v2.txt"))

