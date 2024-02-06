###### -- statistics on rates of pseudogenization -----------------------------

suppressMessages(library(SynExtend))
suppressMessages(library(data.table)) # for rbindlist() function
# suppressMessages(library(gridExtra))
suppressMessages(library(plotrix))

ColVec1 <- c('#e6194B',
             '#3cb44b',
             '#ffe119',
             '#4363d8',
             '#f58231',
             '#42d4f4',
             '#f032e6',
             '#fabed4',
             '#469990',
             '#dcbeff',
             '#9A6324',
             '#fffac8',
             '#800000',
             '#aaffc3',
             '#000075',
             '#a9a9a9',
             '#ffffff',
             '#000000')
ColVec2 <- paste0(ColVec1,
                  "33")
ColVec3 <- paste0(ColVec1,
                  "15")
ColVec4 <- paste0(ColVec1,
                  "75")

targetdir01 <- "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01"
targetdir02 <- "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/Comparisons"

load(file = paste0(targetdir01,
                   "/ComparisonTables_v2.RData"),
     verbose = TRUE)
load(file = paste0(targetdir01,
                   "/ResultE03.RData"),
     verbose = TRUE)

files01 <- list.files(path = targetdir02)
completedcomparisons <- as.integer(unlist(regmatches(x = files01,
                                                     m = gregexpr(pattern = "[0-9]+",
                                                                  text = files01))))

# we want to target comparisons between a 'true' 1000x assembly and a
# 'test' 25x coverage assembly, for every available comparison, select the replicate
# with the true coverage closest to 25x 
target_comps <- by(data = dat4,
                   INDICES = ~ index,
                   FUN = function(x) {
                     x[which.min(abs(x$sr_cov - 25)), ]
                   })
target_comps <- do.call(rbind,
                        target_comps)

ill_only_targets <- target_comps[grep(pattern = "+ x",
                                      x = rownames(target_comps)), ]
hybrid_targets <- target_comps[grep(pattern = "+ 20",
                                    x = rownames(target_comps)), ]

# we need a series of tables with attached metadata
# for each comparison we need a table of coding sequence annotations
# with their count, the number of times they're pseudogenized as IS and as FS
# for both partners of the comparison,
# and then a third table of annotations that 


PBAR <- nrow(ill_only_targets)
stats_res01 <- stats_res02 <- stats_res03 <- vector(mode = "list",
                                                    length = PBAR)

pBar <- txtProgressBar(style = 1)
tstart <- Sys.time()
for (m1 in seq_len(PBAR)) {
  
  load(file = paste0(targetdir02,
                     "/",
                     files01[completedcomparisons == ill_only_targets$comp_id[m1]]),
       verbose = FALSE)
  
  y <- fs_ids <- is_ids <- vector(mode = "list",
                                  length = length(gc01))
  
  for (m2 in seq_along(gc01)) {
    x <- cbind(gc01[[m2]],
               "is" = grepl(pattern = "internal stop",
                            x = gc01[[m2]]$Note),
               "fs" = grepl(pattern = "frameshifted",
                            x = gc01[[m2]]$Note))
    # drop non-coding annotations
    x <- x[x$Coding, ]
    # convert to a data.frame because 'by()' doesn't like DataFrames
    x <- data.frame("Product" = x$Product,
                    "count" = x$Coding,
                    "is" = x$is,
                    "fs" = x$fs)
    
    props <- by(data = x,
                INDICES = ~ Product,
                FUN = function(x) {
                  c(nrow(x),
                    sum(x$is),
                    sum(x$fs))
                })
    props <- do.call(rbind,
                     props)
    props <- data.frame("annotation" = rownames(props),
                        "count" = props[, 1],
                        "is" = props[, 2],
                        "fs" = props[, 3],
                        stringsAsFactors = FALSE,
                        row.names = NULL)
    y[[m2]] <- props
    
    fs_ids[[m2]] <- paste(gc01[[m2]]$ID[grepl(pattern = "frameshifted",
                                              x = gc01[[m2]]$Note)],
                          LETTERS[m2],
                          sep = "_")
    is_ids[[m2]] <- paste(gc01[[m2]]$ID[grepl(pattern = "internal stop",
                                              x = gc01[[m2]]$Note)],
                          LETTERS[m2],
                          sep = "_")
    
  }
  
  # i want every annotation in the test that either:
  # 1: does not appear at all in Cl1
  # or
  # 2: appears in Cl1 in a partnership without a pseudogene partner from the
  # pin
  
  # partners where the test is a pseudogene and the pin is not
  fs_targets_1 <- Cl1[which(!(Cl1$p1 %in% fs_ids[[1L]]) &
                              Cl1$p2 %in% fs_ids[[2L]]), ]
  if (nrow(fs_targets_1) > 0) {
    fs_targets_1 <- fs_targets_1$p2
  } else {
    fs_targets_1 <- NULL
  }
  
  # pseudogenes that appear in the test that have no partner in the pin
  fs_targets_2 <- fs_ids[[2L]][!(fs_ids[[2L]] %in% Cl1$p2)]
  
  if (length(fs_targets_1) > 0) {
    check_ids_1 <- gsub(x = fs_targets_1,
                        pattern = "_B",
                        replacement = "")
  } else {
    check_ids_1 <- NULL
  }
  
  if (length(fs_targets_2) > 0) {
    check_ids_2 <- gsub(x = fs_targets_2,
                        pattern = "_B",
                        replacement = "")
  } else {
    check_ids_2 <- NULL
  }
  
  if (length(check_ids_1) > 0) {
    fs_shifts_1 <- gc01[[2]]$Product[gc01[[2]]$ID %in% check_ids_1]
  } else {
    fs_shifts_1 <- NULL
  }
  if (length(check_ids_2) > 0) {
    fs_additions_1 <- gc01[[2]]$Product[gc01[[2]]$ID %in% check_ids_2]
  } else {
    fs_additions_1 <- NULL
  }
  
  # partners where the test is a pseudogene and the pin is not
  is_targets_1 <- Cl1[which(!(Cl1$p1 %in% is_ids[[1L]]) &
                            Cl1$p2 %in% is_ids[[2L]]), ]
  if (nrow(is_targets_1) > 0) {
    is_targets_1 <- is_targets_1$p2
  } else {
    is_targets_1 <- NULL
  }
  
  # pseudogenes that appear in the test that have no partner in the pin
  is_targets_2 <- is_ids[[2L]][!(is_ids[[2L]] %in% Cl1$p2)]
  
  if (length(is_targets_1) > 0) {
    check_ids_1 <- gsub(x = is_targets_1,
                        pattern = "_B",
                        replacement = "")
  } else {
    check_ids_1 <- NULL
  }
  
  if (length(is_targets_2) > 0) {
    check_ids_2 <- gsub(x = is_targets_2,
                        pattern = "_B",
                        replacement = "")
  } else {
    check_ids_2 <- NULL
  }
  
  if (length(check_ids_1) > 0) {
    is_shifts_1 <- gc01[[2]]$Product[gc01[[2]]$ID %in% check_ids_1]
  } else {
    is_shifts_1 <- NULL
  }
  if (length(check_ids_2) > 0) {
    is_additions_1 <- gc01[[2]]$Product[gc01[[2]]$ID %in% check_ids_2]
  } else {
    is_additions_1 <- NULL
  }
  
  stats_res01[[m1]] <- list(fs_shifts_1,
                            fs_additions_1,
                            is_shifts_1,
                            is_additions_1)
  stats_res02[[m1]] <- y[[1]]
  stats_res03[[m1]] <- y[[2]]
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)
tend <- Sys.time()
print(tend - tstart)

abs_changes <- vector(mode = "list",
                      length = length(stats_res01))
for (m1 in seq_along(abs_change_fs)) {
  abs_changes[[m1]] <- stats_res02[[m1]]
  if (length(stats_res01[[m1]][[1]]) > 0) {
    # shifts
    ph1 <- table(stats_res01[[m1]][[1]])
    ph2 <- names(ph1) %in% stats_res02[[m1]]$annotation
    ph1 <- ph1[ph2]
    if (length(ph1) > 0) {
      ph3 <- match(x = names(ph1),
                   table = stats_res02[[m1]]$annotation)
      abs_changes[[m1]][ph3, 4] <- abs_changes[[m1]][ph3, 4] + unname(ph1)
    }
  }
  if (length(stats_res01[[m1]][[2]]) > 0) {
    # additions
    ph1 <- table(stats_res01[[m1]][[2]])
    ph2 <- names(ph1) %in% stats_res02[[m1]]$annotation
    ph1 <- ph1[ph2]
    if (length(ph1) > 0) {
      ph3 <- match(x = names(ph1),
                   table = stats_res02[[m1]]$annotation)
      abs_changes[[m1]][ph3, c(2, 4)] <- abs_changes[[m1]][ph3, c(2, 4)] + unname(ph1)
    }
  }
  
  # abs_changes[[m1]] <- stats_res02[[m1]]
  if (length(stats_res01[[m1]][[3]]) > 0) {
    # shifts
    ph1 <- table(stats_res01[[m1]][[3]])
    ph2 <- names(ph1) %in% stats_res02[[m1]]$annotation
    ph1 <- ph1[ph2]
    if (length(ph1) > 0) {
      ph3 <- match(x = names(ph1),
                   table = stats_res02[[m1]]$annotation)
      abs_changes[[m1]][ph3, 3] <- abs_changes[[m1]][ph3, 3] + unname(ph1)
    }
  }
  if (length(stats_res01[[m1]][[4]]) > 0) {
    # additions
    ph1 <- table(stats_res01[[m1]][[4]])
    ph2 <- names(ph1) %in% stats_res02[[m1]]$annotation
    ph1 <- ph1[ph2]
    if (length(ph1) > 0) {
      ph3 <- match(x = names(ph1),
                   table = stats_res02[[m1]]$annotation)
      abs_changes[[m1]][ph3, c(2, 3)] <- abs_changes[[m1]][ph3, c(2, 3)] + unname(ph1)
    }
  }
}

# munge these lists down into a single averaged representative for each available
# scientific name

u1 <- table(ill_only_targets$Genus_Species)
pBar <- txtProgressBar(style = 1)
PBAR <- length(u1)
stats_res04 <- stats_res05 <- tot_changes <- set_a <- set_b <- set_c <-  vector(mode = "list",
                                                                                length = PBAR)

for (m1 in seq_along(u1)) {
  u2 <- which(ill_only_targets$Genus_Species == names(u1)[m1])
  if (length(u2) > 1) {
    u3 <- sample(x = length(u2),
                 size = 1,
                 replace = FALSE)
  }
  u3 <- 1
  
  x <- do.call(rbind,
               stats_res02[u2])
  
  z <- list(tapply(X = x,
                   INDEX = x$annotation,
                   FUN = function(x) {
                     sum(x$count)
                   }),
            tapply(X = x,
                   INDEX = x$annotation,
                   FUN = function(x) {
                     sum(x$is)
                   }),
            tapply(X = x,
                   INDEX = x$annotation,
                   FUN = function(x) {
                     sum(x$fs)
                   }))
  z <- do.call(cbind,
               z)
  z1 <- z / length(u2)
  stats_res04[[m1]] <- data.frame("annotation" = rownames(z),
                                  "count" = z[, 1],
                                  "is" = z[, 2],
                                  "fs" = z[, 3],
                                  "avg_count" = z1[, 1],
                                  "avg_is" = z1[, 2],
                                  "avg_fs" = z1[, 3],
                                  row.names = NULL)
  
  x <- do.call(rbind,
               stats_res03[u2])
  
  z <- list(tapply(X = x,
                   INDEX = x$annotation,
                   FUN = function(x) {
                     sum(x$count)
                   }),
            tapply(X = x,
                   INDEX = x$annotation,
                   FUN = function(x) {
                     sum(x$is)
                   }),
            tapply(X = x,
                   INDEX = x$annotation,
                   FUN = function(x) {
                     sum(x$fs)
                   }))
  z <- do.call(cbind,
               z)
  z1 <- z / length(u2)
  stats_res05[[m1]] <- data.frame("annotation" = rownames(z),
                                  "count" = z[, 1],
                                  "is" = z[, 2],
                                  "fs" = z[, 3],
                                  "avg_count" = z1[, 1],
                                  "avg_is" = z1[, 2],
                                  "avg_fs" = z1[, 3],
                                  row.names = NULL)
  
  x <- do.call(rbind,
               abs_changes[u2])
  
  z <- list(tapply(X = x,
                   INDEX = x$annotation,
                   FUN = function(x) {
                     sum(x$count)
                   }),
            tapply(X = x,
                   INDEX = x$annotation,
                   FUN = function(x) {
                     sum(x$is)
                   }),
            tapply(X = x,
                   INDEX = x$annotation,
                   FUN = function(x) {
                     sum(x$fs)
                   }))
  z <- do.call(cbind,
               z)
  z1 <- z / length(u2)
  tot_changes[[m1]] <- data.frame("annotation" = rownames(z),
                                  "count" = z[, 1],
                                  "is" = z[, 2],
                                  "fs" = z[, 3],
                                  "avg_count" = z1[, 1],
                                  "avg_is" = z1[, 2],
                                  "avg_fs" = z1[, 3],
                                  row.names = NULL)
  
  set_a[[m1]] <- do.call(rbind,
                         stats_res02[u2[u3]])
  set_b[[m1]] <- do.call(rbind,
                         stats_res03[u2[u3]])
  set_c[[m1]] <- do.call(rbind,
                         abs_changes[u2[u3]])
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)

stats_res06 <- do.call(rbind,
                       stats_res04)
stats_res06 <- tapply(X = stats_res06,
                      INDEX = stats_res06$annotation,
                      FUN = function(x) {
                        c(c(sum(x$avg_count),
                            sum(x$avg_is),
                            sum(x$avg_fs)) / length(stats_res04),
                          c(sum(x$count),
                            sum(x$is),
                            sum(x$fs)))
                      })
stats_res06 <- do.call(rbind,
                       stats_res06)

stats_res07 <- do.call(rbind,
                       stats_res05)
stats_res07 <- tapply(X = stats_res07,
                      INDEX = stats_res07$annotation,
                      FUN = function(x) {
                        c(c(sum(x$avg_count),
                            sum(x$avg_is),
                            sum(x$avg_fs)) / length(stats_res05),
                          c(sum(x$count),
                            sum(x$is),
                            sum(x$fs)))
                      })
stats_res07 <- do.call(rbind,
                       stats_res07)

tot_changes_a <- do.call(rbind,
                         tot_changes)
tot_changes_a <- tapply(X = tot_changes_a,
                        INDEX = tot_changes_a$annotation,
                        FUN = function(x) {
                          c(c(sum(x$avg_count),
                              sum(x$avg_is),
                              sum(x$avg_fs)) / length(tot_changes),
                            c(sum(x$count),
                              sum(x$is),
                              sum(x$fs)))
                        })
tot_changes_a <- do.call(rbind,
                         tot_changes_a)

set_a1 <- do.call(rbind,
                  set_a)
set_a1 <- tapply(X = set_a1,
                INDEX = set_a1$annotation,
                FUN = function(x) {
                  c(sum(x$count),
                    sum(x$is),
                    sum(x$fs))
                })
set_a1 <- do.call(rbind,
                  set_a1)
set_a1 <- data.frame("annotation" = rownames(set_a1),
                     "count" = set_a1[, 1],
                     "is" = set_a1[, 2],
                     "fs" = set_a1[, 3],
                     row.names = NULL)

set_b1 <- do.call(rbind,
                  set_b)
set_b1 <- tapply(X = set_b1,
                 INDEX = set_b1$annotation,
                 FUN = function(x) {
                   c(sum(x$count),
                     sum(x$is),
                     sum(x$fs))
                 })
set_b1 <- do.call(rbind,
                  set_b1)
set_b1 <- data.frame("annotation" = rownames(set_b1),
                     "count" = set_b1[, 1],
                     "is" = set_b1[, 2],
                     "fs" = set_b1[, 3],
                     row.names = NULL)

set_c1 <- do.call(rbind,
                  set_c)
set_c1 <- tapply(X = set_c1,
                 INDEX = set_c1$annotation,
                 FUN = function(x) {
                   c(sum(x$count),
                     sum(x$is),
                     sum(x$fs))
                 })
set_c1 <- do.call(rbind,
                  set_c1)
set_c1 <- data.frame("annotation" = rownames(set_c1),
                     "count" = set_c1[, 1],
                     "is" = set_c1[, 2],
                     "fs" = set_c1[, 3],
                     row.names = NULL)


u2 <- intersect(rownames(stats_res06),
                rownames(stats_res07))
mat1 <- match(table = rownames(stats_res06),
              x = u2)
mat2 <- match(table = rownames(stats_res07),
              x = u2)

stats_res08 <- data.frame("avg_count_1" = stats_res06[mat1, 1],
                          "avg_is_1" = stats_res06[mat1, 2],
                          "avg_fs_1" = stats_res06[mat1, 3],
                          "tot_count_1" = stats_res06[mat1, 4],
                          "tot_is_1" = stats_res06[mat1, 5],
                          "tot_fs_1" = stats_res06[mat1, 6],
                          "avg_count_2" = stats_res07[mat2, 1],
                          "avg_is_2" = stats_res07[mat2, 2],
                          "avg_fs_2" = stats_res07[mat2, 3],
                          "tot_count_2" = stats_res07[mat2, 4],
                          "tot_is_2" = stats_res07[mat2, 5],
                          "tot_fs_2" = stats_res07[mat2, 6],
                          row.names = u2)

colnames(stats_res06) <- colnames(stats_res07) <- colnames(tot_changes_a) <- c("norm_counts",
                                                                               "norm_is",
                                                                               "norm_fs",
                                                                               "tot_counts",
                                                                               "tot_is",
                                                                               "tot_fs")

# prob_a is the probability of pseudogenization by the weighted mean
# for internal stops FROM THE PIN ASSEMBLIES AT 1000-fold coverage
prob_a <- weighted.mean(x = stats_res06[, "norm_is"] / stats_res06[, "norm_counts"],
                        w = stats_res06[, "norm_counts"])
# prob_b is the probability of pseudogenization by weighted mean
# for frameshifts FROM THE PIN ASSEMBLIES AT 1000-fold coverage
prob_b <- weighted.mean(x = stats_res06[, "norm_fs"] / stats_res06[, "norm_counts"],
                        w = stats_res06[, "norm_counts"])

stats_res09 <- stats_res08[stats_res08$tot_is_2 > 0, ]
stats_res10 <- stats_res08[stats_res08$tot_fs_2 > 0, ]

save(stats_res08,
     stats_res09,
     stats_res10,
     prob_a,
     prob_b,
     file = paste0(targetdir01,
                   "/ResultF01.RData"),
     compress = "xz")

# expected events given the weights
bound_a <- qbinom(p = 0.001 / nrow(tot_changes_a),
                  size = seq(from = 1,
                             to = nrow(tot_changes_a),
                             by = 1),
                  prob = prob_a,
                  lower.tail = FALSE) / seq(from = 1,
                                            to = nrow(tot_changes_a),
                                            by = 1)

val01 <- set_a1$is
val01[val01 > 0] <- (val01[val01 > 0] - 1L)
val02 <- set_c1$is
val02[val02 > 0] <- (val02[val02 > 0] - 1L)
val03 <- set_a1$fs
val03[val03 > 0] <- (val03[val03 > 0] - 1L)
val04 <- set_c1$fs
val04[val04 > 0] <- (val04[val04 > 0] - 1L)

is_probs_1 <- p.adjust(mapply(SIMPLIFY = TRUE,
                              FUN = function(y, z) {
                                pbinom(q = y,
                                       size = z,
                                       prob = prob_a,
                                       lower.tail = FALSE)
                              },
                              y = val01, # set_a1 'successes' offset for the binomial
                              z = set_a1$count),
                       method = "bonferroni")
is_probs_2 <- p.adjust(mapply(SIMPLIFY = TRUE,
                              FUN = function(y, z) {
                                pbinom(q = y,
                                       size = z,
                                       prob = prob_a,
                                       lower.tail = FALSE)
                              },
                              y = val02, # set_c1 'successes' offset for the binomial
                              z = set_c1$count),
                       method = "bonferroni")

fs_probs_1 <- p.adjust(mapply(SIMPLIFY = TRUE,
                              FUN = function(y, z) {
                                pbinom(q = y,
                                       size = z,
                                       prob = prob_b,
                                       lower.tail = FALSE)
                              },
                              y = val03,
                              z = set_a1$count),
                       method = "bonferroni")
fs_probs_2 <- p.adjust(mapply(SIMPLIFY = TRUE,
                              FUN = function(y, z) {
                                pbinom(q = y,
                                       size = z,
                                       prob = prob_b,
                                       lower.tail = FALSE)
                              },
                              y = val04,
                              z = set_c1$count),
                       method = "bonferroni")

suppressWarnings(plot(is_probs_1,
                      is_probs_2,
                      pch = 46,
                      log = "xy",
                      xlim = c(1e-100, 1),
                      ylim = c(1e-100, 1)))
abline(h = 1e-3)
abline(v = 1e-3)
suppressWarnings(plot(fs_probs_1,
                      fs_probs_2,
                      pch = 46,
                      log = "xy",
                      xlim = c(1e-100, 1),
                      ylim = c(1e-100, 1)))
abline(h = 1e-3)
abline(v = 1e-3)

# build a data.frame of the counts, the frequency, the cex bin, and the annotation
# bins:
# 1e-3 to 1e-10
# 1e-10 to 1e-50
# 1e-50 to 1e-100
# < 1e-100
# all args are default but retained for my sanity
# opening the left is not necessary because that's our significance threshold

# all(set_c1$annotation == set_a1$annotation) # should be TRUE

# set bound
signif_bound <- 1e-3

# subset everything
index01 <- fs_probs_2 < signif_bound
fs_probs_3 <- fs_probs_1[index01]
fs_probs_4 <- fs_probs_2[index01]
fs_ann_for_plot <- set_c1$annotation[index01]
fs_plotvals01 <- set_c1$count[index01]
fs_plotvals02 <- set_c1$fs[index01] / set_c1$count[index01]
# axis break accommodation
fs_plotvals03 <- fs_plotvals01
fs_plotvals03[fs_plotvals03 > 2500] <- 3000
fs_cex <- findInterval(x = -log(fs_probs_4),
                       vec = -log(c(1e-3, 1e-10, 1e-50, 1e-100)),
                       rightmost.closed = FALSE,
                       all.inside = FALSE,
                       left.open = FALSE)

# subset everything
index02 <- is_probs_2 < signif_bound
is_probs_3 <- is_probs_1[index02]
is_probs_4 <- is_probs_2[index02]
is_ann_for_plot <- set_c1$annotation[index02]
is_plotvals01 <- set_c1$count[index02]
is_plotvals02 <- set_c1$is[index02] / set_c1$count[index02]
# axis break accommodation
is_plotvals03 <- is_plotvals01
is_plotvals03[is_plotvals03 > 2500] <- 3000
is_cex <- findInterval(x = -log(is_probs_4),
                       vec = -log(c(1e-3, 1e-10, 1e-50, 1e-100)),
                       rightmost.closed = FALSE,
                       all.inside = FALSE,
                       left.open = FALSE)

plot(x = fs_plotvals03,
     y = fs_plotvals02,
     cex = fs_cex,
     col = ifelse(test = fs_probs_3 >= signif_bound,
                  yes = ColVec4[1],
                  no = ColVec4[2]),
     pch = 16,
     xlim = c(0, 2500),
     ylim = c(0, .75),
     # xaxt = "n",
     yaxt = "n",
     yaxs = "i",
     xlab = "Count",
     ylab = "Pseudogenization Rate (%)")
# axis(side = 1,
#      at = c(0, 500, 1000, 1500, 2000, 3000),
#      labels = c("0", "500", "1000", "1500", "2000", ">40000"))
# axis.break(axis = 1,
#            breakpos = 2500)
axis(side = 2,
     at = c(0, .25, .5, .75),
     labels = c(0, 25, 50, 75))

plot(x = is_plotvals03,
     y = is_plotvals02,
     cex = fs_cex,
     col = ifelse(test = is_probs_3 >= signif_bound,
                  yes = ColVec4[1],
                  no = ColVec4[2]),
     pch = 16,
     xlim = c(0, 1000),
     ylim = c(0, .5),
     # xaxt = "n",
     yaxt = "n",
     yaxs = "i",
     xlab = "Count",
     ylab = "Pseudogenization Rate (%)")
axis(side = 2,
     at = c(0, .25, .5),
     labels = c(0, 25, 50))

pdf(file = "tempplot.pdf",
    height = 3.5,
    width = 7)
layout(mat = matrix(data = 1:2,
                    nrow = 1))
par(mar = c(3,3,2,1),
    mgp = c(1.75, .75, 0))
plot(x = fs_plotvals03,
     y = fs_plotvals02,
     cex = fs_cex,
     col = ifelse(test = fs_probs_3 >= signif_bound,
                  yes = ColVec4[1],
                  no = ColVec4[2]),
     pch = 16,
     xlim = c(0, 2500),
     ylim = c(0, .75),
     # xaxt = "n",
     yaxt = "n",
     yaxs = "i",
     xlab = "Count",
     ylab = "Pseudogenization Rate (%)",
     main = "Frameshifts")
# axis(side = 1,
#      at = c(0, 500, 1000, 1500, 2000, 3000),
#      labels = c("0", "500", "1000", "1500", "2000", ">40000"))
# axis.break(axis = 1,
#            breakpos = 2500)
axis(side = 2,
     at = c(0, .25, .5, .75),
     labels = c(0, 25, 50, 75))

par(mar = c(3,2,2,2),
    mgp = c(1.75, .75, 0))
plot(x = is_plotvals03,
     y = is_plotvals02,
     cex = fs_cex,
     col = ifelse(test = is_probs_3 >= signif_bound,
                  yes = ColVec4[1],
                  no = ColVec4[2]),
     pch = 16,
     xlim = c(0, 1000),
     ylim = c(0, .5),
     # xaxt = "n",
     yaxt = "n",
     yaxs = "i",
     xlab = "Count",
     # ylab = "Pseudogenization Rate (%)",
     ylab = "",
     main = "Internal stops")
axis(side = 2,
     at = c(0, .25, .5),
     labels = c(0, 25, 50))
L <- legend(x = 375,
            y = .475,
            legend = rep(NA, 2 * 4),
            # legend = u_meta_1,
            col = rep(ColVec4[1:2], each = 4),
            pch = rep(16, 8),
            pt.cex = rep(1:4, 2),
            ncol = 2,
            bty = 'n',
            x.intersp = 2,
            y.intersp = 1.5,
            inset = 0.02,
            cex = 0.75)
legend(x = L$rect$left + L$rect$w - 120,
       y = L$rect$top,
       # legend = u_meta_1, # can drop 'Illumina' from names -- it's redundant
       legend = c(">= 1e-10",
                  "1e-10 to 1e-50",
                  "1e-50 to 1e-100",
                  "< 1e-100"),
       ncol = 1,
       x.intersp = 1.5,
       y.intersp = 1.5,
       bg = NA,
       bty = "n",
       cex = 0.75)
# x-pos of first column = 627.5
# x-pos of second column = 712.5
text(x = c(420, 555),
     y = rep(.465, 2),
     labels = c("*",
                "**"),
     cex = 1.5,
     xpd = TRUE)
# abline(v = 627.5)
# abline(v = 712.5)
# abline(v = 750)
# abline(v = 425)
# abline(v = 545)
# abline(v = 560)
# abline(h = .465)
# abline(h = .455)
dev.off()

is_plot_df <- data.frame("signif_initial" = is_probs_3,
                         "signif_post" = is_probs_4,
                         "post_count" = is_plotvals01,
                         "post_rate" = is_plotvals02,
                         "annotation" = is_ann_for_plot)

fs_plot_df <- data.frame("signif_initial" = fs_probs_3,
                         "signif_post" = fs_probs_4,
                         "post_count" = fs_plotvals01,
                         "post_rate" = fs_plotvals02,
                         "annotation" = fs_ann_for_plot)

save(fs_plot_df,
     is_plot_df,
     set_a1,
     set_b1,
     set_c1,
     file = paste0(targetdir01,
                   "/ResultG01.RData"),
     compress = "xz")


nvals01 <- weighted.mean(x = (set_c1$is - set_a1$is) / set_a1$count,
                         w = set_a1$count)
nvals02 <- (set_c1$is - set_a1$is)
nvals03 <- nvals02
nvals03[nvals03 > 0] <- nvals03[nvals03 > 0] - 1L
nvals03 <- p.adjust(mapply(SIMPLIFY = TRUE,
                           FUN = function(y, z) {
                             pbinom(q = y,
                                    size = z,
                                    prob = nvals01,
                                    lower.tail = FALSE)
                           },
                           y = nvals03,
                           z = set_a1$count),
                    method = "bonferroni")

nvals04 <- weighted.mean(x = (set_c1$fs - set_a1$fs) / set_a1$count,
                         w = set_a1$count)
nvals05 <- (set_c1$fs - set_a1$fs)
nvals06 <- nvals05
nvals06[nvals06 > 0] <- nvals06[nvals06 > 0] - 1L
nvals06 <- p.adjust(mapply(SIMPLIFY = TRUE,
                           FUN = function(y, z) {
                             pbinom(q = y,
                                    size = z,
                                    prob = nvals04,
                                    lower.tail = FALSE)
                           },
                           y = nvals06,
                           z = set_a1$count),
                    method = "bonferroni")

nvals07 <- nvals03 < signif_bound
nvals08 <- nvals06 < signif_bound



pdf(file = "tempplot.pdf",
    height = 3.5,
    width = 7)
layout(mat = matrix(data = 1:2,
                    nrow = 1))
plot(x = set_a1$count[nvals07],
     y = ((set_c1$is - set_a1$is) / set_a1$count)[nvals07],
     pch = 16,
     xlim = c(0, 1000),
     ylim = c(0, 0.025),
     xlab = "count",
     ylab = "rate of change",
     main = "internal stops")
text(x = set_a1$count[nvals07] - 20,
     y = ((set_c1$is - set_a1$is) / set_a1$count)[nvals07] - 0.0025,
     labels = set_c1$annotation[nvals07],
     xpd = TRUE,
     cex = 0.5)
plot(x = set_a1$count[nvals08],
     y = ((set_c1$fs - set_a1$fs) / set_a1$count)[nvals08],
     pch = 16,
     xlim = c(0, 1000),
     ylim = c(0, 0.07),
     xlab = "count",
     ylab = "rate of change",
     main = "frameshifts")
text(x = set_a1$count[nvals08] + 5,
     y = ((set_c1$fs - set_a1$fs) / set_a1$count)[nvals08] - 0.005,
     labels = set_c1$annotation[nvals08],
     xpd = TRUE,
     cex = 0.5)
dev.off()

vals01 <- is_probs_1
vals01[vals01 == 0] <- min(vals01[vals01 > 0])
# vals01 <- (log(vals01 / min(vals01))) / 50
vals01 <- -log(vals01) / max(-log(vals01))
vals02 <- is_probs_2
vals02[vals02 == 0] <- min(vals02[vals02 > 0])
# vals02 <- (log(vals02 / min(vals02))) / 50
vals02 <- -log(vals02) / max(-log(vals02))
vals03 <- fs_probs_1
vals03[vals03 == 0] <- min(vals03[vals03 > 0])
# vals03 <- (log(vals03 / min(vals03))) / 50
vals03 <- (-log(vals03) / max(-log(vals03)))
vals04 <- fs_probs_2
vals04[vals04 == 0] <- min(vals04[vals04 > 0])
# vals04 <- (log(vals04 / min(vals04))) / 50
vals04 <- (-log(vals04) / max(-log(vals04)))

plot(x = 0,
     y = 0,
     type = "n",
     xlim = c(0, 2500),
     ylim = c(0, 1),
     xlab = "Counts",
     ylab = "pseudogenizations (%)")
points(x = set_c1$count,
       y = set_c1$is / set_c1$count,
       pch = ifelse(test = is_probs_1 >= is_probs_2,
                    yes = 16,
                    no = 15),
       cex = vals02)
points(bound_a, pch = 46)


plot(x = 0,
     y = 0,
     type = "n",
     xlim = c(0, 2500),
     ylim = c(0, 1),
     xlab = "Counts",
     ylab = "pseudogenizations (%)")
points(x = set_c1$count,
       y = set_c1$fs / set_c1$count,
       pch = ifelse(test = fs_probs_1 >= fs_probs_2,
                    yes = 16,
                    no = 15),
       cex = vals04)

plot(x = 0,
     y = 0,
     type = "n",
     xlim = c(0, 1000),
     ylim = c(0, 1),
     xlab = "total occurances",
     ylab = "spu")


t1 <- set_a1$annotation[is_probs_1 > 1e-3 & is_probs_2 <= 1e-3]
t2 <- set_a1$annotation[fs_probs_1 > 1e-3 & fs_probs_2 <= 1e-3]
t3 <- sort(unique(c(t1, t2)))

stats_res06[rownames(stats_res06) %in% t3, ]
range01 <- c(0, max(c(tot_changes_a[rownames(tot_changes_a) %in% t3, 2],
                      tot_changes_a[rownames(tot_changes_a) %in% t3, 3])))

plot(x = 0,
     y = 0,
     type = "n",
     ylim = range01,
     xlim = c(1, length(t3) + 1L),
     xaxt = "n",
     xaxs = "i",
     yaxs = "i",
     ylab = "Frequency",
     xlab = "")
vals01 <- tot_changes_a[rownames(tot_changes_a) %in% t3, 2]
vals02 <- tot_changes_a[rownames(tot_changes_a) %in% t3, 3]
for (m1 in seq_along(t3)) {
  rect(xleft = m1,
       xright = m1 + 0.45,
       ybottom = 0,
       ytop = vals01[m1],
       col = ColVec1[1])
  rect(xleft = m1 + 0.46,
       xright = m1 + 0.9,
       ybottom = 0,
       ytop = vals02[m1],
       col = ColVec1[2])
}

