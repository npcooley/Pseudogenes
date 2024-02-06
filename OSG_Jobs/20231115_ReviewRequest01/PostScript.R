###### -- build an initial table of the parsed annotation comparisons ---------
# each comparison needs to fill in the table rows of:
# # total contigs
# # total nucleotides
# # source contigs
# # source nucleotides
# actual evaluated coverage of the assembly (weighted mean)
# # short read depth
# # long read depth if applicable
# # normalized n50 i.e (n50 / total nucleotides)
# # mean ani
# # mean af
# # source is count
# # source fs count
# # matched is count
# # matched fs count
# # new is count
# # new fs count
# # missing is count
# # missing fs count

# first assembly is the 'pin' assembly, while the second is the 'test' assembly

suppressMessages(library(SynExtend))

ColVec1 <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000')
ColVec2 <- paste0(ColVec1,
                  "33")
ColVec3 <- paste0(ColVec1,
                  "15")

targetdir01 <- "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01"
targetdir02 <- "~/Repos/PseudogenePlots/InputData"
files01 <- list.files(path = targetdir01,
                      pattern = "ResultD")
completedjobs <- as.integer(unlist(regmatches(x = files01,
                                              m = gregexpr(pattern = "[0-9]+",
                                                           text = files01))))
initialmap <- read.table(file = paste0(targetdir01,
                                       "/JobMapD.txt"))
expectedjobs <- initialmap$V1
expectedpinAssembly <- as.integer(unlist(regmatches(x = initialmap$V2,
                                                    m = gregexpr(pattern = "[0-9]+",
                                                                 text = initialmap$V2))))
expectedtestAssembly <- as.integer(unlist(regmatches(x = initialmap$V4,
                                                     m = gregexpr(pattern = "[0-9]+",
                                                                  text = initialmap$V4))))
pBar <- txtProgressBar(style = 1)
PBAR <- length(files01)
initialpin <- 0
res <- vector(mode = "list",
              length = PBAR)

for (m1 in seq_along(files01)) {
  load(file = paste0(targetdir01,
                     "/",
                     files01[m1]),
       verbose = FALSE)
  # load(file = paste0(targetdir01,
  #                    "/ResultC0000702.RData"),
  #      verbose = TRUE)
  w1 <- completedjobs[m1]
  w2 <- expectedjobs[expectedjobs == w1]
  # pinsource <- expectedpinAssembly[w2]
  testsource <- expectedtestAssembly[w2]
  # load in the data for the test assembly for the coverage
  load(file = paste0(targetdir01,
                     "/ResultC",
                     formatC(x = testsource,
                             width = 7,
                             flag = 0),
                     ".RData"),
       verbose = FALSE)
  
  
  res01 <- lengths(diag(syn))
  res02 <- sapply(X = diag(syn),
                  FUN = function(x) {
                    sum(x)
                  })
  z1 <- rle(dep3$V1)
  res03 <- weighted.mean(x = dep3$V3,
                        w = rep(x = (z1$lengths / nrow(dep3)),
                                times = z1$lengths))
  if (!is.null(dep4)) {
    z1 <- rle(dep4$V1)
    res04 <- weighted.mean(x = dep4$V3,
                           w = rep(x = (z1$lengths / nrow(dep4)),
                                   times = z1$lengths))
  } else {
    res04 <- -1
  }
  
  res05 <- unlist(N50List) / res02
  res06 <- strsplit(x = ani_res,
                    split = "\t",
                    fixed = TRUE)[[2]][3:6]
  res06 <- c(mean(as.numeric(res06[1:2])),
             mean(as.numeric(res06[3:4])))
  
  # source frameshifts and internal stops
  res07 <- sum(grepl(pattern = "frameshifted",
                     x = gc01[[1]]$Note))
  res08 <- sum(grepl(pattern = "internal stop",
                     x = gc01[[1]]$Note))
  # matched is # found reference
  res09 <- length(is3)
  # unmatched test # newly generated!
  res10 <- length(is2)
  # unmatched source # non-pseudo in test
  res11 <- length(is1)
  
  # matched fs # found reference
  res12 <- length(fs3)
  # unmatched test # newly generated!
  res13 <- length(fs2)
  # unmatched source # non-pseudo in test
  res14 <- length(fs1)
  
  # test frameshifts and internal stops
  res15 <- sum(grepl(pattern = "frameshifted",
                     x = gc01[[2]]$Note))
  res16 <- sum(grepl(pattern = "internal stop",
                     x = gc01[[2]]$Note))
  # test total genes
  res17 <- sum(gc01[[1]]$Coding)
  res18 <- sum(gc01[[2]]$Coding)
  
  res[[m1]] <- c(res01,
                 res02,
                 res03,
                 res04,
                 res05,
                 res06,
                 res07,
                 res08,
                 res09,
                 res10,
                 res11,
                 res12,
                 res13,
                 res14,
                 res15,
                 res16,
                 res17,
                 res18)
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}

res2 <- do.call(rbind,
                res)

dat <- data.frame("comp_id" = completedjobs,
                  "pin_contigs" = as.integer(res2[, 1]),
                  "test_contigs" = as.integer(res2[, 2]),
                  "pin_tn" = as.integer(res2[, 3]),
                  "test_tn" = as.integer(res2[, 4]),
                  "sr_cov" = res2[, 5],
                  "lr_cov" = res2[, 6],
                  "pin_norm_n50" = res2[, 7],
                  "test_norm_n50" = res2[, 8],
                  "ani" = res2[, 9],
                  "af" = res2[, 10],
                  "pin_fs_count" = as.integer(res2[, 11]),
                  "pin_is_count" = as.integer(res2[, 12]),
                  "found_is" = as.integer(res2[, 13]),
                  "new_is" = as.integer(res2[, 14]),
                  "revert_is" = as.integer(res2[, 15]),
                  "found_fs" = as.integer(res2[, 16]),
                  "new_fs" = as.integer(res2[, 17]),
                  "revert_fs" = as.integer(res2[, 18]),
                  "test_fs_count" = as.integer(res2[, 19]),
                  "test_is_count" = as.integer(res2[, 20]),
                  "source_coding_count" = as.integer(res2[, 21]),
                  "test_coding_count" = as.integer(res2[, 22]))


save(dat,
     file = paste0(targetdir01,
                   "/ResultE01.RData"),
     compress = "xz")

load(file = paste0(targetdir01,
                   "/JobMap_A_B.RData"),
     verbose = TRUE)
load(file = paste0(targetdir02,
                   "/GenBank_SRA_Summary.RData"),
     verbose = TRUE)

map01 <- expectedpinAssembly[expectedjobs %in% completedjobs]
map02 <- expectedtestAssembly[expectedjobs %in% completedjobs]
map03 <- as.integer(rownames(res04))
# row position for sra run, sequencing platforms, and models
map04 <- match(x = map01,
               table = map03)
w1 <- match(x = res04$sr_r,
            table = SRA_Meta$Run)
map05 <- SRA_Meta$ScientificName[w1]

g1 <- grepl(pattern = "^candidatus",
            ignore.case = TRUE,
            x = map05)
if (any(g1)) {
  map05[g1] <- gsub(pattern = "^candidatus ",
                    x = map05[g1],
                    replacement = "",
                    ignore.case = TRUE)
}
map05 <- gsub(pattern = "'",
              replacement = "",
              x = map05)
map05 <- gsub(pattern = '"',
              replacement = "",
              x = map05)
map05 <- unlist(regmatches(m = gregexpr(pattern = "^[^ ]+ [^ ]+",
                                        text = map05),
                           x = map05))

map06 <- as.integer(rownames(res03))
map07 <- match(x = map02,
               table = map06)

# append the genus + species identifier, the short and long read models and platforms, 
apdat <- data.frame("Genus_Species" = map05[map04],
                    "SR_Run" = res04$sr_r[map04],
                    "LR_Run" = res04$lr_r[map04],
                    "SR_Plat" = res04$sr_p[map04],
                    "LR_Plat" = res04$lr_p[map04],
                    "SR_Mod" = res04$sr_mod[map04],
                    "LR_Mod" = res04$lr_mod[map04],
                    "Pin_LR_Cov_Target" = res04$lr_cov_target[map04],
                    "Pin_SR_Cov_Target" = res04$sr_cov_target[map04],
                    "Test_SR_Cov_Target" = res03$sr_cov_target[map07])

dat1 <- cbind(dat,
              apdat)

dat2 <- dat1[dat1$test_tn >= (dat1$pin_tn * 0.8), ]
w1 <- paste(dat2$SR_Run,
            dat2$Pin_LR_Cov_Target,
            sep = " + ")
w2 <- table(w1)
w2 <- names(w2[w2 == 21])
dat3 <- dat2[w1 %in% w2, ]

save(dat,
     dat1,
     dat2,
     dat3,
     file = paste0(targetdir01,
                   "/ResultE02.RData"),
     compress = "xz")

dat4 <- cbind(dat2,
              "index" = w1,
              "cov_log" = log10(dat2$sr_cov / 50))
res_by_fs <- by(data = dat4,
                INDICES = ~ index,
                FUN = function(x) {
                  summary(glm(formula = (abs(test_fs_count - found_fs) / test_coding_count) ~ cov_log,
                              data = x,
                              family = binomial(link = "logit"),
                              weights = test_coding_count))$coefficients
                })
res_by_is <- by(data = dat4,
                INDICES = ~ index,
                FUN = function(x) {
                  summary(glm(formula = (abs(test_is_count - found_is) / test_coding_count) ~ cov_log,
                              data = x,
                              family = binomial(link = "logit"),
                              weights = test_coding_count))$coefficients
                })


# given coefficients, predict |delta pseudogenes| per Mbp
predictResponse <- function(coef,
                            COV,
                            TOT_GENE,
                            TOT_NT) {
  Y <- coef[1] + coef[2]*COV
  Y <- exp(Y)
  Y <- Y/(1 + Y) # |delta pseudogenes| / (total genes)
  Y*TOT_GENE/TOT_NT*1e6 # |delta pseudogenes| / Mbp
}
convertCoverage <- function(cov) {
  50*(10^cov)
}
cov_set <- seq(from = -1,
               to = 1,
               by = 0.01)

# grab the responses
# build out a metadata table
resp_fs <- resp_is <- model_meta <- vector(mode = "list",
                                           length = length(res_by_fs))
for (m1 in seq_along(res_by_fs)) {
  w1 <- which(dat4$index == names(res_by_fs)[m1])
  w2 <- w1[1]
  w3 <- dat4$source_coding_count[w2]
  w4 <- dat4$pin_tn[w2]
  model_meta[[m1]] <- data.frame("Name" = dat4$Genus_Species[w2],
                                 "SR_Model" = dat4$SR_Mod[w2],
                                 "LR_Platform" = dat4$LR_Plat[w2],
                                 "LR_Model" = dat4$LR_Mod[w2],
                                 "Run" = dat4$SR_Run[w2],
                                 "LR_Cov" = dat4$Pin_LR_Cov_Target[w2])
  resp_fs[[m1]] <- predictResponse(coef = res_by_fs[[m1]][, 1L],
                                   COV = cov_set,
                                   TOT_GENE = w3,
                                   TOT_NT = w4)
  resp_is[[m1]] <- predictResponse(coef = res_by_is[[m1]][, 1L],
                                   COV = cov_set,
                                   TOT_GENE = w3,
                                   TOT_NT = w4)
}
model_meta <- do.call(rbind,
                      model_meta)
model_meta <- data.frame("Name" = model_meta[, 1],
                         "SR_Model" = model_meta[, 2],
                         "LR_Platform" = model_meta[, 3],
                         "LR_Model" = model_meta[, 4],
                         "Run" = model_meta[, 5],
                         "LR_Cov" = model_meta[, 6])

save(dat4,
     res_by_fs,
     res_by_is,
     resp_fs,
     resp_is,
     cov_set,
     file = paste0(targetdir01,
                   "/ResultE03.RData"),
     compress = "xz")


# initial plots -- how does the data look in aggregate?
plot(x = 0,
     y = 0,
     xlim = c(5, 500),
     ylim = c(0, 20),
     log = "x",
     xaxs = "i",
     yaxs = "i",
     xlab = "Coverage",
     ylab = "Delta FS / Mbp")
for (m1 in seq_along(resp_fs)) {
  lines(y = resp_fs[[m1]],
        x = convertCoverage(cov = cov_set),
        col = "#11221133")
}

plot(x = 0,
     y = 0,
     xlim = c(5, 500),
     ylim = c(0, 20),
     log = "x",
     xaxs = "i",
     yaxs = "i",
     xlab = "Coverage",
     ylab = "Delta IS / Mbp")
for (m1 in seq_along(resp_is)) {
  lines(y = resp_is[[m1]],
        x = convertCoverage(cov = cov_set),
        col = "#11221133")
}

# how do we want to plot this? what parts of the dataset are intact?
# create a metadata table for the resp objects

# subset attempt 1, look at only LR coverage that was targetted to be close
# to our simulated coverages -- though, pacbio coverage stats don't necessarily match
# the target as well as the ont coverages
s1 <- head(sort(table(model_meta$Name[model_meta$LR_Cov != "50" &
                                        model_meta$LR_Cov != "80"]),
                decreasing = TRUE))
s2 <- names(s1)

u_meta_1 <- unique(model_meta$SR_Model)
u_meta_2 <- unique(model_meta$LR_Model)

pdf(file = paste0(targetdir01,
                  "/InitialPlotting001.pdf"),
    height = 7,
    width = 7)

layout(mat = matrix(data = 1:6,
                    nrow = 3,
                    byrow = TRUE))
par(mar = c(3,3,2,2),
    mgp = c(1.75, .75, 0))
for (m1 in seq_along(s2)) {
  w1 <- which(model_meta$Name == s2[m1] &
                model_meta$LR_Cov != "50" &
                model_meta$LR_Cov != "80")
  val1 <- match(x = model_meta$SR_Model[w1],
                table = u_meta_1)
  val2 <- match(x = model_meta$LR_Model[w1],
                table = u_meta_2)
  
  suppressWarnings(plot(x = 0,
                        y = 0,
                        type = "n",
                        xlim = c(5, 500),
                        ylim = c(0, 10),
                        log = "x",
                        xaxs = "i",
                        yaxs = "i",
                        xlab = "Coverage",
                        ylab = expression("|"*Delta*" Internal stops per Mbp|"),
                        main = s2[m1]))
  for (m2 in seq_along(w1)) {
    lines(y = resp_is[[w1[m2]]],
          x = convertCoverage(cov = cov_set),
          col = ColVec1[val1[m2]],
          lty = val2[m2])
  }
  if (m1 == 1) {
    legend("topright",
           legend = u_meta_1,
           col = ColVec1[seq(length(u_meta_1))],
           lty = 1,
           # cex = 0.5, # plotting to the regular graphics device
           cex = 0.75,
           bty = "n",
           bg = NA)
  } else if (m1 == 2) {
    legend("topright",
           legend = u_meta_2,
           col = "black",
           lty = seq(length(u_meta_2)),
           # cex = 0.5, # plotting to the regular graphics device
           cex = 0.75,
           bty = "n",
           bg = NA)
  }
  
}
dev.off()

pdf(file = paste0(targetdir01,
                  "/InitialPlotting002.pdf"),
    height = 7,
    width = 7)

layout(mat = matrix(data = 1:6,
                    nrow = 3,
                    byrow = TRUE))
par(mar = c(3,3,2,2),
    mgp = c(1.75, .75, 0))
for (m1 in seq_along(s2)) {
  w1 <- which(model_meta$Name == s2[m1] &
                model_meta$LR_Cov != "50" &
                model_meta$LR_Cov != "80")
  val1 <- match(x = model_meta$SR_Model[w1],
                table = u_meta_1)
  val2 <- match(x = model_meta$LR_Model[w1],
                table = u_meta_2)
  
  suppressWarnings(plot(x = 0,
                        y = 0,
                        type = "n",
                        xlim = c(5, 500),
                        ylim = c(0, 10),
                        log = "x",
                        xaxs = "i",
                        yaxs = "i",
                        xlab = "Coverage",
                        ylab = expression("|"*Delta*" Frameshifts per Mbp|"),
                        main = s2[m1]))
  for (m2 in seq_along(w1)) {
    lines(y = resp_fs[[w1[m2]]],
          x = convertCoverage(cov = cov_set),
          col = ColVec1[val1[m2]],
          lty = val2[m2])
  }
  if (m1 == 1) {
    legend("topright",
           legend = u_meta_1,
           col = ColVec1[seq(length(u_meta_1))],
           lty = 1,
           # cex = 0.5, # plotting to the regular graphics device
           cex = 0.75,
           bty = "n",
           bg = NA)
  } else if (m1 == 2) {
    legend("topright",
           legend = u_meta_2,
           col = "black",
           lty = seq(length(u_meta_2)),
           # cex = 0.5, # plotting to the regular graphics device
           cex = 0.75,
           bty = "n",
           bg = NA)
  }
  
}
dev.off()

# everyone scientific name that has at least 2 biosamples remaining in the 
# experimental data
# s3 <- table(model_meta$Name[model_meta$LR_Cov != "50" &
#                               model_meta$LR_Cov != "80"])
# s4 <- names(s3)[s3 > 1]
# 
# for (m1 in seq_along(s4)) {
#   pdf(file = paste0("tmp_plots/InitialPlots",
#                     formatC(x = m1,
#                             width = 3,
#                             flag = 0,
#                             format = "d"),
#                     ".pdf"),
#       height = 3.5,
#       width = 7)
#   
#   layout(mat = matrix(data = 1:2,
#                       nrow = 1))
#   par(mar = c(3,3,2,2),
#       mgp = c(1.75, .75, 0))
#   
#   w1 <- which(model_meta$Name == s4[m1] &
#                 model_meta$LR_Cov != "50" &
#                 model_meta$LR_Cov != "80")
#   val1 <- match(x = model_meta$SR_Model[w1],
#                 table = u_meta_1)
#   val2 <- match(x = model_meta$LR_Model[w1],
#                 table = u_meta_2)
#   
#   suppressWarnings(plot(x = 0,
#                         y = 0,
#                         type = "n",
#                         xlim = c(5, 500),
#                         ylim = c(0, 10),
#                         log = "x",
#                         xaxs = "i",
#                         yaxs = "i",
#                         xlab = "Coverage",
#                         ylab = expression("|"*Delta*" Internal stops per Mbp|"),
#                         main = s4[m1]))
#   for (m2 in seq_along(w1)) {
#     lines(y = resp_is[[w1[m2]]],
#           x = convertCoverage(cov = cov_set),
#           col = ColVec1[val1[m2]],
#           lty = val2[m2])
#   }
#   legend("topright",
#          legend = u_meta_1,
#          col = ColVec1[seq(length(u_meta_1))],
#          lty = 1,
#          # cex = 0.5, # plotting to the regular graphics device
#          cex = 0.75,
#          bty = "n",
#          bg = NA)
#   
#   suppressWarnings(plot(x = 0,
#                         y = 0,
#                         type = "n",
#                         xlim = c(5, 500),
#                         ylim = c(0, 10),
#                         log = "x",
#                         xaxs = "i",
#                         yaxs = "i",
#                         xlab = "Coverage",
#                         ylab = expression("|"*Delta*" Frameshifts per Mbp|"),
#                         main = s4[m1]))
#   for (m2 in seq_along(w1)) {
#     lines(y = resp_fs[[w1[m2]]],
#           x = convertCoverage(cov = cov_set),
#           col = ColVec1[val1[m2]],
#           lty = val2[m2])
#   }
#   legend("topright",
#          legend = u_meta_2,
#          col = "black",
#          lty = seq(length(u_meta_2)),
#          # cex = 0.5, # plotting to the regular graphics device
#          cex = 0.75,
#          bty = "n",
#          bg = NA)
#   dev.off()
# }




