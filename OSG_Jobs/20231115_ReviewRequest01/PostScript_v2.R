###### -- build an initial table of the parsed annotation comparisons v2-------
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
targetdir03 <- "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/Comparisons"
files01 <- list.files(path = targetdir03,
                      pattern = "ResultD")
completedjobs <- as.integer(unlist(regmatches(x = files01,
                                              m = gregexpr(pattern = "[0-9]+",
                                                           text = files01))))
initialmap <- read.table(file = paste0(targetdir01,
                                       "/JobMapD_v2.txt"))
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
  x <- try(load(file = paste0(targetdir03,
                              "/",
                              files01[m1]),
                verbose = FALSE))
  if (is(object = x,
         class2 = "try-error")) {
    next
  }
  
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

# subset to data
dat2 <- dat1[dat1$test_tn >= (dat1$pin_tn * 0.8) &
               dat1$test_tn <= (dat1$pin_tn * 1.2) &
               dat1$sr_cov >= (dat1$Test_SR_Cov_Target * 0.8) &
               dat1$sr_cov <= (dat1$Test_SR_Cov_Target * 1.2) &
               dat1$af >= 0.8, ]
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
                                 "LR_Cov" = dat4$Pin_LR_Cov_Target[w2],
                                 "SR_Cov_005_Reps" = sum(dat4$Test_SR_Cov_Target[w1] == 5),
                                 "SR_Cov_010_Reps" = sum(dat4$Test_SR_Cov_Target[w1] == 10),
                                 "SR_Cov_025_Reps" = sum(dat4$Test_SR_Cov_Target[w1] == 25),
                                 "SR_Cov_050_Reps" = sum(dat4$Test_SR_Cov_Target[w1] == 50),
                                 "SR_Cov_100_Reps" = sum(dat4$Test_SR_Cov_Target[w1] == 100),
                                 "SR_Cov_250_Reps" = sum(dat4$Test_SR_Cov_Target[w1] == 250),
                                 "SR_Cov_500_Reps" = sum(dat4$Test_SR_Cov_Target[w1] == 500))
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

save(dat4,
     model_meta,
     res_by_fs,
     res_by_is,
     resp_fs,
     resp_is,
     cov_set,
     file = paste0(targetdir01,
                   "/ResultE03.RData"),
     compress = "xz")

# subsetting attempt 1, limit to only comparisons where at least one replicate
# was returned for every short read target except 5x, if we didn't assemble (within the previously specified parameter)
# at the lowest target coverage, we just have to live with it

s1 <- head(sort(table(model_meta$Name[model_meta$LR_Cov != "50" &
                                        model_meta$LR_Cov != "80" &
                                        model_meta$SR_Cov_010_Reps >= 1 &
                                        model_meta$SR_Cov_025_Reps >= 1 &
                                        model_meta$SR_Cov_050_Reps >= 1 &
                                        model_meta$SR_Cov_100_Reps >= 1 &
                                        model_meta$SR_Cov_250_Reps >= 1 &
                                        model_meta$SR_Cov_500_Reps >= 1]),
                decreasing = TRUE))
s2 <- names(s1)

u_meta_1 <- unique(model_meta$SR_Model)
u_meta_2 <- unique(model_meta$LR_Model)

pdf(file = paste0(targetdir01,
                  "/InitialPlotting003.pdf"),
    height = 7,
    width = 7)

layout(mat = matrix(data = 1:6,
                    nrow = 3,
                    byrow = TRUE))
par(mar = c(3,3,2,2),
    mgp = c(1.75, .75, 0))
for (m1 in seq_along(s2)) {
  
  # subsetting attempt 1, limit to only comparisons where at least one replicate
  # was returned for every short read target except 5x, if we didn't assemble (within the previously specified parameter)
  # at the lowest target coverage, we just have to live with it
  
  w1 <- which(model_meta$Name == s2[m1] &
                model_meta$LR_Cov != "50" &
                model_meta$LR_Cov != "80" &
                model_meta$SR_Cov_010_Reps >= 1 &
                model_meta$SR_Cov_025_Reps >= 1 &
                model_meta$SR_Cov_050_Reps >= 1 &
                model_meta$SR_Cov_100_Reps >= 1 &
                model_meta$SR_Cov_250_Reps >= 1 &
                model_meta$SR_Cov_500_Reps >= 1)
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
                  "/InitialPlotting004.pdf"),
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
                model_meta$LR_Cov != "80" &
                model_meta$SR_Cov_010_Reps >= 1 &
                model_meta$SR_Cov_025_Reps >= 1 &
                model_meta$SR_Cov_050_Reps >= 1 &
                model_meta$SR_Cov_100_Reps >= 1 &
                model_meta$SR_Cov_250_Reps >= 1 &
                model_meta$SR_Cov_500_Reps >= 1)
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
s3 <- table(model_meta$Name[model_meta$LR_Cov != "50" &
                              model_meta$LR_Cov != "80" &
                              model_meta$SR_Cov_010_Reps >= 1 &
                              model_meta$SR_Cov_025_Reps >= 1 &
                              model_meta$SR_Cov_050_Reps >= 1 &
                              model_meta$SR_Cov_100_Reps >= 1 &
                              model_meta$SR_Cov_250_Reps >= 1 &
                              model_meta$SR_Cov_500_Reps >= 1])
s4 <- names(s3)[s3 > 1]

# # plot out every genus / species name pair that has more than two identifiers
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
#                 model_meta$LR_Cov != "80" &
#                 model_meta$SR_Cov_010_Reps >= 1 &
#                 model_meta$SR_Cov_025_Reps >= 1 &
#                 model_meta$SR_Cov_050_Reps >= 1 &
#                 model_meta$SR_Cov_100_Reps >= 1 &
#                 model_meta$SR_Cov_250_Reps >= 1 &
#                 model_meta$SR_Cov_500_Reps >= 1)
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

# visualise reassemblies based on whether or not the coverage outcome was modestly
# close to the target
plot(dat4$sr_cov - dat4$Test_SR_Cov_Target,
     pch = 16,
     cex = 0.5,
     col = ifelse(test = dat4$sr_cov >= (dat4$Test_SR_Cov_Target * 0.8) & dat4$sr_cov <= (dat4$Test_SR_Cov_Target * 1.2),
                  yes = "#11229933",
                  no = "#99221133"),
     ylim = c(-100, 100),
     xlab = "index",
     ylab = "offset from target coverage")

# erik's plot request:
# beeswarm style plot
# x-axis is bin labels by genus / species name
# y-axis is (additional) pseudogenes (by type) per Mbp at 50x coverage
# and
# minimum coverage to achieve 1 (additional) pseudogene per Mbp (from the 1000x pin / benchmark)

# set the bins:
s1 <- sort(table(model_meta$Name[model_meta$LR_Cov != "50" &
                                   model_meta$LR_Cov != "80" &
                                   model_meta$SR_Cov_010_Reps >= 1 &
                                   model_meta$SR_Cov_025_Reps >= 1 &
                                   model_meta$SR_Cov_050_Reps >= 1 &
                                   model_meta$SR_Cov_100_Reps >= 1 &
                                   model_meta$SR_Cov_250_Reps >= 1 &
                                   model_meta$SR_Cov_500_Reps >= 1]),
           decreasing = TRUE)
s2 <- list(names(s1)[1],
           names(s1)[2],
           names(s1)[3],
           names(s1)[4],
           names(s1)[5],
           names(s1)[6],
           names(s1)[7:length(s1)])

plot(x = 0,
     y = 0,
     type = "n",
     xlim = c(0.25, 7.75),
     ylim = c(0, 10),
     xlab = "",
     ylab = "FS per Mbp @ 50x Coverage",
     xaxs = "i",
     yaxs = "i",
     xaxt = "n")
axis(side = 1,
     at = 1:7,
     labels = NA)
labs01 <- c("C. jejuni",
            "L. monocytogenes",
            "S. aureus",
            "C. coli",
            "E. faecium",
            "E. coli",
            "Other")
for (m1 in seq_along(labs01)) {
  text(labels = bquote(italic(.(labs01[m1]))),
       xpd = TRUE,
       x = m1 - 0.1,
       y = -.6,
       srt = -45,
       adj = 0)
  w1 <- which(model_meta$Name %in% s2[[m1]] &
                model_meta$LR_Cov != "50" &
                model_meta$LR_Cov != "80" &
                model_meta$SR_Cov_010_Reps >= 1 &
                model_meta$SR_Cov_025_Reps >= 1 &
                model_meta$SR_Cov_050_Reps >= 1 &
                model_meta$SR_Cov_100_Reps >= 1 &
                model_meta$SR_Cov_250_Reps >= 1 &
                model_meta$SR_Cov_500_Reps >= 1)
  val1 <- match(x = model_meta$SR_Model[w1],
                table = u_meta_1)
  val2 <- match(x = model_meta$LR_Model[w1],
                table = u_meta_2)
  
  jitterpoints <- rnorm(n = length(w1),
                        mean = m1,
                        sd = 0.10)
  points(y = sapply(X = resp_fs[w1],
                    FUN = function(x) {
                      x[101] # 50 fold cov
                    }),
         x = jitterpoints,
         col = ColVec1[val1],
         pch = val2)
  # for (m2 in seq_along(w1)) {
  #   lines(y = resp_fs[[w1[m2]]],
  #         x = convertCoverage(cov = cov_set),
  #         col = ColVec1[val1[m2]],
  #         lty = val2[m2])
  # }
  
}

plot(x = 0,
     y = 0,
     type = "n",
     xlim = c(0.25, 7.75),
     ylim = c(0, 10),
     xlab = "",
     ylab = "IS per Mbp @ 50x Coverage",
     xaxs = "i",
     yaxs = "i",
     xaxt = "n")
axis(side = 1,
     at = 1:7,
     labels = NA)
labs01 <- c("C. jejuni",
            "L. monocytogenes",
            "S. aureus",
            "C. coli",
            "E. faecium",
            "E. coli",
            "Other")
for (m1 in seq_along(labs01)) {
  text(labels = bquote(italic(.(labs01[m1]))),
       xpd = TRUE,
       x = m1 - 0.1,
       y = -.6,
       srt = -45,
       adj = 0)
  w1 <- which(model_meta$Name %in% s2[[m1]] &
                model_meta$LR_Cov != "50" &
                model_meta$LR_Cov != "80" &
                model_meta$SR_Cov_010_Reps >= 1 &
                model_meta$SR_Cov_025_Reps >= 1 &
                model_meta$SR_Cov_050_Reps >= 1 &
                model_meta$SR_Cov_100_Reps >= 1 &
                model_meta$SR_Cov_250_Reps >= 1 &
                model_meta$SR_Cov_500_Reps >= 1)
  val1 <- match(x = model_meta$SR_Model[w1],
                table = u_meta_1)
  val2 <- match(x = model_meta$LR_Model[w1],
                table = u_meta_2)
  
  jitterpoints <- rnorm(n = length(w1),
                        mean = m1,
                        sd = 0.10)
  points(y = sapply(X = resp_is[w1],
                    FUN = function(x) {
                      x[101] # 50 fold cov
                    }),
         x = jitterpoints,
         col = ColVec1[val1],
         pch = val2)
  # for (m2 in seq_along(w1)) {
  #   lines(y = resp_fs[[w1[m2]]],
  #         x = convertCoverage(cov = cov_set),
  #         col = ColVec1[val1[m2]],
  #         lty = val2[m2])
  # }
  
}

plot(x = 0,
     y = 0,
     type = "n",
     xlim = c(0.25, 7.75),
     ylim = c(0, 500),
     xlab = "",
     ylab = "Min. Cov. for 1 additional FS per Mbp",
     xaxs = "i",
     yaxs = "i",
     xaxt = "n")
axis(side = 1,
     at = 1:7,
     labels = NA)
labs01 <- c("C. jejuni",
            "L. monocytogenes",
            "S. aureus",
            "C. coli",
            "E. faecium",
            "E. coli",
            "Other")
for (m1 in seq_along(labs01)) {
  text(labels = bquote(italic(.(labs01[m1]))),
       xpd = TRUE,
       x = m1 - 0.1,
       y = -.6,
       srt = -45,
       adj = 0)
  w1 <- which(model_meta$Name %in% s2[[m1]] &
                model_meta$LR_Cov != "50" &
                model_meta$LR_Cov != "80" &
                model_meta$SR_Cov_010_Reps >= 1 &
                model_meta$SR_Cov_025_Reps >= 1 &
                model_meta$SR_Cov_050_Reps >= 1 &
                model_meta$SR_Cov_100_Reps >= 1 &
                model_meta$SR_Cov_250_Reps >= 1 &
                model_meta$SR_Cov_500_Reps >= 1)
  val1 <- match(x = model_meta$SR_Model[w1],
                table = u_meta_1)
  val2 <- match(x = model_meta$LR_Model[w1],
                table = u_meta_2)
  
  jitterpoints <- rnorm(n = length(w1),
                        mean = m1,
                        sd = 0.10)
  points(y = sapply(X = resp_fs[w1],
                    FUN = function(x) {
                      # x[101] # 50 fold cov
                      convertCoverage(cov = cov_set[which.min(abs(x - 1))])
                    }),
         x = jitterpoints,
         col = ColVec1[val1],
         pch = val2)
  
}

plot(x = 0,
     y = 0,
     type = "n",
     xlim = c(0.25, 7.75),
     ylim = c(0, 500),
     xlab = "",
     ylab = "Min. Cov. for 1 additional IS per Mbp",
     xaxs = "i",
     yaxs = "i",
     xaxt = "n")
axis(side = 1,
     at = 1:7,
     labels = NA)
labs01 <- c("C. jejuni",
            "L. monocytogenes",
            "S. aureus",
            "C. coli",
            "E. faecium",
            "E. coli",
            "Other")
for (m1 in seq_along(labs01)) {
  text(labels = bquote(italic(.(labs01[m1]))),
       xpd = TRUE,
       x = m1 - 0.1,
       y = -.6,
       srt = -45,
       adj = 0)
  w1 <- which(model_meta$Name %in% s2[[m1]] &
                model_meta$LR_Cov != "50" &
                model_meta$LR_Cov != "80" &
                model_meta$SR_Cov_010_Reps >= 1 &
                model_meta$SR_Cov_025_Reps >= 1 &
                model_meta$SR_Cov_050_Reps >= 1 &
                model_meta$SR_Cov_100_Reps >= 1 &
                model_meta$SR_Cov_250_Reps >= 1 &
                model_meta$SR_Cov_500_Reps >= 1)
  val1 <- match(x = model_meta$SR_Model[w1],
                table = u_meta_1)
  val2 <- match(x = model_meta$LR_Model[w1],
                table = u_meta_2)
  
  jitterpoints <- rnorm(n = length(w1),
                        mean = m1,
                        sd = 0.10)
  points(y = sapply(X = resp_is[w1],
                    FUN = function(x) {
                      # x[101] # 50 fold cov
                      convertCoverage(cov = cov_set[which.min(abs(x - 1))])
                    }),
         x = jitterpoints,
         col = ColVec1[val1],
         pch = val2)
  
}

save(model_meta,
     resp_fs,
     resp_is,
     file = paste0(targetdir02,
                   "/ReviewerRequest01Result_v01.RData"),
     compress = "xz")

###### -- manuscript figure initial -------------------------------------------
# beeswarm style plot
# x-axis is bin labels by genus / species name
# y-axis is (additional) pseudogenes (by type) per Mbp at 50x coverage
# and
# minimum coverage to achieve 1 (additional) pseudogene per Mbp (from the 1000x pin / benchmark)

# set the bins:
s1 <- sort(table(model_meta$Name[model_meta$LR_Cov != "50" &
                                   model_meta$LR_Cov != "80" &
                                   model_meta$SR_Cov_010_Reps >= 1 &
                                   model_meta$SR_Cov_025_Reps >= 1 &
                                   model_meta$SR_Cov_050_Reps >= 1 &
                                   model_meta$SR_Cov_100_Reps >= 1 &
                                   model_meta$SR_Cov_250_Reps >= 1 &
                                   model_meta$SR_Cov_500_Reps >= 1]),
           decreasing = TRUE)
s2 <- list(names(s1)[1],
           names(s1)[2],
           names(s1)[3],
           names(s1)[4],
           names(s1)[5],
           names(s1)[6],
           names(s1)[7],
           names(s1)[8],
           names(s1)[9],
           names(s1)[10:length(s1)])
labs01 <- paste0(c("C. jejuni",
                  "L. monocytogenes",
                  "S. aureus",
                  "C. coli",
                  "E. faecium",
                  "E. coli",
                  "H. pylori",
                  "T. pallidum",
                  "M. bovis",
                  "Other"),
                ", n = ",
                c(unname(s1)[1:9],
                  sum(s1[10:length(s1)])))

pdf(file = "tempplot.pdf",
    width = 7,
    height = 7)
layout(mat = matrix(data = 1:4,
                    byrow = TRUE,
                    ncol = 2))

###### -- top left ------------------------------------------------------------
# column title
# y-axis label
par(mar = c(2.5,3,1.5,0.75),
    mgp = c(1.85, 0.65, 0),
    cex.lab = 0.8,
    cex.axis = 1,
    cex.main = 1,
    cex.sub = 1,
    bg = "white")

plot(x = 0,
     y = 0,
     type = "n",
     xlim = c(0.25, 10.75),
     ylim = c(0, 10),
     xlab = "",
     # ylab = "Additional PGs per Mbp @ 50x Cov.",
     ylab = expression("|"*Delta*" pseudogene| per Mbp at 50-fold coverage"),
     xaxs = "i",
     yaxs = "i",
     xaxt = "n",
     main = "Frameshifts")
axis(side = 1,
     at = seq(length(labs01)),
     labels = NA)
L <- legend(x = 0.5,
            y = 8.25,
            legend = rep(NA, 7 * 5),
            # legend = u_meta_1,
            col = rep(ColVec1[1:7], 5),
            pch = rep(1:5, each = 7),
            ncol = 5,
            bty = 'n',
            x.intersp = 0.65,
            inset = 0.02,
            cex = 0.75)
legend(x = L$rect$left + L$rect$w - 0.45,
       y = L$rect$top,
       # legend = u_meta_1, # can drop 'Illumina' from names -- it's redundant
       legend = c("HiSeq 2500",
                  "HiSeq 2000",
                  "NovaSeq 6000",
                  "MiSeq",
                  "HiSeq 4000",
                  "NextSeq 500",
                  "HiSeq X Ten"),
       ncol = 1,
       x.intersp = 0.65,
       bg = NA,
       bty = "n",
       cex = 0.75)
# x-pos of first column = 0.75
text(x = 0.75 + (0.55 * 0:4),
     y = rep(8.2, 5),
     srt = 45,
     labels = c("Illumina Only",
                "+ PacBio RS",
                "+ MinION",
                "+ GridION",
                "+ PromethION"),
     cex = 0.7,
     xpd = TRUE,
     adj = 0)


for (m1 in seq_along(labs01)) {
  # text(labels = bquote(italic(.(labs01[m1]))),
  #      xpd = TRUE,
  #      x = m1 - 0.1,
  #      y = -.6,
  #      srt = -45,
  #      adj = 0)
  w1 <- which(model_meta$Name %in% s2[[m1]] &
                model_meta$LR_Cov != "50" &
                model_meta$LR_Cov != "80" &
                model_meta$SR_Cov_010_Reps >= 1 &
                model_meta$SR_Cov_025_Reps >= 1 &
                model_meta$SR_Cov_050_Reps >= 1 &
                model_meta$SR_Cov_100_Reps >= 1 &
                model_meta$SR_Cov_250_Reps >= 1 &
                model_meta$SR_Cov_500_Reps >= 1)
  val1 <- match(x = model_meta$SR_Model[w1],
                table = u_meta_1)
  val2 <- match(x = model_meta$LR_Model[w1],
                table = u_meta_2)
  
  jitterpoints <- rnorm(n = length(w1),
                        mean = m1,
                        sd = 0.10)
  if (any(jitterpoints >= m1 + .2)) {
    jitterpoints[jitterpoints >= m1 + .2] <- m1
  }
  if (any(jitterpoints <= m1 - .2)) {
    jitterpoints[jitterpoints <= m1 - .2] <- m1
  }
  points(y = sapply(X = resp_fs[w1],
                    FUN = function(x) {
                      x[101] # 50 fold cov
                    }),
         x = jitterpoints,
         col = ColVec1[val1],
         pch = val2)
}

###### -- top right -----------------------------------------------------------
# column title only
par(mar = c(2.5,3,1.5,0.75),
    mgp = c(1.85, 0.65, 0),
    cex.lab = 0.8,
    cex.axis = 1,
    cex.main = 1,
    cex.sub = 1,
    bg = "white")

plot(x = 0,
     y = 0,
     type = "n",
     xlim = c(0.25, 10.75),
     ylim = c(0, 10),
     xlab = "",
     ylab = "",
     xaxs = "i",
     yaxs = "i",
     xaxt = "n",
     main = "Internal Stops")
axis(side = 1,
     at = seq(length(labs01)),
     labels = NA)
for (m1 in seq_along(labs01)) {
  # text(labels = bquote(italic(.(labs01[m1]))),
  #      xpd = TRUE,
  #      x = m1 - 0.1,
  #      y = -.6,
  #      srt = -45,
  #      adj = 0)
  w1 <- which(model_meta$Name %in% s2[[m1]] &
                model_meta$LR_Cov != "50" &
                model_meta$LR_Cov != "80" &
                model_meta$SR_Cov_010_Reps >= 1 &
                model_meta$SR_Cov_025_Reps >= 1 &
                model_meta$SR_Cov_050_Reps >= 1 &
                model_meta$SR_Cov_100_Reps >= 1 &
                model_meta$SR_Cov_250_Reps >= 1 &
                model_meta$SR_Cov_500_Reps >= 1)
  val1 <- match(x = model_meta$SR_Model[w1],
                table = u_meta_1)
  val2 <- match(x = model_meta$LR_Model[w1],
                table = u_meta_2)
  
  jitterpoints <- rnorm(n = length(w1),
                        mean = m1,
                        sd = 0.10)
  if (any(jitterpoints >= m1 + .2)) {
    jitterpoints[jitterpoints >= m1 + .2] <- m1
  }
  if (any(jitterpoints <= m1 - .2)) {
    jitterpoints[jitterpoints <= m1 - .2] <- m1
  }
  
  points(y = sapply(X = resp_is[w1],
                    FUN = function(x) {
                      x[101] # 50 fold cov
                    }),
         x = jitterpoints,
         col = ColVec1[val1],
         pch = val2)
}

###### -- bottom left ---------------------------------------------------------
# y-axis label
# x-axis labels
par(mar = c(8,3,0,0.75),
    mgp = c(1.85, 0.65, 0),
    cex.lab = 0.8,
    cex.axis = 1,
    cex.main = 1,
    cex.sub = 1,
    bg = "white")
plot(x = 0,
     y = 0,
     type = "n",
     xlim = c(0.25, 10.75),
     ylim = c(0, 500),
     xlab = "",
     # ylab = "Coverage for one pseudogene difference per Mbp",
     ylab = expression("Coverage for 1 "*Delta*" pseudogene per Mbp"),
     xaxs = "i",
     yaxs = "i",
     xaxt = "n",
     yaxt = "n")
axis(side = 1,
     at = seq(length(labs01)),
     labels = NA)
axis(side = 2,
     at = seq(from = 0,
              to = 500,
              by = 100),
     labels = c("0",
                "100",
                "200",
                "300",
                "400",
                "> 400"))
axis.break(axis = 2,
           breakpos = 450)
for (m1 in seq_along(labs01)) {
  text(labels = bquote(italic(.(labs01[m1]))),
       xpd = TRUE,
       x = m1 - 0.3,
       y = -25,
       srt = -60,
       adj = 0,
       cex = 0.8)
  w1 <- which(model_meta$Name %in% s2[[m1]] &
                model_meta$LR_Cov != "50" &
                model_meta$LR_Cov != "80" &
                model_meta$SR_Cov_010_Reps >= 1 &
                model_meta$SR_Cov_025_Reps >= 1 &
                model_meta$SR_Cov_050_Reps >= 1 &
                model_meta$SR_Cov_100_Reps >= 1 &
                model_meta$SR_Cov_250_Reps >= 1 &
                model_meta$SR_Cov_500_Reps >= 1)
  val1 <- match(x = model_meta$SR_Model[w1],
                table = u_meta_1)
  val2 <- match(x = model_meta$LR_Model[w1],
                table = u_meta_2)
  
  jitterpoints <- rnorm(n = length(w1),
                        mean = m1,
                        sd = 0.10)
  if (any(jitterpoints >= m1 + .2)) {
    jitterpoints[jitterpoints >= m1 + .2] <- m1
  }
  if (any(jitterpoints <= m1 - .2)) {
    jitterpoints[jitterpoints <= m1 - .2] <- m1
  }
  
  transf_y_vals <- sapply(X = resp_fs[w1],
                          FUN = function(x) {
                            # some of these fits behave in weird ways, this is within
                            # expectation from the simulated data, but was mostly observed
                            # there with miseq, and more dramatically with MEGAHIT
                            # than with Unicycler
                            convertCoverage(cov = cov_set[which.min(abs(x - 1))])
                          })
  transf_y_vals[transf_y_vals > 400] <- 500
  points(y = transf_y_vals,
         x = jitterpoints,
         col = ColVec1[val1],
         pch = val2)
  
}

###### -- bottom right --------------------------------------------------------
# x axis labels only
par(mar = c(8,3,0,0.75),
    mgp = c(1.85, 0.65, 0),
    cex.lab = 0.8,
    cex.axis = 1,
    cex.main = 1,
    cex.sub = 1,
    bg = "white")
plot(x = 0,
     y = 0,
     type = "n",
     xlim = c(0.25, 10.75),
     ylim = c(0, 500),
     xlab = "",
     ylab = "",
     xaxs = "i",
     yaxs = "i",
     xaxt = "n",
     yaxt = "n")
axis(side = 1,
     at = seq(length(labs01)),
     labels = NA)
axis(side = 2,
     at = seq(from = 0,
              to = 500,
              by = 100),
     labels = c("0",
                "100",
                "200",
                "300",
                "400",
                "> 400"))
axis.break(axis = 2,
           breakpos = 450)
for (m1 in seq_along(labs01)) {
  text(labels = bquote(italic(.(labs01[m1]))),
       xpd = TRUE,
       x = m1 - 0.3,
       y = -25,
       srt = -60,
       adj = 0,
       cex = 0.8)
  w1 <- which(model_meta$Name %in% s2[[m1]] &
                model_meta$LR_Cov != "50" &
                model_meta$LR_Cov != "80" &
                model_meta$SR_Cov_010_Reps >= 1 &
                model_meta$SR_Cov_025_Reps >= 1 &
                model_meta$SR_Cov_050_Reps >= 1 &
                model_meta$SR_Cov_100_Reps >= 1 &
                model_meta$SR_Cov_250_Reps >= 1 &
                model_meta$SR_Cov_500_Reps >= 1)
  val1 <- match(x = model_meta$SR_Model[w1],
                table = u_meta_1)
  val2 <- match(x = model_meta$LR_Model[w1],
                table = u_meta_2)
  
  jitterpoints <- rnorm(n = length(w1),
                        mean = m1,
                        sd = 0.10)
  if (any(jitterpoints >= m1 + .2)) {
    jitterpoints[jitterpoints >= m1 + .2] <- m1
  }
  if (any(jitterpoints <= m1 - .2)) {
    jitterpoints[jitterpoints <= m1 - .2] <- m1
  }
  
  transf_y_vals <- sapply(X = resp_is[w1],
                          FUN = function(x) {
                            # some of these fits behave in weird ways, this is within
                            # expectation from the simulated data, but was mostly observed
                            # there with miseq, and more dramatically with MEGAHIT
                            # than with Unicycler
                            convertCoverage(cov = cov_set[which.min(abs(x - 1))])
                          })
  transf_y_vals[transf_y_vals > 400] <- 500
  
  points(y = transf_y_vals,
         x = jitterpoints,
         col = ColVec1[val1],
         pch = val2)
  
}
dev.off()

# double check a few things

# w1 <- which(model_meta$Name %in% s2[[4]] &
#               model_meta$LR_Cov != "50" &
#               model_meta$LR_Cov != "80" &
#               model_meta$SR_Cov_010_Reps >= 1 &
#               model_meta$SR_Cov_025_Reps >= 1 &
#               model_meta$SR_Cov_050_Reps >= 1 &
#               model_meta$SR_Cov_100_Reps >= 1 &
#               model_meta$SR_Cov_250_Reps >= 1 &
#               model_meta$SR_Cov_500_Reps >= 1)
# 
# suppressWarnings(plot(x = 0,
#                       y = 0,
#                       type = "n",
#                       xlim = c(5, 500),
#                       ylim = c(0, 10),
#                       log = "x",
#                       xaxs = "i",
#                       yaxs = "i",
#                       xlab = "Coverage",
#                       ylab = expression("|"*Delta*" Frameshifts per Mbp|"),
#                       main = s2[[4]]))
# for (m2 in seq_along(w1)) {
#   lines(y = resp_fs[[w1[m2]]],
#         x = convertCoverage(cov = cov_set),
#         col = ColVec1[val1[m2]],
#         lty = val2[m2])
# }


