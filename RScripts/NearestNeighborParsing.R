###### -- load in and reshape data to make access easier and cleaner ----------

# for rollmean function
library(zoo)

# we have ANI / AFR for all within-genera pairs for genera with enough tips
# and we have distances for putative orthologs for full length pairs for:
# Pseudo-Pseudo pairs
# Gene-Gene pairs
# and Pseudo-Gene / Gene-Pseudo pairs -- directionality likely can matter,
# though is not guaranteed to ...

ColVec1 <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000')
ColVec2 <- paste0(ColVec1,
                  "33")

FILES01 <- list.files(path = "~/Repos/Pseudogenes/GeneraSynopsis/",
                      full.names = TRUE)

###### -- Approximate AUC -----------------------------------------------------
# source: https://stackoverflow.com/questions/4954507/calculate-the-area-under-a-curve
# "...effectively calculates the area using the trapezoidal method..."

# ApproximateAUC <- function(x, y) {
#   
#   if (length(x) != length(y)) {
#     stop("Vector lengths do not match")
#   }
#   o <- order(x)
#   if (x[o][1L] != 0L) {
#     x <- c(0, x[o])
#     y <- c(y[o][1], y[o])
#   } else {
#     x <- x[o]
#     y <- y[o]
#   }
#   
#   AUC <- sum(diff(x) * (head(y,
#                              -1) + tail(y,
#                                         -1))) / 2
#   return(AUC)
# }

# don't need this
# load(file = "~/Repos/Pseudogenes/SearchResults.RData",
#      verbose = TRUE)
load(file = "~/Repos/Pseudogenes/SearchResults2.RData",
     verbose = TRUE)

# metadata, GG dists, PP dists, PG dists, GP, dists
RES1 <- RES2 <- RES3 <- RES4 <- RES5 <- vector(mode = "list",
                                               length = length(FILES01))

pBar <- txtProgressBar(style = 1L)
PBAR <- length(FILES01)
for (m1 in seq_along(FILES01)) {
  load(FILES01[m1],
       verbose = FALSE)
  
  RES1[[m1]] <- data.frame("p1" = Res6[, 1L],
                           "p2" = Res6[, 2L],
                           "ANI" = Res1,
                           "AFR" = Res2,
                           "IS" = Res3[, 3L] * 1000, # abs of relative count diff in MB
                           "FS" = Res3[, 1L] * 1000) # abs of relative count diff in MB
  # GG
  RES2[[m1]] <- lapply(X = Res5,
                       FUN = function(x) {
                         x[[5]]
                       })
  # PP
  RES3[[m1]] <- lapply(X = Res5,
                       FUN = function(x) {
                         x[[8]]
                       })
  # GP
  RES4[[m1]] <- lapply(X = Res5,
                       FUN = function(x) {
                         x[[6]]
                       })
  # PG
  RES5[[m1]] <- lapply(X = Res5,
                       FUN = function(x) {
                         x[[7]]
                       })
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)

RES1 <- do.call(rbind,
                RES1)

RES2 <- unlist(RES2, recursive = FALSE)
RES3 <- unlist(RES3, recursive = FALSE)
RES4 <- unlist(RES4, recursive = FALSE)
RES5 <- unlist(RES5, recursive = FALSE)

# subset to high ANI and any recorded pairs with a pseudogene
# w1 <- RES1$ANI > 0.9 &
#   (lengths(RES3) > 0 | lengths(RES4) > 0 | lengths(RES5) > 0)

# RES1 <- RES1[w1, ]
# RES2 <- RES2[w1]
# RES3 <- RES3[w1]
# RES4 <- RES4[w1]
# RES5 <- RES5[w1]

# subset again to nearest neighbor only
UIDs <- unique(unlist(RES1[, 1:2]))
NN <- vector(mode = "integer",
             length = length(UIDs))
pBar <- txtProgressBar(style = 1L)
PBAR <- length(NN)
for (m1 in seq_along(NN)) {
  w2 <- RES1$p1 == UIDs[m1] | RES1$p2 == UIDs[m1]
  w3 <- which.max(RES1$ANI[w2])
  w4 <- which(w2)[w3]
  NN[m1] <- w4
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)

UNN <- sort(unique(NN))

RES1 <- RES1[UNN, ]
RES2 <- RES2[UNN]
RES3 <- RES3[UNN]
RES4 <- RES4[UNN]
RES5 <- RES5[UNN]

# GG vs PP, GP vs PG, respectively for nearest neighbors
RES6 <- RES7 <- RES12 <- vector(mode = "list",
                                length = length(UNN))
# AUCs for both curves
RES8 <- RES9 <- RES10 <- RES11 <- RES13 <- RES14 <- vector(mode = "numeric",
                                                           length = length(UNN))

AHEAD1 <- lengths(RES2)
AHEAD2 <- lengths(RES3)
AHEAD3 <- lengths(RES4)
AHEAD4 <- lengths(RES5)

pBar <- txtProgressBar(style = 1)
PBAR <- length(UNN)
for (m1 in seq_along(UNN)) {
  
  if (AHEAD1[m1] > 0 & AHEAD2[m1] > 0) {
    z1 <- suppressWarnings(ks.test(x = RES2[[m1]], y = RES3[[m1]]))
    if (is(object = z1,
           class2 = "ks.test")) {
      RES6[[m1]] <- c(unname(z1$statistic), unname(z1$p.value))
    } else {
      RES6[[m1]] <- c(-1, -1)
    }
    
    x <- sort(c(RES2[[m1]], 1))
    y <- seq_along(x) / length(x)
    id <- order(x)
    RES8[m1] <- sum(diff(x[id])*rollmean(y[id],2))
    
    x <- sort(c(RES3[[m1]], 1))
    y <- seq_along(x) / length(x)
    id <- order(x)
    RES9[m1] <- sum(diff(x[id])*rollmean(y[id],2))
  } else {
    RES6[[m1]] <- c(-1, -1)
  }
  
  if (AHEAD3[m1] > 0 & AHEAD4[m1] > 0) {
    z2 <- suppressWarnings(ks.test(x = RES4[[m1]], y = RES5[[m1]]))
    if (is(object = z2,
           class2 = "ks.test")) {
      RES7[[m1]] <- c(unname(z2$statistic), unname(z2$p.value))
    } else {
      RES7[[m1]] <- c(-1, -1)
    }
    
    x <- sort(c(RES4[[m1]], 1))
    y <- seq_along(x) / length(x)
    id <- order(x)
    RES10[m1] <- sum(diff(x[id])*rollmean(y[id],2))
    
    x <- sort(c(RES5[[m1]], 1))
    y <- seq_along(x) / length(x)
    id <- order(x)
    RES11[m1] <- sum(diff(x[id])*rollmean(y[id],2))
  } else {
    RES7[[m1]] <- c(-1, -1)
  }
  
  if (AHEAD1[m1] > 0 & (AHEAD3[m1] > 0 | AHEAD4[m1] > 0)) {
    z3 <- suppressWarnings(ks.test(x = RES2[[m1]], y = c(RES5[[m1]], RES4[[m1]])))
    if (is(object = z3,
           class2 = "ks.test")) {
      RES12[[m1]] <- c(unname(z3$statistic), unname(z3$p.value))
    } else {
      RES12[[m1]] <- c(-1, -1)
    }
    
    x <- sort(c(RES2[[m1]], 1))
    y <- seq_along(x) / length(x)
    id <- order(x)
    RES14[m1] <- sum(diff(x[id])*rollmean(y[id],2))
    
    x <- sort(c(RES5[[m1]], RES4[[m1]], 1))
    y <- seq_along(x) / length(x)
    id <- order(x)
    RES13[m1] <- sum(diff(x[id])*rollmean(y[id],2))
  } else {
    RES12[[m1]] <- c(-1, -1)
  }
  
  
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)

z <- do.call(rbind,
             RES6)
z <- data.frame("D" = z[, 1L],
                "p_val" = z[, 2L],
                "GG_AUC" = RES8,
                "PP_AUC" = RES9)
w5 <- which(z$D < 0)
z <- z[-w5, ]
RES6 <- z


z <- do.call(rbind,
             RES7)
z <- data.frame("D" = z[, 1L],
                "p_val" = z[, 2L],
                "GP_AUC" = RES10,
                "PG_AUC" = RES11)
w5 <- which(z$D < 0)
z <- z[-w5, ]
RES7 <- z

z <- do.call(rbind,
             RES12)
z <- data.frame("D" = z[, 1L],
                "p_val" = z[, 2L],
                "GG_AUC" = RES14,
                "PG_GP_AUC" = RES13)
w5 <- which(z$D < 0)
z <- z[-w5, ]
RES12 <- z

plot(x = ifelse(test = RES6$GG_AUC > RES6$PP_AUC,
                yes = RES6$D,
                no = RES6$D * -1),
     y = -log(RES6$p_val),
     pch = 46)
plot(x = ifelse(test = RES12$GG_AUC > RES12$PG_GP_AUC,
                yes = RES12$D,
                no = RES12$D * -1),
     y = -log(RES12$p_val),
     pch = 46)

save(RES1,
     RES6,
     RES7,
     RES12,
     file = "~/Repos/PseudogenePlots/InputData/Counts_Orthos_v01.RData",
     compress = "xz")


# as ANI approaches 100 do pseudogene counts converge to zero
layout(mat = matrix(data = c(1,2),
                    nrow = 1))
par(mar = c(4,3,2,0.75),
    mgp = c(2, 0.75, 0))
plot(x = RES1$ANI,
     y = RES1$IS,
     pch = 46,
     col = ColVec2[1L],
     xlab = "ANI",
     ylab = "\u0394 IS per MB",
     ylim = c(0, 10))
spfit1 <- smooth.spline(x = RES1$ANI,
                        y = RES1$IS,
                        df = 5)
spfit2 <- predict(spfit1)
lines(x = spfit2$x,
      y = spfit2$y,
      col = "black")
plot(x = RES1$ANI,
     y = RES1$FS,
     pch = 46,
     col = ColVec2[2L],
     xlab = "ANI",
     ylab = "\u0394 IS per MB",
     ylim = c(0, 30))
spfit1 <- smooth.spline(x = RES1$ANI,
                        y = RES1$FS,
                        df = 5)
spfit2 <- predict(spfit1)
lines(x = spfit2$x,
      y = spfit2$y,
      col = "black")

# plot ANI vs the p-value of GG-PP distributions


layout(mat = matrix(data = c(1,2),
                    nrow = 1))
par(mar = c(4,3,2,0.75),
    mgp = c(2, 0.75, 0))
plot(x = ifelse(test = RES6$GG_AUC > RES6$PP_AUC,
                yes = RES6$D,
                no = RES6$D * -1),
     y = -log(RES6$p_val),
     pch = 46,
     col = ColVec1[1],
     xlab = "Statistic",
     ylab = "-(log(p-val))")
par(mar = c(4,2.75,2,1))
plot(x = ifelse(test = RES12$GG_AUC > RES12$PG_GP_AUC,
                yes = RES12$D,
                no = RES12$D * -1),
     y = -log(RES12$p_val),
     pch = 46,
     col = ColVec1[2],
     xlab = "Statistic",
     ylab = "-(log(p-val))")

###### -- something else ------------------------------------------------------

# se1 <- names(head(sort(table(EntrezResults$SpeciesName[grepl(pattern = "Xanthomonas",
#                                                              x = EntrezResults$SpeciesName)]), decreasing = TRUE)))
# 
# plot(x = ((PseudosCondensed[, 1L] / EntrezResults$Total_Length) * 1000000)[EntrezResults$SpeciesName %in% se1],
#      y = ((PseudosCondensed[, 4L] / EntrezResults$Total_Length) * 1000000)[EntrezResults$SpeciesName %in% se1],
#      pch = 20,
#      col = ColVec2[match(x = EntrezResults$SpeciesName[EntrezResults$SpeciesName %in% se1],
#                          table = se1)])








