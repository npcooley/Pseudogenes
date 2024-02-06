###### -- build a jobmap ------------------------------------------------------
# a lot is going on here, comments are sprinkled throughout,
# but the gist is that this experiment is a set of re-assemblies being performed
# to address a reviewer critique
# this experiment was sort of originally performed wholesale with neisseria in a way that
# didn't really give us access to what was asked for by the reviewer,
# i.e. does scaling down coverage for real reads on the SRA provide an increase
# in pseudogenes the way that it does with the simulated reads experiment?
# with our original neisseria experiment, we didn't really map out that downscaling
# in a way that would work for that comparison, here we're taking pains to do that
# though we are also subsampling hypothetically appropriate read sets down to a managable
# subset because this experiment will balloon pretty dramatically, i.e. for each read set
# we're going to assemble at 10, 25 50, 100, 250, 500 and 1000 fold coverage, in triplicate

suppressMessages(library(SynExtend))

load(file = "~/Repos/PseudogenePlots/InputData/GenBank_Assembly_Summary.RData",
     verbose = TRUE)
load(file = "~/Repos/PseudogenePlots/InputData/GenBank_SRA_Summary.RData",
     verbose = TRUE)

# for every biosample that has at least one illumina read set with greater than
# approx. 100x coverage build out a set of jobs in triplicate that reassemble
# those reads at 5x 10x 20x 50x 100x 250x 500x and 1000x coverage, if the biosample has
# an accompanying set of long reads that are less than 1000x approximate coverage
# include hybrid reassemblies at 20x 50x and 80x long read coverage

# if there are greater than 100 candidates for short reads alone, subsample to
# 100, same for ONT candidates and PB candidates -- short reads are subsampled by reported sequencer model

# subset to the intersect
u1 <- intersect(x = unique(SRA_Meta$BioSample),
                y = unique(GenBankAccessions$Biosample))

SRA_SubSet <- SRA_Meta[SRA_Meta$BioSample %in% u1, ]
GBA_SubSet <- GenBankAccessions[GenBankAccessions$Biosample %in% u1, ]

# append assembly size to sra table and then append
# approx coverage
mat1 <- match(x = SRA_SubSet$BioSample,
              table = GBA_SubSet$Biosample)
mat2 <- GBA_SubSet$Total_Length[mat1]
mat3 <- GBA_SubSet$FTP[mat1]
SRA_SubSet <- cbind(SRA_SubSet,
                  "assembly_size" = mat2,
                  "assembly_ftp" = mat3)
# when things are smooth, you can expect 90-95 of reads to map cleanly with bowtie
# with 5% or so mapping twice or more
# when things aren't smooth ... shrug?
w1 <- (SRA_SubSet$bases * 0.90) / SRA_SubSet$assembly_size
SRA_SubSet <- cbind(SRA_SubSet,
                    "apprx_cov" = w1)

# drop non-wgs runs and where approx. coverage is outside the range that we're interested in
# step down short read only assemblies from 500 to 5
# when long reads are available, perform hybrid assemblies
# at 20x long read cov, 50x long read cov, and 80x long read cov
# the SRA toolkit can be a bit slow to pull data so we're dropping
# VERY LARGE long read datasets
SRA_SubSet <- SRA_SubSet[SRA_SubSet$LibraryStrategy == "WGS" &
                           SRA_SubSet$LibrarySource == "GENOMIC" &
                           SRA_SubSet$Platform %in% c("ILLUMINA",
                                                      "PACBIO_SMRT",
                                                      "OXFORD_NANOPORE"), ]
w1 <- ifelse(test = SRA_SubSet$Platform == "ILLUMINA",
             yes = ifelse(test = SRA_SubSet$apprx_cov >= 1000 &
                            SRA_SubSet$apprx_cov <= 10000,
                          yes = TRUE,
                          no = FALSE),
             no = ifelse(test = SRA_SubSet$apprx_cov <= 1000 &
                           SRA_SubSet$apprx_cov >= 81,
                         yes = TRUE,
                         no = FALSE))

SRA_SubSet <- SRA_SubSet[w1, ]

# split out to short read only, short read + ont, and short read + pb
w1 <- tapply(X = SRA_SubSet$Platform,
             INDEX = SRA_SubSet$BioSample,
             FUN = function(x) {
               paste(sort(unique(x)), collapse = " + ")
             })

# there are other combinations, but for our sanity and experimental simplicity
# we'll focus on these

target_reads <- list("sr_only" = SRA_SubSet[SRA_SubSet$BioSample %in% names(w1[w1 == "ILLUMINA"]), ],
                     "sr_ont" = SRA_SubSet[SRA_SubSet$BioSample %in% names(w1[w1 == "ILLUMINA + OXFORD_NANOPORE"]), ],
                     "sr_pb" = SRA_SubSet[SRA_SubSet$BioSample %in% names(w1[w1 == "ILLUMINA + PACBIO_SMRT"]), ])

# for each set, for each biosample, grab the largest read set for each available technology

for (m1 in seq_along(target_reads)) {
  
  w1 <- target_reads[[m1]]$BioSample
  t1 <- target_reads[[m1]]$Platform
  w2 <- unique(w1)
  keepval <- vector(mode = "logical",
                    length = length(w1))
  for (m2 in seq_along(w2)) {
    # there should always be illumina reads
    w3 <- which(w1 == w2[m2] &
                  t1 == "ILLUMINA")
    if (length(w3) == 1) {
      keepval[w3] <- TRUE
    } else {
      w4 <- which.max(target_reads[[m1]]$apprx_cov[w3])
      keepval[w3[w4]] <- TRUE
    } # end if else
    
    # check non-illumina
    w3 <- which(w1 == w2[m2] &
                  t1 != "ILLUMINA")
    if (length(w3) == 0) {
      # do nothing
    } else if (length(w3) == 1) {
      keepval[w3] <- TRUE
    } else {
      w4 <- which.max(target_reads[[m1]]$apprx_cov[w3])
      keepval[w3[w4]] <- TRUE
    }
  } # end m2 loop
  
  target_reads[[m1]] <- target_reads[[m1]][keepval, ]
}

# drop everyone who has a name that might not fit with PGAP

for (m1 in seq_along(target_reads)) {
  # deal with species names
  n1 <- target_reads[[m1]]$ScientificName
  g1 <- grepl(pattern = "^candidatus",
              ignore.case = TRUE,
              x = n1)
  if (any(g1)) {
    n1[g1] <- gsub(pattern = "^candidatus ",
                   x = n1[g1],
                   replacement = "",
                   ignore.case = TRUE)
  }
  n1 <- gsub(pattern = "'",
             replacement = "",
             x = n1)
  n1 <- gsub(pattern = '"',
             replacement = "",
             x = n1)
  n1 <- regmatches(m = gregexpr(pattern = "^[^ ]+ [^ ]+",
                                       text = n1),
                   x = n1)
  w1 <- lengths(n1) == 1
  
  
  target_reads[[m1]] <- cbind(target_reads[[m1]][w1, ],
                              "PGAP_name" = unlist(n1[w1]))
}

# subsample to at max 100 reps of the read source type
# for sr only, we can do 100 of each model type that has more than 100 possible targets
# don't do anyone with less
# we play a little fast and loose here with selecting things by model, if the mixed read datasets
# grow in size in an analysis that occurs with a later collection of the source data tables,
# this loop may need to be written a little more strictly

res01 <- vector(mode = "list",
                length = length(target_reads))

for (m1 in seq_along(target_reads)) {
  set.seed(1986)
  if (any(target_reads[[m1]]$Platform != "ILLUMINA")) {
    # non-illumina platforms present, ask for the model table for the non-illumina reads
    t1 <- table(target_reads[[m1]]$Model[target_reads[[m1]]$Platform != "ILLUMINA"])
  } else {
    # only illumina as a platform, just grab the models for everything
    # and call it good
    t1 <- table(target_reads[[m1]]$Model)
  }
  if (any(t1 >= 100)) {
    # subsample to 100 representatives for each model that appears at least 100 times
    w1 <- names(t1[t1 >= 100])
    res01[[m1]] <- vector(mode = "list",
                          length = length(w1))
    for (m2 in seq_along(w1)) {
      w2 <- sample(x = which(target_reads[[m1]]$Model == w1[m2]),
                   size = 100,
                   replace = FALSE)
      w3 <- target_reads[[m1]]$BioSample[w2]
      curr_res <- target_reads[[m1]][target_reads[[m1]]$BioSample %in% w3, ]
      res01[[m1]][[m2]] <- curr_res
     }
  } else {
    # subsample the total target biosamples down to 100 representatives
    res01[[m1]] <- vector(mode = "list",
                          length = 1L)
    w1 <- unique(target_reads[[m1]]$BioSample)
    if (length(w1) >= 100) {
      w2 <- sample(x = w1,
                   size = 100,
                   replace = FALSE)
      curr_res <- target_reads[[m1]][target_reads[[m1]]$BioSample %in% w2, ]
    } else {
      curr_res <- target_reads[[m1]]
    }
    res01[[m1]][[1]] <- curr_res
  }
}

res01 <- unlist(res01,
                recursive = FALSE)
# loop through res01, and build out the jobmap
# we need:
# LR run accession
# SR run accession
# LR cov target
# SR cov target
# LR model
# SR model
# run approx cov

# we play a little fast and loose here, short reads are always
# illumina given the subsetting we've performed
# long reads are either PB or ONT,
# there is no long read only analysis
res02 <- vector(mode = "list",
                length = length(res01))

for (m1 in seq_along(res01)) {
  t1 <- table(res01[[m1]]$Platform)
  if (length(t1) == 1L) {
    # only illumina present
    int_res <- expand.grid(list("biosample" = res01[[m1]]$BioSample,
                                "sr_cov_target" = c(5, 10, 25, 50, 100, 250, 500, 1000),
                                "replicate" = 1:3))
    mat1 <- match(x = int_res$biosample,
                  table = res01[[m1]]$BioSample)
    int_res <- data.frame("assembly_size" = res01[[m1]]$assembly_size[mat1],
                          "sr_r" = res01[[m1]]$Run[mat1],
                          "lr_r" = rep("x",
                                       nrow(int_res)),
                          "sr_p" = res01[[m1]]$Platform[mat1],
                          "lr_p" = rep("x",
                                       nrow(int_res)),
                          "sr_cov_target" = int_res$sr_cov_target,
                          "lr_cov_target" = rep("x",
                                                nrow(int_res)),
                          "sr_approx_cov" = res01[[m1]]$apprx_cov[mat1],
                          "lr_approx_cov" = rep(-1,
                                                nrow(int_res)),
                          "sr_mod" = res01[[m1]]$Model[mat1],
                          "lr_mod" = rep("x",
                                         nrow(int_res)),
                          "replicate" = int_res$replicate)
  } else {
    # illumina + long reads
    ph1 <- res01[[m1]][res01[[m1]]$Platform == "ILLUMINA", ]
    ph2 <- res01[[m1]][res01[[m1]]$Platform != "ILLUMINA", ]
    w1 <- unique(res01[[m1]]$BioSample)
    int_res <- expand.grid(list("biosample" = w1,
                                "sr_cov_target" = c(5, 10, 25, 50, 100, 250, 500, 1000),
                                "lr_cov_target" = c(20, 50, 80),
                                "replicate" = 1:3))
    mat1 <- match(x = int_res$biosample,
                  table = ph1$BioSample)
    mat2 <- match(x = int_res$biosample,
                  table = ph2$BioSample)
    
    int_res <- data.frame("assembly_size" = ph1$assembly_size[mat1],
                          "sr_r" = ph1$Run[mat1],
                          "lr_r" = ph2$Run[mat2],
                          "sr_p" = ph1$Platform[mat1],
                          "lr_p" = ph2$Platform[mat2],
                          "sr_cov_target" = int_res$sr_cov_target,
                          "lr_cov_target" = int_res$lr_cov_target,
                          "sr_approx_cov" = ph1$apprx_cov[mat1],
                          "lr_approx_cov" = ph2$apprx_cov[mat2],
                          "sr_mod" = ph1$Model[mat1],
                          "lr_mod" = ph2$Model[mat2],
                          "replicate" = int_res$replicate)
  }
  res02[[m1]] <- int_res
}
res02 <- do.call(rbind,
                 res02)

res03 <- res02[res02$sr_cov_target != 1000, ]
res04 <- res02[res02$sr_cov_target == 1000 &
                 res02$replicate == 1, ]

save(res03,
     res04,
     file = "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/JobMap_A_B.RData",
     compress = "xz")

set.seed(1986)
write.table(x = res04[, 1:9],
            row.names = TRUE,
            col.names = FALSE,
            append = FALSE,
            quote = FALSE,
            file = "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/JobMapA.txt")

write.table(x = res04[sample(x = nrow(res04),
                             size = 10,
                             replace = FALSE), 1:9],
            row.names = TRUE,
            col.names = FALSE,
            append = FALSE,
            quote = FALSE,
            file = "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/TestMapA.txt")

write.table(x = res03[, 1:9],
            row.names = TRUE,
            col.names = FALSE,
            append = FALSE,
            quote = FALSE,
            file = "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/JobMapB.txt")

write.table(x = res03[sample(x = nrow(res03),
                             size = 10,
                             replace = FALSE), 1:9],
            row.names = TRUE,
            col.names = FALSE,
            append = FALSE,
            quote = FALSE,
            file = "~/Repos/PseudogenePlots/OSG_Jobs/20231115_ReviewRequest01/TestMapB.txt")

