###### -- Mapping Sanger reads ------------------------------------------------

suppressMessages(library(SynExtend))
suppressMessages(library(Rsamtools))

# R requires commandline access to:
# bowtie2 -- in path maybe was installed without conda?
# minimap2 -- managed by conda locally
# samtools -- managed by conda locally
# and their dependencies
# i.e. they must be in R's PATH `Sys.getenv("PATH")` - conda managed things required the conda path stuff

# will generate an RData file of the mappings named `ReadMappingResults.RData` in R's working directory
# the naming scheme in the xlsx file wasn't conveniently parsable, so it was not used and names and IDs were hard-coded.
# things were the second assembly were a little weird so it was not evaluated.

# ARGS can be re-hard coded here based on your file locations, or commented out and fed to RScript as arguments in the command line
# ARGS <- commandArgs(trailingOnly = TRUE)
ARGS <- c("~/Downloads/Primers and Plate Layouts.xlsx", # never actually called on because it's not entirely parseable
          "~/Downloads/Sequence results", # folder containing sequencing results
          "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/593/565/GCF_001593565.1_ASM159356v1/", # nr 51487 on spreadsheet
          "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/350/925/GCF_000350925.1_Esch_coli_KTE181_V1/", # nr 32771 on spreadsheet
          "~/Misc/RFunctions/20220407_BAMPairwise.R", # convert BAM file to a list of DNAStringSets
          "~/Misc/RFunctions/20220406_Minimap2Mapping.R", # map reads to assembly with Minimap
          "~/Misc/RFunctions/20220407_Bowtie2Pipeline.R", # map reads to assembly with Bowtie2
          "~/Misc/RFunctions/20220410_Map_Reads_w_FindSynteny.R") # map reads to assembly / features with FindSynteny

###### -- import pipeline scripts and things ----------------------------------

source(file = ARGS[5],
       echo = FALSE)
source(file = ARGS[6],
       echo = FALSE)
source(file = ARGS[7],
       echo = FALSE)
source(file = ARGS[8],
       echo = FALSE)

###### -- import sanger reads -------------------------------------------------

FILES01 <- list.files(path = ARGS[2L],
                      full.names = TRUE,
                      pattern = ".phd.1")

NS <- unlist(regmatches(x = FILES01,
                        m = gregexpr(pattern = "(?<=results/)([0-9]+)(?=_PREMIX)",
                                     text = FILES01,
                                     perl = TRUE)))


strs <- qual <- vector(mode = "character",
                       length = length(NS))

for (m1 in seq_along(NS)) {
  r <- readLines(FILES01[m1])
  r <- r[(which(r=="BEGIN_DNA")[1] + 1):(which(r=="END_DNA")[1] - 1)]
  r <- strsplit(x = r,
                split = " ",
                fixed = TRUE)
  # strs[m1] <- paste(sapply(r, `[`, 1L), collapse="")
  strs[m1] <- paste(sapply(X = r,
                           FUN = `[`,
                           1L),
                    collapse = "")
  qual[m1] <- paste(rawToChar(as.raw(as.numeric(sapply(r, `[`, 2L)) + 33)), collapse="")
}

strs <- DNAStringSet(strs)
names(strs) <- NS
qual2 <- PhredQuality(qual)

# TrimDNA(strs, "", "", "sequence", qual2)

# Seqs1 <- TrimDNA(myDNAStringSet = strs,
#                  leftPatterns = "",
#                  rightPatterns = "",
#                  type = "sequence",
#                  quality = qual2)
# 
# Seqs2 <- TrimDNA(myDNAStringSet = strs,
#                  leftPatterns = "",
#                  rightPatterns = "",
#                  type = "both",
#                  quality = qual2)
# 
# q1 <- Seqs2[[1]]
# q2 <- qual2
# q3 <- subseq(x = qual2,
#              start = q1@start,
#              end = (q1@start + q1@width - 1))
# q3 <- q3[q1@width > 0]

x <- QualityScaledDNAStringSet(x = strs,
                               quality = qual2)

Seqs3 <- TrimDNA(myDNAStringSet = x,
                 leftPatterns = "",
                 rightPatterns = "",
                 type = "sequence",
                 quality = x@quality)
Seqs3 <- QualityScaledDNAStringSet(x = Seqs2[[2]],
                                   quality = q3)

###### -- import assemblies and gene calls ------------------------------------

ASSEMBLY1 <- strsplit(x = ARGS[3],
                      split = "/",
                      fixed = TRUE)
GENECALLS1 <- paste0(ARGS[3],
                     ASSEMBLY1[[1]][10],
                     "_genomic.gff.gz")
ASSEMBLY1 <- paste0(ARGS[3],
                    ASSEMBLY1[[1]][10],
                    "_genomic.fna.gz")
ASSEMBLY2 <- strsplit(x = ARGS[4],
                      split = "/",
                      fixed = TRUE)
GENECALLS2 <- paste0(ARGS[4],
                     ASSEMBLY2[[1]][10],
                     "_genomic.gff.gz")
ASSEMBLY2 <- paste0(ARGS[4],
                    ASSEMBLY2[[1]][10],
                    "_genomic.fna.gz")

ASSEMBLY1 <- readDNAStringSet(filepath = ASSEMBLY1)
ASSEMBLY2 <- readDNAStringSet(filepath = ASSEMBLY2)
GENECALLS1 <- gffToDataFrame(GFF = GENECALLS1,
                             Verbose = TRUE)
GENECALLS2 <- gffToDataFrame(GFF = GENECALLS2,
                             Verbose = TRUE)

###### -- who belongs where ---------------------------------------------------

# Primer names targeted at assembly 1
Ref1ReadTargets <- c("37",
                     "294",
                     "318",
                     "322",
                     "362",
                     "630",
                     "707",
                     "1305",
                     "1307",
                     "1667",
                     "1756",
                     "2419",
                     "2492",
                     "2826",
                     "3328",
                     "3452",
                     "3363",
                     "3929",
                     "3948",
                     "4126",
                     "4193",
                     "4253",
                     "4569",
                     "4705",
                     "594",
                     "908",
                     "3295",
                     "3762",
                     "3452",
                     "2272")

# Primary names targeted at assembly 2
Ref2ReadTargets <- c("314",
                     "548", # appears twice for some reason
                     "1425", # appears twice for some reason
                     "798",
                     "805",
                     "1374",
                     "377",
                     "2554")

###### -- call pipelines and break out pairs ----------------------------------

res1 <- Minimap2Pipeline(Reference = ASSEMBLY1,
                         Reads = Seqs3[names(Seqs3) %in% Ref1ReadTargets])
res1a <- scanBam(res1)
res1b <- PairwiseBAMImport(BAMImport = res1a[[1]],
                           Reference = ASSEMBLY1)

res2 <- Bowtie2Pipeline(Reference = ASSEMBLY1,
                        Reads = Seqs3[names(Seqs3) %in% Ref1ReadTargets])
res2a <- scanBam(res2)
res2b <- PairwiseBAMImport(BAMImport = res2a[[1]],
                           Reference = ASSEMBLY1)

res3 <- Minimap2Pipeline(Reference = ASSEMBLY2,
                         Reads = Seqs1[names(Seqs1) %in% Ref2ReadTargets])
res3a <- scanBam(res3)
res3b <- PairwiseBAMImport(BAMImport = res3a[[1]],
                           Reference = ASSEMBLY2)

res4 <- Bowtie2Pipeline(Reference = ASSEMBLY2,
                        Reads = Seqs3[(names(Seqs3) %in% Ref2ReadTargets)])
res4a <- scanBam(res4)
res4b <- PairwiseBAMImport(BAMImport = res4a[[1]],
                           Reference = ASSEMBLY2)

# res2b and res4b are the bowtie stringent mappings for the matched sanger seqs
# the targets for assembly 2 are not pseudogenes so we may not find them to be interesting


# sum(names(Seqs1) %in% Ref1ReadTargets)
# [1] 10

# sum(res2a[[1]]$qname %in% Ref1ReadTargets)
# [1] 9

# res2a[[1]]$cigar
# [1] "85M"       "120M1I20M" "104M1I64M" "100M1I73M" "222M"      "151M1S"    "92M"       "90M1I89M"  "57M1I115M"


res5 <- SynMappingPipeline(Reads = Seqs3[names(Seqs3) %in% Ref1ReadTargets],
                           Reference = ASSEMBLY1,
                           GeneCalls = GENECALLS1)
res6 <- SynMappingPipeline(Reads = Seqs3[names(Seqs3) %in% Ref2ReadTargets],
                           Reference = ASSEMBLY2,
                           GeneCalls = GENECALLS2)

# visual inspections,
# repeat for all 10 alignments
# BrowseSeqs(AlignSeqs(res5[[2]][[1]]),
#            highlight = 0)
# GENECALLS1$Type[GENECALLS1$ID == names(res5[[2]][[1]][1])]

# 1: pseudogene - no disagreement
# 2: pseudogene - mismatch
# 3: pseudogene - deleted G in poly G
# 4: pseudogene - no disagreement
# 5: pseudogene - no disagreement
# 6: pseudogene - multiple mismatches + deleted G in poly G
# 7: pseudogene - deleted C in poly C
# 8: pseudogene - deleted C in poly C
# 9: pseudogene - deleted C in poly C
# 10: pseudogene - deleted G in poly G - run size is only 3

save(res1a, # minimap BAM file import -- assembly 1
     res1b, # minimap paired seqs -- assembly 1
     res2a, # bowtie2 BAM file import -- assembly 1
     res2b, # bowtie2 paired seqs -- assembly 1
     res3a, # minimap BAM file import -- assembly 2
     res3b, # minimap paired seqs -- assembly 2
     res4a, # bowtie2 BAM file import -- assembly 2
     res4b, # bowtie2 paired seqs -- assembly 2
     res5, # synteny mappings with reads and stuff
     file = "~/ReadMappingResults.RData",
     compress = "xz")
