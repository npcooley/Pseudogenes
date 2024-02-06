###### -- grab SRA accessions for assemblies ----------------------------------

# probably unnecessary as all we're doing is grabbing stuff from entrez
suppressMessages(library(SynExtend))

ARGS <- commandArgs(trailingOnly = TRUE)

setwd(ARGS[1L])

# form the specified directory, grab a data file that should exist
load(file = paste0(getwd(),
                   "/SearchResults.RData"),
     verbose = TRUE)

#unique biosample accession
ubs <- unique(EntrezResults$Biosample)


TEMP01 <- tempfile()
TOTALLINES <- length(ubs)
# search in batches of 500
INITIALLINE <- 1L
BATCHSIZE <- 500L
FINALLINE <- 500L
TOTALSEARCHES <- ceiling(TOTALLINES / BATCHSIZE)
TEMP02 <- tempfile()
RES <- vector(mode = "list",
              length = TOTALSEARCHES)
pBar <- txtProgressBar(style = 1L)
CONTINUE <- FALSE

# for (m1 in 238:TOTALSEARCHES) {
for (m1 in seq_len(TOTALSEARCHES)) {
  
  TRY <- TRUE
  COUNT <- 0L
  while (TRY) {
    writeLines(text = ubs[INITIALLINE:FINALLINE],
               con = TEMP01)
    ENTREZQUERY <- paste0("epost -db biosample -input",
                          " ",
                          TEMP01,
                          " ",
                          "-format acc",
                          " | ",
                          "elink -target sra",
                          " | ",
                          "efetch -db sra -format runinfo -mode xml",
                          " | ",
                          'xtract -pattern Row -def "NA" -element',
                          " ",
                          "Run spots bases spots_with_mates avglength ",
                          "download_path Experiment LibraryStrategy LibrarySelection ",
                          "LibrarySource Platform Model SRAStudy BioProject ProjectID ",
                          "Sample BioSample SampleType TaxID ScientificName SampleName ",
                          "CenterName Submission Consent",
                          " > ",
                          TEMP02)
    x <- try(system(command = ENTREZQUERY,
                    intern = FALSE,
                    # ignore.stderr = TRUE,
                    # ignore.stdout = FALSE,
                    timeout = 7200))
    print(c(m1, x))
    
    # temp02 gets re-written by default so this should be fine
    x <- try(read.table(TEMP02,
                        sep = "\t",
                        comment.char = ""))
    COUNT <- COUNT + 1L
    # exit the while loop if the search does not return a try error
    # or exit if you've tried the query 5 times and it hasn't succeeded
    if (!is(object = x,
            class2 = "try-error") |
        COUNT > 5L) {
      print(COUNT)
      TRY <- FALSE
    }
  }
  RES[[m1]] <- x
  
  # edit search lines
  INITIALLINE <- FINALLINE + 1L
  FINALLINE <- FINALLINE + BATCHSIZE
  if (FINALLINE > TOTALLINES) {
    FINALLINE <- TOTALLINES
  }
  
}

SRAResults <- do.call(rbind,
                      RES[sapply(X = RES,
                                 FUN = function(x) {
                                   is(object = x,
                                      class2 = "data.frame")
                                 },
                                 simplify = TRUE)])

colnames(SRAResults) <- c("Run",
                          "spots",
                          "bases",
                          "Spots_with_mates",
                          "avglength",
                          "download_path",
                          "Experiment",
                          "LibraryStrategy",
                          "LibrarySelection",
                          "LibrarySource",
                          "Platform",
                          "Model",
                          "SRAStudy",
                          "BioProject",
                          "ProjectID",
                          "Sample",
                          "BioSample",
                          "SampleType",
                          "TaxID",
                          "ScientificName",
                          "SampleName",
                          "CenterName",
                          "Submission",
                          "Consent")

save(SRAResults,
     file = paste0(getwd(),
                   "/SearchResults2.RData"),
     compress = "xz")

