###### -- consolidate data across assembly and SRA ----------------------------

suppressMessages(library(SynExtend))

query1 <- paste("esearch",
                "-db assembly",
                "-query '",
                '("Bacteria"[Organism] OR "Archaea"[Organism])',
                # '"Streptomyces"[Organism]',
                # '"Archaea"[Organism]',
                'AND "latest genbank"[properties]',
                "'",
                "|",
                "esummary",
                "|",
                "xtract -pattern DocumentSummary",
                "-element",
                'FtpPath_GenBank',
                'BioSampleAccn',
                'AssemblyStatus',
                'SubmitterOrganization',
                'SubmissionDate',
                'Organism',
                'Taxid',
                'SpeciesName',
                'ContigN50',
                'ScaffoldN50',
                'Coverage',
                '-block Stat -if "@category" -equals total_length -element Stat')

SEARCHTIMESTART <- Sys.time()
replies01 <- system(command = query1,
                    intern = TRUE,
                    timeout = 60000L)
SEARCHTIMEEND <- Sys.time()
GENBANKSEARCHTIME <- SEARCHTIMEEND - SEARCHTIMESTART
print(GENBANKSEARCHTIME)
# ~2 hrs

x01 <- strsplit(x = replies01,
                split = "\t",
                fixed = TRUE)
x02 <- x01[lengths(x01) == max(lengths(x01))]

x03 <- do.call(rbind,
               x02)

GenBankAccessions <- data.frame("FTP" = x03[, 1L],
                                "Biosample" = x03[, 2L],
                                "Assembly_Status" = x03[, 3L],
                                "Submitter_Org" = x03[, 4L],
                                "Submission_Date" = x03[, 5L],
                                "Organism" = x03[, 6L],
                                "TaxID" = x03[, 7L],
                                "SpeciesName" = x03[, 8L],
                                "ContigN50" = as.integer(x03[, 9L]),
                                "ScaffoldN50" = as.integer(x03[, 10L]),
                                "Coverage" = as.numeric(gsub(pattern = "([^0-9]*)([0-9.]+)([^0-9]*)",
                                                             replacement = "\\2",
                                                             x = x03[, 11L])),
                                "Total_Length" = as.integer(x03[, 12L]))

save(GenBankAccessions,
     file = "~/Repos/PseudogenePlots/InputData/GenBank_Assembly_Summary.RData",
     compress = "xz")

# search the biosample db for sra ids
# w1 <- GenBankAccessions$Biosample[1:1000]

# target ~ 1000 ids to search at a time to live within epost's confines
w1 <- split(x = GenBankAccessions$Biosample,
            f = 1:1600)
w2 <- vector(mode = "list",
             length = length(w1))
for (m1 in seq_along(w1)) {
  TEMP01 <- tempfile()
  writeLines(text = w1[[m1]],
             con = TEMP01)
  
  # search the biosample DB to return identifier sets
  # extract SRA identifiers
  IDQUERY <- paste0("epost -db biosample -input",
                    " ",
                    TEMP01,
                    " ",
                    "-format acc",
                    " | ",
                    "esummary",
                    " | ",
                    "xtract -pattern DocumentSummary",
                    " ",
                    "-element Identifiers")
  
  SRATIMESTART <- Sys.time()
  replies02 <- system(command = IDQUERY,
                      intern = TRUE,
                      timeout = 6000L)
  SRATIMEEND <- Sys.time()
  SRATOTALTIME <- SRATIMEEND - SRATIMESTART
  print(SRATOTALTIME)
  
  w2[[m1]] <- replies02
  print(m1)
}

w3 <- unlist(w2)

# library(stringr)
# table(str_count(string = w3, pattern = "SRA:"))
# 
# 0       1       2       4 
# 348514 1277251     619       1 
# ~1.27m have a single associated SRA ID
# more genbank assemblies lack an associated SRA ID than there are total RefSeq assemblies -proks-

SRAIDs <- unlist(regmatches(x = w3,
                            m = gregexpr(pattern = "(?<=SRA: )([^ ;]+)",
                                         text = w3,
                                         perl = TRUE)))
save(w3,
     SRAIDs,
     file = "~/Repos/PseudogenePlots/InputData/SRAIDList.RData",
     compress = "xz")

# grab the SRA metadata based on the SRA IDs

# target about 500 ids to search at a time to live within epost's confines
# the more ids we try and search for at a single time, the more server errors
# we encounter
w4 <- split(x = SRAIDs,
            f = 1:2600)
w5 <- vector(mode = "list",
             length = length(w4))

for (m1 in seq_along(w5)) {
  TEMP02 <- tempfile()
  writeLines(text = w4[[m1]],
             con = TEMP02)
  
  SRAQUERY <- paste0("epost -db sra -input",
                     " ",
                     TEMP02,
                     " ",
                     "-format acc",
                     " | ",
                     "esummary -format runinfo -mode xml",
                     " | ",
                     "xtract -pattern Row",
                     " ",
                     '-def "NA"',
                     " ",
                     "-element ",
                     "Run spots bases spots_with_mates avglength ",
                     "download_path Experiment LibraryStrategy LibrarySelection ",
                     "LibrarySource Platform Model SRAStudy BioProject ProjectID ",
                     "Sample BioSample SampleType TaxID ScientificName SampleName ",
                     "CenterName Submission Consent")
  
  SRATIMESTART <- Sys.time()
  replies03 <- system(command = SRAQUERY,
                      intern = TRUE,
                      timeout = 6000L)
  SRATIMEEND <- Sys.time()
  SRATOTALTIME <- SRATIMEEND - SRATIMESTART
  print(SRATOTALTIME)
  w5[[m1]] <- replies03
  
  print(m1)
  unlink(TEMP02)
}

save(w5,
     file = "~/tempsrares.RData",
     compress = "xz")
w6 <- unlist(w5)
SRA_Meta <- strsplit(x = w6,
                     split = "\t",
                     fixed = TRUE)
SRA_Meta <- SRA_Meta[lengths(SRA_Meta) == 24]
SRA_Meta <- do.call(rbind,
                    SRA_Meta)

SRA_Meta <- data.frame("Run" = SRA_Meta[, 1],
                       "spots" = as.numeric(SRA_Meta[, 2]),
                       "bases" = as.numeric(SRA_Meta[, 3]),
                       "Spots_with_mates" = as.numeric(SRA_Meta[, 4]),
                       "avglength" = as.numeric(SRA_Meta[, 5]),
                       "download_path" = SRA_Meta[, 6],
                       "Experiment" = SRA_Meta[, 7],
                       "LibraryStrategy" = SRA_Meta[, 8],
                       "LibrarySelection" = SRA_Meta[, 9],
                       "LibrarySource" = SRA_Meta[, 10],
                       "Platform" = SRA_Meta[, 11],
                       "Model" = SRA_Meta[, 12],
                       "SRAStudy" = SRA_Meta[, 13],
                       "BioProject" = SRA_Meta[, 14],
                       "ProjectID" = as.integer(SRA_Meta[, 15]),
                       "Sample" = SRA_Meta[, 16],
                       "BioSample" = SRA_Meta[, 17],
                       "SampleType" = SRA_Meta[, 18],
                       "TaxID" = as.integer(SRA_Meta[, 19]),
                       "ScientificName" = SRA_Meta[, 20],
                       "SampleName" = SRA_Meta[, 21],
                       "CenterName" = SRA_Meta[, 22],
                       "Submission" = SRA_Meta[, 23],
                       "Consent" = SRA_Meta[, 24])

save(SRA_Meta,
     file = "~/Repos/PseudogenePlots/InputData/GenBank_SRA_Summary.RData",
     compress = "xz")


