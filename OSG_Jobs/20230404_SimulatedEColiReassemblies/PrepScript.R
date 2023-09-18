###### -- build a jobmap for read simulations ---------------------------------

suppressMessages(library(SynExtend))
suppressMessages(library(yaml))
suppressMessages(library(stringr))

# pbsim can turn an accuracy nob, and an error distribution nob, but we'll only
# turn the accuracy nob -- we still need to set error distribution for PacBio though

s1 <- expand.grid(list(c(-10, -5, -2, 0, 2, 5, 10),
                       1:3,
                       c(5, 10, 50, 100, 500, 1000),
                       c(100, 150),
                       c("HS25")),
                  stringsAsFactors = FALSE)
s2 <- expand.grid(list(c(-10, -5, -2, 0, 2, 5, 10),
                       1:3,
                       c(5, 10, 50, 100, 500, 1000),
                       c(100, 150),
                       c("HSXt")),
                  stringsAsFactors = FALSE)
s3 <- expand.grid(list(c(-10, -5, -2, 0, 2, 5, 10),
                       1:3,
                       c(5, 10, 50, 100, 500, 1000),
                       c(100, 150, 200, 250),
                       c("MSv3")),
                  stringsAsFactors = FALSE)

JobMap <- do.call(rbind,
                  list(s1, s2, s3))
JobMap <- JobMap[, -2L]
colnames(JobMap) <- c("Qual", "Depth", "Length", "Technology")
JobMap <- cbind("PersistentID" = seq(nrow(JobMap)),
                JobMap)

write.table(x = JobMap,
            file = "~/Repos/20230404_SimulatedEColiReassemblies/JobMapA.txt",
            quote = FALSE,
            append = FALSE,
            row.names = FALSE,
            col.names = FALSE)

TestMap <- JobMap[sample(x = seq(nrow(JobMap)),
                         size = 10,
                         replace = FALSE), ]
write.table(x = TestMap,
            file = "~/Repos/20230404_SimulatedEColiReassemblies/TestMap.txt",
            quote = FALSE,
            append = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# only a single YAML controller and submol file needed:

Y1 <- as.yaml(list("fasta" = list("class" = "File",
                                  "location" = "assembly.fna"),
                   "submol" = list("class" = "File",
                                   "location" = "/srv/submol.yaml"),
                   "supplemental_data" = "{ class: Directory, location: /srv/input-2022-04-14.build6021 }",
                   "report_usage" = "true",
                   "ignore_all_errors" = "true"))
Y1 <- str_replace_all(string = Y1,
                      pattern = "'",
                      replacement = "")

cat(Y1,
    file = "~/Repos/20230404_SimulatedEColiReassemblies/Controller.yaml")

Y2 <- as.yaml(list("topology" = "linear",
                   "organism" = list("genus_species" = "Escherichia coli"),
                   "authors" = list(list("author" = list("last_name" = "Cooley",
                                                         "first_name" = "Nicholas",
                                                         "middle_initial" = "P"))),
                   "contact_info" = list("first_name" = "Nicholas",
                                         "last_name" = "Cooley",
                                         "email" = "npc19@pitt.edu",
                                         "organization" = "University of Pittsburgh",
                                         "department" = "Department of Biomedical Informatics",
                                         "phone" = "412-624-5100",
                                         "street" = "100 Technology Dr.",
                                         "city" = "Pittsburgh",
                                         "postal_code" = "15219",
                                         "country" = "USA")))
Y2 <- gsub(pattern = "(?<=: )([^'\n]+)",
           replacement = "'\\1'",
           x = Y2,
           perl = TRUE)

cat(Y2,
    file = "~/Repos/20230404_SimulatedEColiReassemblies/submol.yaml")



# # e coli k12 reference genome
# REF1 <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
# # camphylobacter jejuni
# # REF1 <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/085/GCF_000009085.1_ASM908v1/GCF_000009085.1_ASM908v1_genomic.fna.gz"
# REF2 <- strsplit(x = REF1,
#                  split = "/",
#                  fixed = TRUE)[[1]][11]
# REF3 <- gsub(pattern = "\\.gz",
#              replacement = "",
#              x = REF2)
# 
# # grab the genome and unzip it
# system(command = paste("wget",
#                        REF1,
#                        "&& gunzip",
#                        REF2))

# local construction
# DNA <- "ASSEMBLY01/assembly.fasta"
# pBar <- txtProgressBar(style = 1L)
# PBAR <- nrow(JobMap)
# for (m1 in seq_len(PBAR)) {
#   TECHNOLOGY <- JobMap[m1, 4L]
#   DEPTH <- JobMap[m1, 3L]
#   QUAL <- JobMap[m1, 2L]
#   if (TECHNOLOGY == "MSv3") {
#     # paired end
#     x <- system(paste("art_illumina",
#                       "-ss MSv3",
#                       "-i",
#                       REF3,
#                       "-l 150",
#                       "-m 300",
#                       "-s 10",
#                       "-o PAIRED_v01",
#                       "-f",
#                       DEPTH,
#                       "-qs",
#                       QUAL,
#                       "-qs2",
#                       QUAL),
#                 intern = TRUE,
#                 ignore.stderr = TRUE)
#     y <- system(paste("unicycler",
#                       "-1 PAIRED_v011.fq",
#                       "-2 PAIRED_v012.fq",
#                       "-o ASSEMBLY01",
#                       "--min_fasta_length 200"),
#                 intern = TRUE,
#                 ignore.stderr = TRUE)
#     if (file.exists(DNA)) {
#       dna <- readDNAStringSet(filepath = DNA)
#       writeXStringSet(x = dna,
#                       filepath = paste0("/UserData/Assembly",
#                                         formatC(m1,
#                                                 width = 4,
#                                                 format = "d",
#                                                 flag = "0"),
#                                         ".fna"))
#     } else {
#       stop(m1)
#     }
#     system(command = "rm -rf ASSEMBLY01 PAIRED_*")
#   } else if (TECHNOLOGY == "ONT") {
#     # ont
#     x <- system(command = paste("pbsim",
#                                 "--strategy wgs",
#                                 "--method qshmm",
#                                 "--qshmm /pbsim3/data/QSHMM-ONT.model",
#                                 "--depth",
#                                 DEPTH,
#                                 "--accuracy-mean",
#                                 QUAL,
#                                 "--genome",
#                                 REF3),
#                 intern = TRUE,
#                 ignore.stderr = TRUE)
#     # reassemble with unicycler
#     y <- system(paste("unicycler",
#                       "-l sd_0001.fastq",
#                       "-o ASSEMBLY01",
#                       "--min_fasta_length 200"),
#                 intern = TRUE,
#                 ignore.stdout = TRUE)
#     if (file.exists(DNA)) {
#       dna <- readDNAStringSet(filepath = DNA)
#       writeXStringSet(x = dna,
#                       filepath = paste0("/UserData/Assembly",
#                                         formatC(m1,
#                                                 width = 4,
#                                                 format = "d",
#                                                 flag = "0"),
#                                         ".fna"))
#     } else {
#       stop(m1)
#     }
#   } else if (TECHNOLOGY == "PB") {
#     # pb
#     x <- system(command = paste("pbsim",
#                                 "--strategy wgs",
#                                 "--method qshmm",
#                                 "--qshmm /pbsim3/data/QSHMM-RSII.model",
#                                 "--depth",
#                                 DEPTH,
#                                 "--accuracy-mean",
#                                 QUAL,
#                                 "--seed",
#                                 ID,
#                                 "--genome",
#                                 REF3),
#                 intern = TRUE,
#                 ignore.stderr = TRUE)
#     # reassemble with unicycler
#     y <- system(paste("unicycler",
#                       "-l sd_0001.fastq",
#                       "-o ASSEMBLY01",
#                       "--min_fasta_length 200"),
#                 intern = TRUE,
#                 ignore.stderr = TRUE)
#     if (file.exists(DNA)) {
#       dna <- readDNAStringSet(filepath = DNA)
#       writeXStringSet(x = dna,
#                       filepath = paste0("/UserData/Assembly",
#                                         formatC(m1,
#                                                 width = 4,
#                                                 format = "d",
#                                                 flag = "0"),
#                                         ".fna"))
#     } else {
#       stop(m1)
#     }
#     system(command = "rm -rf ASSEMBLY01 sd_*")
#   }
#   setTxtProgressBar(pb = pBar,
#                     value = m1 / PBAR)
# }






