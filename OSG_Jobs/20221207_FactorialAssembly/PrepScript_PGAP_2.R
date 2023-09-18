###### -- build jobmap and YAMLs for PGAP -------------------------------------
# match resulting fnas from planned jobs to lines from the original job map

library(SynExtend)
library(yaml)
library(stringr)

InitialMap1 <- "~/Repos/20221207_FactorialAssembly/JobMap_Assembly_2.txt"
InitialMap2 <- "~/Repos/20221207_FactorialAssembly/JobMap_Assembly_3.txt"

GFFLoc <- "~/Data/20230111_Factorial_Annotations/"
FNALoc <- "~/Data/20230111_Factorial_Assemblies/"

DestinationDir <- "~/Repos/20221207_FactorialAssembly/"
dir.create(path = "~/Repos/20221207_FactorialAssembly/YAML")

CompletedAssemblies <- list.files(FNALoc)
CompletedAnnotations <- list.files(GFFLoc)
InitialMap <- rbind(read.table(InitialMap1), read.table(InitialMap2))


# integer in for the job map
w1 <- as.integer(unlist(regmatches(x = CompletedAssemblies,
                                   m = gregexpr(pattern = "[0-9]+",
                                                text = CompletedAssemblies))))
# uncompressed fna characters for the YAML files
w2 <- unlist(regmatches(x = CompletedAssemblies,
                        m = gregexpr(pattern = "(.+)(?=\\.gz)",
                                     perl = TRUE,
                                     text = CompletedAssemblies)))
w3 <- unlist(regmatches(x = CompletedAssemblies,
                        m = gregexpr(pattern = "(?<=Assembly_)(.*)(?=\\.fna)",
                                     perl = TRUE,
                                     text = CompletedAssemblies)))

w4 <- unlist(regmatches(x = CompletedAnnotations,
                        m = gregexpr(pattern = "(?<=Annotation_)(.*)(?=\\.gff)",
                                     text = CompletedAnnotations,
                                     perl = TRUE)))

# slightly different for this round
ActualAssemblies <- rep(FALSE,
                        length(CompletedAssemblies))

pBar <- txtProgressBar(style = 1L)
PBAR <- length(ActualAssemblies)
for (m1 in seq_along(CompletedAssemblies)) {
  x <- readDNAStringSet(filepath = paste0(FNALoc,
                                          "/",
                                          CompletedAssemblies[m1]))
  
  # is the assemblies' GFF in the completed folder
  z <- w3[m1] %in% w4
  if (length(x) > 0 & !z) {
    ActualAssemblies[m1] <- TRUE
    w5 <- InitialMap$V6[which(InitialMap$V5 == w1[m1])]
    w5 <- str_replace(string = w5,
                      pattern = "_",
                      replacement = " ")
    Y1 <- as.yaml(list("fasta" = list("class" = "File",
                                      "location" = paste0("/srv/",
                                                          w2[m1])),
                       "submol" = list("class" = "File",
                                       "location" = paste0("/srv/",
                                                           "submol_",
                                                           w3[m1],
                                                           ".yaml")),
                       "supplemental_data" = "{ class: Directory, location: /srv/input-2022-04-14.build6021 }",
                       "report_usage" = "true",
                       "ignore_all_errors" = "true"))
    Y1 <- str_replace_all(string = Y1,
                          pattern = "'",
                          replacement = "")
    cat(Y1,
        file = paste0("~/Repos/20221207_FactorialAssembly/YAML/Controller_",
                      w3[m1],
                      ".yaml"))

    Y2 <- as.yaml(list("topology" = "linear",
                       "organism" = list("genus_species" = w5),
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
        file = paste0("~/Repos/20221207_FactorialAssembly/YAML/submol_",
                      w3[m1],
                      ".yaml"))
  } # else do nothing
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}

# int to pin job
# controller YAML
# submol YAML
# annotation name
# transferred fasta
# pgap input fasta
JobMap <- data.frame("JobPin" = seq(sum(ActualAssemblies)),
                     "controller" = paste0("Controller_",
                                           w3[ActualAssemblies],
                                           ".yaml"),
                     "submol" = paste0("submol_",
                                       w3[ActualAssemblies],
                                       ".yaml"),
                     "annotation_name" = paste0("Annotation_",
                                                w3[ActualAssemblies],
                                                ".gff"),
                     "compressed" = CompletedAssemblies[ActualAssemblies],
                     "input" = w2[ActualAssemblies])

write.table(x = JobMap,
            file = paste0(DestinationDir,
                          "JobMap_PGAP_2.txt"),
            append = FALSE,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

set.seed(1986)
TestMap <- JobMap[sample(x = nrow(JobMap),
                         size = 10L,
                         replace = FALSE), ]

write.table(x = TestMap,
            file = paste0(DestinationDir,
                          "TestMap.txt"),
            append = FALSE,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

