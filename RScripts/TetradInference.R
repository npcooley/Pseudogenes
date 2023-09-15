
library(SynExtend)
library(rjson)
library(igraph)

ColVec1 <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000')
ColVec2 <- paste0(ColVec1,
                  "33")

InputData <- "~/Repos/Pseudogenes/SearchResults.RData"
CCDDir <- "~/Downloads/causal-cmd-1.3.0/"

load(file = InputData,
     verbose = TRUE)
# load(file = "~/Repos/Pseudogenes/SearchResults2.RData",
#      verbose = TRUE)


# wrangle genus names
g1 <- unlist(regmatches(x = EntrezResults$SpeciesName,
                        m = gregexpr(pattern = "^[^ ]+",
                                     text = EntrezResults$SpeciesName)))
g2 <- grep(pattern = "candidatus",
           ignore.case = TRUE,
           x = g1)
g1[g2] <- unlist(regmatches(x = EntrezResults$SpeciesName[g2],
                            m = gregexpr(pattern = "(?<=Candidatus )([^ ]+)",
                                         text = EntrezResults$SpeciesName[g2],
                                         perl = TRUE)))
g2 <- grep(pattern = "[^A-Za-z]",
           x = g1)
g1[g2] <- unlist(regmatches(x = g1[g2],
                            m = gregexpr(pattern = "[A-Za-z]+",
                                         text = g1[g2])))

# encode year
y1 <- EntrezResults$Submission_Date
y1 <- subseq(x = y1, start = 1, end = 4)
y1 <- as.integer(y1)

# encode assemblers somehow
a1 <- apply(X = AssemblersLogical,
            MARGIN = 1,
            FUN = function(x) {
              which(x)
            })
a1 <- unname(sapply(X = a1,
                    FUN = function(x) {
                      paste0(names(x),
                             collapse = "::")
                    },
                    simplify = TRUE,
                    USE.NAMES = FALSE))
a2 <- unique(a1)
a3 <- match(x = a1,
            table = a2)

# encode technology
t1 <- apply(X = TechnologyLogical,
            MARGIN = 1,
            FUN = function(x) {
              which(x)
            })
t1 <- unname(sapply(X = t1,
                    FUN = function(x) {
                      paste0(names(x),
                             collapse = "::")
                    },
                    simplify = TRUE,
                    USE.NAMES = FALSE))

t2 <- unique(t1)
t3 <- match(x = t1,
            table = t2)
j1 <- unique(EntrezResults$Assembly_Status)
j2 <- match(x = EntrezResults$Assembly_Status,
            table = j1)

# pseudos by mb
frameshiftperkb <- PseudosCondensed[, 1L] / (EntrezResults$Total_Length / 1000000)
internalstopperkb <- PseudosCondensed[, 4L] / (EntrezResults$Total_Length / 1000000)
partialperkb <- PseudosCondensed[, 3L] / (EntrezResults$Total_Length / 1000000)

tetradtable_v01 <- data.frame("genus" = g1,
                              "cov" = EntrezResults$Coverage,
                              "submit_year" = y1,
                              "assembler" = a1,
                              "technology" = t1,
                              "frameshifts" = frameshiftperkb,
                              "stops" = internalstopperkb,
                              "partials" = partialperkb,
                              "assembly_status" = j2,
                              "fragment" = EntrezResults$ContigN50 / EntrezResults$Total_Length)

tetradtable_v02 <- tetradtable_v01[tetradtable_v01$assembler != "Unknown_Other_Uncaptured" |
                                     tetradtable_v01$technology != "Unknown_Other_Uncaptured", ]

w1 <- sort(table(tetradtable_v02$genus), decreasing = TRUE)
w2 <- sort(table(tetradtable_v02$submit_year), decreasing = TRUE)
w3 <- sort(table(tetradtable_v02$assembler), decreasing = TRUE)
w4 <- sort(table(tetradtable_v02$technology), decreasing = TRUE)

w5 <- tetradtable_v02$genus %in% names(w1[w1 >= 1000])
w6 <- tetradtable_v02$submit_year %in% names(w2[w2 >= 1000])
w7 <- tetradtable_v02$assembler %in% names(w3[w3 >= 1000])
w8 <- tetradtable_v02$technology %in% names(w4[w4 >= 1000])
w9 <- tetradtable_v02$cov > 0

# final table, drop low count data, drop unsubmitted metadata
tetradtable_v03 <- tetradtable_v02[w5 & w6 & w7 & w8 & w9, ]

tetradtable_v04 <- data.frame("genus" = match(x = tetradtable_v03$genus,
                                              table = unique(tetradtable_v03$genus)),
                              "cov" = tetradtable_v03$cov,
                              "submit_year" = match(x = tetradtable_v03$submit_year,
                                                    table = unique(tetradtable_v03$submit_year)),
                              "assembler" = match(x = tetradtable_v03$assembler,
                                                  table = unique(tetradtable_v03$assembler)),
                              "technology" = match(x = tetradtable_v03$technology,
                                                   table = unique(tetradtable_v03$technology)),
                              "frameshifts" = tetradtable_v03$frameshifts,
                              "stops" = tetradtable_v03$stops,
                              "partials" = tetradtable_v03$partials,
                              "assembly_status" = tetradtable_v03$assembly_status,
                              "fragment" = tetradtable_v03$fragment)


# run tetrad

setwd(CCDDir)

write.table(x = tetradtable_v04,
            file = "df01.txt",
            quote = FALSE,
            row.names = FALSE,
            append = FALSE,
            sep = "\t")

# generate metadata file to clarify discrete and continuous variables:

tetrad_meta <- list("domains" = list(list("name" = "genus",
                                          "discrete" = TRUE),
                                     list("name" = "cov",
                                          "discrete" = FALSE),
                                     list("name" = "submit_year",
                                          "discrete" = TRUE),
                                     list("name" = "assembler",
                                          "discrete" = TRUE),
                                     list("name" = "technology",
                                          "discrete" = TRUE),
                                     list("name" = "frameshifts",
                                          "discrete" = FALSE),
                                     list("name" = "stops",
                                          "discrete" = FALSE),
                                     list("name" = "partials",
                                          "discrete" = FALSE),
                                     list("name" = "fragment",
                                          "discrete" = FALSE),
                                     list("name" = "assembly_status",
                                          "discrete" = TRUE)))

write(toJSON(tetrad_meta), "metadata_v01.txt")

# generate knowledge file for tiers and forbidden associations

tetrad_knowledge <- c("/knowledge",
                      "",
                      "addtemporal",
                      "",
                      "forbiddirect",
                      "genus assembler", "genus submit_year", "genus technology",
                      "assembler technology", "assembler submit_year", "assembler cov", "assembler genus",
                      "technology genus", "technology submit_year",
                      "submit_year partials", "submit_year frameshifts", "submit_year stops", "submit_year cov", "submit_year assembly_status", "submit_year fragment", 
                      "partials genus", "partials cov", "partials assembler", "partials technology", "partials submit_year", "partials stops", "partials frameshifts", "partials assembly_status", "partials fragment",
                      "frameshifts genus", "frameshifts cov", "frameshifts assembler", "frameshifts technology", "frameshifts submit_year", "frameshifts stops", "frameshifts partials", "frameshifts assembly_status", "frameshifts fragment",
                      "stops partials", "stops frameshifts", "stops genus", "stops cov", "stops assembler", "stops technology", "stops submit_year", "stops assembly_status", "stops fragment",
                      "cov genus", "cov submit_year", "cov technology", "cov assembler",
                      "fragment genus", "fragment cov", "fragment submit_year", "fragment assembler", "fragment technology", "fragment assembly_status",
                      "assembly_status fragment", "assembly_status submit_year", "assembly_status cov", "assembly_status genus", "assembly_status assembler", "assembly_status technology",
                      "",
                      "requiredirect",
                      "")

writeLines(text = tetrad_knowledge,
           con = "knowledge_v01.txt")

# java -jar ./causal-cmd-1.3.0-jar-with-dependencies.jar --algorithm fci --data-type discrete --dataset /Users/sixfoot10/Desktop/Sequencing/Assembly/other/CausalModeling/pseudos_tetrad_v04.txt --delimiter tab --skip-latest --metadata /Users/sixfoot10/Desktop/Sequencing/Assembly/other/CausalModeling/metadata_v1.txt --knowledge /Users/sixfoot10/Desktop/Sequencing/Assembly/other/CausalModeling/knowledge_v1.txt --test cg-lr-test --numberResampling 100
tetrad_command <- paste0("java -jar ./causal-cmd-1.3.0-jar-with-dependencies.jar",
                         " --algorithm fci",
                         " --data-type discrete",
                         " --dataset ", paste0(getwd(), "/df01.txt"),
                         " --delimiter tab",
                         " --skip-latest",
                         " --metadata ", paste0(getwd(), "/metadata_v01.txt"),
                         " --knowledge ", paste0(getwd(), "/knowledge_v01.txt"),
                         " --test cg-lr-test",
                         " --numberResampling 100")

system(command = tetrad_command)

# A --> B
# present
# A is a cause of B. It may be a direct or indirect cause that may include other measured variables. Also, there may be an unmeasured confounder of A and B.
# absent
# B is not a cause of A.
# A <-> B
# present
# There is an unmeasured variable (call it L) that is a cause of A and B. There may be measured variables along the causal pathway from L to A or from L to B.
# absent
# A is not a cause of B. B is not a cause of A.
# A o-> B
# present
# Either A is a cause of B, or there is an unmeasured variable that is a cause of A and B, or both.
# absent
# B is not a cause of A.
# A o-o B
# Exactly one of the following holds: (a) A is a cause of B, or (b) B is a cause of A, or (c) there is an unmeasured variable that is a cause of A and B, or (d) both a and c, or (e) both b and c.

tetresult <- list.files()
e0 <- grep(pattern = "^fci",
           x = tetresult)
tetresult <- tetresult[e0[length(e0)]]
tetresult <- readLines(tetresult)
# tetrad edges
e1 <- grepl(pattern = "^[0-9]+\\.",
            x = tetresult)
edges1 <- tetresult[e1]
# extract edges, edge types and weights
edges2 <- unlist(regmatches(x = edges1,
                            m = gregexpr(pattern = "(?<= )([a-z_]+ [->o<]{3} [a-z_]+)(?= )",
                                         text = edges1,
                                         perl = TRUE)))
edges3 <- regmatches(x = edges1,
                     m = gregexpr(pattern = "(?<=:)([^;]+)(?=;)",
                                  text = edges1,
                                  perl = TRUE))
edges3 <- as.numeric(sapply(X = edges3,
                            FUN = function(x) {
                              x[1]
                            }))
edges4 <- strsplit(x = edges2,
                   split = " ",
                   fixed = TRUE)
df <- data.frame("p1" = sapply(X = edges4,
                               FUN = function(x) {
                                 x[1]
                               }),
                 "p2" = sapply(X = edges4,
                               FUN = function(x) {
                                 x[3]
                               }),
                 "weight" = edges3,
                 "type" = sapply(X = edges4,
                                 FUN = function(x) {
                                   x[2]
                                 }))
df$p1 <- gsub(pattern = "_",
              replacement = " ",
              x = df$p1)
df$p2 <- gsub(pattern = "_",
              replacement = " ",
              x = df$p2)

g <- graph_from_data_frame(d = df,
                           directed = TRUE)
E(g)$weight <- df$weight
# l <- layout_as_tree(g)
# set.seed(1986)
# plot(g,
#      layout=l,
#      edge.arrow.size = 0.1)
# 
# set.seed(1986)
# plot(g,
#      layout = layout.fruchterman.reingold,
#      edge.arrow.size = 0.1)


# available nodes:
avl_nodes <- colnames(tetradtable_v04)
avl_nodes <- gsub(x = avl_nodes,
                  pattern = "_",
                  replacement = " ")
# submitter controlled: 4
# sequencing/assembly outcome, non-pseudo related: 5
# pseudos: 6
node_cols <- ColVec1[c(4,5,4,4,4,6,6,6,5,5)]
# node_cols <- ColVec1[c(1,2,1,1,1,3,3,3,2,2)]
pres_nodes <- attr(V(g), "names")

ColVec3 <- ColVec1[match(x = df$type,
                         table = unique(df$type))]
weights2 <- df$weight * 100
weights2[weights2 == 100] <- 99
ColVec3 <- paste0(ColVec3,
                  formatC(x = weights2,
                          width = 2,
                          flag = 0,
                          format = "d"))

# l <- layout_with_dh(g)
pdf(file = "test.pdf",
    width = 3.5,
    height = 3.5)
l <- layout_with_sugiyama(g)
set.seed(1986)
par(mgp = c(2,1,0),
    mar = c(2,2,1,1))
plot(g,
     layout=l,
     # vertex.shape = "rectangle",
     # vertex.size = (strwidth(attr(V(g), "names")) + strwidth("oo")) * 25,
     # vertex.size2 = strheight("I") * 2 * 15,
     vertex.shape = "circle",
     # vertex.size = max((strwidth(attr(V(g), "names")) + strwidth("o")) * 20),
     # vertex.size = seq_along(attr(V(g), "names")) * 5,
     vertex.size = 40,
     vertex.label.cex = 0.3,
     vertex.color = node_cols[match(x = pres_nodes,
                                    table = avl_nodes)],
     # edge.arrow.size = 1,
     edge.color = ColVec3,
     # edge.width = 2,
     vertex.frame.width = 0,
     vertex.label.color = "black")
legend("bottomleft",
       # legend = unique(df$type),
       legend = c("causal", "causal*"),
       lty = 1,
       lwd = 2,
       col = ColVec1[seq_along(unique(df$type))],
       cex = 0.5)
dev.off()

save(df,
     tetradtable_v04,
     file = "~/Repos/PseudogenePlots/InputData/TetradEdges.RData",
     compress = "xz")

