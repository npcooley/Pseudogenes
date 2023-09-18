###### -- find missing jobs ---------------------------------------------------

suppressMessages(library(SynExtend))


ColVec1 <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000')
ColVec2 <- paste0(ColVec1,
                  "33")
ColVec3 <- paste0(ColVec1,
                  "15")

# expects arguments in the format of 
# folder location of set comparison results
# file location of job map
# file location of original entrez search results
# file location for the output
ARGS <- commandArgs(trailingOnly = TRUE)

if (length(ARGS) == 0L) {
  stop("when run from the terminal this script expects trailing arguments")
}

ARGS <- c("~/Data/20230902_SetComparison",
          "~/Repos/20230814_SetComparison/JobMap01.txt",
          "~/Repos/20230814_SetComparison/SearchResults.RData",
          "~/Repos/20230814_SetComparison/JobMap03.txt")

load(file = ARGS[3L],
     verbose = TRUE)
jobmap <- read.table(file = ARGS[2L])

path01 <- ARGS[1L]
files01 <- list.files(path = path01)
completedjobs <- as.integer(unlist(regmatches(x = files01,
                                              m = gregexpr(pattern = "[0-9]+",
                                                           text = files01))))
# jobs don't always complete neatly on the grid, accountting for that,
# build out a table of the newer refseq gene counts for the genomes that are present in the completed set
completedjobmap <- jobmap[completedjobs, ]
missingjobs <- jobmap[!(seq(nrow(jobmap)) %in% completedjobs), ]

write.table(x = missingjobs,
            file = ARGS[4],
            quote = FALSE,
            append = FALSE,
            row.names = FALSE,
            col.names = FALSE)


