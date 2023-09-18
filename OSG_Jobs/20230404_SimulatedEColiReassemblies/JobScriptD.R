###### -- grab ANI between the source genome and reassemblies -----------------

suppressMessages(library(SynExtend))

ARGS <- commandArgs(trailingOnly = TRUE)


RefAssembly <- readDNAStringSet(filepath = "ref.fna.gz")
RefGC <- gffToDataFrame(GFF = "annot.gff",
                        Verbose = TRUE)

seqs1 <- readDNAStringSet(filepath = ARGS[2L])
gc2 <- gffToDataFrame(GFF = ARGS[3L],
                      Verbose = TRUE)

genes1 <- ExtractBy(x = RefGC[(RefGC$Type != "pseudogene") & RefGC$Coding, ],
                    y = RefAssembly)
writeXStringSet(x = genes1,
                filepath = "seqs1.fna")
genes2 <- ExtractBy(x = gc2[(gc2$Type != "pseudogene") & gc2$Coding, ],
                    y = seqs1)
writeXStringSet(x = genes2,
                filepath = "seqs2.fna")

# system(command = "chmod +x ANIcalculator_v1/*")

ANICall <- "ANIcalculator -genome1fna seqs1.fna -genome2fna seqs2.fna -outfile RES.txt -outdir res"

# paste0("ANIcalculator",
#        " -genome1fna ",
#        SEQ1,
#        " -genome2fna ",
#        SEQ2,
#        " -outfile ",
#        OUTFILE,
#        " -outdir ",
#        OUTDIR)

system(command = ANICall,
       intern = TRUE,
       ignore.stderr = TRUE)

x <- readLines("res/RES.txt")

save(x,
     file = paste0("ANIRes",
                   formatC(x = as.integer(ARGS[1]),
                           width = 6,
                           flag = 0,
                           format = "d"),
                   ".RData"),
     compress = "xz")

