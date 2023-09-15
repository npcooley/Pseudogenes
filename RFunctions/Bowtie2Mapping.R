###### -- Bowtie2 pipeline ----------------------------------------------------

# mimic Breseq initially
Bowtie2Pipeline <- function(Reference,
                            Reads) {
  # write out assembly to map
  tmp01 <- tempfile()
  writeXStringSet(x = Reference,
                  filepath = tmp01,
                  append = FALSE,
                  format = "fasta")
  # print(paste0("Assembly written as: ",
  #              tmp01))
  # write out reads
  tmp02 <- tempfile()
  writeQualityScaledXStringSet(x = Reads,
                               filepath = tmp02,
                               append = FALSE)
  # print(paste0("Reads written as: ",
  #              tmp02))
  # create an index with samtools? not sure why?
  # breseq does this so shrug?
  FAIDX <- paste0("samtools faidx ",
                  tmp01)
  system(command = FAIDX,
         intern = FALSE)
  
  # create a directory for the bowtie index
  tmp05 <- tempfile()
  BOWTIE_INDEX <- paste0("bowtie2-build -q ",
                         '"',
                         tmp01,
                         '"',
                         " ",
                         '"',
                         tmp05,
                         '"')
  # print(paste0("Assembly indexed as: ",
  #              tmp01))
  system(command = BOWTIE_INDEX,
         intern = FALSE)
  
  # output for sam file
  tmp03 <- tempfile()
  # print(paste0("Samfile pointed to: ",
  #              tmp03,
  #              ".sam"))
  # run bowtie 2
  # currently hard coded with breseq's initial defaults
  BOWTIE2 <- paste0("bowtie2",
                    " ",
                    "-t", # print time
                    " ",
                    "--no-unal", # suppress records for unaligned sequences
                    " ",
                    "-p 12", # processors
                    " ",
                    "-L 31", # length of seed substrings 31 < & > 3
                    # breseq's second less stringent call is:
                    # "-L 23"
                    " ",
                    "--ma 1", # match bonus
                    " ",
                    "--mp 3", # max mismatch penalty
                    " ",
                    "--np 0", # penalty for ambiguity characters
                    " ",
                    "--rdg 2,3", # read gap open
                    " ",
                    "--rfg 2,3", # reference gap open
                    " ",
                    "--ignore-quals",
                    " ",
                    "--local", # short read alignments
                    " ",
                    "-i S,1,0.25", # interval between seed substrings
                    # per Bowtie default is 'S,1,1.15'
                    " ",
                    "--score-min L,1,0.9", # minimum acceptable alignment score
                    # per Bowtie
                    # '--score-min G,20,8' is default for local mode
                    # '--score-min L,-0.6,0.6' is default for end-to-end
                    # Breseq's second less stringent call is:
                    # "--score-min L,6,0.2"
                    " ",
                    "-k 2000", # number of reports per READ
                    " ",
                    "--reorder", # force sam output to match order of input reads
                    " ",
                    "-x",
                    " ",
                    tmp05,
                    " ",
                    "-U",
                    " ",
                    tmp02,
                    " ",
                    "-S",
                    " ",
                    tmp03)
  
  system(command = BOWTIE2,
         intern = FALSE)
  # breseq sends unmatched sequences to a second command with less stringent match
  # requirements and a smaller substring seed
  # likely we can skip this for now
  
  # output for BAM file
  tmp04 <- tempfile()
  # print(paste0("Bamfile pointed to: ",
  #              tmp04))
  # run samtools to convert the samfile to a bamfile
  SAMTOOLS <- paste0("samtools view -S -b ",
                     tmp03,
                     " > ",
                     tmp04)
  system(command = SAMTOOLS,
         intern = FALSE)
  
  Res <- BamFile(tmp04)
  # return an r-object that can be targetted by scanBam()
  return(Res)
}



