###### -- map short reads to a reference with minimap2 ------------------------

# given sequences that are in R currently!

Minimap2Pipeline <- function(Reference,
                             Reads,
                             type = "sr") {
  # initialize the tempfiles we need
  TMP1 <- tempfile()
  TMP2 <- tempfile()
  TMP3 <- tempfile()
  TMP4 <- tempfile()
  # write out reference and reads to tempfiles
  writeXStringSet(x = Reference,
                  filepath = TMP1,
                  append = FALSE)
  writeXStringSet(x = Reads,
                  filepath = TMP2,
                  append = FALSE)
  if (type == "sr") {
    # run minimap2 to map short reads
    MINIMAP2 <- paste0("minimap2",
                       " -ax sr ",
                       TMP1,
                       " ",
                       TMP2,
                       " > ",
                       TMP3)
  }
  
  system(command = MINIMAP2,
         intern = FALSE)
  # run samtools to convert the samfile to a bamfile
  SAMTOOLS <- paste0("samtools view -S -b ",
                     TMP3,
                     " > ",
                     TMP4)
  system(command = SAMTOOLS,
         intern = FALSE)
  
  Res <- BamFile(TMP4)
  # return an r-object that can be targetted by scanBam()
  return(Res)
}



