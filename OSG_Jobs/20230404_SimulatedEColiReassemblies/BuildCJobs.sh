#!/bin/bash
# create a jobmap for the B job node

VAL=1
# collect current completed jobs
shopt -s nullglob
for file in Assembly_*.fna.gz; do
  # test the that file is a regular file that is not zero bytes
  [[ -f $file && -s $file ]] && \
  # extract the original integer assignment with zero padding
  foo=$(echo "$file" | sed -e 's/[^0-9]//g') && \
  # write out the new line with the input file and the output file name
  printf '%s\n' "$VAL $file Annotation_$foo.gff" && \
  # increment the integer
  ((VAL++))
done > JobMapB.txt



