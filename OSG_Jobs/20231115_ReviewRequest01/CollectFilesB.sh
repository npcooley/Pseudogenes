#!/bin/bash

shopt -s nullglob
for file in ResultC*.fna.gz; do
  [[ -f $file && -s $file ]] && printf '%s\n' "$file"
done > CompleteAssemblyFiles.txt

