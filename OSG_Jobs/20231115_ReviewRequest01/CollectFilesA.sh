#!/bin/bash

shopt -s nullglob
for file in ResultC*.RData; do
  [[ -f $file && -s $file ]] && printf '%s\n' "$file"
done > CompleteRDataFiles.txt

