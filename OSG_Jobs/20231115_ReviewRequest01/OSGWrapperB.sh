#!/bin/bash

Rscript JobScriptB.R ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10}

if [ -e ResultB*.fna.gz ]
then
  exit 0
else
  exit 1
fi