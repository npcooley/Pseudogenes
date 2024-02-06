#!/bin/bash

Rscript JobScriptD.R ${1} ${2} ${3} ${4} ${5}

if [ -e Result*.RData ]
then
  exit 0
else
  exit 1
fi
