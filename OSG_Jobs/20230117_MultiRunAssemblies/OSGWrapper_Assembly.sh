#!/bin/bash

# running inside a singularity container, ENV commands in the dockefile aren't run, so need to export the PATH
# This path needs to match where the executable was installed in the dockerfile
# export PATH=/blast/ncbi-blast-2.9.0+/bin:$PATH
# export PATH=/hmmer/hmmer-3.3.1/bin:$PATH

mkdir ncbi
mkdir fastq
mkdir sra

export TMPDIR=$(pwd)
export NCBI_SETTINGS=user-settings.mkfg

Rscript JobScript_Assembly.R ${1} ${2} ${3} ${4}

# if [ -e Assembly_Unicycler*.fna.gz ]
# then
#   exit 0
# else
#   exit 1
# fi
